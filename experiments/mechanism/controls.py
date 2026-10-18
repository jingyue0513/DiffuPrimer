"""Mechanism controls; all random choices have explicit, non-clock seeds."""
from __future__ import annotations

import time
import numpy as np
import torch

from diffuprimer.cnn_model import MASK, encode_conditions, encode_dna
from diffuprimer.decode import interval_table, grammar_candidates, sample_grammar
from diffuprimer.sequence import PrimerPair, labels_to_pair
from diffuprimer.constraints import assess_pair


@torch.inference_mode()
def iterative(model, template, conditions, allowed, temperatures, seed, steps,
              frozen=False, budget=20, temperature=4.0):
    """Identical state/reveal/beam code for full and frozen-first-score sampling."""
    model.eval()
    rng = np.random.default_rng(seed)
    dna = torch.tensor(encode_dna(template))[None].expand(budget, -1)
    c = torch.tensor(encode_conditions(conditions))[None].expand(budget, -1)
    noisy = np.full((budget, len(template)), MASK, dtype=np.int64)
    previous, failures = [None] * budget, [False] * budget
    logits, calls = None, 0
    for t in range(steps, 0, -1):
        if logits is None or not frozen:
            logits = model(dna, torch.tensor(noisy), torch.full((budget,), t / steps), c).numpy()
            calls += 1
        for i in range(budget):
            if failures[i]:
                continue
            pair = sample_grammar(logits[i], conditions, rng, temperature=temperature,
                                  known=noisy[i], previous=previous[i],
                                  span_allowed=allowed, span_temperatures=temperatures)
            if pair is None:
                failures[i] = True
                continue
            previous[i] = pair
            clean = pair.labels(len(template))
            reveal = (noisy[i] == MASK) & (rng.random(len(template)) < 1 / t)
            noisy[i, reveal] = clean[reveal]
    return [None if failures[i] else labels_to_pair(noisy[i]) for i in range(budget)], {
        "network_calls": calls, "network_examples": calls * budget,
        "decoder_calls": steps * budget, "steps": steps, "frozen_scores": frozen}


@torch.inference_mode()
def first_scores(model, template, conditions):
    model.eval()
    dna = torch.tensor(encode_dna(template))[None]
    c = torch.tensor(encode_conditions(conditions))[None]
    return model(dna, torch.full_like(dna, MASK), torch.ones(1), c)[0].numpy()


def all_support(logits, conditions, allowed, temperatures):
    """Enumerate EVERY geometry/single-oligo/joint-Tm-compatible interval pair.

    Reverse intervals are indexed by their end to avoid a quadratic transcript
    scan. No beam, score cutoff, junction annotation or specificity is used.
    Heterodimers remain in the shared final pair assessment, as in the original.
    """
    intervals = interval_table(len(logits), conditions, 1)
    starts, ends = intervals.T
    pf = np.r_[0., np.cumsum(np.asarray(logits[:, 1] - logits[:, 0], dtype=np.float64))]
    pr = np.r_[0., np.cumsum(np.asarray(logits[:, 2] - logits[:, 0], dtype=np.float64))]
    sf, sr = pf[ends] - pf[starts], pr[ends] - pr[starts]
    rights = np.flatnonzero(allowed[1])
    rights = rights[np.argsort(ends[rights], kind="stable")]
    rend = ends[rights]
    coordinates, scores = [], []
    for left in np.flatnonzero(allowed[0]):
        lo = np.searchsorted(rend, starts[left] + conditions.amplicon_min, side="left")
        hi = np.searchsorted(rend, starts[left] + conditions.amplicon_max, side="right")
        r = rights[lo:hi]
        ok = starts[r] >= ends[left]
        ok &= np.abs(temperatures[1, r] - temperatures[0, left]) <= conditions.max_tm_difference
        r = r[ok]
        if len(r):
            coordinates.append(np.column_stack((np.full(len(r), starts[left]),
                np.full(len(r), ends[left]), starts[r], ends[r])))
            scores.append(sf[left] + sr[r])
    if not coordinates:
        return np.empty((0, 4), dtype=np.int32), np.empty(0)
    return np.concatenate(coordinates).astype(np.int32), np.concatenate(scores)


def weighted_order(scores, seed, count, temperature=4.0):
    """Gumbel top-k: Plackett-Luce sampling without replacement, no underflow."""
    if not len(scores):
        return np.empty(0, dtype=int)
    keys = (np.asarray(scores) - np.max(scores)) / temperature
    keys = keys + np.random.default_rng(seed).gumbel(size=len(scores))
    count = min(count, len(scores))
    chosen = np.argpartition(keys, len(keys) - count)[-count:] if count < len(keys) else np.arange(len(keys))
    return chosen[np.argsort(-keys[chosen], kind="stable")]


def pairs_from_order(coordinates, order, template, count):
    pairs, seen = [], set()
    for i in order:
        pair = PrimerPair(*map(int, coordinates[i]))
        sequence_key = pair.sequences(template)
        if sequence_key not in seen:
            pairs.append(pair)
            seen.add(sequence_key)
        if len(pairs) == count:
            break
    return pairs


def diverse_indices(pairs, assessments, count=20):
    """Feasible-only farthest-first amplicon selection; Tm then pool order ties.

    First candidate has the smallest shared Tm penalty. All subsequent choices
    maximize minimum Jaccard distance to the selected set. No reference screen.
    """
    eligible = np.array([i for i, a in enumerate(assessments) if a["computational_feasible"]], dtype=int)
    if not len(eligible):
        return []
    fs = np.array([pairs[i].forward_start for i in eligible])
    re = np.array([pairs[i].reverse_end for i in eligible])
    penalties = np.array([assessments[i]["shared_temperature_penalty"] for i in eligible])
    selected, available, mindist = [], np.ones(len(eligible), bool), np.ones(len(eligible))
    for _ in range(min(count, len(eligible))):
        candidates = np.flatnonzero(available)
        rank = np.lexsort((eligible[candidates], penalties[candidates], -mindist[candidates]))
        chosen = candidates[rank[0]]
        selected.append(int(eligible[chosen]))
        available[chosen] = False
        intersection = np.maximum(0, np.minimum(re, re[chosen]) - np.maximum(fs, fs[chosen]))
        union = re - fs + re[chosen] - fs[chosen] - intersection
        mindist = np.minimum(mindist, 1 - intersection / union)
    return selected


def strong(model, template, conditions, allowed, temperatures, seed,
           pool_size=20, diverse=False):
    logits = first_scores(model, template, conditions)
    coordinates, scores = all_support(logits, conditions, allowed, temperatures)
    # Draw extra coordinates only to allow synthesis-sequence deduplication.
    order = weighted_order(scores, seed, len(scores))
    pairs = pairs_from_order(coordinates, order, template, pool_size)
    assessments = [assess_pair(template, p, conditions) for p in pairs]
    chosen = diverse_indices(pairs, assessments) if diverse else list(range(min(20, len(pairs))))
    return [pairs[i] for i in chosen], [assessments[i] for i in chosen], {
        "network_calls": 1, "network_examples": 1, "support_size": len(scores),
        "pool_size": len(pairs), "selection": "farthest_first_tm_tiebreak" if diverse else "weighted_without_replacement",
        "all_support": True}


def timed_direct(model, template, conditions, allowed, temperatures, seed,
                 cpu_budget, max_pool=4096, clock=time.process_time):
    """Anytime direct search: only checkpoints completed by the CPU deadline count.

    The reference budget is full diffusion's network+search+assessment CPU cost.
    Shared span preparation is outside BOTH budgets. Full-reference screens are
    never used for selection. Record actual overshoot separately; do not pretend
    a Python operation can be preempted. Candidate sequences are saved for replay.
    """
    start = clock()
    trace, accepted_pairs, accepted_assessments = [], [], []
    logits = first_scores(model, template, conditions)
    seen, pool, assessments = set(), [], []

    def checkpoint(stage):
        nonlocal accepted_pairs, accepted_assessments
        chosen = diverse_indices(pool, assessments)
        elapsed = clock() - start
        on_time = elapsed <= cpu_budget
        trace.append({"stage": stage, "cpu_seconds": elapsed, "pool_size": len(pool),
                      "selected": len(chosen), "accepted_before_deadline": on_time})
        if on_time:
            accepted_pairs = [pool[i] for i in chosen]
            accepted_assessments = [assessments[i] for i in chosen]
        return on_time

    # A cheap beam result is available before the exhaustive support is built.
    bp, bs = grammar_candidates(logits, conditions, span_allowed=allowed, span_temperatures=temperatures)
    for i in weighted_order(bs, seed, len(bs)):
        pair = bp[i]
        key = pair.sequences(template)
        if key not in seen:
            seen.add(key); pool.append(pair); assessments.append(assess_pair(template, pair, conditions))
        if len(pool) >= 20 or clock() - start >= cpu_budget:
            break
    checkpoint("beam20")
    support_size = None
    if clock() - start < cpu_budget:
        coordinates, scores = all_support(logits, conditions, allowed, temperatures)
        support_size = len(scores)
        if clock() - start < cpu_budget:
            order = weighted_order(scores, seed, len(scores))
            thresholds = set([40, 80, 160, 320, 640, 1280, 2560, max_pool])
            for i in order:
                if clock() - start >= cpu_budget:
                    break
                pair = PrimerPair(*map(int, coordinates[i]))
                key = pair.sequences(template)
                if key in seen:
                    continue
                seen.add(key); pool.append(pair); assessments.append(assess_pair(template, pair, conditions))
                if len(pool) in thresholds and not checkpoint("full_support"):
                    break
                if len(pool) >= max_pool:
                    break
            if clock() - start < cpu_budget:
                checkpoint("final")
    elapsed = clock() - start
    accepted_times = [r["cpu_seconds"] for r in trace if r["accepted_before_deadline"]]
    return accepted_pairs, accepted_assessments, {
        "network_calls": 1, "network_examples": 1, "support_size": support_size,
        "pool_size": len(pool), "cpu_budget_seconds": cpu_budget, "trace": trace,
        "accepted_checkpoint_cpu_seconds": max(accepted_times, default=0.),
        "actual_cpu_seconds": elapsed, "overshoot_cpu_seconds": max(0., elapsed - cpu_budget),
        "selection": "deadline_committed_farthest_first_tm_tiebreak", "max_pool": max_pool}
