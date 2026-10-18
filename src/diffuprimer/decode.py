from __future__ import annotations

import numpy as np

from .constraints import Conditions
from .sequence import PrimerPair


def interval_table(length, conditions, label):
    starts, ends = [], []
    for size in range(conditions.primer_min_length, conditions.primer_max_length + 1):
        start = np.arange(max(0, length - size + 1), dtype=np.int32)
        starts.append(start)
        ends.append(start + size)
    if not starts:
        return np.empty((0, 2), dtype=np.int32)
    return np.stack([np.concatenate(starts), np.concatenate(ends)], axis=1)


def _prefix(values):
    return np.concatenate([[0], np.cumsum(values)])


def grammar_candidates(logits, conditions: Conditions, known=None, fixed_forward=None,
                       fixed_reverse=None, beam=64, reverse_per_forward=8, previous=None,
                       span_allowed=None):
    """Beam approximation to a normalized two-interval x0 distribution.

    Candidate scores are sums of log p(label)-log p(background) over F/R.
    Every pair obeys geometry and all revealed labels. Thermodynamics is evaluated
    AFTER a fixed number of pair proposals, never used here to refill a budget.
    """
    logits = np.asarray(logits, dtype=np.float64)
    length = len(logits)
    if logits.shape != (length, 3) or not np.all(np.isfinite(logits)):
        raise ValueError("Expected finite length x 3 logits")
    known = np.full(length, 3, dtype=np.int64) if known is None else np.asarray(known)
    if known.shape != (length,) or np.any((known < 0) | (known > 3)):
        raise ValueError("Invalid revealed labels")
    intervals = interval_table(length, conditions, 1)
    if len(intervals) == 0:
        return [], np.empty(0)
    starts, ends = intervals.T
    energy_f = _prefix(logits[:, 1] - logits[:, 0])
    energy_r = _prefix(logits[:, 2] - logits[:, 0])
    scores_f = energy_f[ends] - energy_f[starts]
    scores_r = energy_r[ends] - energy_r[starts]
    prefix_known = [_prefix(known == label) for label in range(3)]
    totals = [int(np.sum(known == label)) for label in range(3)]

    def known_counts(start, end):
        return np.array([p[end] - p[start] for p in prefix_known])

    f_ok = np.ones(len(intervals), dtype=bool)
    r_ok = np.ones(len(intervals), dtype=bool)
    if span_allowed is not None:
        f_ok &= span_allowed[0]
        r_ok &= span_allowed[1]
    for label, flags in [(1, f_ok), (2, r_ok)]:
        # A span cannot cover revealed 0 or opposite-direction foreground.
        for other in {0, 1, 2} - {label}:
            flags &= (prefix_known[other][ends] - prefix_known[other][starts]) == 0
        # Since a complete assignment contains exactly one span per direction,
        # all known sites of that direction must lie inside it.
        flags &= (prefix_known[label][ends] - prefix_known[label][starts]) == totals[label]
    if fixed_forward is not None:
        f_ok &= (starts == fixed_forward[0]) & (ends == fixed_forward[1])
    if fixed_reverse is not None:
        r_ok &= (starts == fixed_reverse[0]) & (ends == fixed_reverse[1])
    f_indices = np.flatnonzero(f_ok)
    if len(f_indices) > beam:
        selected = np.argpartition(scores_f[f_indices], -beam)[-beam:]
        f_indices = f_indices[selected]
    pairs, scores = [], []
    for index in f_indices:
        fs, fe = map(int, intervals[index])
        r_indices = np.flatnonzero(r_ok & (starts >= fe)
            & (ends - fs >= conditions.amplicon_min) & (ends - fs <= conditions.amplicon_max))
        if len(r_indices) > reverse_per_forward:
            selected = np.argpartition(scores_r[r_indices], -reverse_per_forward)[-reverse_per_forward:]
            r_indices = r_indices[selected]
        for right in r_indices:
            rs, re = map(int, intervals[right])
            pairs.append(PrimerPair(fs, fe, rs, re))
            scores.append(float(scores_f[index] + scores_r[right]))
    if previous is not None:
        previous.validate("N" * length)
        clean = previous.labels(length)
        compatible = np.all((known == 3) | (known == clean))
        in_range = (conditions.amplicon_min <= previous.amplicon_length <= conditions.amplicon_max
                    and conditions.primer_min_length <= previous.forward_end - previous.forward_start <= conditions.primer_max_length
                    and conditions.primer_min_length <= previous.reverse_end - previous.reverse_start <= conditions.primer_max_length)
        anchors_ok = ((fixed_forward is None or (previous.forward_start, previous.forward_end) == tuple(fixed_forward))
                      and (fixed_reverse is None or (previous.reverse_start, previous.reverse_end) == tuple(fixed_reverse)))
        if compatible and in_range and anchors_ok and previous not in pairs:
            pairs.append(previous)
            scores.append(float(energy_f[previous.forward_end] - energy_f[previous.forward_start]
                                + energy_r[previous.reverse_end] - energy_r[previous.reverse_start]))
    return pairs, np.asarray(scores)


def sample_grammar(logits, conditions, rng, temperature=1.0, **kwargs):
    if temperature <= 0:
        raise ValueError("Temperature must be positive")
    pairs, scores = grammar_candidates(logits, conditions, **kwargs)
    if not pairs:
        return None
    probabilities = np.exp((scores - scores.max()) / temperature)
    probabilities /= probabilities.sum()
    return pairs[int(rng.choice(len(pairs), p=probabilities))]
