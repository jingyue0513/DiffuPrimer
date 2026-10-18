from __future__ import annotations

import re

import numpy as np
import primer3

from .constraints import oligo_properties
from .decode import interval_table
from .sequence import PrimerPair, reverse_complement


def feasible_oligo_spans(template, conditions):
    """Sequence-only hard single-oligo constraints, shared by learned and rule samplers."""
    intervals = interval_table(len(template), conditions, 1)
    allowed = [np.zeros(len(intervals), dtype=bool), np.zeros(len(intervals), dtype=bool)]
    for index, (start, end) in enumerate(intervals):
        sequence = template[start:end]
        if set(sequence) - set("ACGT"):
            continue
        gc = 100 * (sequence.count("G") + sequence.count("C")) / len(sequence)
        if not conditions.gc_min <= gc <= conditions.gc_max:
            continue
        if max(len(run) for run in re.findall(r"A+|C+|G+|T+", sequence)) > conditions.max_homopolymer:
            continue
        tm = primer3.calc_tm(sequence, **conditions.thermo_kwargs())
        if not conditions.tm_min <= tm <= conditions.tm_max:
            continue
        for side, primer in enumerate((sequence, reverse_complement(sequence))):
            properties = oligo_properties(primer, conditions)
            allowed[side][index] = (conditions.tm_min <= properties["tm"] <= conditions.tm_max
                and properties["hairpin_tm"] <= conditions.max_structure_tm
                and properties["homodimer_tm"] <= conditions.max_structure_tm)
    return intervals, allowed


def compatible_span_temperatures(template, conditions, spans):
    """Temperatures of allowed spans, in the same interval-table ordering."""
    intervals, allowed = spans
    temperatures = np.full((2, len(intervals)), np.nan)
    for side in [0, 1]:
        for index in np.flatnonzero(allowed[side]):
            start, end = intervals[index]
            sequence = template[start:end]
            if side == 1:
                sequence = reverse_complement(sequence)
            temperatures[side, index] = oligo_properties(sequence, conditions)["tm"]
    return temperatures


def rule_span_sampler(template, conditions, budget, seed=0, spans=None, fixed_forward=None, fixed_reverse=None, span_temperatures=None):
    """Uniform F span followed by uniform compatible R span; not uniform over all pairs."""
    rng = np.random.default_rng(seed)
    intervals, allowed = feasible_oligo_spans(template, conditions) if spans is None else spans
    starts, ends = intervals.T if len(intervals) else (np.empty(0), np.empty(0))
    f_ok, r_ok = allowed[0].copy(), allowed[1].copy()
    if fixed_forward is not None:
        f_ok &= (starts == fixed_forward[0]) & (ends == fixed_forward[1])
    if fixed_reverse is not None:
        r_ok &= (starts == fixed_reverse[0]) & (ends == fixed_reverse[1])
    forward = np.flatnonzero(f_ok)
    pairs = []
    for _ in range(budget):
        if not len(forward):
            pairs.append(None)
            continue
        left = int(rng.choice(forward))
        fs, fe = map(int, intervals[left])
        compatible = r_ok & (starts >= fe) & (ends - fs >= conditions.amplicon_min) & (ends - fs <= conditions.amplicon_max)
        if span_temperatures is not None:
            compatible &= np.abs(span_temperatures[1] - span_temperatures[0, left]) <= conditions.max_tm_difference
        reverse = np.flatnonzero(compatible)
        if not len(reverse):
            pairs.append(None)
            continue
        rs, re = map(int, intervals[int(rng.choice(reverse))])
        pairs.append(PrimerPair(fs, fe, rs, re))
    return pairs
