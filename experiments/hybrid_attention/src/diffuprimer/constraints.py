from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import primer3

from .sequence import PrimerPair


@dataclass(frozen=True)
class Conditions:
    primer_min_length: int = 18
    primer_max_length: int = 25
    amplicon_min: int = 80
    amplicon_max: int = 200
    tm_min: float = 57.0
    tm_max: float = 63.0
    tm_opt: float = 60.0
    max_tm_difference: float = 2.0
    gc_min: float = 35.0
    gc_max: float = 65.0
    max_homopolymer: int = 5
    max_structure_tm: float = 47.0
    monovalent_mM: float = 50.0
    divalent_mM: float = 1.5
    dntp_mM: float = 0.6
    primer_nM: float = 50.0

    def __post_init__(self):
        if not (0 < self.primer_min_length <= self.primer_max_length <= 60):
            raise ValueError("Invalid primer length range")
        if not (2 * self.primer_min_length <= self.amplicon_min <= self.amplicon_max):
            raise ValueError("Invalid amplicon range")
        if not self.tm_min <= self.tm_opt <= self.tm_max:
            raise ValueError("Invalid melting temperature range")

    def thermo_kwargs(self):
        return dict(mv_conc=self.monovalent_mM, dv_conc=self.divalent_mM,
                    dntp_conc=self.dntp_mM, dna_conc=self.primer_nM)


@lru_cache(maxsize=100000)
def oligo_properties(sequence, conditions):
    if not sequence or set(sequence) - set("ACGT"):
        return {"valid_sequence": False}
    kwargs = conditions.thermo_kwargs()
    longest = max(len(run) for run in __import__("re").findall(r"A+|C+|G+|T+", sequence))
    return {
        "valid_sequence": True, "length": len(sequence),
        "gc_percent": 100 * (sequence.count("G") + sequence.count("C")) / len(sequence),
        "tm": float(primer3.calc_tm(sequence, **kwargs)), "max_homopolymer": longest,
        "hairpin_tm": float(primer3.calc_hairpin(sequence, **kwargs).tm),
        "homodimer_tm": float(primer3.calc_homodimer(sequence, **kwargs).tm),
    }


def assess_pair(template, pair: PrimerPair, conditions: Conditions):
    failures = []
    try:
        pair.validate(template)
    except ValueError:
        return {"geometry_pass": False, "thermodynamic_pass": False,
                "computational_feasible": False, "failures": ["invalid_geometry"],
                "specificity_status": "not_assessed", "experimental_status": "not_tested"}
    forward, reverse = pair.sequences(template)
    if not conditions.amplicon_min <= pair.amplicon_length <= conditions.amplicon_max:
        failures.append("amplicon_length")
    properties = [oligo_properties(primer, conditions) for primer in (forward, reverse)]
    geometry_pass = not failures and all(conditions.primer_min_length <= len(p) <= conditions.primer_max_length for p in (forward, reverse))
    for name, prop in zip(("forward", "reverse"), properties):
        if not prop["valid_sequence"]:
            failures.append(f"{name}_ambiguous_sequence")
            continue
        for key, low, high in [("length", conditions.primer_min_length, conditions.primer_max_length),
                               ("tm", conditions.tm_min, conditions.tm_max),
                               ("gc_percent", conditions.gc_min, conditions.gc_max)]:
            if not low <= prop[key] <= high:
                failures.append(f"{name}_{key}")
        for key, limit in [("max_homopolymer", conditions.max_homopolymer),
                           ("hairpin_tm", conditions.max_structure_tm),
                           ("homodimer_tm", conditions.max_structure_tm)]:
            if prop[key] > limit:
                failures.append(f"{name}_{key}")
    delta = heterodimer = score = None
    if all(p["valid_sequence"] for p in properties):
        delta = abs(properties[0]["tm"] - properties[1]["tm"])
        heterodimer = float(primer3.calc_heterodimer(forward, reverse, **conditions.thermo_kwargs()).tm)
        if delta > conditions.max_tm_difference:
            failures.append("tm_difference")
        if heterodimer > conditions.max_structure_tm:
            failures.append("heterodimer_tm")
        # A shared transparent temperature score; not an experimental success probability.
        score = sum(abs(p["tm"] - conditions.tm_opt) for p in properties) + delta
    return {"geometry_pass": geometry_pass, "thermodynamic_pass": not failures,
            "computational_feasible": not failures, "failures": failures,
            "forward_properties": properties[0], "reverse_properties": properties[1],
            "tm_difference": delta, "heterodimer_tm": heterodimer,
            "shared_temperature_penalty": score,
            "specificity_status": "not_assessed", "experimental_status": "not_tested"}
