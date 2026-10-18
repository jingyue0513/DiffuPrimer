from __future__ import annotations

import primer3

from .constraints import Conditions
from .sequence import PrimerPair


def primer3_candidates(template, conditions: Conditions, budget=20, fixed_forward=None, fixed_reverse=None):
    args = {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": template}
    if fixed_forward is not None:
        args["SEQUENCE_PRIMER"] = fixed_forward
    if fixed_reverse is not None:
        args["SEQUENCE_PRIMER_REVCOMP"] = fixed_reverse
    settings = {
        "PRIMER_TASK": "generic", "PRIMER_NUM_RETURN": budget,
        "PRIMER_MIN_SIZE": conditions.primer_min_length,
        "PRIMER_OPT_SIZE": 21, "PRIMER_MAX_SIZE": conditions.primer_max_length,
        "PRIMER_MIN_TM": conditions.tm_min, "PRIMER_OPT_TM": conditions.tm_opt,
        "PRIMER_MAX_TM": conditions.tm_max, "PRIMER_PAIR_MAX_DIFF_TM": conditions.max_tm_difference,
        "PRIMER_MIN_GC": conditions.gc_min, "PRIMER_MAX_GC": conditions.gc_max,
        "PRIMER_MAX_POLY_X": conditions.max_homopolymer,
        "PRIMER_PRODUCT_SIZE_RANGE": [[conditions.amplicon_min, conditions.amplicon_max]],
        "PRIMER_SALT_MONOVALENT": conditions.monovalent_mM,
        "PRIMER_SALT_DIVALENT": conditions.divalent_mM,
        "PRIMER_DNTP_CONC": conditions.dntp_mM, "PRIMER_DNA_CONC": conditions.primer_nM,
        "PRIMER_MAX_HAIRPIN_TH": conditions.max_structure_tm,
        "PRIMER_MAX_SELF_ANY_TH": conditions.max_structure_tm,
        "PRIMER_MAX_SELF_END_TH": conditions.max_structure_tm,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": conditions.max_structure_tm,
        "PRIMER_PAIR_MAX_COMPL_END_TH": conditions.max_structure_tm,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1, "PRIMER_EXPLAIN_FLAG": 1,
    }
    result = primer3.bindings.design_primers(args, settings)
    candidates = []
    for i in range(result.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        left, llen = result[f"PRIMER_LEFT_{i}"]
        # Primer3 reports the RIGHT primer's 5' template coordinate and length.
        right, rlen = result[f"PRIMER_RIGHT_{i}"]
        pair = PrimerPair(left, left + llen, right - rlen + 1, right + 1)
        assert pair.sequences(template) == (result[f"PRIMER_LEFT_{i}_SEQUENCE"], result[f"PRIMER_RIGHT_{i}_SEQUENCE"])
        candidates.append({**pair.as_dict(template), "primer3_penalty": float(result[f"PRIMER_PAIR_{i}_PENALTY"])})
    explain = {k: v for k, v in result.items() if k.endswith("EXPLAIN") or k in ["PRIMER_ERROR", "PRIMER_WARNING"]}
    return candidates, explain
