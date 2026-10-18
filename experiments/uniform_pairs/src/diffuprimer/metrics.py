from __future__ import annotations

import numpy as np


def paired_bootstrap_difference(left, right, seed=0, repeats=10000):
    """Resample target genes as paired units, not correlated primer pairs."""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.shape != right.shape or left.ndim != 1 or not len(left):
        raise ValueError("Require nonempty equally sized paired target arrays")
    finite = np.isfinite(left) & np.isfinite(right)
    differences = left[finite] - right[finite]
    if not len(differences):
        return {"n": 0, "mean_difference": None, "ci95": None}
    rng = np.random.default_rng(seed)
    resampled = []
    for offset in range(0, repeats, 1000):
        count = min(1000, repeats - offset)
        indices = rng.integers(0, len(differences), (count, len(differences)))
        resampled.extend(differences[indices].mean(1))
    return {"n": len(differences), "mean_difference": float(differences.mean()),
            "ci95": np.quantile(resampled, [0.025, 0.975]).tolist(), "bootstrap_repeats": repeats,
            "unit": "target_gene", "excluded_nonfinite_pairs": int((~finite).sum())}


def summarize_candidates(candidates, budget):
    if len(candidates) > budget or budget < 1:
        raise ValueError("Invalid candidate count or budget")
    feasible = [row for row in candidates if row.get("assessment", {}).get("computational_feasible", False)]
    unique = {(row["forward_primer"], row["reverse_primer"]) for row in feasible}
    positions = {(row["forward_start"], row["reverse_end"]) for row in feasible}
    scores = [row["assessment"]["shared_temperature_penalty"] for row in feasible]
    scores = [score for score in scores if score is not None]
    intervals = sorted(positions)
    distances = []
    for index, (left, right) in enumerate(intervals):
        for other_left, other_right in intervals[:index]:
            intersection = max(0, min(right, other_right) - max(left, other_left))
            union = (right - left) + (other_right - other_left) - intersection
            distances.append(1 - intersection / union)
    return {
        "candidate_budget": budget, "proposals_returned": len(candidates),
        "structural_valid_count": sum(row.get("structural_valid", True) for row in candidates),
        "computational_feasible_count": len(feasible), "unique_feasible_pairs": len(unique),
        "distinct_feasible_amplicons": len(positions), "target_covered": bool(unique),
        "mean_feasible_amplicon_jaccard_distance": float(np.mean(distances)) if distances else None,
        # Missing proposals count as failures; duplicates are not replaced.
        "budget_feasibility_rate": len(feasible) / budget,
        "returned_feasibility_rate": len(feasible) / len(candidates) if candidates else None,
        "best_shared_temperature_penalty": min(scores) if scores else None,
        "specificity_status": "not_assessed", "experimental_status": "not_tested"}
