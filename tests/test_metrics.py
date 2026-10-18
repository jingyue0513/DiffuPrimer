from diffuprimer.metrics import paired_bootstrap_difference, summarize_candidates


def test_empty_and_duplicate_proposals_have_budget_cost():
    row = {"forward_primer": "ACGT", "reverse_primer": "TGCA", "forward_start": 0,
           "reverse_end": 80, "assessment": {"computational_feasible": True, "shared_temperature_penalty": 1.5}}
    report = summarize_candidates([row, dict(row)], 5)
    assert report["unique_feasible_pairs"] == 1
    assert report["budget_feasibility_rate"] == 2 / 5
    assert report["proposals_returned"] == 2


def test_bootstrap_resamples_paired_genes():
    result = paired_bootstrap_difference([1, 2, 3], [0, 1, 2], repeats=1000)
    assert result["ci95"] == [1, 1]
    assert result["mean_difference"] == 1
    assert result["unit"] == "target_gene"
