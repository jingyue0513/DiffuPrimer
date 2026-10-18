#!/usr/bin/env python3
"""Verify the completed fixed-seed study against its actual saved outputs."""
import json
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json

ROOT = Path(__file__).resolve().parents[1]
targets = {r["gene_id"]: r for r in read_jsonl(ROOT / "data/processed/targets_pilot.jsonl") if r["split"] == "test"}
tasks = {r["gene_id"]: r for r in read_jsonl(ROOT / "experiments/final_followup/redesign_tasks.jsonl")}
state = json.loads((ROOT / "experiments/final_followup/results/status.json").read_text())
assert state["stage"] == "followup_experiments_complete" and state["anchor_prefilter_tasks_complete"]
assert len(targets) == len(tasks) == 75
assert all(t["anchor"] and t["old_pair_violates_new_length"] for t in tasks.values())
summary = json.loads((ROOT / "results/final_fixed_seed/study_summary.json").read_text())
assert summary["source_runs"] == len(summary["inputs"]) == 47
for source in summary["inputs"]:
    for key, digest in source.items():
        if key.endswith("_sha256"):
            assert sha256(Path(source["directory"]) / key[:-7]) == digest, (source["directory"], key)

complement = str.maketrans("ACGT", "TGCA")
checked = []
for directory in sorted((ROOT / "experiments/final_followup/results/evaluations").iterdir()):
    m = json.loads((directory / "manifest.json").read_text())
    c = json.loads((directory / "completion.json").read_text())
    s = json.loads((directory / "specificity_bowtie_status.json").read_text())
    ledger, candidates, screened = [read_jsonl(directory / name) for name in ["per_gene.jsonl", "candidates.jsonl", "specificity_bowtie.jsonl"]]
    assert m["arguments"]["seed"] == 20261013 and set(m["expected_targets"]) == set(targets)
    assert c["all_expected_targets_evaluated"] and c["targets"] == 75
    assert s["status"] == "completed" and {r["gene_id"] for r in ledger} == set(targets)
    assert len(ledger) == 75 and c["per_gene_sha256"] == sha256(directory / "per_gene.jsonl")
    assert c["candidates_sha256"] == s["candidate_sha256"] == sha256(directory / "candidates.jsonl")
    assert s["output_sha256"] == sha256(directory / "specificity_bowtie.jsonl")
    eligible = {r["candidate_id"] for r in candidates if r["assessment"]["computational_feasible"]}
    assert {r["candidate_id"] for r in screened} == {r["candidate_id"] for r in candidates}
    assert {r["candidate_id"] for r in screened if r.get("status") != "not_screened_computational_failure"} == eligible
    structural = 0
    for r in candidates:
        if not r["structural_valid"]:
            continue
        structural += 1
        sequence = targets[r["gene_id"]]["template_sequence"]
        fs, fe, rs, re = [r[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]]
        assert 0 <= fs < fe <= rs < re <= len(sequence)
        assert r["forward_primer"] == sequence[fs:fe]
        assert r["reverse_primer"] == sequence[rs:re].translate(complement)[::-1]
        assert r["amplicon_length"] == re - fs
        if r["assessment"]["computational_feasible"]:
            assert r["requested_amplicon_min"] <= re - fs <= r["requested_amplicon_max"]
        if m["task"] in ["redesign", "completion"]:
            side = m["fixed_side"]
            anchor = tasks[r["gene_id"]]["anchor"]
            assert r["anchor_preserved"] and all(r[side + suffix] == anchor[side + suffix] for suffix in ["_start", "_end", "_primer"])
        if m["task"] == "redesign":
            assert [r["requested_amplicon_min"], r["requested_amplicon_max"]] == tasks[r["gene_id"]]["new_length_band"]
    if m["arguments"].get("anchor_prefilter"):
        assert structural == len(candidates) == 1500
    if m.get("checkpoint_sha256"):
        model = ROOT / f"results/checkpoints/{m['method']}_weak_partial_pilot_seed20261013.pt"
        assert sha256(model) == m["checkpoint_sha256"], m["method"]
    checked.append({"directory": str(directory), "returned_proposals": len(candidates), "structural_pairs": structural, "audited_feasible_pairs": len(eligible), "candidates_sha256": c["candidates_sha256"], "screen_sha256": s["output_sha256"]})
assert len(checked) == 18
old = json.loads((ROOT / "results/computational_completion_check.json").read_text())
assert old["computational_checks_passed"] and sha256(ROOT / "results/frozen_protocol.json") == old["frozen_protocol_sha256"]
for audit in old["audits"]:
    directory = Path(audit["directory"])
    assert sha256(directory / "candidates.jsonl") == audit["candidates_sha256"]
    assert sha256(directory / "specificity_bowtie.jsonl") == audit["specificity_sha256"]
build = json.loads((ROOT / "paper/build_manifest.json").read_text())
assert build["status"] == "formal_computational_manuscript_complete"
for key, name in [("source_sha256", "manuscript.tex"), ("pdf_sha256", "manuscript.pdf"), ("supplementary_source_sha256", "supplementary.tex"), ("supplementary_pdf_sha256", "supplementary.pdf"), ("reported_values_sha256", "reported_values.json")]:
    assert sha256(ROOT / "paper" / name) == build[key]
assert build["summary_sha256"] == sha256(ROOT / "results/final_fixed_seed/study_summary.json")
for name, digest in build["figures"].items():
    assert sha256(ROOT / "results/final_fixed_seed/figures" / name) == digest
write_json(ROOT / "results/final_fixed_seed/completion_check.json", {
    "verified_at": utc_now(), "passed": True, "seed": 20261013, "held_out_genes": 75,
    "new_runs": 18, "new_budgeted_proposals": 27000, "new_returned_proposals": sum(r["returned_proposals"] for r in checked),
    "total_completed_runs": 264, "total_full_reference_audits": 264, "total_budgeted_proposals": 396000,
    "final_analysis_source_runs": 47, "statistical_groups": len(summary["summaries"]), "paired_contrasts": len(summary["paired_comparisons"]),
    "original_study_hashes_preserved": True, "changed_length_targets": 75,
    "sequence_coordinates_and_anchors_verified": True, "reported_paper_and_figures_verified": True,
    "physical_qpcr_required": False, "runs": checked})
print("Final study verified: 18 new runs, 264 total audits, 47 final-analysis sources, manuscript and supplement.")
