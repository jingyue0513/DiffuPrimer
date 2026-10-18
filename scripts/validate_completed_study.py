#!/usr/bin/env python3
"""Verify final computational artifacts against their recorded source hashes."""
import json
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]


def load(relative):
    return json.loads((root / relative).read_text())


if load("results/computational_artifact_status.json")["stage"] != "computational_artifacts_complete_independent_assays_pending":
    raise ValueError("Final artifact driver has not completed")
frozen = load("results/frozen_protocol.json")
for relative, digest in frozen["sha256"].items():
    if sha256(root / relative) != digest:
        raise ValueError("Frozen protocol changed: " + relative)
reference = load("data/processed/reference_manifest.json")
for record in reference["files"]:
    if sha256(root / "data/raw/refseq" / record["name"]) != record["sha256"]:
        raise ValueError("Reference sequence or annotation changed")
spliced = load("results/spliced_reference_verification.json")
if spliced["targets"] != 500 or spliced["total_differing_bases"] != 0 or spliced["reference_manifest_sha256"] != sha256(root / "data/processed/blastdb/genome.manifest.json"):
    raise ValueError("Spliced reference verification differs")
if spliced["target_sha256"] != sha256(root / "data/processed/targets_pilot.jsonl") or spliced["exon_context_sha256"] != sha256(root / "results/frozen_target_exon_context.json"):
    raise ValueError("Verified target or exon geometry changed")
integrity = load("results/checks/candidate_integrity/check.json")
if integrity["completed_runs"] != 246 or not all(integrity[k] for k in ["candidate_sequence_geometry_checks_passed",
        "complete_frozen_cohorts", "candidate_budget_accounting_passed", "fixed_primer_sequence_preservation_passed"]):
    raise ValueError("Incomplete candidate integrity checks")
audits = []
for run in integrity["runs"]:
    directory = Path(run["directory"])
    if sha256(directory / "candidates.jsonl") != run["candidates_sha256"] or sha256(directory / "completion.json") != run["completion_sha256"]:
        raise ValueError("Checked candidate output changed")
    status = json.loads((directory / "specificity_bowtie_status.json").read_text())
    screen = directory / "specificity_bowtie.jsonl"
    if status["status"] != "completed" or status["candidate_sha256"] != run["candidates_sha256"] or status["output_sha256"] != sha256(screen):
        raise ValueError("Incomplete full reference audit: " + str(directory))
    audits.append({"directory": str(directory), "candidates_sha256": run["candidates_sha256"],
                   "specificity_sha256": status["output_sha256"]})

analyses = ["results/temperature_tuning/selected_test_summary.json",
    "results/temperature_tuning/selected_conditions_summary.json", "results/guide_ablation_summary.json"]
for relative, expected in zip(analyses, [13, 87, 19]):
    analysis = load(relative)
    if analysis.get("RNA_only") or analysis["screen"] != "bowtie" or len(analysis["inputs"]) != expected:
        raise ValueError("Incomplete final comparison matrix")
    for item in analysis["inputs"]:
        directory = Path(item["directory"])
        for filename, key in [("completion.json", "completion_sha256"), ("specificity_bowtie.jsonl", "specificity_sha256")]:
            if sha256(directory / filename) != item[key]:
                raise ValueError("Statistical source changed")
    for row in read_jsonl((root / relative).with_suffix(".per_gene.jsonl")):
        expected_seeds = [20261013] if row["method"] == "primer3" else [20261013, 20261014, 20261015]
        if row["seeds"] != expected_seeds:
            raise ValueError("Unexpected final comparison seed denominator")

duration = load("experiments/epoch100/results/duration_reference_comparison.json")
if len(duration["summaries"]) != 8 or len(duration["paired_comparisons"]) != 40:
    raise ValueError("Incomplete duration comparison")
for item in duration["inputs"]:
    path = Path(item["summary"])
    if sha256(path) != item["summary_sha256"] or sha256(path.with_suffix(".per_gene.jsonl")) != item["per_gene_sha256"]:
        raise ValueError("Duration comparison source changed")
training = load("results/figures/training/training_checks.json")
if len(training["runs"]) != 6 or not training["losses_finite"] or not training["hundred_epoch_prefix_equality"] or any(r["completed_epochs"] != 100 for r in training["runs"]):
    raise ValueError("Incomplete long training trajectories")
selection = load("results/temperature_tuning/full_audit_consistency.json")
if not selection["exact_RNA_decision_and_subject_role_agreement"] or not selection["exact_selection_endpoint_agreement"]:
    raise ValueError("Model-selection reference audits disagree")
if selection["selection_sha256"] != sha256(root / "results/temperature_tuning/selection.json"):
    raise ValueError("Selected temperature changed")

runtime = load("results/runtime_benchmark/summary.json")
timings = read_jsonl(root / "results/runtime_benchmark/per_gene.jsonl")
if runtime["targets"] != 75 or runtime["repetitions"] != 3 or not runtime["all_candidate_coordinates_match"] or len(timings) != 1125:
    raise ValueError("Incomplete measured runtime experiment")
if len({(r["method"], r["gene_id"], r["repetition"]) for r in timings}) != 1125 or not all(r["candidate_coordinates_exactly_match_recorded_run"] for r in timings):
    raise ValueError("Runtime cohort or candidate reproduction differs")
figures = load("results/figures/selected_study/provenance.json")
for path, digest in {**figures["source_summaries"], **figures["source_records"]}.items():
    if sha256(Path(path)) != digest:
        raise ValueError("Figure source changed")
figure_names = ["primary_RNA_endpoints", "genomic_product_screen", "exon_geometry_partition", "proposal_budget", "condition_and_completion",
                "denoising_steps", "guidance_mechanism", "training_duration_endpoints", "generation_runtime"]
for name in figure_names:
    for suffix in ["pdf", "png"]:
        if not (root / f"results/figures/selected_study/{name}.{suffix}").is_file():
            raise ValueError("Missing scientific figure")
tables = load("results/tables/provenance.json")
for path, digest in tables["source_summaries"].items():
    if sha256(Path(path)) != digest:
        raise ValueError("Statistical table source changed")
for filename, record in tables["artifacts"].items():
    if sha256(root / "results/tables" / filename) != record["sha256"]:
        raise ValueError("Statistical CSV table changed")

lab = load("results/lab_plan/status.json")
for filename, digest in lab["panel_artifacts_sha256"].items():
    if sha256(root / "results/lab_plan" / filename) != digest:
        raise ValueError("Laboratory panel source changed")
scope_path = root / "results/current_research_scope.json"
scope = json.loads(scope_path.read_text()) if scope_path.exists() else {}
public_data_scope = (scope.get("study_scope") == "public_data_computational_research"
                     and scope.get("physical_qpcr_required") is False)
paused = (root / "results/paper_writing_paused.flag").exists()
if lab["experiments_conducted"] or (not public_data_scope and not paused):
    raise ValueError("Study status requires explicit update before manuscript work")
output = {"verified_at": utc_now(), "computational_checks_passed": True,
    "completed_sampling_runs": 246, "completed_full_reference_audits": len(audits),
    "recorded_proposals": sum(r["proposals"] for r in integrity["runs"]),
    "selected_condition_runs": 87, "held_out_genes": 75, "neural_model_seeds": 3,
    "long_training_trajectories": 6, "duration_comparison_groups": 8,
    "runtime_measurements": len(timings), "scientific_figure_sets": len(figure_names), "statistical_tables": len(tables["artifacts"]),
    "frozen_protocol_sha256": sha256(root / "results/frozen_protocol.json"),
    "reference_assembly": reference["reference_assembly"], "spliced_reference_exactly_verified_targets": 500, "audits": audits,
    "laboratory_assays_conducted": False, "laboratory_panel_designed_assays": lab["designed_assays"],
    "study_scope": scope.get("study_scope", "original_scope_with_independent_qpcr"),
    "physical_qpcr_required": not public_data_scope,
    "manuscript_work_paused": paused,
    "remaining_required_work": "Chinese review webpage and formal computational manuscript" if public_data_scope else "Independent qPCR measurements and their analysis, followed by the formal manuscript"}
write_json(root / "results/computational_completion_check.json", output)
print(json.dumps({k: v for k, v in output.items() if k != "audits"}, indent=2))
