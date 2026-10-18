#!/usr/bin/env python3
"""Build a portable Chinese research report from recorded computational results."""
import base64
import csv
import json
import shutil
from collections import Counter
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.metrics import paired_bootstrap_difference

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "reports/chinese_research"
OUTPUT.mkdir(parents=True, exist_ok=True)
EVIDENCE = OUTPUT / "evidence"
EVIDENCE.mkdir(exist_ok=True)
inputs = {}


def load(relative, lines=False):
    path = ROOT / relative
    inputs[relative] = sha256(path)
    return read_jsonl(path) if lines else json.loads(path.read_text())


completion = load("results/computational_completion_check.json")
scope = load("results/current_research_scope.json")
if not completion["computational_checks_passed"] or scope["physical_qpcr_required"]:
    raise ValueError("Report requires the verified computational study and current scope")
primary = load("results/temperature_tuning/selected_test_summary.json")
conditions = load("results/temperature_tuning/selected_conditions_summary.json")
training = load("results/figures/training/training_checks.json")
selection = load("results/temperature_tuning/selection.json")
duration = load("experiments/epoch100/results/duration_reference_comparison.json")
guidance = load("results/guide_ablation_summary.json")
decoder = load("results/complete_secondary_summary.json")
exon = load("results/exon_specificity_association.json")
runtime = load("results/runtime_benchmark/summary.json")
dataset = load("data/processed/dataset_partial_pilot.manifest.json")
reference = load("data/processed/reference_manifest.json")
partition = load("results/partition_audit/summary.json")
snapshot = load("results/experiment_snapshot.json")
integrity = load("results/checks/candidate_integrity/check.json")
acquisition = load("results/acquisition_status.json")
config = load("configs/pilot.json")
targets = {row["gene_id"]: row for row in load("data/processed/targets_pilot.jsonl", True)}
per_gene = load("results/temperature_tuning/selected_test_summary.per_gene.jsonl", True)
gene_ids = sorted({r["gene_id"] for r in per_gene}, key=lambda g: targets[g]["gene_symbol"].lower())
if len(gene_ids) != 75 or sum(s["sampling_complete_runs"] for s in snapshot["summaries"]) != 246:
    raise ValueError("Report cohort or experiment inventory differs")
for summary in [primary, conditions]:
    if summary["RNA_only"] or summary["validation_selection_sha256"] != sha256(ROOT / "results/temperature_tuning/selection.json"):
        raise ValueError("Report must use selected complete RNA/genome audits")
    for row in summary["inputs"]:
        directory = Path(row["directory"])
        if sha256(directory / "completion.json") != row["completion_sha256"] or sha256(directory / "specificity_bowtie.jsonl") != row["specificity_sha256"]:
            raise ValueError("Statistical inputs changed after analysis")

# All five methods at the same master seed are available in the candidate browser.
# Aggregate tables retain all three seeds, except deterministic Primer3.
candidates = []
for item in primary["inputs"]:
    directory = Path(item["directory"])
    if "seed20261013_" not in directory.name:
        continue
    rows = load(str((directory / "candidates.jsonl").relative_to(ROOT)), True)
    screens = {r["candidate_id"]: r for r in load(str((directory / "specificity_bowtie.jsonl").relative_to(ROOT)), True)}
    for row in rows:
        a = row["assessment"]
        s = screens.get(row["candidate_id"], {})
        candidates.append({
            "id": row["candidate_id"], "gene": row["gene_id"], "method": row["method"],
            "index": row["proposal_index"] + 1, "f": row.get("forward_primer"), "r": row.get("reverse_primer"),
            "fs": row.get("forward_start"), "fe": row.get("forward_end"), "rs": row.get("reverse_start"), "re": row.get("reverse_end"),
            "length": row.get("amplicon_length"), "feasible": a["computational_feasible"],
            "rna": bool(a["computational_feasible"] and s.get("gene_specificity_screen_pass")),
            "gflag": s.get("genomic_DNA_product_detected"), "failures": a.get("failures", []),
            "ftm": a.get("forward_properties", {}).get("tm"), "rtm": a.get("reverse_properties", {}).get("tm"),
            "penalty": a.get("shared_temperature_penalty"),
        })
if len(candidates) != 7500 or set(Counter(c["method"] for c in candidates).values()) != {1500}:
    raise ValueError("Candidate browser must contain 75 × 20 × 5 recorded proposals")

figure_names = ["primary_RNA_endpoints", "genomic_product_screen", "exon_geometry_partition", "proposal_budget",
                "condition_and_completion", "denoising_steps", "guidance_mechanism", "training_duration_endpoints", "generation_runtime"]
final_figures = ROOT / "results/final_fixed_seed/figures"
if (final_figures / "provenance.json").exists():
    figure_names = ["00_method_overview", "01_candidate_quality", "02_budget_and_diversity", "03_length_condition_control",
                    "04_completion_and_redesign", "05_components_and_denoising", "06_exon_structure", "training_duration_endpoints", "generation_runtime"]
figures = {}
for name in figure_names + ["training_curves"]:
    relative = (f"results/figures/training/{name}.png" if name == "training_curves" else
                f"results/final_fixed_seed/figures/{name}.png" if name[:2].isdigit() else f"results/figures/selected_study/{name}.png")
    path = ROOT / relative
    inputs[relative] = sha256(path)
    figures[name] = "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()
csvpath = ROOT / "results/figures/selected_study/figure_values.csv"
inputs[str(csvpath.relative_to(ROOT))] = sha256(csvpath)
with csvpath.open() as handle:
    figure_values = list(csv.DictReader(handle))
for row in figure_values:
    for key in ["budget", "mean", "ci95_lower", "ci95_upper", "genes"]:
        row[key] = float(row[key])

fixed_source = "results/final_fixed_seed/study_summary.json" if (ROOT / "results/final_fixed_seed/study_summary.json").exists() else "results/final_fixed_seed/existing_quality.json"
fixed = load(fixed_source) if (ROOT / fixed_source).exists() else None
followup_status = load("experiments/final_followup/results/status.json") if (ROOT / "experiments/final_followup/results/status.json").exists() else None
if fixed is not None:
    fixed_rows = load(str(Path(fixed_source).with_suffix(".per_gene.jsonl")), True)
    rename = {"one_no_genome_pair": "one_RNA_pass_pair_without_genomic_flag", "two_no_genome_pairs": "two_RNA_pass_pairs_without_genomic_flag"}
    def normalized(row):
        result = {rename.get(k, k): v for k, v in row.items()}
        result.update(study="temperature_selected", genes=75, seeds=[20261013])
        result["variant"] = "primary" if row["variant"] == "full" else row["variant"]
        for k, v in row.items():
            if k.endswith("_n"):
                result[rename.get(k[:-2], k[:-2]) + "_gene_denominator"] = v
        return result
    groups = [normalized(r) for r in fixed["summaries"]]
    comparisons = []
    for row in fixed["paired_comparisons"]:
        r = dict(row)
        r.update(endpoint=rename.get(r["endpoint"], r["endpoint"]), left=r["left"].replace("_full", "_primary"),
                 right=r["right"].replace("_full", "_primary"))
        comparisons.append(r)
    primary = {"summaries": [r for r in groups if r["context"] == "amp80-200" and r["variant"] == "primary"],
               "paired_comparisons": comparisons, "seed": 20261013}
    conditions = {"summaries": groups, "paired_comparisons": comparisons, "seed": 20261013}
    per_gene = [normalized(r) for r in fixed_rows if r["context"] == "amp80-200" and r["variant"] == "full"]
    figure_values = []
    metric_keys = ["unique_RNA_pass_pairs", "two_RNA_pass_pairs", "unique_RNA_pass_pairs_without_genomic_flag", "two_RNA_pass_pairs_without_genomic_flag",
                   "RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty", "RNA_no_genome_fraction"]
    normalized_records = [normalized(r) for r in fixed_rows]
    for group in groups:
        own = [r for r in normalized_records if all(r[k] == group[k] for k in ["method", "context", "variant", "budget"])]
        for key in metric_keys:
            values = [r[key] for r in own if r.get(key) is not None]
            if values:
                estimate = paired_bootstrap_difference(values, [0] * len(values), seed=20261013)
                figure_values.append({"study": "temperature_selected", "method": group["method"], "context": group["context"],
                    "variant": group["variant"], "budget": group["budget"], "endpoint": key, "mean": estimate["mean_difference"],
                    "ci95_lower": estimate["ci95"][0], "ci95_upper": estimate["ci95"][1], "genes": estimate["n"]})

evidence_paths = [
    "results/computational_completion_check.json", "results/current_research_scope.json",
    "data/processed/dataset_partial_pilot.manifest.json", "data/processed/reference_manifest.json",
    "results/acquisition_status.json", "results/partition_audit/summary.json",
    "results/temperature_tuning/selection.json", "results/temperature_tuning/selected_test_summary.json",
    "results/temperature_tuning/selected_test_summary.per_gene.jsonl", "results/temperature_tuning/selected_conditions_summary.json",
    "results/figures/training/training_checks.json", "results/figures/training/training_history.csv",
    "experiments/epoch100/results/duration_reference_comparison.json", "results/runtime_benchmark/summary.json",
    "results/checks/candidate_integrity/check.json", "results/guide_ablation_summary.json",
    "results/complete_secondary_summary.json", "results/exon_specificity_association.json",
    "results/spliced_reference_verification.json", "results/frozen_protocol.json",
    "results/figures/selected_study/provenance.json", "results/tables/provenance.json",
]
evidence_paths += [str(p.relative_to(ROOT)) for p in sorted((ROOT / "results/tables").glob("*.csv"))]
if fixed is not None:
    evidence_paths += [fixed_source, str(Path(fixed_source).with_suffix(".per_gene.jsonl")), "experiments/final_followup/protocol.json",
                       "experiments/final_followup/redesign_tasks.jsonl", "results/final_fixed_seed/paper_plan.json"]
    if (ROOT / "results/final_fixed_seed/study_summary.json").exists():
        evidence_paths += ["results/final_fixed_seed/completion_check.json", "results/final_fixed_seed/study_summary.groups.csv", "results/final_fixed_seed/study_summary.contrasts.csv",
                           "experiments/final_followup/anchor_prefilter_protocol.json"]
evidence_files = []
for number, relative in enumerate(evidence_paths, 1):
    path = ROOT / relative
    name = f"{number:02d}_{path.name}"
    shutil.copyfile(path, EVIDENCE / name)
    inputs[relative] = sha256(path)
    evidence_files.append({"path": relative, "href": "evidence/" + name, "sha256": inputs[relative], "bytes": path.stat().st_size})
shutil.copyfile(ROOT / "reports/assets/FONT_LICENSE.txt", OUTPUT / "FONT_LICENSE.txt")
data = {
    "built_at": utc_now(), "completion": {k: v for k, v in completion.items() if k != "audits"}, "scope": scope,
    "primary": primary, "conditions": conditions, "training": training, "selection": selection,
    "duration": duration, "guidance": guidance, "decoder": decoder, "exon": exon, "runtime": runtime,
    "dataset": {k: dataset[k] for k in ["rows", "by_split", "by_evidence", "exclusions"]},
    "reference": reference, "partition": partition, "inventory": snapshot["summaries"],
    "integrity": {k: v for k, v in integrity.items() if k != "runs"}, "acquisition": acquisition, "config": config,
    "genes": [{"id": g, "symbol": targets[g]["gene_symbol"], "accession": targets[g]["transcript_accession_version"], "length": len(targets[g]["template_sequence"])} for g in gene_ids],
    "per_gene": [r for r in per_gene if r["budget"] == 20], "candidates": candidates,
    "figures": figures, "figure_values": figure_values, "evidence": evidence_files,
    "fixed_seed": 20261013 if fixed is not None else None, "followup_status": followup_status,
}
if (ROOT / "paper/build_manifest.json").exists():
    build = json.loads((ROOT / "paper/build_manifest.json").read_text())
    data["manuscript_complete"] = build.get("status") == "formal_computational_manuscript_complete"
    if data["manuscript_complete"]:
        for filename in ["manuscript.pdf", "supplementary.pdf"]:
            shutil.copyfile(ROOT / "paper" / filename, OUTPUT / filename)
            inputs["paper/" + filename] = sha256(ROOT / "paper" / filename)
if followup_status and followup_status.get("anchor_prefilter_tasks_complete"):
    data["completion"]["completed_sampling_runs"] = 246 + followup_status["new_runs"]
    data["completion"]["completed_full_reference_audits"] = 246 + followup_status["new_runs"]
    data["completion"]["budgeted_proposals"] = 369000 + followup_status["new_runs"] * 1500
template = ROOT / "reports/templates/chinese_research.html"
inputs[str(template.relative_to(ROOT))] = sha256(template)
serialized = json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
html = template.read_text().replace("__STUDY_DATA__", serialized)
font_path = ROOT / "reports/assets/report_chinese_subset.otf"
inputs[str(font_path.relative_to(ROOT))] = sha256(font_path)
font_data = base64.b64encode(font_path.read_bytes()).decode()
font_css = ("<style>@font-face{font-family:DiffuPrimerCN;src:url(data:font/otf;base64," + font_data
            + ") format('opentype');font-display:swap;unicode-range:U+2000-206F,U+2190-21FF,U+3000-9FFF,U+FF00-FFEF;}"
            + "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI','PingFang SC','Microsoft YaHei',DiffuPrimerCN,sans-serif}"
            + ".hero h1{text-wrap:balance}</style>")
html = html.replace("</head>", font_css + "\n</head>")
(OUTPUT / "index.html").write_text(html)
write_json(OUTPUT / "build_manifest.json", {"built_at": data["built_at"], "scope": scope["study_scope"],
    "page_sha256": sha256(OUTPUT / "index.html"), "builder_sha256": sha256(Path(__file__)), "inputs": inputs,
    "recorded_master_seed_candidates": len(candidates), "figure_count": len(figures),
    "held_out_genes": len(gene_ids), "external_runtime_dependencies": [], "manuscript_complete": data.get("manuscript_complete", False)})
print(str(OUTPUT / "index.html"))
