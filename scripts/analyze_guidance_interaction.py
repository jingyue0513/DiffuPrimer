#!/usr/bin/env python3
"""Estimate paired differences in response to the joint-Tm intervention."""
import json
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference


root = Path(__file__).resolve().parents[1]
source = root / "results/guide_ablation_summary.json"
summary = json.loads(source.read_text())
records = read_jsonl(source.with_suffix(".per_gene.jsonl"))
expected = {r["gene_id"] for r in read_jsonl(root / "data/processed/targets_pilot.jsonl") if r["split"] == "test"}
groups = {}
for study in ["original", "joint_tm"]:
    for method in ["direct", "diffusion", "rule_sampler"]:
        own = [r for r in records if r["study"] == study and r["method"] == method and
               r["context"] == "amp80-200" and r["variant"] == "primary" and r["budget"] == 20]
        if len(own) != 75 or {r["gene_id"] for r in own} != expected or any(r["seeds"] != [20261013, 20261014, 20261015] for r in own):
            raise ValueError("Incomplete paired guidance cohort")
        groups[(study, method)] = {r["gene_id"]: r for r in own}
for item in summary["inputs"]:
    directory = Path(item["directory"])
    if sha256(directory / "completion.json") != item["completion_sha256"] or sha256(directory / "specificity_bowtie.jsonl") != item["specificity_sha256"]:
        raise ValueError("Guidance analysis inputs changed")
fields = ["unique_feasible_pairs", "budget_feasibility_rate", "unique_RNA_pass_pairs", "one_RNA_pass_pair",
          "two_RNA_pass_pairs", "unique_RNA_pass_pairs_without_genomic_flag", "two_RNA_pass_pairs_without_genomic_flag",
          "RNA_pass_amplicon_distance", "distinct_RNA_pass_amplicons", "best_RNA_pass_temperature_penalty"]
rows, contrasts = [], []
for comparator in ["direct", "rule_sampler"]:
    for field in fields:
        own = []
        for gene in sorted(expected):
            values = {(study, method): groups[(study, method)][gene][field]
                      for study in ["original", "joint_tm"] for method in ["diffusion", comparator]}
            if any(v is None for v in values.values()):
                continue
            diffusion_gain = values[("joint_tm", "diffusion")] - values[("original", "diffusion")]
            comparator_gain = values[("joint_tm", comparator)] - values[("original", comparator)]
            row = {"gene_id": gene, "comparator": comparator, "endpoint": field,
                   "diffusion_guidance_response": diffusion_gain, "comparator_guidance_response": comparator_gain,
                   "response_difference": diffusion_gain - comparator_gain}
            own.append(row)
        if own:
            estimate = paired_bootstrap_difference([r["diffusion_guidance_response"] for r in own],
                                                  [r["comparator_guidance_response"] for r in own], seed=20261013)
            contrasts.append({"comparator": comparator, "endpoint": field, "total_target_genes": 75,
                              "contrast": "(joint minus single guidance in diffusion) minus (joint minus single guidance in comparator)", **estimate})
            rows.extend(own)
output = root / "results/guidance_interaction.json"
write_jsonl(output.with_suffix(".per_gene.jsonl"), rows)
write_json(output, {"at": utc_now(), "guide_summary_sha256": sha256(source), "comparisons": contrasts,
                   "unit": "Target gene after averaging three model seeds", "budget": 20,
                   "scope": "Mechanism analysis at fixed temperature 1; complete RNA/genome audit; available-gene denominators recorded for quality and diversity"})
print(json.dumps([r for r in contrasts if r["endpoint"] == "unique_RNA_pass_pairs"], indent=2))
