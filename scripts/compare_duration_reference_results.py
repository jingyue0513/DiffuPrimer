#!/usr/bin/env python3
"""Compare actual 30/100-epoch endpoints under identical generation settings."""
import json
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.metrics import paired_bootstrap_difference


root = Path(__file__).resolve().parents[1]
experiment = root / "experiments/epoch100"
fields = ["unique_feasible_pairs", "budget_feasibility_rate", "unique_RNA_pass_pairs", "one_RNA_pass_pair",
          "two_RNA_pass_pairs", "unique_RNA_pass_pairs_without_genomic_flag", "two_RNA_pass_pairs_without_genomic_flag",
          "RNA_pass_amplicon_distance", "distinct_RNA_pass_amplicons", "best_RNA_pass_temperature_penalty"]
summaries, comparisons, inputs = [], [], []
for split, short in [("validation", root / "results/bowtie_audit/validation_summary.json"),
                     ("test", root / "results/primary_bowtie_summary.json")]:
    paths = {30: short, 100: experiment / f"results/{split}_bowtie_summary.json"}
    cohorts = {}
    for epochs, path in paths.items():
        summary = json.loads(path.read_text())
        if summary["split"] != split or summary.get("RNA_only") or summary["screen"] != "bowtie":
            raise ValueError("Duration comparison needs the complete, matching reference audit")
        records = read_jsonl(path.with_suffix(".per_gene.jsonl"))
        own = [r for r in records if r["study"] == "original" and r["context"] == "amp80-200" and
               r["variant"] == "primary" and r["budget"] == 20 and r["method"] in ["direct", "diffusion"]]
        if len(own) != 150 or any(r["seeds"] != [20261013, 20261014, 20261015] for r in own):
            raise ValueError("Duration cohort or seed denominator differs")
        cohorts[epochs] = {(r["method"], r["gene_id"]): r for r in own}
        if len(cohorts[epochs]) != 150:
            raise ValueError("Repeated duration cohort gene")
        for group in summary["summaries"]:
            if group["study"] == "original" and group["context"] == "amp80-200" and group["variant"] == "primary" and group["budget"] == 20 and group["method"] in ["direct", "diffusion"]:
                summaries.append({"split": split, "epochs": epochs, **group})
        for item in summary["inputs"]:
            directory = Path(item["directory"])
            if sha256(directory / "completion.json") != item["completion_sha256"] or sha256(directory / "specificity_bowtie.jsonl") != item["specificity_sha256"]:
                raise ValueError("Duration analysis source changed")
        inputs.append({"split": split, "epochs": epochs, "summary_sha256": sha256(path), "summary": str(path),
                       "per_gene_sha256": sha256(path.with_suffix(".per_gene.jsonl"))})
    if cohorts[30].keys() != cohorts[100].keys():
        raise ValueError("Different target genes in duration groups")
    for method in ["direct", "diffusion"]:
        keys = sorted(k for k in cohorts[30] if k[0] == method)
        for field in fields:
            paired = [k for k in keys if cohorts[30][k][field] is not None and cohorts[100][k][field] is not None]
            if not paired:
                continue
            estimate = paired_bootstrap_difference([cohorts[100][k][field] for k in paired],
                                                  [cohorts[30][k][field] for k in paired], seed=20261013)
            comparisons.append({"split": split, "method": method, "endpoint": field,
                                "comparison": "100 minus 30 epochs", "total_target_genes": 75, **estimate})
write_json(experiment / "results/duration_reference_comparison.json", {"at": utc_now(), "inputs": inputs,
           "summaries": summaries, "paired_comparisons": comparisons,
           "scope": "Same seeds, whole-transcript inference, single-oligo guidance, temperature 1, 20 proposals, full RNA/genome audit; original 30-epoch and selected-temperature studies retained"})
print(json.dumps({"groups": len(summaries), "paired_endpoints": len(comparisons)}), flush=True)
