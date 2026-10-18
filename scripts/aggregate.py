#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference


parser = argparse.ArgumentParser()
parser.add_argument("directories", nargs="+")
parser.add_argument("--output", default="results/aggregate.json")
parser.add_argument("--calibration", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
runs, per_gene_rows = [], []
for raw in args.directories:
    directory = root / raw
    manifest = json.loads((directory / "manifest.json").read_text())
    completion = json.loads((directory / "completion.json").read_text())
    if not completion["all_expected_targets_evaluated"]:
        raise ValueError("Incomplete evaluation cannot be summarized as complete")
    if completion["per_gene_sha256"] != sha256(directory / "per_gene.jsonl") or completion["candidates_sha256"] != sha256(directory / "candidates.jsonl"):
        raise ValueError("Evaluation files changed after completion")
    if manifest["stage"] == "calibration" and not args.calibration:
        raise ValueError("Calibration results must remain explicitly labeled")
    rows = read_jsonl(directory / "per_gene.jsonl")
    candidates = read_jsonl(directory / "candidates.jsonl")
    specificity_path = directory / "specificity.jsonl"
    specificity = {row["candidate_id"]: row for row in read_jsonl(specificity_path)} if specificity_path.exists() else {}
    if specificity:
        status = json.loads((directory / "specificity_status.json").read_text())
        if status["candidate_sha256"] != sha256(directory / "candidates.jsonl"):
            raise ValueError("Specificity screen uses different candidates")
    label = directory.name
    per_gene = []
    for row in rows:
        own = [c for c in candidates if c["gene_id"] == row["gene_id"]]
        eligible = [c for c in own if c["assessment"]["computational_feasible"] and specificity.get(c["candidate_id"], {}).get("gene_specificity_screen_pass")]
        unique_specific = len({(c["forward_primer"], c["reverse_primer"]) for c in eligible}) if specificity else None
        record = {**row, "run": label, "unique_gene_screen_pass_pairs": unique_specific,
                  "two_gene_screen_pass_pairs": unique_specific >= 2 if unique_specific is not None else None}
        per_gene.append(record)
        per_gene_rows.append(record)
    mean = lambda key: float(np.mean([row[key] for row in per_gene]))
    summary = {"run": label, "method": manifest["method"], "stage": manifest["stage"], "genes": len(per_gene),
        "candidate_budget": manifest["arguments"]["budget"], "seed": manifest["arguments"]["seed"],
        "mean_unique_feasible_pairs": mean("unique_feasible_pairs"),
        "two_feasible_pairs_coverage": float(np.mean([x["unique_feasible_pairs"] >= 2 for x in per_gene])),
        "target_coverage": mean("target_covered"), "mean_budget_feasibility_rate": mean("budget_feasibility_rate"),
        "median_runtime_seconds": float(np.median([x["elapsed_seconds"] for x in per_gene])),
        "mean_unique_gene_screen_pass_pairs": mean("unique_gene_screen_pass_pairs") if specificity else None,
        "two_gene_screen_pass_pairs_coverage": mean("two_gene_screen_pass_pairs") if specificity else None,
        "source": str(directory), "rows": per_gene}
    runs.append(summary)
comparisons = []
for index, run in enumerate(runs):
    for baseline in runs[:index]:
        left = {row["gene_id"]: row for row in run["rows"]}
        right = {row["gene_id"]: row for row in baseline["rows"]}
        if set(left) != set(right):
            continue
        keys = ["unique_feasible_pairs", "budget_feasibility_rate"]
        if all(row["unique_gene_screen_pass_pairs"] is not None for row in run["rows"] + baseline["rows"]):
            keys += ["unique_gene_screen_pass_pairs", "two_gene_screen_pass_pairs"]
        for key in keys:
            ordered = sorted(left)
            comparisons.append({"left": run["run"], "right": baseline["run"], "endpoint": key,
                                **paired_bootstrap_difference([left[g][key] for g in ordered], [right[g][key] for g in ordered], seed=20261013)})
output = {"generated_at": utc_now(), "stage": "calibration" if args.calibration else "pilot",
          "runs": [{k: v for k, v in run.items() if k != "rows"} for run in runs], "paired_comparisons": comparisons,
          "interpretation": "No assay success rates are present; specificity endpoints refer only to the declared in-silico rule."}
write_json(root / args.output, output)
write_jsonl((root / args.output).with_suffix(".per_gene.jsonl"), per_gene_rows)
print(json.dumps(output["runs"], indent=2))
