#!/usr/bin/env python3
"""Compare complete validation experiments and isolate each implementation change."""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference


root = Path(__file__).resolve().parents[1]
expected = {r["gene_id"] for r in read_jsonl(root / "data/processed/targets_pilot.jsonl") if r["split"] == "validation"}
metrics = ["unique_feasible_pairs", "budget_feasibility_rate", "distinct_feasible_amplicons",
    "mean_feasible_amplicon_jaccard_distance", "best_shared_temperature_penalty", "elapsed_seconds"]
runs, inputs = [], []


def collect(study, paths):
    for path in paths:
        manifest = json.loads((path / "manifest.json").read_text())
        completion = json.loads((path / "completion.json").read_text())
        if not completion["all_expected_targets_evaluated"]:
            raise ValueError("Incomplete validation experiment")
        if completion["per_gene_sha256"] != sha256(path / "per_gene.jsonl"):
            raise ValueError("Validation output changed")
        rows = read_jsonl(path / "per_gene.jsonl")
        if {r["gene_id"] for r in rows} != expected:
            raise ValueError("Validation comparison must include all 75 frozen genes")
        for row in rows:
            runs.append({**row, "study": study, "seed": manifest["arguments"]["seed"], "source": str(path)})
        inputs.append({"study": study, "directory": str(path), "manifest_sha256": sha256(path / "manifest.json"),
            "completion_sha256": sha256(path / "completion.json")})


primary = [root / f"results/evaluations/{mode}_{'pilot' if mode == 'rule_sampler' else 'weak_partial_pilot'}_validation_seed{seed}_amp80-200_oligo_guided"
    for seed in [20261013, 20261014, 20261015] for mode in ["direct", "diffusion", "rule_sampler"]]
baseline = root / "results/evaluations/primer3_pilot_validation_seed20261013_amp80-200"
if not all((p / "completion.json").exists() for p in primary + [baseline]):
    raise ValueError("Primary full validation incomplete")
collect("primary", primary + [baseline])
duration_seeds = [20261013, 20261014, 20261015] if (root / "experiments/epoch100/three_seed_protocol.json").exists() else [20261013]
long = [root / f"experiments/epoch100/results/evaluations/{mode}_weak_partial_pilot_validation_seed{seed}_amp80-200_oligo_guided" for seed in duration_seeds for mode in ["direct", "diffusion"]]
if all((p / "completion.json").exists() for p in long):
    collect("epoch100", long)
joint = [root / f"experiments/joint_tm/results/evaluations/{mode}_{'pilot' if mode == 'rule_sampler' else 'weak_partial_pilot'}_validation_seed{seed}_amp80-200_oligo_guided_joint_tm"
    for seed in [20261013, 20261014, 20261015] for mode in ["direct", "diffusion", "rule_sampler"]]
if all((p / "completion.json").exists() for p in joint):
    collect("joint_tm", joint)
grouped = defaultdict(list)
for row in runs:
    grouped[(row["study"], row["method"], row["gene_id"])].append(row)
averaged = []
for (study, method, gene), rows in sorted(grouped.items()):
    record = {"study": study, "method": method, "gene_id": gene, "seeds": [r["seed"] for r in rows]}
    for metric in metrics:
        own = [r[metric] for r in rows if r[metric] is not None]
        record[metric] = float(np.mean(own)) if own else None
        record[metric + "_seed_denominator"] = len(own)
    averaged.append(record)
summaries = []
for study, method in sorted({(r["study"], r["method"]) for r in averaged}):
    own = [r for r in averaged if r["study"] == study and r["method"] == method]
    record = {"study": study, "method": method, "genes": len(own), "seeds": sorted({s for r in own for s in r["seeds"]})}
    for metric in metrics:
        values = [r[metric] for r in own if r[metric] is not None]
        record[metric] = float(np.mean(values)) if values else None
        record[metric + "_gene_denominator"] = len(values)
    summaries.append(record)
comparisons = []


def compare(left_study, left_method, right_study, right_method, purpose, same_master_seed=False):
    if same_master_seed:
        left = {r["gene_id"]: r for r in runs if r["study"] == left_study and r["method"] == left_method and r["seed"] == 20261013}
        right = {r["gene_id"]: r for r in runs if r["study"] == right_study and r["method"] == right_method and r["seed"] == 20261013}
    else:
        left = {r["gene_id"]: r for r in averaged if r["study"] == left_study and r["method"] == left_method}
        right = {r["gene_id"]: r for r in averaged if r["study"] == right_study and r["method"] == right_method}
    if not left or not right:
        return
    for metric in metrics[:-1]:
        genes = sorted(g for g in set(left) & set(right) if left[g][metric] is not None and right[g][metric] is not None)
        if not genes:
            continue
        estimate = paired_bootstrap_difference([left[g][metric] for g in genes], [right[g][metric] for g in genes], seed=20261013)
        comparisons.append({"left_study": left_study, "left_method": left_method, "right_study": right_study,
            "right_method": right_method, "purpose": purpose, "endpoint": metric, "master_seed_only": same_master_seed, **estimate})


for method in ["direct", "rule_sampler", "primer3"]:
    compare("primary", "diffusion", "primary", method, "Primary method comparison")
for method in ["direct", "diffusion"]:
    compare("epoch100", method, "primary", method, "Training-duration effect under matched seeds", same_master_seed=len(duration_seeds) == 1)
for method in ["direct", "diffusion", "rule_sampler"]:
    compare("joint_tm", method, "primary", method, "Joint-pair-guidance effect")
for method in ["direct", "rule_sampler"]:
    compare("joint_tm", "diffusion", "joint_tm", method, "Method difference under matched joint-pair guidance")
compare("joint_tm", "diffusion", "primary", "primer3", "Joint-guided diffusion versus Primer3")
write_jsonl(root / "results/validation_experiment_comparison.per_gene.jsonl", averaged)
write_json(root / "results/validation_experiment_comparison.json", {"generated_at": utc_now(), "inputs": inputs,
    "summaries": summaries, "paired_comparisons": comparisons,
    "comparison_policy": "Training-duration comparisons use the same declared seeds in both durations. Three-seed comparisons average within gene before paired gene bootstrap; single-seed duration checks are identified explicitly. Coverage-restricted score/diversity comparisons state their denominators."})
print(json.dumps(summaries, indent=2))
