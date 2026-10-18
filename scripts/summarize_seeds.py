#!/usr/bin/env python3
"""Average seed-level target outcomes before paired gene bootstrap comparisons."""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference


root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--analysis-only", action="store_true")
args = parser.parse_args()
aggregate_path = root / "results/primary_aggregate.json"
aggregate = json.loads(aggregate_path.read_text())
rows = read_jsonl(aggregate_path.with_suffix(".per_gene.jsonl"))
run_meta = {row["run"]: row for row in aggregate["runs"]}
expected = {row["gene_id"] for row in read_jsonl(root / "data/processed/targets_pilot.jsonl") if row["split"] == "test"}
if len(expected) != 75:
    raise ValueError("Unexpected frozen test gene count")
groups = defaultdict(list)
for row in rows:
    if row["unique_gene_screen_pass_pairs"] is None:
        raise ValueError("Primary specificity screen incomplete")
    groups[(row["method"], row["gene_id"])].append(row)
methods = sorted({key[0] for key in groups})
averaged = []
metrics = ["unique_feasible_pairs", "budget_feasibility_rate", "unique_gene_screen_pass_pairs", "two_gene_screen_pass_pairs"]
for method in methods:
    own_genes = {gene for model, gene in groups if model == method}
    if own_genes != expected:
        raise ValueError("Method does not retain all frozen test genes")
    for gene in sorted(expected, key=int):
        values = groups[(method, gene)]
        count = 1 if method == "primer3" else 3
        if len(values) != count:
            raise ValueError("Missing or duplicate seed/target rows")
        record = {"method": method, "gene_id": gene, "seeds": len(values)}
        for key in metrics:
            record[key] = float(np.mean([row[key] for row in values]))
        distances = [row.get("mean_feasible_amplicon_jaccard_distance") for row in values]
        distances = [value for value in distances if value is not None]
        record["mean_feasible_amplicon_jaccard_distance"] = float(np.mean(distances)) if distances else None
        averaged.append(record)
summaries = []
for method in methods:
    own = [row for row in averaged if row["method"] == method]
    summary = {"method": method, "genes": len(own)}
    for key in metrics:
        summary[key] = float(np.mean([row[key] for row in own]))
    distances = [row["mean_feasible_amplicon_jaccard_distance"] for row in own if row["mean_feasible_amplicon_jaccard_distance"] is not None]
    summary["mean_feasible_amplicon_jaccard_distance"] = float(np.mean(distances)) if distances else None
    summary["jaccard_distance_gene_denominator"] = len(distances)
    seed_runs = [run for run in aggregate["runs"] if run["method"] == method]
    summary["seed_range_mean_unique_gene_screen_pass_pairs"] = [min(run["mean_unique_gene_screen_pass_pairs"] for run in seed_runs), max(run["mean_unique_gene_screen_pass_pairs"] for run in seed_runs)]
    summaries.append(summary)
comparisons = []
lookup = {(row["method"], row["gene_id"]): row for row in averaged}
for method in methods:
    if method == "diffusion":
        continue
    for key in ["unique_gene_screen_pass_pairs", "two_gene_screen_pass_pairs"]:
        genes = sorted(expected, key=int)
        estimate = paired_bootstrap_difference([lookup[("diffusion", gene)][key] for gene in genes],
            [lookup[(method, gene)][key] for gene in genes], seed=20261013)
        comparisons.append({"left": "diffusion", "right": method, "endpoint": key, **estimate})
write_jsonl(root / "results/primary_seed_averaged.per_gene.jsonl", averaged)
write_json(root / "results/primary_seed_averaged.json", {"generated_at": utc_now(),
    "input_sha256": sha256(aggregate_path), "per_gene_input_sha256": sha256(aggregate_path.with_suffix(".per_gene.jsonl")),
    "summaries": summaries, "paired_comparisons": comparisons,
    "interpretation": "Seeds averaged within gene; paired genes bootstrapped. Confidence intervals are conditional on the observed trained seeds. No laboratory outcomes exist."})
text = "Seed outcomes were averaged within each target before gene-level paired comparisons. The analysis retained all 75 test genes.\\par\n"
for comparison in comparisons:
    if comparison["endpoint"] != "unique_gene_screen_pass_pairs":
        continue
    low, high = comparison["ci95"]
    right = comparison["right"].replace("_", "\\_")
    text += (f"The difference in mean unique RNA-screen-pass pairs per target, diffusion minus {right}, was "
             f"{comparison['mean_difference']:.2f} (paired gene-bootstrap 95\\% interval {low:.2f} to {high:.2f}).\\par\n")
text += "These intervals condition on the observed three learned/sampling seeds and measure computational endpoints, not experimental success.\n"
if not args.analysis_only:
    (root / "paper/tables/paired_estimates.tex").write_text(text)
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.5), constrained_layout=True)
labels = [row["method"] for row in summaries]
axes[0].bar(labels, [row["unique_gene_screen_pass_pairs"] for row in summaries], color="#167f9d")
axes[0].set_ylabel("Mean unique RNA-screen-pass pairs / gene")
axes[1].bar(labels, [row["two_gene_screen_pass_pairs"] for row in summaries], color="#819c59")
axes[1].set_ylabel("Fraction of targets with at least two pairs")
axes[1].set_ylim(0, 1)
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
figure_dir = root / ("results/figures" if args.analysis_only else "paper/figures")
figure_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(figure_dir / "primary_seed_averaged.pdf")
fig.savefig(figure_dir / "primary_seed_averaged.png", dpi=220)
plt.close(fig)
print(json.dumps({"summaries": summaries, "comparisons": comparisons}, indent=2))
