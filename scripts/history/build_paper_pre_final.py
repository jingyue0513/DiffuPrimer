#!/usr/bin/env python3
import argparse
import json
import subprocess
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json


parser = argparse.ArgumentParser()
parser.add_argument("--aggregate", help="Completed non-calibration computational results")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
if (root / "results/paper_writing_paused.flag").exists():
    print("Manuscript work deferred until computational experiments are complete")
    raise SystemExit(0)
paper = root / "paper"
figures, tables = paper / "figures", paper / "tables"
figures.mkdir(exist_ok=True)
tables.mkdir(exist_ok=True)
source = root / "data/processed/primerbank_pilot.jsonl"
rows = read_jsonl(source) if source.exists() else []
status = json.loads((root / "results/acquisition_status.json").read_text())
split_manifest = json.loads((root / "data/splits/pilot.manifest.json").read_text())
mapping = Counter(row["mapping_status"] for row in rows)
evidence = Counter(row["evidence_type"] for row in rows)
dataset_manifest = json.loads((root / "data/processed/dataset_partial_pilot.manifest.json").read_text())
complete = status["completed_genes"] == status["planned_genes"]
data_text = (f"At the recorded acquisition snapshot, {status['completed_genes']} of {status['planned_genes']} planned genes had been queried. "
             f"The ledger contained {len(rows)} primer-pair records, including {mapping['unique_exact_pair_remap']} uniquely remapped pairs. "
             f"{sum(row.get('experimental_record_available', False) for row in rows)} records linked an experimental artifact. "
             f"The frozen reference-based partition contained {split_manifest['counts']['train']} training, "
             f"{split_manifest['counts']['validation']} validation, and {split_manifest['counts']['test']} test genes. "
             f"Combining deduplicated source positions with computational teacher pairs yielded {dataset_manifest['rows']:,} position-labelled records. "
             + r"\begin{center}\begin{tabular}{lr}\toprule Positional dataset & Records \\ \midrule" + "\n"
             + "\n".join(key.replace("_", "\\_") + f" & {value:,} \\\\" for key, value in sorted(dataset_manifest["by_evidence"].items()))
             + "\n\\midrule\n"
             + "\n".join(key.title() + f" & {value:,} \\\\" for key, value in sorted(dataset_manifest["by_split"].items()))
             + "\n\\bottomrule\\end{tabular}\\end{center}\n")
(tables / "data_status.tex").write_text(data_text)
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
fig, axes = plt.subplots(1, 3, figsize=(10, 3.1), constrained_layout=True)
axes[0].bar(["exact remap", "excluded"], [mapping["unique_exact_pair_remap"], len(rows) - mapping["unique_exact_pair_remap"]], color=["#167b9f", "#d08b54"])
axes[0].set_ylabel("Source pairs")
axes[0].set_title("Reference mapping")
lengths = [row["amplicon_length"] for row in rows]
if lengths:
    axes[1].hist(lengths, bins=30, color="#167b9f", edgecolor="white")
axes[1].axvspan(80, 200, alpha=0.13, color="#ca7000", label="primary condition")
axes[1].set_xlabel("Source product length (bp)")
axes[1].set_title("Observed source distribution")
axes[1].legend(frameon=False, fontsize=7)
axes[2].bar(["train", "validation", "test"], [split_manifest["counts"][key] for key in ["train", "validation", "test"]], color=["#167b9f", "#819c59", "#d08b54"])
axes[2].set_ylabel("Frozen genes")
axes[2].set_title("Gene/component partition")
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
fig.savefig(figures / "data_audit.pdf")
fig.savefig(figures / "data_audit.png", dpi=220)
plt.close(fig)
if args.aggregate:
    aggregate_path = root / args.aggregate
    aggregate = json.loads(aggregate_path.read_text())
    if aggregate["stage"] == "calibration" or any(run["stage"] == "calibration" for run in aggregate["runs"]):
        raise ValueError("Calibration runs cannot populate the manuscript's comparative results")
    runs = aggregate["runs"]
    text = "All 75 frozen test genes were evaluated with a budget of 20 proposals per target.\\par\n"
    text += r"\begin{center}\small\begin{tabular}{lrrrr}\toprule Method/seed & Genes & Unique feasible & Feasibility & RNA-screen pairs \\ \midrule" + "\n"
    for run in runs:
        label = run["method"].replace("_", "\\_") + "/" + str(run["seed"])[-2:]
        specificity = "pending" if run["mean_unique_gene_screen_pass_pairs"] is None else f"{run['mean_unique_gene_screen_pass_pairs']:.2f}"
        text += f"{label} & {run['genes']} & {run['mean_unique_feasible_pairs']:.2f} & {run['mean_budget_feasibility_rate']:.3f} & {specificity} \\\\\n"
    text += "\\bottomrule\\end{tabular}\\end{center}\nUnique-pair counts are means per target. Feasibility uses the full proposal budget as denominator. RNA-screen pairs satisfy the common thermodynamic and gene-specificity screen.\n"
    (tables / "computational_results.tex").write_text(text)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
    labels = [run["method"] + "\n" + str(run["seed"])[-2:] for run in runs]
    axes[0].bar(labels, [run["mean_unique_feasible_pairs"] for run in runs], color="#167b9f")
    axes[0].set_ylabel("Mean unique feasible pairs / target")
    axes[1].bar(labels, [run["mean_budget_feasibility_rate"] for run in runs], color="#819c59")
    axes[1].set_ylabel("Feasibility / fixed proposal budget")
    axes[1].set_ylim(0, 1)
    for ax in axes:
        ax.tick_params(axis="x", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(figures / "computational_comparison.pdf")
    fig.savefig(figures / "computational_comparison.png", dpi=220)
    plt.close(fig)
for _ in range(2):
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "manuscript.tex"], cwd=paper,
                   check=True, stdout=(paper / "build.stdout").open("w"), stderr=subprocess.STDOUT)
write_json(paper / "build_manifest.json", {"built_at": utc_now(), "manuscript_sha256": sha256(paper / "manuscript.tex"),
    "pdf_sha256": sha256(paper / "manuscript.pdf"), "source_data_sha256": sha256(source),
    "aggregate_sha256": sha256(root / args.aggregate) if args.aggregate else None,
    "status": "working draft; independent assays and final submission requirements pending"})
print(str(paper / "manuscript.pdf"))
