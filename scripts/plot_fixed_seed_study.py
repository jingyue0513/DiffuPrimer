#!/usr/bin/env python3
"""Six scientific result figures and an overview for the completed fixed-seed study."""
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.metrics import paired_bootstrap_difference

ROOT = Path(__file__).resolve().parents[1]
base = ROOT / "results/final_fixed_seed"
summary = json.loads((base / "study_summary.json").read_text())
rows = read_jsonl(base / "study_summary.per_gene.jsonl")
out = base / "figures"
out.mkdir(exist_ok=True)
methods = ["primer3", "uniform_pair_sampler", "rule_sampler", "direct", "diffusion"]
labels = {"primer3": "Primer3", "uniform_pair_sampler": "Uniform pairs", "rule_sampler": "Forward-first", "direct": "Direct", "diffusion": "DiffuPrimer"}
colors = {"primer3": "#db9c27", "uniform_pair_sampler": "#37926c", "rule_sampler": "#a4a7a9", "direct": "#3d7aa8", "diffusion": "#855ab0"}
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
values = []


def estimate(method, endpoint, context="amp80-200", variant="full", budget=20):
    own = [r for r in rows if r["method"] == method and r["context"] == context and r["variant"] == variant and r["budget"] == budget]
    data = [r[endpoint] for r in own if r[endpoint] is not None]
    if not data:
        return None
    ci = paired_bootstrap_difference(data, np.zeros(len(data)), seed=20261013)
    values.append({"method": method, "endpoint": endpoint, "context": context, "variant": variant,
                   "budget": budget, "mean": ci["mean_difference"], "ci95": ci["ci95"], "n": ci["n"]})
    return ci["mean_difference"], *ci["ci95"], ci["n"]


def bars(ax, endpoint, context="amp80-200", subset=methods, title=None):
    for i, method in enumerate(subset):
        e = estimate(method, endpoint, context)
        if e:
            mean, lo, hi, n = e
            ax.bar(i, mean, color=colors[method], width=.65)
            ax.vlines(i, lo, hi, color="black", lw=.8)
            ax.hlines([lo, hi], i-.06, i+.06, color="black", lw=.8)
            if n != 75:
                ax.annotate(f"n={n}", (i, hi), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=7)
    ax.set_xticks(range(len(subset)), [labels[m] for m in subset], rotation=20, ha="right")
    ax.set_ylim(bottom=0)
    if title:
        ax.set_title(title)


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(out / (name + "." + ext), bbox_inches="tight", dpi=220)
    plt.close(fig)


fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
for ax, key, title in zip(axes.flat,
    ["unique_RNA_pass_pairs_without_genomic_flag", "RNA_no_genome_fraction", "unique_RNA_pass_pairs", "two_RNA_pass_pairs"],
    ["RNA-specific pairs without genomic product flag", "Fraction retained after genomic-product screening", "Unique RNA-specific pairs", "Targets with at least two RNA-specific pairs"]):
    bars(ax, key, title=title)
save(fig, "01_candidate_quality")

fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
for ax, key, title in zip(axes.flat,
    ["unique_RNA_pass_pairs_without_genomic_flag", "unique_RNA_pass_pairs", "RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty"],
    ["RNA-specific pairs without genomic flag", "Unique RNA-specific pairs", "Matched neural models: amplicon diversity", "Matched neural models: best Tm penalty"]):
    subset = ["direct", "diffusion"] if key in ["RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty"] else methods
    for method in subset:
        es = [estimate(method, key, budget=k) for k in [1, 2, 5, 10, 20]]
        ks = [k for k, e in zip([1, 2, 5, 10, 20], es) if e]
        es = [e for e in es if e]
        ax.plot(ks, [e[0] for e in es], "o-", label=labels[method], color=colors[method], lw=1.3, ms=3)
        ax.fill_between(ks, [e[1] for e in es], [e[2] for e in es], color=colors[method], alpha=.12)
    ax.set(title=title, xlabel="Proposal budget", xticks=[1, 2, 5, 10, 20], ylim=(0, None))
axes[0, 0].legend(frameon=False, fontsize=7)
axes[1, 0].legend(frameon=False, fontsize=7)
save(fig, "02_budget_and_diversity")

fig, axes = plt.subplots(2, 3, figsize=(12, 6.6), layout="constrained")
for column, context in enumerate(["amp80-120", "amp121-160", "amp161-200"]):
    bars(axes[0, column], "unique_RNA_pass_pairs_without_genomic_flag", context, title=context.replace("amp", "") + " bp")
    bars(axes[1, column], "RNA_no_genome_fraction", context)
    axes[0, column].set_ylabel("Pairs without genomic product flag")
    axes[1, column].set_ylabel("Fraction among RNA-specific pairs")
save(fig, "03_length_condition_control")

fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
for ax, context, title in zip(axes.flat,
    ["amp80-200_fix_forward", "amp80-200_fix_reverse", "redesign_forward", "redesign_reverse"],
    ["Forward-primer completion", "Reverse-primer completion", "Changed-length redesign: retain forward", "Changed-length redesign: retain reverse"]):
    bars(ax, "unique_RNA_pass_pairs", context, subset=["direct", "diffusion"], title=title)
    ax.set_ylabel("Unique RNA-specific pairs")
save(fig, "04_completion_and_redesign")

fig, axes = plt.subplots(2, 3, figsize=(12, 6.8), layout="constrained")
variants = ["raw", "geometry_only", "single_oligo", "full"]
variant_labels = ["Raw labels", "Geometry", "Single-oligo", "Joint Tm"]
for ax, key, title in zip(axes[0], ["geometry_valid_rate", "unique_RNA_pass_pairs", "best_RNA_pass_temperature_penalty"],
                          ["Structural validity", "RNA-specific pairs", "Best paired Tm penalty"]):
    for i, variant in enumerate(variants):
        e = estimate("diffusion", key, variant=variant)
        if e:
            ax.bar(i, e[0], color=colors["diffusion"], alpha=.4+i*.2, width=.65)
            ax.vlines(i, e[1], e[2], color="black", lw=.8)
        else:
            ax.text(i, 0, "No RNA\npairs", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(4), variant_labels, rotation=20, ha="right")
    ax.set(title=title, ylim=(0, None))
for ax, key, title in zip(axes[1], ["unique_RNA_pass_pairs", "RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty"],
                          ["Denoising: RNA-specific pairs", "Denoising: amplicon diversity", "Denoising: best Tm penalty"]):
    es = [estimate("diffusion", key, variant=v) for v in ["steps1", "steps4", "steps8", "full"]]
    ax.plot([1, 4, 8, 16], [e[0] for e in es], "o-", color=colors["diffusion"])
    ax.fill_between([1, 4, 8, 16], [e[1] for e in es], [e[2] for e in es], color=colors["diffusion"], alpha=.13)
    ax.set(title=title, xlabel="Denoising steps", xticks=[1, 4, 8, 16], ylim=(0, None))
save(fig, "05_components_and_denoising")

# Structure attribution uses exactly the master-seed candidates in the main comparison.
main_inputs = [i for i in summary["inputs"] if any(r["run"] == i["directory"] and r["context"] == "amp80-200" and r["variant"] == "full" for r in rows)]
annotation_ids = {}
for line in (ROOT / "results/selected_exon_context_audit.per_candidate.jsonl").open():
    a = json.loads(line)
    if any(a["run_directory"] == i["directory"] for i in main_inputs):
        annotation_ids[a["candidate_id"]] = a
categories = ["primer_junction", "amplicon_junction", "same_exon", "polyA", "unmapped"]
category_rows = []
for item in main_inputs:
    directory = Path(item["directory"])
    candidates = read_jsonl(directory / "candidates.jsonl")
    screens = {r["candidate_id"]: r for r in read_jsonl(directory / "specificity_bowtie.jsonl")}
    seen = set()
    counts = defaultdict(lambda: defaultdict(int))
    for c in candidates:
        if not c["assessment"]["computational_feasible"]:
            continue
        s = screens[c["candidate_id"]]
        if not s.get("gene_specificity_screen_pass") or s["genomic_DNA_product_detected"]:
            continue
        key = (c["gene_id"], c["forward_primer"], c["reverse_primer"])
        if key in seen:
            continue
        seen.add(key)
        a = annotation_ids[c["candidate_id"]]
        category = "polyA" if a["forward_overlaps_terminal_polyA"] or a["reverse_overlaps_terminal_polyA"] else "primer_junction" if a["forward_spans_junction"] or a["reverse_spans_junction"] else "amplicon_junction" if a["amplicon_crosses_annotated_junction"] else "same_exon" if a["entire_amplicon_within_one_exon"] else "unmapped"
        counts[c["gene_id"]][category] += 1
    method = json.loads((directory / "manifest.json").read_text())["method"]
    for gene in {r["gene_id"] for r in rows}:
        category_rows.append({"method": method, "gene_id": gene, **{c: counts[gene][c] for c in categories}})
write_json(base / "exon_partition.json", {"at": utc_now(), "seed": 20261013, "rows": category_rows,
    "annotation_sha256": sha256(ROOT / "results/selected_exon_context_audit.per_candidate.jsonl")})
fig, ax = plt.subplots(figsize=(8.3, 4.2), layout="constrained")
bottom = np.zeros(5)
for key, color, label in [("primer_junction", "#3d7aa8", "Primer spans junction"),
                           ("amplicon_junction", "#37926c", "Amplicon spans junction only"), ("polyA", "#db9c27", "Terminal poly(A)")]:
    data = [np.mean([r[key] for r in category_rows if r["method"] == m]) for m in methods]
    ax.bar(range(5), data, bottom=bottom, color=color, label=label, width=.65)
    bottom += np.asarray(data)
ax.set_xticks(range(5), [labels[m] for m in methods], rotation=20, ha="right")
ax.set_ylabel("RNA-specific pairs without genomic product flag")
ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.02), ncol=1, frameon=False, fontsize=8)
save(fig, "06_exon_structure")

fig, ax = plt.subplots(figsize=(11, 4.4))
ax.set(xlim=(0, 11), ylim=(0, 4.4))
ax.axis("off")
boxes = [(0.15, 1.5, "Fixed transcript\n+ design conditions\n+ optional primer anchor", "#e7f2ee"),
         (3.1, 1.5, "Conditional denoiser\n0 / F / R position labels\nabsorbing-mask diffusion", "#e9eef6"),
         (6.05, 1.5, "Interval decoder\n+ joint Tm guidance\n+ anchor preservation", "#eee8f4"),
         (9.0, 1.5, "Candidate pairs\nRNA / genome\nassessment", "#f6eedf")]
for x, y, label, color in boxes:
    width = 2.5 if x < 9 else 1.8
    ax.add_patch(FancyBboxPatch((x, y), width, 1.3, boxstyle="round,pad=.12", facecolor=color, edgecolor="#7d9595"))
    ax.text(x+width/2, y+.65, label, ha="center", va="center", fontsize=10)
    if x < 9:
        ax.annotate("", xy=(x+2.9, y+.65), xytext=(x+width+.15, y+.65), arrowprops={"arrowstyle": "->", "color": "#506d73"})
ax.text(5.5, 3.7, "DiffuPrimer: conditional position generation on a fixed template", ha="center", fontsize=14, weight="bold")
ax.text(5.5, .65, "De novo pairs  |  Multiple alternatives  |  Length conditioning  |  Fixed-primer completion  |  Changed-condition redesign", ha="center", fontsize=10)
save(fig, "00_method_overview")
write_json(out / "provenance.json", {"at": utc_now(), "seed": 20261013,
    "summary_sha256": sha256(base / "study_summary.json"), "records_sha256": sha256(base / "study_summary.per_gene.jsonl"),
    "plotter_sha256": sha256(Path(__file__)), "values": values,
    "figures": {p.name: sha256(p) for p in sorted(out.glob("*.pdf"))}})
print(out)
