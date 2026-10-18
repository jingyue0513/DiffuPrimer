#!/usr/bin/env python3
"""Publication-sized figures from the completed, source-verified attention study."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from diffuprimer.io import sha256, utc_now, write_json

ORDER = ["primer3", "qprimer", "primercast", "uniform_pair_sampler", "rule_sampler", "cnn_direct",
         "cnn_diffusion", "attention_rope_direct", "attention_rope_diffusion",
         "attention_direct", "attention_diffusion", "attention_hybrid_direct", "attention_hybrid_diffusion"]
LABELS = {"primer3": "Primer3", "qprimer": "qPrimer", "primercast": "PrimerCast", "uniform_pair_sampler": "Uniform pairs",
          "rule_sampler": "Forward-first", "cnn_direct": "CNN direct", "cnn_diffusion": "CNN diffusion",
          "attention_direct": "Attention direct", "attention_diffusion": "DiffuPrimer (attention)"}
LABELS.update({"attention_rope_direct":"RoPE direct", "attention_rope_diffusion":"RoPE diffusion",
               "attention_hybrid_direct":"Hybrid direct", "attention_hybrid_diffusion":"DiffuPrimer (hybrid)"})
COLORS = {"primer3": "#C29B48", "qprimer": "#DA7F64", "primercast": "#B78094", "uniform_pair_sampler": "#ACB4BE",
          "rule_sampler": "#6D8294", "cnn_direct": "#7CB2C9", "cnn_diffusion": "#7581B8",
          "attention_direct": "#80BBB5", "attention_diffusion": "#21877F"}
COLORS.update({"attention_rope_direct":"#A8A1C3", "attention_rope_diffusion":"#8B76A8",
               "attention_hybrid_direct":"#80BBB5", "attention_hybrid_diffusion":"#21877F"})
ENDPOINT_LABELS = {
    "unique_RNA_pass_pairs_without_genomic_flag": "RNA-specific pairs without genomic products",
    "RNA_no_genome_fraction": "Retained fraction of RNA-specific pairs (%)",
    "unique_RNA_pass_pairs": "Unique RNA-specific primer pairs",
    "two_no_genome_pairs": "Targets with ≥2 retained pairs (%)",
    "RNA_pass_amplicon_distance": "Amplicon Jaccard distance",
    "best_RNA_pass_temperature_penalty": "Best paired Tm deviation (°C)",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default="experiments/attention/results/study_summary.json")
    parser.add_argument("--output", default="experiments/attention/results/figures")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--primary-method", help="Primary diffusion label; defaults to declared protocol/summary primary_methods")
    parser.add_argument("--direct-method", help="Matched direct predictor label")
    parser.add_argument("--reference-summary", default="experiments/attention/results/study_summary.json", help="For external_test only: separate original-cohort summary for generalization panels")
    args = parser.parse_args()
    summary_path = ROOT / args.summary
    summary_bytes = summary_path.read_bytes()
    summary = json.loads(summary_bytes)
    summary_sha = hashlib.sha256(summary_bytes).hexdigest()
    if summary["status"] != "complete" and not args.allow_partial:
        parser.error("Study is incomplete; --allow-partial renders completed sources with explicit status")
    out = ROOT / args.output
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "axes.titlesize": 8,
                         "axes.labelsize": 7, "xtick.labelsize": 6.8, "ytick.labelsize": 7,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": .6, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.facecolor": "white", "figure.facecolor": "white"})
    groups = {(s["context"], s["method"], s["variant"], s["budget"]): s for s in summary["summaries"]}
    comparisons = {(c["context"], c["left_method"], c["right_method"], c["right_variant"], c["budget"], c["endpoint"]): c
                   for c in summary["paired_comparisons"]}
    observed = {g["method"] for g in summary["summaries"]}
    primary = args.primary_method or summary.get("primary_methods", {}).get("diffusion") or next(
        (m for m in ("attention_hybrid_diffusion", "attention_rope_diffusion", "attention_diffusion") if m in observed), "attention_diffusion")
    direct = args.direct_method or summary.get("primary_methods", {}).get("direct") or primary.removesuffix("diffusion") + "direct"
    values, generated, absent, reference_meta = [], [], [], None
    partial = summary["status"] != "complete"

    def get(method, key, context="amp80-200", variant="full", budget=20):
        g = groups.get((context, method, variant, budget))
        if g is None or g[key] is None:
            return None
        scale = 100 if key in ("RNA_no_genome_fraction", "two_no_genome_pairs", "geometry_valid_rate") else 1
        value = (g[key] * scale, np.asarray(g[key + "_ci95"]) * scale, g[key + "_n"])
        values.append({"method": method, "endpoint": key, "context": context, "variant": variant,
                       "budget": budget, "mean": value[0], "ci95": value[1].tolist(), "n": value[2], "scale": scale})
        return value

    def dots(ax, key, methods, context="amp80-200", labels=True):
        available = [(m, get(m, key, context)) for m in methods]
        available = [(m, value) for m, value in available if value is not None]
        for i, (m, (mean, ci, n)) in enumerate(available):
            ax.plot(ci, [i, i], color=COLORS[m], lw=1.5, solid_capstyle="round")
            ax.scatter(mean, i, c=COLORS.get(m,"#647B87"), s=23 if m == primary else 16, zorder=4,
                       edgecolors="white", linewidths=.4)
        ax.set_yticks(range(len(available)), [LABELS[m] if labels else "" for m, _ in available])
        ax.invert_yaxis()
        ax.set_xlim(left=0)
        ax.grid(axis="x", color="#E4E8EB", lw=.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        return available

    def save(fig, name, caption, alt, fixed_width=False):
        if partial:
            fig.suptitle("ONGOING STUDY · completed sources only", color="#9B6536", fontsize=8, y=1.07)
        files = []
        for extension in ("pdf", "svg", "png"):
            path = out / (name + "." + extension)
            save_options = {} if fixed_width else {"bbox_inches": "tight", "pad_inches": .04}
            fig.savefig(path, dpi=350, **save_options)
            files.append({"path": str(path.relative_to(ROOT)), "sha256": sha256(path)})
        plt.close(fig)
        generated.append({"name": name, "files": files, "caption": caption, "alt_text": alt})

    # Three compact quantitative panels, sized for a 178 mm double column.
    method_order = ORDER + sorted(observed - set(ORDER))
    for m in observed:
        LABELS.setdefault(m, m.replace("_", " "))
        COLORS.setdefault(m,"#647B87")
    present = [m for m in method_order if ("amp80-200", m, "full", 20) in groups]
    if present:
        fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.8), layout="constrained")
        keys = ["unique_RNA_pass_pairs_without_genomic_flag", "RNA_no_genome_fraction", "unique_RNA_pass_pairs"]
        titles = ["A  Retained candidate count", "B  Retained fraction", "C  RNA-specific supply"]
        for index, (ax, key, title) in enumerate(zip(axes, keys, titles)):
            dots(ax, key, present, labels=index == 0)
            ax.set_title(title, loc="left", fontsize=7.5)
            ax.set_xlabel("Pairs per gene" if index != 1 else "Retained fraction (%)")
        save(fig, "01_candidate_quality", "Full-reference primer-pair quality at K=20. Points indicate target-gene means; lines show 95% gene-bootstrap intervals. Count endpoints retain all test genes. Retained fraction is defined only for genes with RNA-specific pairs; denominators are supplied in figure_values.json.",
             "Three horizontal comparison plots show retained primer-pair counts, retained fraction, and RNA-specific counts for the methods with completed test sources.")

    neural = [m for m in method_order if m.endswith(("_direct", "_diffusion")) and ("amp80-200", m, "full", 20) in groups]
    if primary in neural:
        fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.4), layout="constrained")
        keys = ["unique_RNA_pass_pairs_without_genomic_flag", "RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty"]
        for ax, key, title in zip(axes[0], keys, ["A  Retained candidates", "B  Amplicon diversity", "C  Paired Tm deviation"]):
            dots(ax, key, neural)
            ax.set_title(title, fontsize=7.1, loc="left")
            ax.set_xlabel("Pairs per gene" if key == keys[0] else "Distance" if key == keys[1] else "°C")
        primary_label = "Hybrid" if "hybrid" in primary else "RoPE" if "rope" in primary else "Attention"
        paired = [(primary, "cnn_diffusion", f"{primary_label} − CNN\n(diffusion)"),
                  (direct, "cnn_direct", f"{primary_label} − CNN\n(direct)"),
                  (primary, direct, f"Diffusion − direct\n({primary_label.lower()})")]
        if primary == "attention_hybrid_diffusion" and "attention_rope_diffusion" in neural:
            paired.append((primary,"attention_rope_diffusion","Hybrid − RoPE\n(diffusion)"))
        for ax, key, letter in zip(axes[1], keys, "DEF"):
            available_pairs = [p for p in paired if ("amp80-200", p[0], p[1], "full", 20, key) in comparisons]
            for i, (left, right, label) in enumerate(available_pairs):
                c = comparisons.get(("amp80-200", left, right, "full", 20, key))
                if c:
                    ax.plot(c["ci95"], [i, i], color=COLORS[left], lw=1.4)
                    ax.scatter(c["mean_difference"], i, c=COLORS[left], s=18, zorder=4)
            ax.axvline(0, color="#A9B0B6", lw=.8, ls="--")
            ax.set_yticks(range(len(available_pairs)), [p[2] for p in available_pairs])
            ax.invert_yaxis()
            ax.set_title(f"{letter}  Paired difference", loc="left")
            if not available_pairs:
                ax.text(.5,.5,"Paired sources pending",transform=ax.transAxes,ha="center",va="center",color="#7B858B",fontsize=7)
            ax.set_xlabel("Paired mean difference")
            ax.grid(axis="x", color="#E4E8EB", lw=.5)
            ax.tick_params(axis="y", length=0, labelsize=6.5)
            ax.spines["left"].set_visible(False)
        save(fig, "02_architecture_comparison", "Matched neural architectures and generation mechanisms at K=20. Upper panels show gene-level means with marginal intervals. Lower panels show paired mean differences with gene-paired 95% bootstrap intervals. Diversity and temperature contrasts use jointly defined genes and their reported denominators. Positive diversity and negative temperature differences indicate improvements.",
             "Candidate retention, amplicon diversity and temperature deviation are compared for CNN and attention direct/diffusion models, followed by three paired architecture and generation-mechanism contrasts.")
    else:
        absent.append({"name": "02_architecture_comparison", "reason": "attention diffusion test source has not completed"})

    contexts = [("amp80-120", "A  Length: 80–120 bp"), ("amp121-160", "B  Length: 121–160 bp"),
                ("amp161-200", "C  Length: 161–200 bp"), ("amp80-200_fix_forward", "D  Retain F: completion"),
                ("amp80-200_fix_reverse", "E  Retain R: completion"), ("redesign_forward", "F  Retain F: redesign"),
                ("redesign_reverse", "G  Retain R: redesign")]
    if any((ctx, primary, "full", 20) in groups for ctx, _ in contexts):
        subset = ["primer3", "qprimer", "cnn_diffusion", direct, primary]
        # Shared method labels once per row preserve printable plot widths.
        # Explicit margins avoid constrained_layout reserving four copies of long labels.
        fig, axes = plt.subplots(2, 4, figsize=(178 / 25.4, 4.3))
        fig.subplots_adjust(left=.22, right=.99, bottom=.10, top=.90, wspace=.35, hspace=.74)
        destinations = [axes[0, 0], axes[0, 1], axes[0, 2], *axes[1]]
        titles = ["A  Length\n80–120 bp", "B  Length\n121–160 bp", "C  Length\n161–200 bp",
                  "D  Completion\nretain F", "E  Completion\nretain R",
                  "F  Redesign\nretain F", "G  Redesign\nretain R"]
        for index, (ax, (context, _title), title) in enumerate(zip(destinations, contexts, titles)):
            dots(ax, "unique_RNA_pass_pairs_without_genomic_flag", subset, context,
                 labels=index in (0, 3))
            ax.set_title(title, loc="left", fontsize=7.8, pad=7)
            ax.tick_params(axis="y", labelsize=7.2)
            ax.set_xlabel("Retained pairs / gene", fontsize=7.2)
        axes[0, 3].axis("off")
        axes[0, 3].text(.02, .95, "One model\nFlexible modes", va="top", fontsize=9, weight="bold", color="#254955")
        axes[0, 3].text(.02, .55, "Same genes\nSame budget\nSame screen", va="top", fontsize=7.5, linespacing=1.5, color="#526778")
        save(fig, "03_flexible_design", "Length-conditioned design, fixed-primer completion and changed-length redesign at K=20. Retained pairs are unique RNA-specific pairs without predicted genomic products. Counts include all held-out genes, including targets with no returned candidates. F and R refer to the preserved primer side. All attention tasks reuse the same trained weights.",
             "Seven panels compare reference-screened primer-pair counts across three product-length bands, completion with either primer retained, and changed-length redesign with either primer retained.", fixed_width=True)
    else:
        absent.append({"name": "03_flexible_design", "reason": "attention flexible-task test sources have not completed"})

    if primary in present:
        fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.4), layout="constrained")
        keys = ["unique_RNA_pass_pairs_without_genomic_flag", "RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty"]
        for ax, key, title in zip(axes, keys, ["A  Retained candidates", "B  Amplicon diversity", "C  Paired Tm deviation"]):
            for method in [m for m in present if m in ("primer3", "qprimer", "primercast", "cnn_diffusion", direct, primary)]:
                points = [(budget, get(method, key, budget=budget)) for budget in (1, 2, 5, 10, 20)]
                points = [(b, p) for b, p in points if p is not None]
                if not points:
                    continue
                xs = [p[0] for p in points]
                ys = [p[1][0] for p in points]
                ax.plot(xs, ys, "o-", color=COLORS[method], label=LABELS[method], markersize=2.7, lw=1.1)
                ax.fill_between(xs, [p[1][1][0] for p in points], [p[1][1][1] for p in points], color=COLORS[method], alpha=.1)
            ax.set_title(title, fontsize=7.1, loc="left")
            ax.set(xlabel="Proposal budget K", xticks=[1, 5, 10, 20], ylim=(0, None))
            ax.grid(axis="y", color="#E4E8EB", lw=.5)
        axes[0].legend(frameon=False, fontsize=5.8, loc="upper left")
        save(fig, "04_budget_and_diversity", "Recorded candidate prefixes at K=1, 2, 5, 10 and 20. Lines show gene means and shading shows 95% gene-bootstrap intervals. Invalid proposals, duplicates and shortages remain in K. Diversity requires at least two distinct amplicon intervals and temperature endpoints require RNA-specific pairs; denominators are supplied for every point.",
             "Three budget curves compare retained primer-pair supply, positional diversity and best paired temperature deviation for completed methods.")
    else:
        absent.append({"name": "04_budget_and_diversity", "reason": "attention diffusion test source has not completed"})

    # Independent cohorts have separate means/intervals; never pair different genes.
    if summary.get("evaluation_split") == "external_test":
        reference_path = ROOT / args.reference_summary
        if present and reference_path.is_file():
            reference_bytes = reference_path.read_bytes()
            reference = json.loads(reference_bytes)
            original = {(g["context"],g["method"],g["variant"],g["budget"]):g for g in reference["summaries"]}
            common = [m for m in present if ("amp80-200",m,"full",20) in original]
            if common:
                fig, axes = plt.subplots(2,2,figsize=(178 / 25.4,4.5))
                fig.subplots_adjust(left=.22,right=.98,bottom=.16,top=.91,wspace=.22,hspace=.46)
                for ax,key,title in zip(axes.flat,
                    ["unique_RNA_pass_pairs_without_genomic_flag","RNA_no_genome_fraction","RNA_pass_amplicon_distance","best_RNA_pass_temperature_penalty"],
                    ["A  Retained candidates","B  Retained fraction","C  Amplicon diversity","D  Paired Tm deviation"]):
                    for i,m in enumerate(common):
                        for offset,current,label,marker in [(-.13,original[("amp80-200",m,"full",20)],"Original held-out genes","s"),(.13,groups[("amp80-200",m,"full",20)],"New held-out genes","o")]:
                            if current[key] is None:
                                continue
                            scale=100 if key=="RNA_no_genome_fraction" else 1
                            mean=current[key]*scale
                            ci=np.asarray(current[key+"_ci95"])*scale
                            y=i+offset
                            ax.plot(ci,[y,y],color=COLORS[m],lw=1.1,alpha=.5 if offset<0 else 1)
                            ax.scatter(mean,y,edgecolors=COLORS[m],facecolors="white" if offset<0 else COLORS[m],marker=marker,s=15,zorder=4,linewidths=.7)
                            values.append({"cohort":label,"method":m,"endpoint":key,"budget":20,"mean":mean,"ci95":ci.tolist(),"n":current[key+"_n"],"total_genes":current["total_genes"]})
                    ax.set_yticks(range(len(common)),[LABELS[m] for m in common] if ax in (axes[0,0],axes[1,0]) else [])
                    ax.invert_yaxis();ax.set_xlim(left=0)
                    ax.set_title(title,loc="left")
                    ax.set_xlabel("Pairs per gene" if key=="unique_RNA_pass_pairs_without_genomic_flag" else "%" if key=="RNA_no_genome_fraction" else "Jaccard distance" if key=="RNA_pass_amplicon_distance" else "°C")
                    ax.grid(axis="x",color="#E4E8EB",lw=.5)
                    ax.spines["left"].set_visible(False);ax.tick_params(axis="y",length=0,labelsize=6.8)
                from matplotlib.lines import Line2D
                fig.legend(handles=[Line2D([],[],marker="s",linestyle="none",markerfacecolor="white",markeredgecolor="#647B87",label=f"Original cohort (N={reference['total_test_genes']})"),Line2D([],[],marker="o",linestyle="none",color="#647B87",label=f"New cohort (N={summary['total_test_genes']})")],frameon=False,fontsize=7.4,loc="lower center",bbox_to_anchor=(.57,.015),ncol=2)
                save(fig,"05_external_generalization",f"Separate original and new held-out gene cohorts (N={reference['total_test_genes']} and N={summary['total_test_genes']}). New genes use the same mouse RefSeq source and were frozen before candidate/model outcomes; they do not constitute a new species or independent assay dataset. Hollow squares and filled circles indicate cohort-specific means and gene-bootstrap intervals. No genes are pooled or paired across cohorts. Each method's within-cohort paired contrasts remain in its own summary; model/checkpoint identifiers are recorded in source manifests.",
                     "Four panels compare reference-screened candidate counts, retained fractions, positional diversity and temperature deviation on the original and new gene cohorts using separate means and confidence intervals.",fixed_width=True)
                reference_meta={"path":str(reference_path),"sha256":hashlib.sha256(reference_bytes).hexdigest(),"total_genes":reference["total_test_genes"],"comparison":"Separate cohort means; no cross-cohort gene pairing or pooling"}
            else:
                absent.append({"name":"05_external_generalization","reason":"No completed matching original/external method groups"})
        else:
            absent.append({"name":"05_external_generalization","reason":"External completed sources or original-cohort summary not yet available"})

    write_json(out / "figure_values.json", {"at": utc_now(), "status": summary["status"], "values": values})
    write_json(out / "provenance.json", {"at": utc_now(), "status": summary["status"],
                                         "summary_sha256": summary_sha, "source_runs": summary["source_runs"],
                                         "plot_script_sha256":sha256(Path(__file__)),
                                         "planned_new_runs": summary["planned_new_runs"],
                                         "pending_planned_runs": summary["pending_planned_runs"],
                                         "primary_method":primary,"direct_method":direct,
                                         "evaluation_split":summary.get("evaluation_split","test"),"total_genes":summary["total_test_genes"],"reference_cohort":reference_meta,
                                         "generated": generated, "not_generated": absent,
                                         "style": "White background, restrained method colors, 178 mm nominal width, editable vector PDF/SVG, gene-bootstrap intervals"})
    print(json.dumps({"status": summary["status"], "figures": [f["name"] for f in generated], "not_generated": absent, "output": str(out)}))


if __name__ == "__main__":
    main()
