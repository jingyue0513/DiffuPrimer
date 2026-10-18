#!/usr/bin/env python3
"""Recorded-run secondary endpoints; absent screens remain absent, never passes."""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference


def overlap_distance(candidates):
    intervals = sorted({(c["forward_start"], c["reverse_end"]) for c in candidates})
    distances = []
    for i, (start, end) in enumerate(intervals):
        for other_start, other_end in intervals[:i]:
            shared = max(0, min(end, other_end) - max(start, other_start))
            union = end - start + other_end - other_start - shared
            distances.append(1 - shared / union)
    return float(np.mean(distances)) if distances else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+")
    parser.add_argument("--output", default="results/secondary_summary.json")
    parser.add_argument("--screen", choices=["blast", "bowtie"], default="blast")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    expected = {r["gene_id"] for r in read_jsonl(root / "data/processed/targets_pilot.jsonl") if r["split"] == "test"}
    rows, inputs = [], []
    for raw in args.directories:
        directory = root / raw
        manifest = json.loads((directory / "manifest.json").read_text())
        completion = json.loads((directory / "completion.json").read_text())
        if manifest["stage"] == "calibration" or manifest["arguments"]["split"] != "test":
            raise ValueError("Secondary final summary requires held-out noncalibration runs")
        if not completion["all_expected_targets_evaluated"]:
            raise ValueError("Incomplete evaluation")
        for key in ["per_gene", "candidates"]:
            if sha256(directory / (key + ".jsonl")) != completion[key + "_sha256"]:
                raise ValueError("Evaluation changed after completion")
        per_gene = read_jsonl(directory / "per_gene.jsonl")
        if {r["gene_id"] for r in per_gene} != expected:
            raise ValueError("All frozen test genes must remain in the denominator")
        candidate_groups = defaultdict(list)
        for c in read_jsonl(directory / "candidates.jsonl"):
            candidate_groups[c["gene_id"]].append(c)
        suffix = "_bowtie" if args.screen == "bowtie" else ""
        screen_path = directory / f"specificity{suffix}.jsonl"
        screens = {r["candidate_id"]: r for r in read_jsonl(screen_path)} if screen_path.exists() else {}
        if screen_path.exists():
            screen_status = json.loads((directory / f"specificity{suffix}_status.json").read_text())
            if screen_status["candidate_sha256"] != completion["candidates_sha256"]:
                raise ValueError("Specificity screen is stale")
            if "output_sha256" in screen_status and sha256(screen_path) != screen_status["output_sha256"]:
                raise ValueError("Specificity output changed after screening")
            if set(screens) != {c["candidate_id"] for cs in candidate_groups.values() for c in cs}:
                raise ValueError("Specificity screen candidate accounting incomplete")
        arguments = manifest["arguments"]
        context = f"amp{arguments['amplicon_min']}-{arguments['amplicon_max']}"
        if arguments["repair"]:
            context += "_fix_" + arguments["repair"]
        if arguments["raw"]:
            variant = "raw"
        elif arguments["steps"] == 1:
            variant = "one_step"
        elif arguments["steps"] is not None:
            variant = "steps" + str(arguments["steps"])
        elif manifest["method"] == "diffusion" and not arguments["oligo_guidance"]:
            variant = "geometry_only"
        else:
            variant = "primary"
        for record in per_gene:
            own = candidate_groups[record["gene_id"]]
            feasible = [c for c in own if c["assessment"]["computational_feasible"]]
            specific = [c for c in feasible if screens.get(c["candidate_id"], {}).get("gene_specificity_screen_pass")]
            scores = [c["assessment"]["shared_temperature_penalty"] for c in specific]
            genome_screened = bool(specific) and all("genomic_DNA_product_detected" in screens[c["candidate_id"]] for c in specific)
            row = {**record, "run": directory.name, "seed": arguments["seed"], "context": context,
                "variant": variant, "geometry_per_budget": sum(c["assessment"].get("geometry_pass", False) for c in own) / arguments["budget"],
                "structure_per_budget": record["structural_valid_count"] / arguments["budget"],
                "specificity_assessed": screen_path.exists(), "specificity_screen": args.screen,
                "specificity_status": screen_status["status"] if screen_path.exists() else "not_assessed",
                "best_RNA_pass_temperature_penalty": min(scores) if scores else None,
                "RNA_pass_amplicon_distance": overlap_distance(specific) if screen_path.exists() else None,
                "unique_RNA_pass_pairs": len({(c["forward_primer"], c["reverse_primer"]) for c in specific}) if screen_path.exists() else None,
                "unique_RNA_pass_distinct_amplicons": len({(c["forward_start"], c["reverse_end"]) for c in specific}) if screen_path.exists() else None,
                "any_RNA_pass_genomic_flag": any(screens[c["candidate_id"]].get("genomic_DNA_product_detected", False) for c in specific) if genome_screened else None,
                "any_RNA_pass_without_genomic_flag": any(not screens[c["candidate_id"]]["genomic_DNA_product_detected"] for c in specific) if genome_screened else (False if screen_path.exists() and not specific else None)}
            rows.append(row)
        inputs.append({"directory": str(directory), "manifest_sha256": sha256(directory / "manifest.json"),
            "completion_sha256": sha256(directory / "completion.json"), "specificity_sha256": sha256(screen_path) if screen_path.exists() else None})
    groups = defaultdict(list)
    for row in rows:
        groups[(row["context"], row["method"], row["variant"], row["gene_id"])].append(row)
    metrics = ["geometry_per_budget", "structure_per_budget", "budget_feasibility_rate", "unique_feasible_pairs",
        "distinct_feasible_amplicons", "mean_feasible_amplicon_jaccard_distance", "best_shared_temperature_penalty",
        "elapsed_seconds", "best_RNA_pass_temperature_penalty", "RNA_pass_amplicon_distance", "unique_RNA_pass_pairs",
        "unique_RNA_pass_distinct_amplicons", "any_RNA_pass_genomic_flag", "any_RNA_pass_without_genomic_flag"]
    averaged = []
    for (context, method, variant, gene), records in sorted(groups.items()):
        if len({r["seed"] for r in records}) != len(records):
            raise ValueError("Duplicate run/seed for one context")
        item = {"context": context, "method": method, "variant": variant, "gene_id": gene,
            "observed_seeds": [r["seed"] for r in records], "seed_count": len(records)}
        for key in metrics:
            values = [r[key] for r in records if r.get(key) is not None]
            item[key] = float(np.mean(values)) if values else None
            item[key + "_seed_denominator"] = len(values)
        averaged.append(item)
    grouped = defaultdict(list)
    for row in averaged:
        grouped[(row["context"], row["method"], row["variant"])].append(row)
    summaries = []
    for (context, method, variant), records in sorted(grouped.items()):
        summary = {"context": context, "method": method, "variant": variant, "genes": len(records),
            "seed_count_range": [min(r["seed_count"] for r in records), max(r["seed_count"] for r in records)]}
        for key in metrics:
            values = [r[key] for r in records if r[key] is not None]
            summary[key] = float(np.mean(values)) if values else None
            summary[key + "_gene_denominator"] = len(values)
        summaries.append(summary)
    comparisons = []
    for (context, method, variant), left_rows in grouped.items():
        if method != "diffusion" or variant != "primary":
            continue
        for (other_context, other_method, other_variant), right_rows in grouped.items():
            if context != other_context or (other_method == method and other_variant == variant):
                continue
            left = {r["gene_id"]: r for r in left_rows}
            right = {r["gene_id"]: r for r in right_rows}
            comparison_seeds = None
            if other_method == method and other_variant != variant:
                left_seeds, right_seeds = set(left_rows[0]["observed_seeds"]), set(right_rows[0]["observed_seeds"])
                if left_seeds != right_seeds:
                    comparison_seeds = sorted(left_seeds & right_seeds)
                    if not comparison_seeds:
                        raise ValueError("Ablation has no matching primary model seed")
                    for current, current_variant in [(left, variant), (right, other_variant)]:
                        for gene in current:
                            own = groups[(context, method, current_variant, gene)]
                            current[gene] = {key: float(np.mean(values)) if values else None for key in metrics
                                             for values in [[r[key] for r in own if r["seed"] in comparison_seeds and r.get(key) is not None]]}
            for key in ["best_shared_temperature_penalty", "best_RNA_pass_temperature_penalty",
                    "mean_feasible_amplicon_jaccard_distance", "RNA_pass_amplicon_distance", "unique_feasible_pairs"]:
                joint = sorted(g for g in set(left) & set(right) if left[g][key] is not None and right[g][key] is not None)
                if not joint:
                    continue
                estimate = paired_bootstrap_difference([left[g][key] for g in joint], [right[g][key] for g in joint], seed=20261013)
                comparisons.append({"context": context, "left": "diffusion_primary", "right": other_method + "_" + other_variant,
                    "endpoint": key, "joint_coverage_genes": len(joint), "total_test_genes": len(expected),
                    "matched_ablation_seeds": comparison_seeds, **estimate})
    output = root / args.output
    write_jsonl(output.with_suffix(".per_gene.jsonl"), averaged)
    write_json(output, {"generated_at": utc_now(), "specificity_screen": args.screen,
        "inputs": inputs, "summaries": summaries, "paired_comparisons": comparisons,
        "interpretation": "Exploratory secondary endpoints; seeds averaged within gene. Unequal-seed same-model ablations use the shared seeds for paired contrasts; pooled group summaries retain their stated seed counts. Missing diversity/score values have explicit restricted denominators. Genomic flags describe in-silico products, not contamination outcomes."})
    print(json.dumps({"summaries": len(summaries), "paired_comparisons": len(comparisons), "output": str(output)}))


if __name__ == "__main__":
    main()
