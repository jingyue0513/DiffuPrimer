#!/usr/bin/env python3
"""Matched gene-level endpoints for original and joint-guidance studies."""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference
from summarize_secondary import overlap_distance


parser = argparse.ArgumentParser()
parser.add_argument("directories", nargs="+")
parser.add_argument("--screen", choices=["bowtie", "blast"], default="bowtie")
parser.add_argument("--split", choices=["test", "validation"], default="test")
parser.add_argument("--output", default="results/bowtie_audit/summary.json")
parser.add_argument("--selection", help="Frozen validation-only temperature selection for a matched final comparison")
parser.add_argument("--rna-only", action="store_true", help="Use separately completed RNA screens; genomic endpoints remain missing")
parser.add_argument("--selected-conditions", action="store_true", help="Include the frozen condition/completion/step tasks at selected temperatures")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
expected = {r["gene_id"] for r in read_jsonl(root / "data/processed/targets_pilot.jsonl") if r["split"] == args.split}
selection = json.loads((root / args.selection).read_text()) if args.selection else None
if selection and args.split != "test":
    raise ValueError("Selected parameter comparison requires the test partition")
if args.selected_conditions and not selection:
    raise ValueError("Selected condition tasks require frozen validation-only selection")
rows, inputs = [], []
suffix = "_bowtie" if args.screen == "bowtie" else ""
if args.rna_only:
    if args.screen != "bowtie":
        raise ValueError("RNA-only completion is defined for the independent Bowtie audit")
    suffix = "_bowtie_rna"
for raw in args.directories:
    directory = root / raw
    manifest = json.loads((directory / "manifest.json").read_text())
    completion = json.loads((directory / "completion.json").read_text())
    arguments = manifest["arguments"]
    if manifest["stage"] == "calibration" or arguments["split"] != args.split or not completion["all_expected_targets_evaluated"]:
        raise ValueError("Incomplete or wrong-partition final evaluation")
    for key in ["per_gene", "candidates"]:
        if sha256(directory / (key + ".jsonl")) != completion[key + "_sha256"]:
            raise ValueError("Evaluation changed")
    records = read_jsonl(directory / "per_gene.jsonl")
    if {r["gene_id"] for r in records} != expected or len(records) != len(expected):
        raise ValueError("Frozen gene denominator changed")
    screen_path = directory / f"specificity{suffix}.jsonl"
    screen_status = json.loads((directory / f"specificity{suffix}_status.json").read_text())
    if screen_status["candidate_sha256"] != completion["candidates_sha256"]:
        raise ValueError("Stale specificity output")
    if "output_sha256" in screen_status and sha256(screen_path) != screen_status["output_sha256"]:
        raise ValueError("Specificity output changed")
    screens = {r["candidate_id"]: r for r in read_jsonl(screen_path)}
    candidates = read_jsonl(directory / "candidates.jsonl")
    if set(screens) != {r["candidate_id"] for r in candidates}:
        raise ValueError("Incomplete candidate screening ledger")
    by_gene = defaultdict(list)
    for row in candidates:
        by_gene[row["gene_id"]].append(row)
    study = "joint_tm" if arguments.get("joint_tm_guidance", False) else "original"
    if arguments.get("temperature", 1.0) != 1.0:
        study += "_temperature" + str(arguments["temperature"])
    if selection:
        method = manifest["method"]
        if method in ["direct", "diffusion"]:
            if arguments["temperature"] != selection["selected"][method]["temperature"]:
                raise ValueError("Run temperature differs from validation-only selection")
            if not arguments.get("joint_tm_guidance"):
                raise ValueError("Selected learned run must use matched joint guidance")
        elif method in ["rule_sampler", "uniform_pair_sampler"] and not arguments.get("joint_tm_guidance"):
            raise ValueError("Selected rule baseline must use matched joint guidance")
        if arguments["raw"] or arguments["budget"] != 20:
            raise ValueError("Selected comparisons require the frozen guided K20 setting")
        if args.selected_conditions:
            if (arguments["amplicon_min"],arguments["amplicon_max"]) not in [(80,200),(80,120),(121,160),(161,200)]:
                raise ValueError("Undeclared selected condition range")
            if arguments["repair"] and (arguments["amplicon_min"],arguments["amplicon_max"]) != (80,200):
                raise ValueError("Undeclared selected repair condition")
            if arguments["steps"] and (method != "diffusion" or arguments["steps"] not in [1,4,8] or arguments["repair"]):
                raise ValueError("Undeclared selected step setting")
        elif arguments["amplicon_min"] != 80 or arguments["amplicon_max"] != 200 or arguments["repair"] or arguments["steps"]:
            raise ValueError("Selected primary endpoint requires the frozen broad-band K20 setting")
        study = "temperature_selected"
    context = f"amp{arguments['amplicon_min']}-{arguments['amplicon_max']}"
    if arguments["repair"]:
        context += "_fix_" + arguments["repair"]
    variant = "raw" if arguments["raw"] else "steps" + str(arguments["steps"]) if arguments["steps"] else "geometry_only" if manifest["method"] == "diffusion" and not arguments["oligo_guidance"] else "primary"
    for record in records:
        own = by_gene[record["gene_id"]]
        for budget in [1, 2, 5, 10, arguments["budget"]]:
            if budget > arguments["budget"]:
                continue
            prefix = own[:budget]
            feasible = [c for c in prefix if c["assessment"]["computational_feasible"]]
            specific = [c for c in feasible if screens[c["candidate_id"]].get("gene_specificity_screen_pass")]
            safe = [] if args.rna_only else [c for c in specific if not screens[c["candidate_id"]]["genomic_DNA_product_detected"]]
            unique = {(c["forward_primer"], c["reverse_primer"]) for c in specific}
            unique_safe = {(c["forward_primer"], c["reverse_primer"]) for c in safe}
            row = {"study": study, "method": manifest["method"], "context": context, "variant": variant,
                "gene_id": record["gene_id"], "seed": arguments["seed"], "budget": budget,
                "unique_feasible_pairs": len({(c["forward_primer"], c["reverse_primer"]) for c in feasible}),
                "budget_feasibility_rate": len(feasible)/budget,
                "unique_RNA_pass_pairs": len(unique), "one_RNA_pass_pair": bool(unique), "two_RNA_pass_pairs": len(unique) >= 2,
                "unique_RNA_pass_pairs_without_genomic_flag": None if args.rna_only else len(unique_safe),
                "two_RNA_pass_pairs_without_genomic_flag": None if args.rna_only else len(unique_safe) >= 2,
                "RNA_pass_amplicon_distance": overlap_distance(specific),
                "distinct_RNA_pass_amplicons": len({(c["forward_start"], c["reverse_end"]) for c in specific}),
                "best_RNA_pass_temperature_penalty": min((c["assessment"]["shared_temperature_penalty"] for c in specific), default=None)}
            rows.append(row)
    inputs.append({"directory": str(directory), "completion_sha256": sha256(directory / "completion.json"),
        "specificity_sha256": sha256(screen_path), "specificity_status_sha256": sha256(directory / f"specificity{suffix}_status.json")})
keys = ["unique_feasible_pairs", "budget_feasibility_rate", "unique_RNA_pass_pairs", "one_RNA_pass_pair", "two_RNA_pass_pairs",
        "unique_RNA_pass_pairs_without_genomic_flag", "two_RNA_pass_pairs_without_genomic_flag",
        "RNA_pass_amplicon_distance", "distinct_RNA_pass_amplicons", "best_RNA_pass_temperature_penalty"]
groups = defaultdict(list)
for row in rows:
    groups[tuple(row[k] for k in ["study", "method", "context", "variant", "budget", "gene_id"])].append(row)
averaged = []
for (study, method, context, variant, budget, gene), records in sorted(groups.items()):
    seeds = [r["seed"] for r in records]
    if len(seeds) != len(set(seeds)):
        raise ValueError("Duplicate study/method/context/budget/gene/seed")
    row = {"study": study, "method": method, "context": context, "variant": variant, "budget": budget,
           "gene_id": gene, "seeds": sorted(seeds)}
    for key in keys:
        values = [r[key] for r in records if r[key] is not None]
        row[key] = float(np.mean(values)) if values else None
        row[key + "_seed_denominator"] = len(values)
    averaged.append(row)
grouped = defaultdict(list)
for row in averaged:
    grouped[tuple(row[k] for k in ["study", "method", "context", "variant", "budget"])].append(row)
summaries = []
for (study, method, context, variant, budget), records in sorted(grouped.items()):
    if {r["gene_id"] for r in records} != expected:
        raise ValueError("Incomplete gene-level endpoint group")
    counts = [len(r["seeds"]) for r in records]
    if min(counts) != max(counts):
        raise ValueError("Incomplete seed accounting")
    summary = {"study": study, "method": method, "context": context, "variant": variant, "budget": budget,
               "genes": len(records), "seeds": records[0]["seeds"]}
    for key in keys:
        values = [r[key] for r in records if r[key] is not None]
        summary[key] = float(np.mean(values)) if values else None
        summary[key + "_gene_denominator"] = len(values)
    summaries.append(summary)
comparisons = []
for (study, method, context, variant, budget), left_rows in grouped.items():
    if method != "diffusion" or variant != "primary":
        continue
    for (other_study, other_method, other_context, other_variant, other_budget), right_rows in grouped.items():
        if (context, budget) != (other_context, other_budget):
            continue
        matched = study == other_study
        shared_baseline = other_method == "primer3"
        guidance_effect = study == "joint_tm" and other_study == "original" and other_method == "diffusion" and other_variant == "primary"
        if not (matched or shared_baseline or guidance_effect) or (study, method, variant) == (other_study, other_method, other_variant):
            continue
        left, right = {r["gene_id"]: r for r in left_rows}, {r["gene_id"]: r for r in right_rows}
        comparison_seeds = None
        if other_method == method and other_variant != variant:
            left_seeds, right_seeds = set(left_rows[0]["seeds"]), set(right_rows[0]["seeds"])
            if left_seeds != right_seeds:
                comparison_seeds = sorted(left_seeds & right_seeds)
                if not comparison_seeds:
                    raise ValueError("Ablation has no matching primary model seed")
                for current,current_study,current_variant in [(left,study,variant),(right,other_study,other_variant)]:
                    for gene in current:
                        own=groups[(current_study,method,context,current_variant,budget,gene)]
                        current[gene]={key:float(np.mean(values)) if values else None for key in keys
                                       for values in [[r[key] for r in own if r['seed'] in comparison_seeds and r[key] is not None]]}
        for key in keys:
            genes = sorted(g for g in expected if left[g][key] is not None and right[g][key] is not None)
            if not genes:
                continue
            estimate = paired_bootstrap_difference([left[g][key] for g in genes], [right[g][key] for g in genes], seed=20261013)
            comparisons.append({"left_study": study, "left": method + "_" + variant, "right_study": other_study,
                "right": other_method + "_" + other_variant, "context": context, "budget": budget, "endpoint": key,
                "total_target_genes": len(expected), "matched_ablation_seeds": comparison_seeds, **estimate})
output = root / args.output
write_jsonl(output.with_suffix(".per_seed.jsonl"), rows)
write_jsonl(output.with_suffix(".per_gene.jsonl"), averaged)
write_json(output, {"generated_at": utc_now(), "split": args.split, "screen": args.screen, "inputs": inputs,
    "RNA_only": args.rna_only,
    "validation_selection_sha256": sha256(root / args.selection) if args.selection else None,
    "summaries": summaries, "paired_comparisons": comparisons,
    "policy": "All frozen targets included. Seeds averaged within gene, 10000 paired gene bootstraps. Diversity and quality explicitly condition on available pairs. Prefix budgets reuse original proposal order."})
print(json.dumps({"output": str(output), "groups": len(summaries), "comparisons": len(comparisons)}), flush=True)
