#!/usr/bin/env python3
"""Gene-level analysis of the fixed-seed final study and candidate quality."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference
from summarize_secondary import overlap_distance

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--existing-only", action="store_true")
args = parser.parse_args()
output = ROOT / "results/final_fixed_seed"
output.mkdir(exist_ok=True)
seed = 20261013
targets = {r["gene_id"]: r for r in read_jsonl(ROOT / "data/processed/targets_pilot.jsonl") if r["split"] == "test"}
sources = json.loads((ROOT / "results/temperature_tuning/selected_conditions_summary.json").read_text())["inputs"]
directories = [Path(r["directory"]) for r in sources
               if json.loads((Path(r["directory"]) / "manifest.json").read_text())["arguments"]["seed"] == seed]
if not args.existing_only:
    stage = json.loads((ROOT / "experiments/final_followup/results/status.json").read_text())["stage"]
    if stage != "followup_experiments_complete":
        raise ValueError("Follow-up experiments must complete before the final analysis")
    directories += sorted((ROOT / "experiments/final_followup/results/evaluations").glob("*"))
    corrected = any("anchor_prefilter" in p.name for p in directories)
    if corrected:
        # Replace the affected recipe, while retaining its first implementation on disk.
        directories = [p for p in directories if not (
            "anchor_prefilter" not in p.name and
            (("_fix_reverse_" in p.name and p.name.startswith(("direct_", "diffusion_")))
             or p.name in [f"{m}_redesign_reverse_seed{seed}" for m in ["direct", "diffusion"]]))]
rows, inputs = [], []
for directory in directories:
    manifest = json.loads((directory / "manifest.json").read_text())
    completion = json.loads((directory / "completion.json").read_text())
    status = json.loads((directory / "specificity_bowtie_status.json").read_text())
    candidates = read_jsonl(directory / "candidates.jsonl")
    screened = read_jsonl(directory / "specificity_bowtie.jsonl")
    ledger = read_jsonl(directory / "per_gene.jsonl")
    screens = {r["candidate_id"]: r for r in screened}
    if not completion["all_expected_targets_evaluated"] or {r["gene_id"] for r in ledger} != set(targets):
        raise ValueError("Incomplete gene cohort")
    if completion["candidates_sha256"] != sha256(directory / "candidates.jsonl") or status["output_sha256"] != sha256(directory / "specificity_bowtie.jsonl"):
        raise ValueError("Result source changed")
    by_gene = defaultdict(list)
    for c in candidates:
        by_gene[c["gene_id"]].append(c)
    own_args = manifest["arguments"]
    new = "task" in manifest
    method = manifest["method"]
    if new:
        context = ("redesign_" + manifest["fixed_side"] if manifest["task"] == "redesign" else
                   "amp80-200_fix_" + manifest["fixed_side"] if manifest["task"] == "completion" else "amp80-200")
        variant = manifest["variant"]
    else:
        context = f"amp{own_args['amplicon_min']}-{own_args['amplicon_max']}"
        if own_args["repair"]:
            context += "_fix_" + own_args["repair"]
        variant = "steps" + str(own_args["steps"]) if own_args["steps"] else "full"
    ledger_by_gene = {r["gene_id"]: r for r in ledger}
    for gene in sorted(targets):
        for budget in [1, 2, 5, 10, 20]:
            own = by_gene[gene][:budget]
            feasible = [c for c in own if c["assessment"]["computational_feasible"]]
            rna = [c for c in feasible if screens[c["candidate_id"]].get("gene_specificity_screen_pass")]
            no_genome = [c for c in rna if not screens[c["candidate_id"]]["genomic_DNA_product_detected"]]
            distinct = lambda cs: len({(c["forward_primer"], c["reverse_primer"]) for c in cs})
            count, retained = distinct(rna), distinct(no_genome)
            scores = [c["assessment"]["shared_temperature_penalty"] for c in rna]
            anchor = [c for c in own if c["structural_valid"]]
            length_pass = [c for c in feasible if c.get("requested_amplicon_min", own_args.get("amplicon_min", 80)) <= c["amplicon_length"] <= c.get("requested_amplicon_max", own_args.get("amplicon_max", 200))]
            record = {"gene_id": gene, "gene_symbol": targets[gene]["gene_symbol"], "method": method,
                "context": context, "variant": variant, "budget": budget, "seed": seed,
                "unique_feasible_pairs": distinct(feasible), "unique_RNA_pass_pairs": count,
                "unique_RNA_pass_pairs_without_genomic_flag": retained,
                "RNA_no_genome_fraction": retained / count if count else None,
                "one_RNA_pass_pair": float(count >= 1), "two_RNA_pass_pairs": float(count >= 2),
                "one_no_genome_pair": float(retained >= 1), "two_no_genome_pairs": float(retained >= 2),
                "RNA_pass_amplicon_distance": overlap_distance(rna),
                "best_RNA_pass_temperature_penalty": min(scores) if scores else None,
                "budget_feasibility_rate": len(feasible) / budget,
                "geometry_valid_rate": sum(c["structural_valid"] for c in own) / budget,
                "length_compliant_feasible_rate": len(length_pass) / budget,
                "anchor_preservation_rate": sum(c.get("anchor_preserved", True) for c in anchor) / len(anchor) if anchor else None,
                "run": str(directory), "original_amplicon_length": ledger_by_gene[gene].get("original_amplicon_length"),
                "conditions": ledger_by_gene[gene].get("conditions", manifest["conditions"])}
            rows.append(record)
    inputs.append({"directory": str(directory), **{name + "_sha256": sha256(directory / name) for name in
        ["manifest.json", "completion.json", "candidates.jsonl", "specificity_bowtie.jsonl", "specificity_bowtie_status.json"]}})
keys = ["unique_feasible_pairs", "unique_RNA_pass_pairs", "unique_RNA_pass_pairs_without_genomic_flag",
        "RNA_no_genome_fraction", "one_RNA_pass_pair", "two_RNA_pass_pairs", "one_no_genome_pair", "two_no_genome_pairs",
        "RNA_pass_amplicon_distance", "best_RNA_pass_temperature_penalty", "budget_feasibility_rate", "geometry_valid_rate",
        "length_compliant_feasible_rate", "anchor_preservation_rate"]
groups = defaultdict(list)
for r in rows:
    groups[r["context"], r["method"], r["variant"], r["budget"]].append(r)
summaries, contrasts = [], []
for (context, method, variant, budget), own in sorted(groups.items()):
    if len(own) != 75:
        raise ValueError("Duplicate or missing gene rows")
    group = {"context": context, "method": method, "variant": variant, "budget": budget, "seed": seed, "total_genes": 75}
    for key in keys:
        values = [r[key] for r in own if r[key] is not None]
        group[key] = float(np.mean(values)) if values else None
        group[key + "_n"] = len(values)
    summaries.append(group)
for (context, method, variant, budget), own in sorted(groups.items()):
    if method != "diffusion" or variant != "full":
        continue
    left = {r["gene_id"]: r for r in own}
    rights = [("direct", "full"), ("primer3", "full"), ("uniform_pair_sampler", "full"), ("rule_sampler", "full")]
    if context == "amp80-200":
        rights += [("diffusion", v) for v in ["steps1", "steps4", "steps8", "single_oligo", "geometry_only", "raw"]]
    for other_method, other_variant in rights:
        right = groups.get((context, other_method, other_variant, budget))
        if not right:
            continue
        other = {r["gene_id"]: r for r in right}
        for key in keys:
            genes = sorted(g for g in targets if left[g][key] is not None and other[g][key] is not None)
            if genes:
                estimate = paired_bootstrap_difference([left[g][key] for g in genes], [other[g][key] for g in genes], seed=seed)
                contrasts.append({"context": context, "left": "diffusion_full", "right": other_method + "_" + other_variant,
                    "budget": budget, "endpoint": key, "total_genes": 75, **estimate})
filename = "existing_quality" if args.existing_only else "study_summary"
write_jsonl(output / (filename + ".per_gene.jsonl"), rows)
write_json(output / (filename + ".json"), {"at": utc_now(), "seed": seed, "source_runs": len(inputs),
    "inputs": inputs, "summaries": summaries, "paired_comparisons": contrasts,
    "analysis_unit": "target gene; one fixed training and generation seed",
    "quality_fraction": "unique RNA-pass pairs without genomic flag / unique RNA-pass pairs; defined only when at least one RNA-pass pair exists",
    "additional_experiments_complete": not args.existing_only})
for label, records in [("groups", summaries), ("contrasts", contrasts)]:
    with (output / (filename + "." + label + ".csv")).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in records for k in r)))
        writer.writeheader()
        writer.writerows([{k: json.dumps(v) if isinstance(v, list) else v for k, v in r.items()} for r in records])
print(json.dumps({"source_runs": len(inputs), "groups": len(summaries), "contrasts": len(contrasts), "output": str(output / (filename + ".json"))}))
