#!/usr/bin/env python3
"""Analyze completed attention/qPrimer runs together with the frozen CNN study.

Each source must contain the complete frozen test cohort and a completed full-
reference screen. --allow-partial permits missing *runs*, never a smaller gene
cohort or treating pending specificity assessments as negative results.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference, summarize_candidates
from summarize_secondary import overlap_distance

BUDGETS = (1, 2, 5, 10, 20)
ENDPOINTS = (
    "unique_feasible_pairs", "unique_RNA_pass_pairs",
    "unique_RNA_pass_pairs_without_genomic_flag", "RNA_no_genome_fraction",
    "one_RNA_pass_pair", "two_RNA_pass_pairs", "one_no_genome_pair",
    "two_no_genome_pairs", "RNA_pass_amplicon_distance",
    "best_RNA_pass_temperature_penalty", "budget_feasibility_rate",
    "geometry_valid_rate", "length_compliant_feasible_rate",
    "anchor_preservation_rate", "distinct_feasible_amplicons",
    "mean_feasible_amplicon_jaccard_distance", "best_shared_temperature_penalty",
)
REQUIRED = ("manifest.json", "completion.json", "candidates.jsonl",
            "per_gene.jsonl", "specificity_bowtie.jsonl",
            "specificity_bowtie_status.json")


class PendingRun(Exception):
    """A source has not finished generation and full-reference assessment."""


def resolve(raw):
    p = Path(raw)
    return p if p.is_absolute() else ROOT / p


def run_identity(manifest, old=False, aliases=None, directory=None):
    a = manifest["arguments"]
    method = manifest["method"]
    if old and method in ("diffusion", "direct"):
        method = "cnn_" + method
    elif method in ("attention_diffusion", "attention_direct"):
        architecture = str(manifest.get("architecture", "")).lower()
        if "hybrid" in architecture:
            method = method.replace("attention_", "attention_hybrid_", 1)
        elif "rope" in architecture or "rotary" in architecture:
            method = method.replace("attention_", "attention_rope_", 1)
    if aliases:
        method = aliases.get(Path(directory).name if directory else "", aliases.get(method, aliases.get(manifest["method"], method)))
    task = manifest.get("task", "design")
    side = manifest.get("fixed_side") or a.get("fixed_side") or a.get("repair")
    context = f"amp{a.get('amplicon_min', 80)}-{a.get('amplicon_max', 200)}"
    if task == "redesign":
        context = "redesign_" + side
    elif task == "completion" or side:
        context += "_fix_" + side
    variant = manifest.get("variant", "full")
    steps = a.get("steps")
    if variant == "full" and method.endswith("diffusion") and steps is not None and int(steps) != 16:
        variant = "steps" + str(steps)
    return context, method, variant


def planned_names(protocol, split="test"):
    planned = protocol.get("expected_runs", protocol.get("evaluation_runs", protocol.get("run_names", [])))
    if isinstance(planned, dict):
        planned = list(planned)
    names = []
    for entry in planned:
        if isinstance(entry, dict):
            entry = entry.get("directory", entry.get("name", entry.get("run")))
        if not entry:
            raise ValueError("Protocol run entries need a directory/name")
        names.append(Path(entry).name)
    if len(names) != len(set(names)):
        raise ValueError("Duplicate planned evaluation names")
    if split == "test" and "qprimer" in str(protocol.get("additional_baseline", "")).lower():
        seed = protocol.get("seed", 20261013)
        names.extend(f"qprimer_design_{lo}-{hi}_test_seed{seed}" for lo, hi in
                     ((80, 200), (80, 120), (121, 160), (161, 200)))
        names.extend(f"qprimer_{task}_{side}_test_seed{seed}" for task in
                     ("completion", "redesign") for side in ("forward", "reverse"))
    return set(names)


def load_source(directory, expected, seed, pinned=None, split="test"):
    missing = [name for name in REQUIRED if not (directory / name).is_file()]
    if missing:
        raise PendingRun("missing " + ", ".join(missing))
    manifest = json.loads((directory / "manifest.json").read_text())
    args = manifest["arguments"]
    if args.get("split", manifest.get("split", "test")) != split:
        raise PendingRun("not the expected frozen split")
    if args.get("seed") != seed:
        raise ValueError(f"Different generation seed in {directory}")
    if args.get("budget", 20) != 20:
        raise ValueError(f"Different proposal budget in {directory}")
    completion = json.loads((directory / "completion.json").read_text())
    status = json.loads((directory / "specificity_bowtie_status.json").read_text())
    if not completion.get("all_expected_targets_evaluated") or status.get("status") != "completed":
        raise PendingRun("generation or full-reference screen is not complete")
    hashes = {name + "_sha256": sha256(directory / name) for name in REQUIRED}
    if pinned:
        for key, value in pinned.items():
            if key.endswith("_sha256") and key in hashes and value != hashes[key]:
                raise ValueError(f"Frozen study source changed: {directory}/{key}")
    if completion.get("candidates_sha256") != hashes["candidates.jsonl_sha256"]:
        raise PendingRun("candidate file changed after recorded completion")
    if completion.get("per_gene_sha256") != hashes["per_gene.jsonl_sha256"]:
        raise PendingRun("gene ledger changed after recorded completion")
    if status.get("candidate_sha256") != hashes["candidates.jsonl_sha256"]:
        raise PendingRun("full-reference screen does not match completed candidates")
    if status.get("output_sha256") != hashes["specificity_bowtie.jsonl_sha256"]:
        raise PendingRun("full-reference output does not match completion hash")
    candidates = read_jsonl(directory / "candidates.jsonl")
    ledger = read_jsonl(directory / "per_gene.jsonl")
    assessed = read_jsonl(directory / "specificity_bowtie.jsonl")
    genes = [r["gene_id"] for r in ledger]
    if len(genes) != len(expected) or set(genes) != set(expected):
        raise ValueError(f"Duplicate, missing or foreign test genes in {directory}")
    declared = manifest.get("expected_targets")
    if declared is not None and (len(declared) != len(expected) or set(declared) != set(expected)):
        raise ValueError(f"Manifest cohort mismatch in {directory}")
    ids = [c["candidate_id"] for c in candidates]
    screen_ids = [s["candidate_id"] for s in assessed]
    if len(ids) != len(set(ids)) or len(screen_ids) != len(set(screen_ids)):
        raise ValueError(f"Duplicate candidate/screen IDs in {directory}")
    if set(ids) != set(screen_ids):
        raise PendingRun("specificity candidate accounting is incomplete")
    screens = {s["candidate_id"]: s for s in assessed}
    for c in candidates:
        if c["gene_id"] not in expected:
            raise ValueError(f"Candidate on a foreign test gene in {directory}")
        if c["assessment"].get("computational_feasible"):
            s = screens[c["candidate_id"]]
            if "gene_specificity_screen_pass" not in s or "genomic_DNA_product_detected" not in s:
                raise PendingRun("eligible candidate has an unfinished RNA/genome screen")
    return manifest, candidates, {r["gene_id"]: r for r in ledger}, screens, hashes


def gene_records(directory, loaded, targets, seed, old=False, aliases=None):
    manifest, candidates, ledger, screens, hashes = loaded
    context, method, variant = run_identity(manifest, old, aliases, directory)
    by_gene = defaultdict(list)
    for c in candidates:
        by_gene[c["gene_id"]].append(c)
    records = []
    for gene in sorted(targets):
        proposals = sorted(by_gene[gene], key=lambda c: c.get("proposal_index", 0))
        if len(proposals) > 20:
            raise ValueError(f"Candidate budget exceeded in {directory}:{gene}")
        for budget in BUDGETS:
            own = proposals[:budget]
            base = summarize_candidates(own, budget)
            eligible = [c for c in own if c["assessment"]["computational_feasible"]]
            specific = [c for c in eligible if screens[c["candidate_id"]]["gene_specificity_screen_pass"]]
            retained = [c for c in specific if not screens[c["candidate_id"]]["genomic_DNA_product_detected"]]
            distinct = lambda cs: len({(c["forward_primer"], c["reverse_primer"]) for c in cs})
            count, kept = distinct(specific), distinct(retained)
            temps = [c["assessment"]["shared_temperature_penalty"] for c in specific]
            temps = [t for t in temps if t is not None]
            valid = [c for c in own if c["structural_valid"]]
            conditions = ledger[gene].get("conditions", manifest["conditions"])
            length_pass = [c for c in eligible if
                           c.get("requested_amplicon_min", conditions["amplicon_min"]) <= c["amplicon_length"] <=
                           c.get("requested_amplicon_max", conditions["amplicon_max"])]
            records.append({
                "gene_id": gene, "gene_symbol": targets[gene]["gene_symbol"],
                "method": method, "context": context, "variant": variant,
                "budget": budget, "seed": seed, "total_genes": len(targets),
                "unique_feasible_pairs": base["unique_feasible_pairs"],
                "unique_RNA_pass_pairs": count,
                "unique_RNA_pass_pairs_without_genomic_flag": kept,
                "RNA_no_genome_fraction": kept / count if count else None,
                "one_RNA_pass_pair": float(count >= 1), "two_RNA_pass_pairs": float(count >= 2),
                "one_no_genome_pair": float(kept >= 1), "two_no_genome_pairs": float(kept >= 2),
                "RNA_pass_amplicon_distance": overlap_distance(specific),
                "best_RNA_pass_temperature_penalty": min(temps) if temps else None,
                "budget_feasibility_rate": base["budget_feasibility_rate"],
                "geometry_valid_rate": len(valid) / budget,
                "length_compliant_feasible_rate": len(length_pass) / budget,
                "anchor_preservation_rate": sum(c.get("anchor_preserved", True) is not False for c in valid) / len(valid) if valid else None,
                "distinct_feasible_amplicons": base["distinct_feasible_amplicons"],
                "mean_feasible_amplicon_jaccard_distance": base["mean_feasible_amplicon_jaccard_distance"],
                "best_shared_temperature_penalty": base["best_shared_temperature_penalty"],
                "run": str(directory), "source_family": "frozen_cnn_study" if old else "attention_comparison",
                "original_amplicon_length": ledger[gene].get("original_amplicon_length"),
                "conditions": conditions,
            })
    return records, {"directory": str(directory), "name": directory.name,
                     "method": method, "context": context, "variant": variant,
                     "architecture": manifest.get("architecture", "dilated_cnn" if old and method.startswith("cnn_") else None),
                     "checkpoint_sha256":manifest.get("checkpoint_sha256"),
                     "source_family": "frozen_cnn_study" if old else "attention_comparison", **hashes}


def summarize(rows, seed, repeats):
    groups = defaultdict(list)
    for r in rows:
        groups[r["context"], r["method"], r["variant"], r["budget"]].append(r)
    summaries, contrasts = [], []
    for (context, method, variant, budget), own in sorted(groups.items()):
        if len(own) != len({r["gene_id"] for r in own}):
            raise ValueError(f"Duplicate runs for {context}/{method}/{variant}/{budget}")
        item = {"context": context, "method": method, "variant": variant, "budget": budget,
                "seed": seed, "total_genes": len(own), "run": own[0]["run"]}
        for key in ENDPOINTS:
            values = [r[key] for r in own if r[key] is not None]
            ci = paired_bootstrap_difference(values, np.zeros(len(values)), seed=seed, repeats=repeats) if values else None
            item[key] = ci["mean_difference"] if ci else None
            item[key + "_n"] = len(values)
            item[key + "_ci95"] = ci["ci95"] if ci else None
        summaries.append(item)
    for (context, method, variant, budget), own in sorted(groups.items()):
        # All defined comparisons are retained, including architecture and task effects.
        if variant != "full" or not (method.endswith(("_diffusion", "_direct")) or method in ("qprimer", "primercast")):
            continue
        left = {r["gene_id"]: r for r in own}
        for (other_context, other_method, other_variant, other_budget), right_rows in sorted(groups.items()):
            if context != other_context or budget != other_budget or (method, variant) == (other_method, other_variant):
                continue
            if method.endswith("_direct") and not (other_method.endswith("_direct") and other_variant == "full"):
                continue
            if method in ("qprimer", "primercast") and (other_method,other_variant) != ("primer3","full"):
                continue
            if method == "cnn_diffusion" and other_method.startswith("attention_"):
                continue
            right = {r["gene_id"]: r for r in right_rows}
            for key in ENDPOINTS:
                genes = sorted(g for g in left.keys() & right.keys() if left[g][key] is not None and right[g][key] is not None)
                if not genes:
                    continue
                estimate = paired_bootstrap_difference([left[g][key] for g in genes], [right[g][key] for g in genes], seed=seed, repeats=repeats)
                contrasts.append({"context": context, "left": method + "_" + variant,
                                  "right": other_method + "_" + other_variant,
                                  "left_method": method, "left_variant": variant,
                                  "right_method": other_method, "right_variant": other_variant,
                                  "budget": budget, "endpoint": key, "total_genes": len(left),
                                  "joint_gene_ids": genes, **estimate})
    return summaries, contrasts


def csv_file(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in records for k in r)))
        writer.writeheader()
        writer.writerows([{k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in r.items()} for r in records])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-partial", action="store_true", help="Include only fully completed/screened runs while other runs are pending")
    parser.add_argument("--protocol", default="experiments/attention/protocol.json")
    parser.add_argument("--old-summary", default="results/final_fixed_seed/study_summary.json")
    parser.add_argument("--evaluations", nargs="+", help="One or more completed evaluation directories; defaults to protocol.evaluation_directories")
    parser.add_argument("--output", default="experiments/attention/results/study_summary.json")
    parser.add_argument("--bootstrap-repeats", type=int, default=10000)
    parser.add_argument("--targets", default="data/processed/targets_pilot.jsonl", help="Frozen target cohort JSONL")
    parser.add_argument("--split", choices=["test", "external_test"], help="Held-out evaluation split; inferred when the target file contains only external_test")
    parser.add_argument("--no-old-sources", action="store_true", help="Analyze only supplied evaluation directories, without adding the original frozen CNN cohort")
    args = parser.parse_args()
    old_path, protocol_path = resolve(args.old_summary), resolve(args.protocol)
    old = json.loads(old_path.read_text()) if not args.no_old_sources else {"inputs":[],"seed":20261013}
    protocol = json.loads(protocol_path.read_text()) if protocol_path.is_file() else {}
    target_path = resolve(args.targets)
    all_targets = read_jsonl(target_path)
    evaluation_split = args.split or ("external_test" if {r["split"] for r in all_targets} == {"external_test"} else "test")
    if evaluation_split == "external_test" and not args.no_old_sources:
        parser.error("External gene analysis must use --no-old-sources; original test genes cannot be pooled implicitly")
    planned = planned_names(protocol,evaluation_split)
    if not args.allow_partial and not planned:
        parser.error("Default analysis requires a protocol with expected_runs; use --allow-partial for ongoing work")
    seed = protocol.get("seed", old["seed"])
    if not args.no_old_sources and seed != old["seed"]:
        raise ValueError("Old/new studies must use the same fixed seed")
    cohorts = {split: {r["gene_id"]: r for r in all_targets if r["split"] == split} for split in (evaluation_split, "validation")}
    targets = cohorts[evaluation_split]
    if not targets:
        raise ValueError("No targets in the requested evaluation split")
    selected_rows = [r for r in all_targets if r["split"] == evaluation_split]
    if len(targets) != len(selected_rows):
        raise ValueError("Frozen target file contains duplicate gene IDs")
    rows, inputs, skipped, validation_inputs = [], [], [], []
    for pin in old["inputs"]:
        directory = resolve(pin["directory"])
        try:
            loaded = load_source(directory, targets, seed, pin)
        except PendingRun as error:
            raise ValueError(f"Frozen completed source unavailable: {directory}: {error}") from error
        records, source = gene_records(directory, loaded, targets, seed, old=True)
        rows.extend(records)
        inputs.append(source)
    evaluation_roots = args.evaluations or protocol.get("evaluation_directories", ["experiments/attention/results/evaluations"])
    if isinstance(evaluation_roots, str):
        evaluation_roots = [evaluation_roots]
    directories = sorted({p.resolve() for root in evaluation_roots for p in resolve(root).glob("*") if p.is_dir()})
    completed_new = set()
    observed = {p.name for p in directories}
    for directory in directories:
        try:
            manifest_path = directory / "manifest.json"
            manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
            split = manifest.get("arguments", {}).get("split", evaluation_split)
            if split not in cohorts:
                raise PendingRun("not a validation or held-out test source")
            loaded = load_source(directory, cohorts[split], seed, split=split)
        except PendingRun as error:
            skipped.append({"directory": str(directory), "name": directory.name, "reason": str(error)})
            continue
        if split == "validation":
            validation_inputs.append({"directory": str(directory), "name": directory.name,
                                      "split": split, "total_genes": len(cohorts[split]), **loaded[-1]})
            completed_new.add(directory.name)
            continue
        records, source = gene_records(directory, loaded, targets, seed, aliases=protocol.get("source_method_aliases"))
        source["evaluation_split"] = split
        if split == "external_test":
            source["source_family"] = "external_heldout_genes"
            for record in records:
                record["source_family"] = "external_heldout_genes"
                record["evaluation_split"] = split
        rows.extend(records)
        source["planned"] = directory.name in planned
        inputs.append(source)
        completed_new.add(directory.name)
    for name in sorted(planned - observed):
        skipped.append({"name": name, "reason": "planned run has not been created"})
    pending = sorted(planned - completed_new)
    if not args.allow_partial and pending:
        raise ValueError(f"Planned runs are unfinished: {', '.join(pending)}")
    summaries, contrasts = summarize(rows, seed, args.bootstrap_repeats)
    output = resolve(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    complete = bool(planned) and not pending
    doc = {"at": utc_now(), "seed": seed, "source_runs": len(inputs), "inputs": inputs,
           "summaries": summaries, "paired_comparisons": contrasts,
           "analysis_unit": "target gene; one fixed training and generation seed",
           "budgets": list(BUDGETS), "endpoints": list(ENDPOINTS), "total_test_genes": len(targets),
           "status": "complete" if complete else "partial_completed_sources_only",
           "additional_experiments_complete": complete,
           "allow_partial": args.allow_partial,
           "planned_new_runs": len(planned), "completed_new_runs": len(completed_new),
           "pending_planned_runs": pending, "skipped_sources": skipped,
           "validation_inputs": validation_inputs,
           "validation_policy": "Validation sources establish protocol completion only; no validation genes enter test summaries or comparisons",
           "evaluation_directories": [str(resolve(p)) for p in evaluation_roots],
           "primary_methods": protocol.get("primary_methods", {}),
           "architecture_labels": "Pure rotary attention and CNN-attention hybrid are separate methods, selected from each recorded manifest",
           "old_summary_sha256": sha256(old_path) if not args.no_old_sources else None,
           "old_sources_included": not args.no_old_sources,
           "evaluation_split":evaluation_split,"targets_path":str(target_path),"targets_sha256":sha256(target_path),
           "cohort_policy":"Separate frozen external genes; no original study genes pooled" if evaluation_split == "external_test" else "Original frozen held-out test genes",
           "protocol_sha256": sha256(protocol_path) if protocol_path.is_file() else None,
           "quality_fraction": "unique RNA-pass pairs without genomic flag / unique RNA-pass pairs; defined only when RNA-pass count > 0",
           "prefix_accounting": "Recorded proposal prefixes; invalid proposals, duplicates and shortages remain in K",
           "conditional_endpoints": "Diversity and temperature endpoints use only defined genes; paired contrasts use their joint intersection and report n",
           "intervals": "Gene-paired percentile bootstrap, 95%; all defined comparisons retained, unadjusted for multiplicity",
           "method_labels": {"cnn_diffusion": "Existing CNN diffusion", "cnn_direct": "Existing CNN direct predictor",
                             "attention_diffusion": "Attention diffusion", "attention_direct": "Attention direct predictor", "qprimer": "qPrimer upstream design module with shared pair-task conditions",
                             "attention_rope_diffusion": "Pure rotary-attention diffusion", "attention_rope_direct": "Pure rotary-attention direct predictor",
                             "attention_hybrid_diffusion": "CNN-attention hybrid diffusion", "attention_hybrid_direct": "CNN-attention hybrid direct predictor",
                             "primercast": "PrimerCast upstream trained design/ranking modules with shared pair-task conditions"}}
    write_jsonl(output.with_suffix(".per_gene.jsonl"), rows)
    write_json(output, doc)
    csv_file(output.with_suffix(".groups.csv"), summaries)
    csv_file(output.with_suffix(".contrasts.csv"), contrasts)
    print(json.dumps({"status": doc["status"], "source_runs": len(inputs), "new_runs": len(completed_new),
                      "pending_planned_runs": pending, "groups": len(summaries), "contrasts": len(contrasts), "output": str(output)}))


if __name__ == "__main__":
    main()
