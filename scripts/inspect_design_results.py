#!/usr/bin/env python3
"""Inspect completed experiments: budget curves, failure modes and gene strata."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+")
    parser.add_argument("--output", default="results/design_diagnostics.json")
    parser.add_argument("--screen", choices=["blast", "bowtie"], default="blast")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target_rows = read_jsonl(root / "data/processed/targets_pilot.jsonl")
    targets = {row["gene_id"]: row for row in target_rows}
    sources = {row["gene_id"] for row in read_jsonl(root / "data/processed/primerbank_pilot.jsonl")}
    records, run_summaries = [], []
    for raw in args.directories:
        directory = root / raw
        manifest = json.loads((directory / "manifest.json").read_text())
        complete = json.loads((directory / "completion.json").read_text())
        if not complete["all_expected_targets_evaluated"] or complete["candidates_sha256"] != sha256(directory / "candidates.jsonl") or complete["per_gene_sha256"] != sha256(directory / "per_gene.jsonl"):
            raise ValueError("Require complete, unchanged experiment outputs")
        candidate_groups = defaultdict(list)
        for row in read_jsonl(directory / "candidates.jsonl"):
            candidate_groups[row["gene_id"]].append(row)
        suffix = "_bowtie" if args.screen == "bowtie" else ""
        screen_path = directory / f"specificity{suffix}.jsonl"
        screens = {row["candidate_id"]: row for row in read_jsonl(screen_path)} if screen_path.exists() else {}
        if screens:
            status = json.loads((directory / f"specificity{suffix}_status.json").read_text())
            if status["candidate_sha256"] != complete["candidates_sha256"]:
                raise ValueError("Stale specificity screen")
            if "output_sha256" in status and sha256(screen_path) != status["output_sha256"]:
                raise ValueError("Specificity output changed")
        failures = Counter()
        forward_lengths, reverse_lengths = Counter(), Counter()
        forward_tm, reverse_tm, tm_differences = [], [], []
        per_gene = read_jsonl(directory / "per_gene.jsonl")
        expected = {row["gene_id"] for row in target_rows if row["split"] == manifest["arguments"]["split"]}
        if {row["gene_id"] for row in per_gene} != expected:
            raise ValueError("Missing frozen target genes")
        budget = manifest["arguments"]["budget"]
        for row in per_gene:
            gene = row["gene_id"]
            own = sorted(candidate_groups[gene], key=lambda item: item["proposal_index"])
            failures.update(failure for item in own for failure in item["assessment"]["failures"])
            for item in own:
                if "forward_primer" in item:
                    forward_lengths[len(item["forward_primer"])] += 1
                    reverse_lengths[len(item["reverse_primer"])] += 1
                assessment = item["assessment"]
                if assessment.get("tm_difference") is not None:
                    forward_tm.append(assessment["forward_properties"]["tm"])
                    reverse_tm.append(assessment["reverse_properties"]["tm"])
                    tm_differences.append(assessment["tm_difference"])
            failures["missing_proposal"] += budget - len(own)
            template = targets[gene]["template_sequence"]
            feasible = [item for item in own if item["assessment"]["computational_feasible"]]
            curve = []
            for k in [1, 2, 5, 10, 20]:
                if k > budget:
                    continue
                prefix = own[:k]
                passed = [item for item in prefix if item["assessment"]["computational_feasible"]]
                unique = {(item["forward_primer"], item["reverse_primer"]) for item in passed}
                specific = {(item["forward_primer"], item["reverse_primer"]) for item in passed
                    if screens.get(item["candidate_id"], {}).get("gene_specificity_screen_pass")}
                curve.append({"prefix_budget": k, "unique_feasible": len(unique), "at_least_two_feasible": len(unique) >= 2,
                    "unique_RNA_pass": len(specific) if screens else None,
                    "at_least_two_RNA_pass": len(specific) >= 2 if screens else None})
            records.append({**row, "run": directory.name, "seed": manifest["arguments"]["seed"],
                "split": manifest["arguments"]["split"], "template_length": len(template),
                "template_gc_percent": 100 * (template.count("G") + template.count("C")) / len(template),
                "length_stratum": "<1000" if len(template) < 1000 else ("1000-2999" if len(template) < 3000 else ">=3000"),
                "source_record_available": gene in sources, "prefix_budget_curve": curve,
                "feasible_duplicate_proposals": len(feasible) - len({(item["forward_primer"], item["reverse_primer"]) for item in feasible}),
                "failed_candidate_ids": [item["candidate_id"] for item in own if not item["assessment"]["computational_feasible"]]})
        own_records = [row for row in records if row["run"] == directory.name]
        strata = []
        for stratum in ["<1000", "1000-2999", ">=3000"]:
            own = [row for row in own_records if row["length_stratum"] == stratum]
            if own:
                strata.append({"template_length": stratum, "genes": len(own),
                    "mean_unique_feasible": float(np.mean([row["unique_feasible_pairs"] for row in own])),
                    "mean_amplicon_distance": float(np.mean([row["mean_feasible_amplicon_jaccard_distance"] for row in own
                        if row["mean_feasible_amplicon_jaccard_distance"] is not None])) if any(row["mean_feasible_amplicon_jaccard_distance"] is not None for row in own) else None,
                    "median_runtime_seconds": float(np.median([row["elapsed_seconds"] for row in own]))})
        curves = []
        for k in [1, 2, 5, 10, 20]:
            points = [point for row in own_records for point in row["prefix_budget_curve"] if point["prefix_budget"] == k]
            if not points:
                continue
            summary = {"prefix_budget": k, "genes": len(points)}
            for key in ["unique_feasible", "at_least_two_feasible", "unique_RNA_pass", "at_least_two_RNA_pass"]:
                values = [point[key] for point in points if point[key] is not None]
                summary[key] = float(np.mean(values)) if values else None
            curves.append(summary)
        run_summaries.append({"run": directory.name, "method": manifest["method"], "seed": manifest["arguments"]["seed"],
            "split": manifest["arguments"]["split"], "genes": len(own_records),
            "failure_counts": dict(failures.most_common()), "failure_rate_denominator": len(own_records) * budget,
            "forward_length_histogram": dict(sorted(forward_lengths.items())), "reverse_length_histogram": dict(sorted(reverse_lengths.items())),
            "mean_forward_tm": float(np.mean(forward_tm)) if forward_tm else None,
            "mean_reverse_tm": float(np.mean(reverse_tm)) if reverse_tm else None,
            "pair_tm_correlation": float(np.corrcoef(forward_tm, reverse_tm)[0, 1]) if len(forward_tm) > 1 and np.std(forward_tm) > 0 and np.std(reverse_tm) > 0 else None,
            "tm_difference_quantile_levels": [0, .25, .5, .75, .95, 1],
            "tm_difference_quantiles": np.quantile(tm_differences, [0, .25, .5, .75, .95, 1]).tolist() if tm_differences else None,
            "length_strata": strata, "prefix_budget_curves": curves,
            "no_feasible_design_genes": [row["gene_id"] for row in own_records if not row["target_covered"]],
            "completion_sha256": sha256(directory / "completion.json")})
    output = root / args.output
    write_jsonl(output.with_suffix(".per_gene.jsonl"), records)
    write_json(output, {"generated_at": utc_now(), "screen": args.screen, "runs": run_summaries,
        "prefix_policy": "Prefixes of the original ordered K proposals; not independent regenerated K-sized batches. Primer3 prefixes retain its returned ranking.",
        "failure_policy": "One candidate can violate several constraints; category counts are not mutually exclusive."})
    print(json.dumps({"runs": len(run_summaries), "output": str(output)}))


if __name__ == "__main__":
    main()
