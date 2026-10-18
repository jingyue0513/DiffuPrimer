#!/usr/bin/env python3
"""Describe exon geometry of the frozen, RNA-passing primary designs."""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import paired_bootstrap_difference


root = Path(__file__).resolve().parents[1]
source = root / "results/temperature_tuning/selected_test_summary.json"
summary = json.loads(source.read_text())
inputs = summary["inputs"]
directories = {r["directory"] for r in inputs}
annotation_path = root / "results/selected_exon_context_audit.per_candidate.jsonl"
annotations = {}
for row in read_jsonl(annotation_path):
    if row["run_directory"] in directories:
        annotations[(row["run_directory"], row["candidate_id"])] = row

categories = ["primer_spans_junction", "amplicon_crosses_junction_only", "same_exon", "terminal_polyA", "unmapped"]
records, sources = [], []
for item in inputs:
    directory = Path(item["directory"])
    for filename, key in [("completion.json", "completion_sha256"), ("specificity_bowtie.jsonl", "specificity_sha256"),
                          ("specificity_bowtie_status.json", "specificity_status_sha256")]:
        if sha256(directory / filename) != item[key]:
            raise ValueError("Primary analysis source changed")
    complete = json.loads((directory / "completion.json").read_text())
    if sha256(directory / "candidates.jsonl") != complete["candidates_sha256"]:
        raise ValueError("Candidate source changed")
    manifest = json.loads((directory / "manifest.json").read_text())
    screens = {r["candidate_id"]: r for r in read_jsonl(directory / "specificity_bowtie.jsonl")}
    by_gene = defaultdict(list)
    for row in read_jsonl(directory / "candidates.jsonl"):
        by_gene[row["gene_id"]].append(row)
    for gene, rows in by_gene.items():
        counts = {c: {"RNA_pass": 0, "RNA_pass_without_genomic_flag": 0} for c in categories}
        seen = set()
        for row in rows:
            if not row["assessment"]["computational_feasible"]:
                continue
            screen = screens[row["candidate_id"]]
            if not screen["gene_specificity_screen_pass"]:
                continue
            pair = (row["forward_primer"], row["reverse_primer"])
            if pair in seen:
                continue
            seen.add(pair)
            context = annotations[(str(directory), row["candidate_id"])]
            if "forward_spans_junction" not in context:
                category = "unmapped"
            elif context["forward_overlaps_terminal_polyA"] or context["reverse_overlaps_terminal_polyA"]:
                category = "terminal_polyA"
            elif context["forward_spans_junction"] or context["reverse_spans_junction"]:
                category = "primer_spans_junction"
            elif context["amplicon_crosses_annotated_junction"]:
                category = "amplicon_crosses_junction_only"
            else:
                category = "same_exon"
            counts[category]["RNA_pass"] += 1
            counts[category]["RNA_pass_without_genomic_flag"] += not screen["genomic_DNA_product_detected"]
        records.append({"method": manifest["method"], "gene_id": gene, "seed": manifest["arguments"]["seed"], "categories": counts})
    sources.append({**item, "candidates_sha256": complete["candidates_sha256"]})

averaged = []
for method in sorted({r["method"] for r in records}):
    for gene in sorted({r["gene_id"] for r in records}):
        own = [r for r in records if r["method"] == method and r["gene_id"] == gene]
        expected = [20261013] if method == "primer3" else [20261013, 20261014, 20261015]
        if sorted(r["seed"] for r in own) != expected:
            raise ValueError("Exon association gene/seed denominator differs")
        averaged.append({"method": method, "gene_id": gene, "seeds": expected,
            "categories": {c: {k: float(np.mean([r["categories"][c][k] for r in own]))
                for k in ["RNA_pass", "RNA_pass_without_genomic_flag"]} for c in categories}})

groups, comparisons = [], []
for method in sorted({r["method"] for r in averaged}):
    own = [r for r in averaged if r["method"] == method]
    if len(own) != 75:
        raise ValueError("Expected all 75 held-out genes")
    for category in categories:
        total = sum(r["categories"][category]["RNA_pass"] for r in own)
        no_flag = sum(r["categories"][category]["RNA_pass_without_genomic_flag"] for r in own)
        groups.append({"method": method, "category": category, "genes": 75,
            "mean_unique_RNA_pass_pairs": total / 75,
            "mean_unique_RNA_pass_pairs_without_genomic_flag": no_flag / 75,
            "pooled_fraction_without_genomic_flag": no_flag / total if total else None})
    if method != "diffusion":
        left = {r["gene_id"]: r for r in averaged if r["method"] == "diffusion"}
        right = {r["gene_id"]: r for r in own}
        genes = sorted(left)
        for category in categories:
            endpoint = "RNA_pass_without_genomic_flag"
            comparisons.append({"comparison": "diffusion minus " + method, "category": category,
                **paired_bootstrap_difference([left[g]["categories"][category][endpoint] for g in genes],
                    [right[g]["categories"][category][endpoint] for g in genes], seed=20261013)})
output = root / "results/exon_specificity_association.json"
write_jsonl(output.with_suffix(".per_gene.jsonl"), averaged)
write_json(output, {"at": utc_now(), "inputs": sources, "exon_records_sha256": sha256(annotation_path),
    "groups": groups, "paired_comparisons": comparisons,
    "interpretation": "Descriptive partition of unique RNA-passing pairs by annotated geometry; gene is the bootstrap unit and model seeds are averaged within gene. Association does not establish causation. Genomic flag refers to potential products under the fixed alignment rule."})
print(json.dumps({"groups": len(groups), "paired_comparisons": len(comparisons)}))
