#!/usr/bin/env python3
"""Compare recorded limited BLAST and independent full-length mismatch screens."""
import json
from collections import Counter, defaultdict
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl


root = Path(__file__).resolve().parents[1]
aggregate = json.loads((root / "results/primary_aggregate.json").read_text())
rows, runs, inputs = [], [], []
for run in aggregate["runs"]:
    directory = root / run["source"]
    source = directory / "candidates.jsonl"
    expected_hash = sha256(source)
    for name in ["specificity_status.json", "specificity_bowtie_status.json"]:
        if json.loads((directory / name).read_text())["candidate_sha256"] != expected_hash:
            raise ValueError("Stale screen")
    old = {r["candidate_id"]: r for r in read_jsonl(directory / "specificity.jsonl")}
    new = {r["candidate_id"]: r for r in read_jsonl(directory / "specificity_bowtie.jsonl")}
    candidates = read_jsonl(source)
    if set(old) != set(new) or set(old) != {r["candidate_id"] for r in candidates}:
        raise ValueError("Candidate ledgers differ")
    counts = Counter()
    strata = defaultdict(Counter)
    for row in candidates:
        if not row["assessment"]["computational_feasible"]:
            continue
        a, b = old[row["candidate_id"]], new[row["candidate_id"]]
        changed_RNA = a["gene_specificity_screen_pass"] != b["gene_specificity_screen_pass"]
        changed_genome = a["genomic_DNA_product_detected"] != b["genomic_DNA_product_detected"]
        key = min(len(row["forward_primer"]), len(row["reverse_primer"]))
        flags = {"feasible_proposals": True, "BLAST_RNA_pass": a["gene_specificity_screen_pass"],
            "Bowtie_RNA_pass": b["gene_specificity_screen_pass"], "RNA_classification_changed": changed_RNA,
            "new_off_target_failure": a["gene_specificity_screen_pass"] and bool(b["off_target_products"]),
            "genomic_classification_changed": changed_genome,
            "new_genomic_product_flag": not a["genomic_DNA_product_detected"] and b["genomic_DNA_product_detected"]}
        for flag, value in flags.items():
            counts[flag] += int(value)
            strata[key][flag] += int(value)
        if changed_RNA or changed_genome:
            rows.append({"run": directory.name, "method": run["method"], "gene_id": row["gene_id"],
                "candidate_id": row["candidate_id"], "forward_primer": row["forward_primer"], "reverse_primer": row["reverse_primer"],
                "BLAST_RNA_pass": a["gene_specificity_screen_pass"], "Bowtie_RNA_pass": b["gene_specificity_screen_pass"],
                "BLAST_genomic_flag": a["genomic_DNA_product_detected"], "Bowtie_genomic_flag": b["genomic_DNA_product_detected"],
                "Bowtie_off_target_products": b["off_target_products"], "Bowtie_unresolved_products": b["unresolved_products"]})
    runs.append({"run": directory.name, "method": run["method"], **dict(counts),
        "by_minimum_primer_length": {str(k): dict(v) for k,v in sorted(strata.items())}})
    inputs.append({"directory": str(directory), "candidate_sha256": expected_hash,
        "BLAST_sha256": sha256(directory / "specificity.jsonl"), "Bowtie_sha256": sha256(directory / "specificity_bowtie.jsonl")})
write_jsonl(root / "results/bowtie_audit/screen_comparison.changed_candidates.jsonl", rows)
write_json(root / "results/bowtie_audit/screen_comparison.json", {"analyzed_at": utc_now(), "inputs": inputs,
    "runs": runs, "changed_candidates": len(rows),
    "comparison": "Identical feasible proposals and pairing kernel. BLAST local full-HSP acceptance versus all-strata full-length ungapped <=2 substitutions; both require perfect original terminal 5 bases."})
print(json.dumps({"runs": len(runs), "changed_candidates": len(rows)}), flush=True)
