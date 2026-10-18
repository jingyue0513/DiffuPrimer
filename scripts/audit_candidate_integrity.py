#!/usr/bin/env python3
"""Verify recorded sequences, cohort accounting and fixed-primer preservation."""
import argparse
import json
from collections import Counter
from pathlib import Path

from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.sequence import reverse_complement


parser = argparse.ArgumentParser()
parser.add_argument("directories", nargs="*")
parser.add_argument("--output", default="results/checks/candidate_integrity/check.json")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
targets_path = root / "data/processed/targets_pilot.jsonl"
targets = {r["gene_id"]: r for r in read_jsonl(targets_path)}
if args.directories:
    directories = [root / path for path in args.directories]
else:
    bases = [root / "results/evaluations"] + sorted((root / "experiments").glob("*/results/evaluations"))
    directories = sorted({p.parent for base in bases for p in base.glob("*/completion.json")})
anchor_cache, records = {}, []
for directory in directories:
    manifest_path, completion_path = directory / "manifest.json", directory / "completion.json"
    manifest, completion = [json.loads(p.read_text()) for p in [manifest_path, completion_path]]
    if manifest["stage"] == "calibration":
        continue
    candidate_path = directory / "candidates.jsonl"
    if sha256(candidate_path) != completion["candidates_sha256"] or manifest["target_sha256"] != sha256(targets_path):
        raise ValueError(f"Changed experiment inputs: {directory}")
    arguments = manifest["arguments"]
    expected = {gene for gene, target in targets.items() if target["split"] == arguments["split"]}
    ledger = read_jsonl(directory / "per_gene.jsonl")
    if set(manifest["expected_targets"]) != expected or len(ledger) != len(expected) or {r["gene_id"] for r in ledger} != expected:
        raise ValueError(f"Incomplete cohort: {directory}")
    candidates = read_jsonl(candidate_path)
    if len({r["candidate_id"] for r in candidates}) != len(candidates):
        raise ValueError(f"Repeated candidate identifier: {directory}")
    counts = Counter(r["gene_id"] for r in candidates)
    if set(counts) - expected or any(counts[r["gene_id"]] != r["proposals_returned"] or
                                   counts[r["gene_id"]] > arguments["budget"] for r in ledger):
        raise ValueError(f"Candidate budget accounting differs: {directory}")
    if not completion["all_expected_targets_evaluated"]:
        raise ValueError(f"Incomplete experiment: {directory}")
    structural, preserved, alternate_positions = 0, 0, 0
    side = arguments["repair"]
    for row in candidates:
        if not row["structural_valid"]:
            if row["assessment"]["computational_feasible"]:
                raise ValueError("An invalid structure was counted as feasible")
            continue
        structural += 1
        sequence = targets[row["gene_id"]]["template_sequence"]
        f0, f1, r0, r1 = [row[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]]
        if not (0 <= f0 < f1 <= r0 < r1 <= len(sequence)):
            raise ValueError(f"Invalid recorded coordinates: {row['candidate_id']}")
        if row["forward_primer"] != sequence[f0:f1] or row["reverse_primer"] != reverse_complement(sequence[r0:r1]) or row["amplicon_length"] != r1 - f0:
            raise ValueError(f"Sequence/coordinate disagreement: {row['candidate_id']}")
        if side:
            key = (row["gene_id"], json.dumps(manifest["conditions"], sort_keys=True))
            if key not in anchor_cache:
                anchor, _ = primer3_candidates(sequence, Conditions(**manifest["conditions"]), budget=1)
                anchor_cache[key] = anchor[0] if anchor else None
            anchor = anchor_cache[key]
            if anchor is None or row[f"{side}_primer"] != anchor[f"{side}_primer"]:
                raise ValueError(f"Fixed primer changed: {row['candidate_id']}")
            same_positions = all(row[f"{side}_{boundary}"] == anchor[f"{side}_{boundary}"] for boundary in ["start", "end"])
            if not same_positions and manifest["method"] != "primer3":
                raise ValueError(f"Fixed binding positions changed: {row['candidate_id']}")
            alternate_positions += int(not same_positions)
            preserved += 1
    records.append({"directory": str(directory), "method": manifest["method"], "split": arguments["split"],
                    "genes": len(expected), "proposals": len(candidates), "structural_pairs": structural,
                    "fixed_side": side, "fixed_primer_sequence_preserved_pairs": preserved,
                    "primer3_identical_sequence_alternate_positions": alternate_positions,
                    "candidates_sha256": sha256(candidate_path), "completion_sha256": sha256(completion_path)})
write_json(root / args.output, {"at": utc_now(), "runs": records, "completed_runs": len(records),
                                "candidate_sequence_geometry_checks_passed": True,
                                "complete_frozen_cohorts": True, "candidate_budget_accounting_passed": True,
                                "fixed_primer_sequence_preservation_passed": True,
                                "targets_sha256": sha256(targets_path)})
print(json.dumps({"completed_runs": len(records), "proposals": sum(r["proposals"] for r in records),
                  "fixed_primer_pairs": sum(r["fixed_primer_sequence_preserved_pairs"] for r in records)}))
