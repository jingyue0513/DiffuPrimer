#!/usr/bin/env python3
"""Draft a blinded laboratory plan without creating or inferring assay outcomes."""
import argparse
import csv
import json
import random
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl


parser = argparse.ArgumentParser()
parser.add_argument("--primer3-directory")
parser.add_argument("--diffusion-directory")
parser.add_argument("--screen", choices=["bowtie", "blast"], default="bowtie")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
directory = root / "results/lab_plan"
directory.mkdir(parents=True, exist_ok=True)
targets_path = root / "data/processed/targets_pilot.jsonl"
targets = sorted([row for row in read_jsonl(targets_path) if row["split"] == "test"], key=lambda row: int(row["gene_id"]))
rng = random.Random(20261013 + 7000)
selection = rng.sample(targets, min(40, len(targets)))
selected = [{key: row[key] for key in ["gene_id", "gene_symbol", "transcript_accession_version", "cluster_id"]} for row in selection]
selection_path = directory / "target_selection.jsonl"
if selection_path.exists() and read_jsonl(selection_path) != selected:
    raise ValueError("Laboratory gene selection already frozen with different contents")
write_jsonl(selection_path, selected)
criteria_path = directory / "success_criteria_pending.json"
if not criteria_path.exists():
    write_json(criteria_path, {
        "status": "pending_laboratory_definition_not_frozen", "technical_replicates": 3,
        "cq_repeatability_threshold": None, "efficiency_minimum": None, "efficiency_maximum": None,
        "product_identity_rule": None, "no_template_rule": None, "no_reverse_transcription_rule": None,
        "sample_materials": None, "instrument": None,
        "instruction": "Laboratory must finalize criteria and material/instrument metadata before candidate unblinding. Empty values are not successful assays."})
if not args.primer3_directory or not args.diffusion_directory:
    write_json(directory / "status.json", {"prepared_at": utc_now(), "status": "target_selection_only_candidate_screen_pending",
        "selected_genes": len(selected), "selection_sha256": sha256(selection_path), "experiments_conducted": False})
    print(str(selection_path))
    raise SystemExit(0)
key, candidate_ledger, reaction_rows, curve_rows = [], [], [], []
source_records = []
for method, raw in [("primer3", args.primer3_directory), ("diffusion", args.diffusion_directory)]:
    source = root / raw
    completion = json.loads((source / "completion.json").read_text())
    manifest = json.loads((source / "manifest.json").read_text())
    suffix = "_bowtie" if args.screen == "bowtie" else ""
    screen_path = source / f"specificity{suffix}.jsonl"
    screen_status_path = source / f"specificity{suffix}_status.json"
    screen_status = json.loads(screen_status_path.read_text())
    if manifest["arguments"]["split"] != "test" or not completion["all_expected_targets_evaluated"]:
        raise ValueError("Laboratory candidates require complete held-out evaluations")
    if sha256(source / "candidates.jsonl") != completion["candidates_sha256"] or screen_status["candidate_sha256"] != completion["candidates_sha256"]:
        raise ValueError("Laboratory candidate or screening inputs changed")
    if "output_sha256" in screen_status and sha256(screen_path) != screen_status["output_sha256"]:
        raise ValueError("Laboratory specificity output changed")
    source_records.append({"method": method, "directory": str(source), "candidates_sha256": completion["candidates_sha256"],
                           "screen_status_sha256": sha256(screen_status_path), "screen_sha256": sha256(screen_path)})
    candidates = read_jsonl(source / "candidates.jsonl")
    screens = {row["candidate_id"]: row for row in read_jsonl(screen_path)}
    if set(screens) != {r["candidate_id"] for r in candidates}:
        raise ValueError("Incomplete laboratory screening ledger")
    for target in selected:
        eligible = [row for row in candidates if row["gene_id"] == target["gene_id"] and row["assessment"]["computational_feasible"]
                    and screens.get(row["candidate_id"], {}).get("gene_specificity_screen_pass", False)]
        eligible.sort(key=lambda row: (row["assessment"]["shared_temperature_penalty"], row["candidate_id"]))
        unique, seen = [], set()
        for row in eligible:
            if (row["forward_primer"], row["reverse_primer"]) not in seen:
                seen.add((row["forward_primer"], row["reverse_primer"]))
                unique.append(row)
        for rank in [1, 2]:
            if len(unique) < rank:
                candidate_ledger.append({"gene_id": target["gene_id"], "method": method, "rank": rank, "status": "no_design"})
                continue
            row = unique[rank - 1]
            blind_id = f"B{rng.randrange(10**8):08d}"
            if any(item["blind_id"] == blind_id for item in key):
                raise ValueError("Blind identifier collision")
            record = {"blind_id": blind_id, "gene_id": row["gene_id"], "method": method, "rank": rank,
                      "candidate_id": row["candidate_id"], "status": "planned_not_tested"}
            key.append(record)
            candidate_ledger.append(record)
            sequence_data = {"blind_id": blind_id, "gene_id": row["gene_id"], "forward_primer": row["forward_primer"],
                "reverse_primer": row["reverse_primer"], "expected_amplicon_bp": row["amplicon_length"],
                "genomic_DNA_product_detected": screens[row["candidate_id"]].get("genomic_DNA_product_detected"),
                "primer_blast_review_status": "pending", "assay_status": "not_conducted"}
            for replicate in [1, 2, 3]:
                reaction_rows.append({**sequence_data, "reaction_type": "sample", "technical_replicate": replicate})
            for control in ["NTC", "noRT"]:
                reaction_rows.append({**sequence_data, "reaction_type": control, "technical_replicate": 1})
            for level in range(1, 6):
                for replicate in [1, 2, 3]:
                    curve_rows.append({"blind_id": blind_id, "gene_id": row["gene_id"], "dilution_level": level,
                        "technical_replicate": replicate, "assay_status": "not_conducted"})
rng.shuffle(reaction_rows)
fields = ["plate", "well", "blind_id", "gene_id", "reaction_type", "technical_replicate", "forward_primer", "reverse_primer",
          "expected_amplicon_bp", "genomic_DNA_product_detected", "primer_blast_review_status", "assay_status",
          "sample_id", "cq", "melting_single_peak", "gel_single_expected_band", "control_outcome", "raw_data_file"]
with (directory / "reaction_record_TEMPLATE.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for index, row in enumerate(reaction_rows):
        plate = index // 96 + 1
        position = index % 96
        writer.writerow({**row, "plate": plate, "well": f"{'ABCDEFGH'[position // 12]}{position % 12 + 1:02d}"})
rng.shuffle(curve_rows)
curve_fields = ["plate", "well", "blind_id", "gene_id", "dilution_level", "technical_replicate", "assay_status",
                "input_quantity", "input_unit", "sample_id", "cq", "raw_data_file"]
with (directory / "standard_curve_record_TEMPLATE.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=curve_fields)
    writer.writeheader()
    for index, row in enumerate(curve_rows):
        position = index % 96
        writer.writerow({**row, "plate": index // 96 + 1, "well": f"{'ABCDEFGH'[position // 12]}{position % 12 + 1:02d}"})
write_jsonl(directory / "blinding_key_restricted.jsonl", key)
write_jsonl(directory / "candidate_selection_ledger.jsonl", candidate_ledger)
oligos = {}
for row in reaction_rows:
    for side in ["forward", "reverse"]:
        sequence = row[f"{side}_primer"]
        oligos.setdefault(sequence, set()).add(f"{row['blind_id']}:{side}")
with (directory / "oligo_catalog.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["oligo_id", "sequence_5_to_3", "length_nt", "blind_assay_roles"])
    writer.writeheader()
    for index, sequence in enumerate(sorted(oligos), 1):
        writer.writerow({"oligo_id": f"O{index:04d}", "sequence_5_to_3": sequence, "length_nt": len(sequence),
                         "blind_assay_roles": ";".join(sorted(oligos[sequence]))})
write_json(directory / "status.json", {"prepared_at": utc_now(), "status": "draft_lab_plan_actual_resources_and_criteria_pending",
    "selected_genes": len(selected), "designed_assays": len(key), "planned_reactions": len(reaction_rows),
    "planned_standard_curve_reactions": len(curve_rows), "standard_curve_levels": 5,
    "unique_oligos": len(oligos),
    "standard_curve_quantity_policy": "Laboratory records actual positive quantities in a consistent unit at five dilution levels; no concentration or assay result is inferred",
    "screen": args.screen, "candidate_sources": source_records,
    "candidate_selection_rule": "Forty frozen random test genes; two unique RNA-screen-passing feasible pairs per method, ranked by common temperature penalty and candidate identifier",
    "panel_artifacts_sha256": {name: sha256(directory / name) for name in [
        "target_selection.jsonl", "candidate_selection_ledger.jsonl", "blinding_key_restricted.jsonl",
        "reaction_record_TEMPLATE.csv", "standard_curve_record_TEMPLATE.csv", "oligo_catalog.csv"]},
    "no_design_slots": sum(row["status"] == "no_design" for row in candidate_ledger),
    "selection_sha256": sha256(selection_path), "experiments_conducted": False,
    "warning": "Reaction records are unfilled templates. No oligos have been ordered and no samples processed.",
    "primer3_candidates_sha256": sha256(root / args.primer3_directory / "candidates.jsonl"),
    "diffusion_candidates_sha256": sha256(root / args.diffusion_directory / "candidates.jsonl")})
print(str(directory / "reaction_record_TEMPLATE.csv"))
