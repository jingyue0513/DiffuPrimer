#!/usr/bin/env python3
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from diffuprimer.acquisition import load_reference_sequences
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.sequence import PrimerPair


parser = argparse.ArgumentParser()
parser.add_argument("--allow-incomplete", action="store_true", help="Calibration only; cannot serve as the frozen final dataset")
parser.add_argument("--freeze-partial", action="store_true", help="Freeze an explicitly incomplete source pilot with recorded source unavailability")
parser.add_argument("--include-teacher", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
selection_path = root / "data/processed/pilot_gene_selection.jsonl"
source_path = root / "data/processed/primerbank_pilot.jsonl"
split_path = root / "data/splits/pilot.jsonl"
selection = read_jsonl(selection_path)
ledger = {x["gene_id"]: x["status"] for x in read_jsonl(root / "data/processed/acquisition_ledger.jsonl")}
incomplete = [x["gene_id"] for x in selection if ledger.get(x["gene_id"]) not in {"complete", "no_records_returned"}]
if incomplete and not args.allow_incomplete and not args.freeze_partial:
    raise ValueError(f"Acquisition incomplete for {len(incomplete)} genes; final dataset is not frozen")
split_rows = read_jsonl(split_path)
splits = {x["gene_id"]: x for x in split_rows}
split_manifest = json.loads(split_path.with_suffix(".manifest.json").read_text())
if split_manifest["split_sha256"] != sha256(split_path) or split_manifest["selection_sha256"] != sha256(selection_path):
    raise ValueError("Frozen split integrity failure")
source_bytes = source_path.read_bytes()
source_hash = hashlib.sha256(source_bytes).hexdigest()
source = [json.loads(line) for line in source_bytes.splitlines() if line.strip()]
snapshot = root / f"data/processed/acquisition_snapshot_{source_hash[:12]}.jsonl"
snapshot.write_bytes(source_bytes)
evidence_path = root / "data/processed/published_evidence_examples.jsonl"
published_evidence = {row["primerbank_id"]: row for row in read_jsonl(evidence_path)} if evidence_path.exists() else {}
dataset, exclusions, seen = [], [], set()
for row in source:
    if row["primerbank_id"] in published_evidence:
        assertion = published_evidence[row["primerbank_id"]]
        if any(row[key] != assertion[key] for key in ["forward_primer", "reverse_primer", "amplicon_length"]):
            raise ValueError("Published outcome does not match retrieved primer sequences")
        row = {**row, **{key: assertion[key] for key in ["evidence_type", "experimental_result", "failure_reason"]},
               "outcome_evidence_source": assertion["source_url"], "outcome_evidence_sha256": assertion["source_pdf_sha256"]}
    if row["mapping_status"] != "unique_exact_pair_remap":
        exclusions.append({"primerbank_id": row["primerbank_id"], "reason": row["mapping_status"]})
        continue
    pair = PrimerPair(*(row[key] for key in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
    pair.validate(row["template_sequence"])
    if pair.sequences(row["template_sequence"]) != (row["forward_primer"], row["reverse_primer"]):
        raise ValueError("Mapped coordinates do not reproduce synthesis sequences")
    key = (row["gene_id"], row["forward_primer"], row["reverse_primer"], row["transcript_accession_version"])
    if key in seen:
        exclusions.append({"primerbank_id": row["primerbank_id"], "reason": "duplicate_pair"})
        continue
    seen.add(key)
    dataset.append({**row, "split": splits[row["gene_id"]]["split"], "cluster_id": splits[row["gene_id"]]["cluster_id"],
                    "training_label_semantics": "position imitation; experimental success only when separately supported"})
if args.freeze_partial and not (root / "results/source_access_block.json").exists():
    raise ValueError("Partial-source freeze requires a documented source-access block")
if args.include_teacher:
    teacher_path = root / "data/processed/computational_candidates.jsonl"
    teacher_manifest = json.loads(teacher_path.with_suffix(".manifest.json").read_text())
    if teacher_manifest["dataset_sha256"] != sha256(teacher_path):
        raise ValueError("Teacher dataset integrity failure")
    for row in read_jsonl(teacher_path):
        if row["split"] == "test":
            raise ValueError("Test gene present in teacher supervision")
        if splits[row["gene_id"]]["split"] != row["split"]:
            raise ValueError("Teacher row violates frozen split")
        key = (row["gene_id"], row["forward_primer"], row["reverse_primer"], row["transcript_accession_version"])
        if key not in seen:
            seen.add(key)
            dataset.append(row)
suffix = "partial_pilot" if args.freeze_partial else ("calibration" if incomplete else "pilot")
output = root / f"data/processed/dataset_{suffix}.jsonl"
write_jsonl(output, dataset)
write_jsonl(root / f"data/processed/exclusions_{suffix}.jsonl", exclusions)
accessions = {x["reference_transcripts"][0] for x in selection}
seqs = load_reference_sequences(root / "data/raw/refseq/GCF_000001635.27_GRCm39_rna.fna.gz", accessions)
targets = [{**x, "transcript_accession_version": x["reference_transcripts"][0],
            "template_sequence": seqs[x["reference_transcripts"][0]], **splits[x["gene_id"]],
            "target_selection_rule": "lexicographically first annotated NM_ accession fixed before observing PrimerBank records"}
           for x in selection if x["reference_transcripts"][0] in seqs]
write_jsonl(root / "data/processed/targets_pilot.jsonl", targets)
manifest = {"created_at": utc_now(), "stage": "frozen partial-source pilot" if args.freeze_partial else ("calibration" if incomplete else "acquisition-complete pilot"),
            "partial_source_frozen": args.freeze_partial,
            "incomplete_genes": incomplete, "source_sha256": source_hash, "source_snapshot": str(snapshot),
            "dataset_sha256": sha256(output), "split_sha256": sha256(split_path),
            "target_sha256": sha256(root / "data/processed/targets_pilot.jsonl"),
            "published_evidence_sha256": sha256(evidence_path) if evidence_path.exists() else None,
            "teacher_sha256": sha256(root / "data/processed/computational_candidates.jsonl") if args.include_teacher else None,
            "rows": len(dataset), "by_split": dict(Counter(x["split"] for x in dataset)),
            "by_evidence": dict(Counter(x["evidence_type"] for x in dataset)), "exclusions": len(exclusions)}
write_json(output.with_suffix(".manifest.json"), manifest)
print(json.dumps(manifest, indent=2))
