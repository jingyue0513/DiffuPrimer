#!/usr/bin/env python3
"""A declared Primer3 teacher dataset. No computed pair is an assay success."""
import json
import math
from pathlib import Path

import numpy as np

from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.sequence import PrimerPair


root = Path(__file__).resolve().parents[1]
config = json.loads((root / "configs/pilot.json").read_text())
conditions = Conditions(**config["conditions"])
target_path = root / "data/processed/targets_pilot.jsonl"
targets = [row for row in read_jsonl(target_path) if row["split"] in {"train", "validation"}]
destination = root / "data/processed/computational_candidates.jsonl"
ledger_path = root / "data/processed/computational_candidates_ledger.jsonl"
rows = read_jsonl(destination) if destination.exists() else []
ledger = read_jsonl(ledger_path) if ledger_path.exists() else []
completed = {row["gene_id"] for row in ledger}
for target in targets:
    if target["gene_id"] in completed:
        continue
    template = target["template_sequence"]
    window = config["model"]["window_length"]
    count = min(5, max(1, math.ceil(len(template) / window)))
    starts = sorted(set(int(value) for value in np.linspace(0, max(0, len(template) - window), count)))
    regions = [(offset, template[offset:offset + window]) for offset in starts]
    if len(template) > window:
        regions.append((0, template))
    candidates, seen = [], set()
    for offset, sequence in regions:
        proposals, _ = primer3_candidates(sequence, conditions, budget=4)
        for proposal in proposals:
            coordinates = [proposal[key] + offset for key in ["forward_start", "forward_end", "reverse_start", "reverse_end"]]
            pair = PrimerPair(*coordinates)
            f, r = pair.sequences(template)
            if (f, r) in seen:
                continue
            seen.add((f, r))
            assessment = assess_pair(template, pair, conditions)
            if not assessment["computational_feasible"]:
                continue
            candidates.append({"species": "Mus musculus", "gene_id": target["gene_id"],
                "transcript_accession_version": target["transcript_accession_version"],
                "template_sequence": template, **pair.as_dict(template),
                "split": target["split"], "cluster_id": target["cluster_id"],
                "evidence_type": "computational_candidate", "experimental_result": "not_tested",
                "experimental_record_available": False, "source_database": "Primer3 local teacher",
                "source_url": "https://primer3.org/manual.html", "primerbank_id": None,
                "reference_assembly": config["reference_assembly"], "annotation_version": config["annotation"],
                "teacher_primer3_penalty": proposal["primer3_penalty"], "teacher_window_offset": offset,
                "teacher_conditions": config["conditions"], "assessment": assessment,
                "training_label_semantics": "Rule-based position imitation; no experimental efficacy evidence"})
    rows.extend(candidates)
    ledger.append({"gene_id": target["gene_id"], "candidates": len(candidates), "at": utc_now()})
    write_jsonl(destination, rows)
    write_jsonl(ledger_path, ledger)
    print(json.dumps(ledger[-1]), flush=True)
write_json(destination.with_suffix(".manifest.json"), {
    "generated_at": utc_now(), "target_sha256": sha256(target_path), "config_sha256": sha256(root / "configs/pilot.json"),
    "dataset_sha256": sha256(destination), "genes_attempted": len(ledger), "candidates": len(rows),
    "target_selection": "All frozen train and validation genes; no test genes used for teacher supervision",
    "teacher": "Primer3 in up to five deterministic 512-base windows plus full transcript, up to four pairs per region",
    "limitations": "Teacher supervision may transmit Primer3 preferences. Comparison with that teacher does not establish independent experimental generalization."})
