#!/usr/bin/env python3
"""Freeze the computational recipe before opening final test comparisons."""
import json
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
files = [root / "configs/pilot.json", root / "configs/analysis_plan.json", root / "requirements.lock"]
files += sorted((root / "src/diffuprimer").glob("*.py"))
files += [root / "scripts" / name for name in ["train.py", "evaluate.py", "screen_specificity.py", "aggregate.py", "generate_weak_supervision.py", "prepare_data.py"]]
files += [root / "data/processed/dataset_partial_pilot.jsonl", root / "data/processed/targets_pilot.jsonl", root / "data/splits/pilot.jsonl"]
hashes = {str(path.relative_to(root)): sha256(path) for path in files}
destination = root / "results/frozen_protocol.json"
if destination.exists():
    existing = json.loads(destination.read_text())
    if existing["sha256"] != hashes:
        raise ValueError("Frozen protocol changed; document a versioned correction before new test runs")
    print("Frozen protocol unchanged")
else:
    write_json(destination, {"frozen_at": utc_now(), "version": "weak-partial-pilot-v1", "sha256": hashes,
                            "source_limitation": "PrimerBank acquisition incomplete after HTTP 403; cached source records plus declared Primer3 teacher supervision",
                            "laboratory_status": "not conducted", "test_targets": 75})
    print(str(destination))
