#!/usr/bin/env python3
"""Finish specificity and paired summaries for the inference-step experiment."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import utc_now, write_json


root = Path(__file__).resolve().parents[1]
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
directory = root / "results/step_sweep"
write_json(directory / "analysis_status.json", {"stage": "waiting_for_completed_runs", "at": utc_now()})
status_path = directory / "status.json"
while not status_path.exists() or json.loads(status_path.read_text())["stage"] != "sampling_complete_specificity_pending":
    if status_path.exists() and json.loads(status_path.read_text())["stage"] == "sampling_error":
        raise RuntimeError("Step-sweep sampling failed")
    time.sleep(5)
aggregate = root / "results/primary_aggregate.json"
while not aggregate.exists():
    if json.loads((root / "results/study_status.json").read_text())["stage"] == "pipeline_error":
        raise RuntimeError("Core study failed")
    time.sleep(5)
all_primary = json.loads(aggregate.read_text())["runs"]
primary_sources = [run["source"] for run in all_primary]
sixteen = [run["source"] for run in all_primary if run["method"] == "diffusion"]
extra = json.loads(status_path.read_text())["directories"]
write_json(directory / "analysis_status.json", {"stage": "waiting_for_full_length_reference_audit", "at": utc_now()})
while not all((Path(path) / "specificity_bowtie_status.json").exists() for path in sixteen + extra):
    if json.loads((root / "results/bowtie_audit/driver_status.json").read_text())["stage"] == "audit_error":
        raise RuntimeError("Full-length inference-step specificity audit failed")
    time.sleep(5)
with (directory / "statistics.log").open("w") as handle:
    subprocess.run([python, "scripts/summarize_secondary.py", *sixteen, *extra, "--screen", "bowtie", "--output", "results/step_sweep/summary.json"],
        cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
write_json(directory / "analysis_status.json", {"stage": "analyzed", "at": utc_now(), "step_levels": [1, 4, 8, 16],
    "seeds_per_level": 3, "specificity_screen": "all_strata_full_length_bowtie"})
