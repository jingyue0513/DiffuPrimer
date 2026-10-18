#!/usr/bin/env python3
"""Freeze the independent assay panel from completed, selected test designs."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
directory = root / "results/lab_plan"
directory.mkdir(parents=True, exist_ok=True)
state = directory / "panel_driver_status.json"
selection_path = root / "results/temperature_tuning/selection.json"
write_json(state, {"stage": "waiting_for_selected_design_screens", "at": utc_now()})
while not selection_path.exists():
    time.sleep(10)
selection_hash = sha256(selection_path)
selection = json.loads(selection_path.read_text())
temperature = selection["selected"]["diffusion"]["temperature"]
experiment = {0.5: "temperature_05", 1.0: "joint_tm", 2.0: "temperature_2", 4.0: "temperature_4"}[temperature]
diffusion = root / f"experiments/{experiment}/results/evaluations/diffusion_weak_partial_pilot_test_seed20261013_amp80-200_oligo_guided_joint_tm"
primer3 = root / "results/evaluations/primer3_pilot_test_seed20261013_amp80-200"
while not all((path / "specificity_bowtie_status.json").exists() for path in [diffusion, primer3]):
    audit_path = root / "results/bowtie_audit/driver_status.json"
    if audit_path.exists() and json.loads(audit_path.read_text())["stage"] == "audit_error":
        write_json(state, {"stage": "upstream_error", "at": utc_now()})
        raise ValueError("Independent screening requires repair")
    time.sleep(10)
if sha256(selection_path) != selection_hash:
    raise ValueError("Model selection changed while preparing the assay panel")
command = [str(root / ".venv/bin/python"), "scripts/prepare_lab_plan.py", "--screen", "bowtie",
           "--primer3-directory", str(primer3), "--diffusion-directory", str(diffusion)]
with (directory / "panel_preparation.log").open("w") as handle:
    subprocess.run(command, cwd=root, env=dict(os.environ, PYTHONPATH=str(root / "src")),
                   check=True, stdout=handle, stderr=subprocess.STDOUT)
write_json(state, {"stage": "panel_and_reaction_templates_prepared", "at": utc_now(),
                  "selection_sha256": selection_hash, "diffusion_temperature": temperature,
                  "model_seed": 20261013, "command": command,
                  "lab_plan_status_sha256": sha256(directory / "status.json"),
                  "experiments_conducted": False})
