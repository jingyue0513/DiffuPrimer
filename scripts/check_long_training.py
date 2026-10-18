#!/usr/bin/env python3
"""Validate the training-duration experiment after both 100-epoch runs finish."""
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import utc_now, write_json


root = Path(__file__).resolve().parents[1]
experiment = root / "experiments/epoch100"
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
status_path = experiment / "results/status.json"
write_json(status_path, {"stage": "waiting_for_training", "at": utc_now()})
for mode in ["direct", "diffusion"]:
    history = experiment / f"results/{mode}_weak_partial_pilot_seed20261013_training.json"
    while not history.exists() or len(json.loads(history.read_text())["history"]) < 100:
        time.sleep(5)


def run(mode):
    command = [python, str(experiment / "scripts/evaluate.py"), "--split", "validation", "--seed", "20261013",
        "--checkpoint", f"results/checkpoints/{mode}_weak_partial_pilot_seed20261013.pt", "--oligo-guidance"]
    with (experiment / "results" / (mode + "_validation.log")).open("w") as handle:
        subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    return mode


write_json(status_path, {"stage": "full_validation", "at": utc_now()})
try:
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(run, mode) for mode in ["direct", "diffusion"]]
        for future in as_completed(futures):
            print(json.dumps({"completed": future.result(), "at": utc_now()}), flush=True)
    directories = [str(experiment / f"results/evaluations/{mode}_weak_partial_pilot_validation_seed20261013_amp80-200_oligo_guided") for mode in ["direct", "diffusion"]]
    subprocess.run([python, "scripts/inspect_design_results.py", *directories, "--output", "experiments/epoch100/results/validation_diagnostics.json"],
        cwd=root, env=env, check=True)
    write_json(status_path, {"stage": "validation_complete", "at": utc_now(), "directories": directories})
except Exception as exc:
    write_json(status_path, {"stage": "validation_error", "at": utc_now(), "error": str(exc)})
    raise
