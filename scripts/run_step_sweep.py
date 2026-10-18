#!/usr/bin/env python3
"""Evaluate the number of reverse steps using the already trained denoisers."""
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
directory = root / "results/step_sweep"
directory.mkdir(parents=True, exist_ok=True)
planned = {"frozen_at": utc_now(), "hypothesis": "Iterative denoising changes candidate uniqueness and positional diversity under a fixed output budget",
    "reverse_steps": [1, 4, 8, 16], "seeds": [20261013, 20261014, 20261015], "test_genes": 75,
    "candidate_budget": 20, "trained_steps": 16, "training_unchanged": True,
    "time_policy": "Coarser inference schedules use normalized times contained in the 16-step training grid",
    "source_protocol_sha256": sha256(root / "results/frozen_protocol.json"), "sixteen_step_runs": "Reuse the three primary diffusion test runs"}
write_json(directory / "protocol.json", planned)
jobs, directories = [], []
for seed in planned["seeds"]:
    for steps in [1, 4, 8]:
        command = [python, "scripts/evaluate.py", "--split", "test", "--seed", str(seed),
            "--checkpoint", f"results/checkpoints/diffusion_weak_partial_pilot_seed{seed}.pt", "--oligo-guidance", "--steps", str(steps)]
        label = f"seed{seed}_steps{steps}"
        jobs.append((label, command))
        directories.append(str(root / f"results/evaluations/diffusion_weak_partial_pilot_test_seed{seed}_amp80-200_oligo_guided_steps{steps}"))
state = {"stage": "sampling", "started_at": utc_now(), "completed": [], "directories": directories}


def run(label, command):
    with (directory / (label + ".log")).open("w") as handle:
        subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    return label


try:
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(run, *job) for job in jobs]
        for future in as_completed(futures):
            state["completed"].append(future.result())
            state["updated_at"] = utc_now()
            write_json(directory / "status.json", state)
            print(json.dumps({"completed": state["completed"][-1], "at": state["updated_at"]}), flush=True)
    state["stage"] = "sampling_complete_specificity_pending"
    write_json(directory / "status.json", state)
except Exception as exc:
    state["stage"], state["error"] = "sampling_error", str(exc)
    write_json(directory / "status.json", state)
    raise
