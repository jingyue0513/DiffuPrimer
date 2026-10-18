#!/usr/bin/env python3
"""Three-seed full-validation experiment for matched joint-pair guidance."""
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


experiment = Path(__file__).resolve().parent
root = experiment.parents[1]
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(experiment / "src"))
protocol = json.loads((experiment / "implemented_protocol.json").read_text())
for relative, checksum in protocol["sha256"].items():
    if sha256(experiment / relative) != checksum:
        raise ValueError("Joint-guidance implementation changed")
check = (experiment / "results/checks_correct_import.log").read_text()
if "23 passed" not in check:
    raise ValueError("Joint-guidance checks have not passed")
jobs = []
directories = []
for seed in [20261013, 20261014, 20261015]:
    for mode in ["direct", "diffusion", "rule_sampler"]:
        command = [python, str(experiment / "scripts/evaluate.py"), "--split", "validation", "--seed", str(seed), "--joint-tm-guidance"]
        if mode == "rule_sampler":
            command += ["--rule-sampler"]
            stage = "pilot"
        else:
            command += ["--checkpoint", str(root / f"results/checkpoints/{mode}_weak_partial_pilot_seed{seed}.pt"), "--oligo-guidance"]
            stage = "weak_partial_pilot"
        label = f"{mode}_{seed}"
        directory = experiment / f"results/evaluations/{mode}_{stage}_validation_seed{seed}_amp80-200_oligo_guided_joint_tm"
        jobs.append((label, command))
        directories.append(str(directory))
state = {"stage": "full_validation", "started_at": utc_now(), "completed": [], "directories": directories}
status_path = experiment / "results/status.json"
write_json(status_path, state)


def run(label, command):
    with (experiment / "results" / (label + ".log")).open("w") as handle:
        subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    return label


try:
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(run, *job) for job in jobs]
        for future in as_completed(futures):
            state["completed"].append(future.result())
            state["updated_at"] = utc_now()
            write_json(status_path, state)
            print(json.dumps({"completed": state["completed"][-1], "at": state["updated_at"]}), flush=True)
    baseline = root / "results/evaluations/primer3_pilot_validation_seed20261013_amp80-200"
    subprocess.run([python, "scripts/inspect_design_results.py", *directories, str(baseline), "--output", "experiments/joint_tm/results/validation_diagnostics.json"],
        cwd=root, env=dict(env, PYTHONPATH=str(root / "src")), check=True)
    state["stage"] = "validation_complete"
    write_json(status_path, state)
except Exception as exc:
    state["stage"], state["error"] = "validation_error", str(exc)
    write_json(status_path, state)
    raise
