#!/usr/bin/env python3
"""Calibrate learned sampling temperature on the complete validation partition."""
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
directory = root / "results/temperature_tuning"
directory.mkdir(parents=True, exist_ok=True)
python = str(root / ".venv/bin/python")
joint = root / "experiments/joint_tm"
env = dict(os.environ, PYTHONPATH=str(joint / "src"))
jobs, directories = [], []
for temperature, label in [(0.5, "05"), (2.0, "2"), (4.0, "4")]:
    experiment = root / ("experiments/temperature_" + label)
    (experiment / "scripts").mkdir(parents=True, exist_ok=True)
    (experiment / "results").mkdir(exist_ok=True)
    for name, target in [("data", root / "data"), ("configs", root / "configs"), ("src", joint / "src")]:
        if not (experiment / name).exists():
            (experiment / name).symlink_to(target, target_is_directory=True)
    script = experiment / "scripts/evaluate.py"
    if script.exists() and sha256(script) != sha256(joint / "scripts/evaluate.py"):
        raise ValueError("Temperature experiment evaluation source differs")
    shutil.copyfile(joint / "scripts/evaluate.py", script)
    for seed in [20261013, 20261014, 20261015]:
        for mode in ["direct", "diffusion"]:
            command = [python, str(root / "scripts/evaluate_cached_spans.py"), "--evaluation-script", str(script),
                "--split", "validation", "--seed", str(seed), "--joint-tm-guidance", "--oligo-guidance",
                "--temperature", str(temperature), "--checkpoint", str(root / f"results/checkpoints/{mode}_weak_partial_pilot_seed{seed}.pt")]
            job_label = f"{mode}_seed{seed}_temperature{label}"
            jobs.append((job_label, command))
            directories.append(str(experiment / f"results/evaluations/{mode}_weak_partial_pilot_validation_seed{seed}_amp80-200_oligo_guided_joint_tm"))
plan = {"frozen_at": utc_now(), "temperature_levels": [0.5, 1.0, 2.0, 4.0], "seeds": [20261013, 20261014, 20261015],
    "genes_per_run": 75, "budget": 20, "new_validation_runs": len(jobs), "default_temperature_policy": "Reuse the six completed joint-Tm learned validation runs at temperature 1",
    "selection_rule": "For each learned method, maximize mean unique RNA-screen-pass pairs across all validation genes and seeds; ties minimize mean best RNA-pass shared-temperature penalty, then lower sampling temperature",
    "selection_screen": "Independent all-strata full-length ungapped mismatch screen",
    "selected_test_policy": "Evaluate each selected learned method on all 75 test genes with all three frozen checkpoints; reuse temperature-1 runs if selected",
    "checkpoint_policy": "Thirty-epoch frozen primary checkpoints unchanged; this calibrates sampling only",
    "joint_implementation_sha256": sha256(joint / "implemented_protocol.json"),
    "driver_sha256": sha256(Path(__file__)), "commands": [command for _, command in jobs]}
plan_path = directory / "protocol.json"
if plan_path.exists():
    if json.loads(plan_path.read_text())["commands"] != plan["commands"]:
        raise ValueError("Temperature validation recipe changed")
else:
    write_json(plan_path, plan)
state = {"stage": "validation_sampling", "started_at": utc_now(), "directories": directories, "completed": []}
status_path = directory / "status.json"
write_json(status_path, state)


def run(label, command):
    with (directory / (label + ".log")).open("a") as handle:
        subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    return label


try:
    with ThreadPoolExecutor(max_workers=3) as executor:
        for future in as_completed([executor.submit(run, *job) for job in jobs]):
            state["completed"].append(future.result())
            state["updated_at"] = utc_now()
            write_json(status_path, state)
            print(json.dumps({"completed": state["completed"][-1], "at": state["updated_at"]}), flush=True)
    state["stage"] = "validation_sampling_complete_specificity_pending"
    write_json(status_path, state)
except Exception as exc:
    state["stage"], state["error"] = "validation_sampling_error", str(exc)
    write_json(status_path, state)
    raise
