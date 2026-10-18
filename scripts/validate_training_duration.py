#!/usr/bin/env python3
"""Evaluate all 100-epoch checkpoints in the original validation setting."""
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
experiment = root / "experiments/epoch100"
state_path = experiment / "results/three_seed_validation_status.json"
write_json(state_path, {"stage": "waiting_for_training", "at": utc_now()})
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
try:
    training_path = experiment / "results/three_seed_status.json"
    while True:
        stage = json.loads(training_path.read_text())["stage"] if training_path.exists() else "not_started"
        if stage == "all_three_seed_training_complete":
            break
        if stage == "training_error":
            raise ValueError("Duration training failed")
        time.sleep(10)
    write_json(state_path, {"stage": "validation_sampling", "at": utc_now()})
    jobs = [(method, seed) for seed in [20261013, 20261014, 20261015] for method in ["direct", "diffusion"]]

    def run(job):
        method, seed = job
        checkpoint = experiment / f"results/checkpoints/{method}_weak_partial_pilot_seed{seed}.pt"
        output = experiment / f"results/evaluations/{method}_weak_partial_pilot_validation_seed{seed}_amp80-200_oligo_guided"
        if not (output / "completion.json").exists():
            command = [python, str(root / "scripts/evaluate_cached_spans.py"), "--evaluation-script",
                       str(experiment / "scripts/evaluate.py"), "--split", "validation", "--seed", str(seed),
                       "--oligo-guidance", "--checkpoint", str(checkpoint)]
            with (experiment / f"results/training_logs/validation_{method}_{seed}.log").open("w") as handle:
                subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
        manifest = json.loads((output / "manifest.json").read_text())
        completion = json.loads((output / "completion.json").read_text())
        if manifest["checkpoint_sha256"] != sha256(checkpoint) or not completion["all_expected_targets_evaluated"]:
            raise ValueError("Duration evaluation provenance differs")
        return str(output)

    with ThreadPoolExecutor(max_workers=2) as executor:
        directories = list(executor.map(run, jobs))
    with (experiment / "results/three_seed_validation_analysis.log").open("w") as handle:
        subprocess.run([python, "scripts/compare_validation_experiments.py"], cwd=root, env=env,
                       check=True, stdout=handle, stderr=subprocess.STDOUT)
        subprocess.run([python, "scripts/plot_training_curves.py"], cwd=root, env=env,
                       check=True, stdout=handle, stderr=subprocess.STDOUT)
    write_json(state_path, {"stage": "validation_comparison_complete", "at": utc_now(), "directories": directories,
                           "summary_sha256": sha256(root / "results/validation_experiment_comparison.json"),
                           "scope": "Training-duration sensitivity on validation genes; thermodynamic endpoints, not new qPCR observations"})
except Exception as exc:
    write_json(state_path, {"stage": "duration_validation_error", "at": utc_now(), "error": str(exc)})
    raise
