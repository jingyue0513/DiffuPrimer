#!/usr/bin/env python3
"""Complete the declared 100-epoch duration check for all model seeds."""
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
experiment = root / "experiments/epoch100"
protocol = json.loads((experiment / "protocol.json").read_text())
state_path = experiment / "results/three_seed_status.json"
script = experiment / "scripts/train.py"
if sha256(script) != protocol["script_hashes"]["train.py"]:
    raise ValueError("Training-duration implementation changed")
extension = experiment / "three_seed_protocol.json"
planned = {
    "seeds": [20261013, 20261014, 20261015], "methods": ["direct", "diffusion"],
    "epochs": 100, "primary_epochs": 30, "reused_completed_models": 2,
    "checkpoint_selection": "Lowest validation objective within each fixed training trajectory",
    "dataset_sha256": protocol["dataset_sha256"], "split_sha256": protocol["split_sha256"],
    "training_script_sha256": sha256(script), "initial_protocol_sha256": sha256(experiment / "protocol.json"),
    "primary_results_retained": True, "test_outcomes_used_for_checkpoint_selection": False,
    "concurrency": 2, "threads_per_model": 8,
}
if extension.exists():
    if json.loads(extension.read_text())["plan"] != planned:
        raise ValueError("Previously recorded duration extension differs")
else:
    write_json(extension, {"frozen_at": utc_now(), "plan": planned})
logs = experiment / "results/training_logs"
logs.mkdir(parents=True, exist_ok=True)
write_json(state_path, {"stage": "training", "at": utc_now(), "planned_models": 6})


def run(job):
    method, seed = job
    name = f"{method}_weak_partial_pilot_seed{seed}"
    history_path = experiment / f"results/{name}_training.json"
    if not history_path.exists() or len(json.loads(history_path.read_text())["history"]) != 100:
        command = [str(root / ".venv/bin/python"), str(script), "--mode", method, "--seed", str(seed),
                   "--epochs", "100", "--dataset", "data/processed/dataset_partial_pilot.jsonl",
                   "--allow-weak", "--partial-source"]
        with (logs / f"{name}.log").open("w") as handle:
            subprocess.run(command, cwd=root, env=dict(os.environ, PYTHONPATH=str(root / "src")),
                           check=True, stdout=handle, stderr=subprocess.STDOUT)
    extended = json.loads(history_path.read_text())
    primary = json.loads((root / f"results/{name}_training.json").read_text())
    if len(extended["history"]) != 100:
        raise ValueError("Incomplete duration trajectory")
    for a, b in zip(primary["history"], extended["history"][:30]):
        if any(a[key] != b[key] for key in ["epoch", "training_loss", "validation_loss"]):
            raise ValueError("First 30 epochs no longer reproduce the primary trajectory")
    first = min(primary["history"], key=lambda r: r["validation_loss"])
    best = min(extended["history"], key=lambda r: r["validation_loss"])
    result = {"method": method, "seed": seed, "completed_epochs": 100,
              "best_primary_epoch": first["epoch"], "best_extended_epoch": best["epoch"],
              "best_primary_validation_loss": first["validation_loss"],
              "best_extended_validation_loss": best["validation_loss"],
              "first_30_epochs_exactly_match": True,
              "checkpoint_sha256": sha256(Path(extended["checkpoint"])),
              "primary_checkpoint_sha256": sha256(Path(primary["checkpoint"])),
              "history_sha256": sha256(history_path)}
    print(json.dumps(result), flush=True)
    return result


try:
    jobs = [(method, seed) for seed in planned["seeds"] for method in planned["methods"]]
    with ThreadPoolExecutor(max_workers=2) as executor:
        models = list(executor.map(run, jobs))
    write_json(state_path, {"stage": "all_three_seed_training_complete", "at": utc_now(), "models": models,
                           "protocol_sha256": sha256(extension), "test_evaluation_pending": True})
except Exception as exc:
    write_json(state_path, {"stage": "training_error", "at": utc_now(), "error": str(exc)})
    raise
