#!/usr/bin/env python3
"""Incrementally audit completed experiments without changing frozen screens."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
directory = root / "results/bowtie_audit"
directory.mkdir(parents=True, exist_ok=True)
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
state_path = directory / "driver_status.json"


def status(stage, **extra):
    write_json(state_path, {"stage": stage, "at": utc_now(), **extra})
    print(json.dumps({"stage": stage, "at": utc_now(), **extra}), flush=True)


def execute(label, command):
    with (directory / (label + ".log")).open("a") as handle:
        subprocess.run([python, *command], cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)


try:
    status("waiting_for_reference_indexes")
    for name in ["transcriptome", "genome"]:
        while not (root / f"data/processed/bowtiedb/{name}.manifest.json").exists():
            time.sleep(5)
    last_evaluation_set = None
    while True:
        evaluations = []
        bases = [root / "results/evaluations", root / "experiments/joint_tm/results/evaluations"]
        bases += [root / "experiments/uniform_pairs/results/evaluations"]
        bases += [path / "results/evaluations" for path in sorted((root / "experiments").glob("temperature_*"))]
        for base in bases:
            for completion_path in base.glob("*/completion.json"):
                path = completion_path.parent
                manifest = json.loads((path / "manifest.json").read_text())
                if manifest["stage"] == "calibration":
                    continue
                evaluations.append(path)
        pending = [path for path in evaluations if not (path / "specificity_bowtie_status.json").exists()]
        if pending:
            status("screening_completed_runs", total_completed_evaluations=len(evaluations), pending_runs=len(pending))
            execute("screening", ["scripts/screen_bowtie.py", *map(str, pending), "--threads", "4"])
        current_evaluation_set = tuple(sorted(map(str, evaluations)))
        changed = current_evaluation_set != last_evaluation_set or bool(pending)
        if changed:
            for split in ["validation", "test"]:
                own = [path for path in evaluations if json.loads((path / "manifest.json").read_text())["arguments"]["split"] == split]
                if own:
                    execute(split + "_statistics", ["scripts/summarize_specificity_audit.py", *map(str, own), "--split", split,
                        "--output", f"results/bowtie_audit/{split}_summary.json"])
            last_evaluation_set = current_evaluation_set
        primary = json.loads((root / "results/study_status.json").read_text())
        joint = json.loads((root / "experiments/joint_tm/results/test_status.json").read_text())
        temperature = json.loads((root / "results/temperature_tuning/selection_status.json").read_text())
        uniform_path = root / "experiments/uniform_pairs/results/status.json"
        uniform = json.loads(uniform_path.read_text()) if uniform_path.exists() else {"stage": "not_started"}
        condition_path = root / "results/temperature_tuning/selected_conditions_status.json"
        condition = json.loads(condition_path.read_text()) if condition_path.exists() else {"stage": "not_started"}
        if primary["stage"] == "pipeline_error" or joint["stage"] == "test_sampling_error":
            raise ValueError("Experiment driver failed")
        if temperature["stage"] == "selection_error":
            raise ValueError("Temperature selection driver failed")
        if uniform["stage"] == "test_sampling_error":
            raise ValueError("Uniform compatible-pair baseline failed")
        if condition["stage"] == "selected_condition_error":
            raise ValueError("Selected-temperature condition study failed")
        finished = primary["stage"] == "computational_pilot_complete_independent_qpcr_pending" and joint["stage"] == "test_sampling_complete_specificity_pending"
        finished &= temperature["stage"] in ["selected_test_sampling_complete_specificity_pending", "selected_test_analyzed"]
        finished &= uniform["stage"] == "test_sampling_complete_specificity_pending"
        finished &= condition["stage"] == "selected_condition_sampling_complete_specificity_pending"
        if finished:
            # Producers may finish new evaluations during a long-running screen.
            # Re-scan only after all producer completion barriers are published.
            latest = []
            for base in bases:
                for completion_path in base.glob("*/completion.json"):
                    candidate = completion_path.parent
                    if json.loads((candidate / "manifest.json").read_text())["stage"] != "calibration":
                        latest.append(candidate)
            if tuple(sorted(map(str, latest))) != current_evaluation_set or any(not (path / "specificity_bowtie_status.json").exists() for path in latest):
                last_evaluation_set = None
                continue
            status("all_completed_experiments_audited", completed_evaluations=len(evaluations), manuscript_work_paused=True)
            break
        if changed:
            status("waiting_for_remaining_experiments", audited_runs=len(evaluations))
        time.sleep(15)
except Exception as exc:
    status("audit_error", error=str(exc))
    raise
