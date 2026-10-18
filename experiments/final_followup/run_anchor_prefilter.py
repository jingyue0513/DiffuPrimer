#!/usr/bin/env python3
"""Reproduce the focused anchored-support correction with the fixed seed."""
import os
import subprocess
from pathlib import Path
from diffuprimer.io import utc_now, write_json

root = Path(__file__).resolve().parents[2]
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
directories = []
for task in ["completion", "redesign"]:
    workers = []
    for method in ["direct", "diffusion"]:
        directory = root / f"experiments/final_followup/results/evaluations/{method}_{task}_reverse_seed20261013_anchor_prefilter"
        directories.append(directory)
        if (directory / "completion.json").exists():
            continue
        log = (root / f"experiments/final_followup/results/{method}_{task}_reverse_anchor_prefilter.log").open("a")
        command = [python, "experiments/final_followup/evaluate.py", "--task", task, "--method", method,
                   "--fixed-side", "reverse", "--anchor-prefilter", "--seed", "20261013"]
        workers.append((subprocess.Popen(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT), log))
    for worker, log in workers:
        code = worker.wait()
        log.close()
        if code:
            raise RuntimeError("Anchored-support generation failed")
pending = [d for d in directories if not (d / "specificity_bowtie_status.json").exists()]
if pending:
    with (root / "experiments/final_followup/results/anchor_prefilter_screening.log").open("a") as log:
        subprocess.run([python, "scripts/screen_bowtie.py", *map(str, pending), "--threads", "4"], cwd=root,
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
write_json(root / "experiments/final_followup/results/status.json", {"at": utc_now(),
    "stage": "followup_experiments_complete", "seed": 20261013, "new_runs": 18, "anchor_prefilter_tasks_complete": True})
