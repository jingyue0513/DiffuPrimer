#!/usr/bin/env python3
"""Run a fixed, matched held-out duration comparison and full reference audit."""
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
experiment = root / "experiments/epoch100"
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
state_path = experiment / "results/followup_status.json"
plan_path = experiment / "three_seed_evaluation_protocol.json"
directories = {split: [str(experiment / f"results/evaluations/{method}_weak_partial_pilot_{split}_seed{seed}_amp80-200_oligo_guided")
    for seed in [20261013, 20261014, 20261015] for method in ["direct", "diffusion"]] for split in ["validation", "test"]}
plan = {"directories": directories, "seeds": [20261013, 20261014, 20261015], "budget": 20,
        "conditions": json.loads((root / "configs/pilot.json").read_text())["conditions"],
        "sampling_temperature": 1.0, "guidance": "Original single-oligo guidance in both duration groups",
        "checkpoints": "Best validation-objective checkpoints from the declared 100-epoch trajectories",
        "primary_results_retained": True, "all_test_genes": 75, "new_validation_sampling_runs": 4,
        "new_test_sampling_runs": 6, "split_sha256": sha256(root / "data/splits/pilot.jsonl"),
        "training_protocol_sha256": sha256(experiment / "three_seed_protocol.json"),
        "evaluate_source_sha256": sha256(experiment / "scripts/evaluate.py"),
        "test_policy": "Run all six test evaluations regardless of the validation-duration effect; fixed model and decoding settings; no test-based retuning",
        "audit_policy": "Reuse complete reference alignments; wait for the primary audit before becoming the only alignment-cache writer"}
if plan_path.exists():
    if json.loads(plan_path.read_text())["plan"] != plan:
        raise ValueError("Duration follow-up plan changed")
else:
    write_json(plan_path, {"frozen_at": utc_now(), "plan": plan})


def state(stage, **extra):
    write_json(state_path, {"stage": stage, "at": utc_now(), **extra})
    print(json.dumps({"stage": stage, **extra}), flush=True)


def wait(path, required, errors):
    while True:
        stage = json.loads(path.read_text())["stage"] if path.exists() else "not_started"
        if stage == required:
            return
        if stage in errors:
            raise ValueError(f"Upstream experiment failed: {path}/{stage}")
        time.sleep(10)


try:
    state("waiting_for_three_seed_validation")
    wait(experiment / "results/three_seed_validation_status.json", "validation_comparison_complete", ["duration_validation_error"])
    state("test_sampling")

    def run(job):
        method, seed = job
        command = [python, str(root / "scripts/evaluate_cached_spans.py"), "--evaluation-script",
                   str(experiment / "scripts/evaluate.py"), "--split", "test", "--seed", str(seed),
                   "--oligo-guidance", "--checkpoint", str(experiment / f"results/checkpoints/{method}_weak_partial_pilot_seed{seed}.pt")]
        with (experiment / f"results/training_logs/test_{method}_{seed}.log").open("a") as handle:
            subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)

    jobs = [(method, seed) for seed in plan["seeds"] for method in ["direct", "diffusion"]]
    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(run, jobs))
    state("test_sampling_complete_waiting_for_primary_audit")
    wait(root / "results/bowtie_audit/driver_status.json", "all_completed_experiments_audited", ["audit_error"])
    state("full_reference_audit")
    with (experiment / "results/followup_reference_audit.log").open("a") as handle:
        subprocess.run([python, "scripts/screen_bowtie.py", *directories["validation"], *directories["test"], "--threads", "4"],
                       cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    for split in ["validation", "test"]:
        with (experiment / f"results/followup_{split}_statistics.log").open("w") as handle:
            subprocess.run([python, "scripts/summarize_specificity_audit.py", *directories[split], "--split", split,
                            "--output", str(experiment / f"results/{split}_bowtie_summary.json")],
                           cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    with (experiment / "results/duration_reference_comparison.log").open("w") as handle:
        subprocess.run([python, "scripts/compare_duration_reference_results.py"], cwd=root, env=env,
                       check=True, stdout=handle, stderr=subprocess.STDOUT)
    state("duration_followup_analyzed", runs=12, protocol_sha256=sha256(plan_path), manuscript_work_paused=True)
except Exception as exc:
    state("duration_followup_error", error=str(exc))
    raise
