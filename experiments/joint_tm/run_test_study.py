#!/usr/bin/env python3
"""Matched test, condition-control, completion and step experiments."""
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
if json.loads((experiment / "results/status.json").read_text())["stage"] != "validation_complete":
    raise ValueError("Full validation must complete before test study")
jobs, directories = [], []
for seed in [20261013, 20261014, 20261015]:
    for mode in ["direct", "diffusion", "rule_sampler"]:
        base = [python, str(experiment / "scripts/evaluate.py"), "--split", "test", "--seed", str(seed), "--joint-tm-guidance"]
        stage = "pilot" if mode == "rule_sampler" else "weak_partial_pilot"
        if mode == "rule_sampler":
            base += ["--rule-sampler"]
        else:
            base += ["--checkpoint", str(root / f"results/checkpoints/{mode}_weak_partial_pilot_seed{seed}.pt"), "--oligo-guidance"]
        conditions = [(80, 200, None, None)] + [(a, b, None, None) for a, b in [(80, 120), (121, 160), (161, 200)]]
        conditions += [(80, 200, role, None) for role in ["forward", "reverse"]]
        if mode == "diffusion":
            conditions += [(80, 200, None, steps) for steps in [1, 4, 8]]
        for lower, upper, role, steps in conditions:
            command = base + ["--amplicon-min", str(lower), "--amplicon-max", str(upper)]
            suffix = f"{mode}_{stage}_test_seed{seed}_amp{lower}-{upper}"
            if role:
                command += ["--repair", role]
                suffix += "_fix_" + role
            suffix += "_oligo_guided"
            if steps:
                command += ["--steps", str(steps)]
                suffix += "_steps" + str(steps)
            suffix += "_joint_tm"
            directory = experiment / "results/evaluations" / suffix
            label = f"{mode}_{seed}_{lower}_{upper}_{role or 'free'}_{steps or 16}"
            jobs.append((label, command, directory))
            directories.append(str(directory))
planned = {"frozen_at": utc_now(), "stage": "post_validation_followup", "test_genes": 75,
    "seeds": [20261013, 20261014, 20261015], "candidate_budget": 20, "runs": len(jobs),
    "implemented_protocol_sha256": sha256(experiment / "implemented_protocol.json"),
    "validation_comparison_sha256": sha256(root / "results/validation_experiment_comparison.json"),
    "selection": "Joint Tm guidance selected from full validation; training and checkpoints unchanged",
    "comparators": "Identical joint guidance for diffusion, direct and rule; frozen Primer3 runs reused",
    "conditions": "Three disjoint amplicon bands; fixed forward and reverse completion; diffusion 1/4/8/16 steps",
    "primary_results_policy": "Retain original single-oligo-guided study separately; this is a follow-up study",
    "commands": [command for _, command, _ in jobs]}
protocol_path = experiment / "results/test_protocol.json"
if protocol_path.exists():
    previous = json.loads(protocol_path.read_text())
    if previous["commands"] != planned["commands"]:
        raise ValueError("Test protocol changed")
else:
    write_json(protocol_path, planned)
state = {"stage": "test_sampling", "started_at": utc_now(), "completed": [], "directories": directories}
status_path = experiment / "results/test_status.json"
write_json(status_path, state)
cache_check = json.loads((root / "results/checks/span_cache/check.json").read_text())
if not cache_check["exact_output_equality"] or cache_check["cache_sha256"] != sha256(root / "scripts/span_cache.py"):
    raise ValueError("Deterministic span-cache check missing or stale")
write_json(experiment / "results/cache_execution.json", {"started_at": utc_now(),
    "original_test_protocol_sha256": sha256(protocol_path), "wrapper_sha256": sha256(root / "scripts/evaluate_cached_spans.py"),
    "cache_sha256": sha256(root / "scripts/span_cache.py"), "equivalence_check_sha256": sha256(root / "results/checks/span_cache/check.json"),
    "policy": "Identical evaluation source, conditions, checkpoints, proposals and RNG. Disk cache only reuses deterministic span arrays.",
    "timing_policy": "Mixed fresh/cached per-gene timings are excluded from runtime comparisons; use a separate common-cache benchmark."})


def run(label, command, directory):
    with (experiment / "results" / (label + ".log")).open("a") as handle:
        actual = [python, str(root / "scripts/evaluate_cached_spans.py"), "--evaluation-script", command[1], *command[2:]]
        subprocess.run(actual, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
    if not (directory / "completion.json").exists():
        raise ValueError("Expected complete output missing: " + str(directory))
    return label


try:
    # Submit broad-band runs first so primary endpoint can be inspected promptly.
    jobs.sort(key=lambda job: not ("80_200_free_16" in job[0]))
    with ThreadPoolExecutor(max_workers=3) as executor:
        for future in as_completed([executor.submit(run, *job) for job in jobs]):
            state["completed"].append(future.result())
            state["updated_at"] = utc_now()
            write_json(status_path, state)
            print(json.dumps({"completed": state["completed"][-1], "at": state["updated_at"]}), flush=True)
    state["stage"] = "test_sampling_complete_specificity_pending"
    write_json(status_path, state)
except Exception as exc:
    state["stage"], state["error"] = "test_sampling_error", str(exc)
    write_json(status_path, state)
    raise
