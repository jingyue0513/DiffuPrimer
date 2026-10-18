#!/usr/bin/env python3
"""Wait for recorded core runs, then screen and analyze every experiment."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import utc_now, write_json


root = Path(__file__).resolve().parents[1]
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "src"))
state_path = root / "results/result_checks_status.json"


def update(stage, **kwargs):
    state = {"stage": stage, "at": utc_now(), **kwargs}
    write_json(state_path, state)
    print(json.dumps(state), flush=True)


def run(label, arguments):
    update(label)
    with (root / "results/study_logs" / (label + ".log")).open("w") as handle:
        subprocess.run([python, *arguments], cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)


def driver_status():
    state = json.loads((root / "results/study_status.json").read_text())
    if state["stage"] == "pipeline_error":
        raise RuntimeError("Core driver failed; inspect study_status.json before continuing")
    return state


def wait_for_audit(directories):
    update("waiting_for_full_length_reference_audit", required_runs=len(directories))
    while not all((Path(path) / "specificity_bowtie_status.json").exists() for path in directories):
        audit = json.loads((root / "results/bowtie_audit/driver_status.json").read_text())
        if audit["stage"] == "audit_error":
            raise RuntimeError("Full-length reference audit failed")
        driver_status()
        time.sleep(5)


try:
    update("waiting_for_primary_results")
    aggregate_path = root / "results/primary_aggregate.json"
    while not aggregate_path.exists():
        driver_status()
        time.sleep(5)
    primary = [row["source"] for row in json.loads(aggregate_path.read_text())["runs"]]
    run("primary_seed_statistics", ["scripts/summarize_seeds.py", "--analysis-only"])
    wait_for_audit(primary)
    run("primary_full_length_reference_statistics", ["scripts/summarize_specificity_audit.py", *primary,
        "--output", "results/primary_bowtie_summary.json"])
    run("primary_secondary_statistics", ["scripts/summarize_secondary.py", *primary, "--screen", "bowtie", "--output", "results/primary_secondary_summary.json"])
    run("primary_failure_diagnostics", ["scripts/inspect_design_results.py", *primary, "--screen", "bowtie", "--output", "results/primary_design_diagnostics.json"])
    update("waiting_for_condition_and_completion_experiments")
    while driver_status()["stage"] != "computational_pilot_complete_independent_qpcr_pending":
        time.sleep(5)
    extended = driver_status()["extended_evaluation_directories"]
    ablation_directories = []
    for path in (root / "results/evaluations").glob("diffusion_weak_partial_pilot_test_seed20261013_amp80-200*"):
        manifest = json.loads((path / "manifest.json").read_text())
        args = manifest["arguments"]
        if args["raw"] or args["steps"] == 1 or (not args["oligo_guidance"] and not args["repair"]):
            ablation_directories.append(str(path))
    extra = extended + ablation_directories
    wait_for_audit(extra)
    run("complete_secondary_statistics", ["scripts/summarize_secondary.py", *primary, *extra, "--screen", "bowtie", "--output", "results/complete_secondary_summary.json"])
    run("complete_failure_diagnostics", ["scripts/inspect_design_results.py", *primary, *extra, "--screen", "bowtie", "--output", "results/complete_design_diagnostics.json"])
    update("core_experiments_analyzed", primary_runs=len(primary), extended_runs=len(extra),
           specificity_screen="all_strata_full_length_bowtie", manuscript_work_paused=True)
except Exception as exc:
    update("analysis_error", error=str(exc))
    raise
