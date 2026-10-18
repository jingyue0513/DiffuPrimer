#!/usr/bin/env python3
"""Run recorded computational work; finish as a pilot, never as a validated paper."""
import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extended", action="store_true", help="Run all condition and completion analyses after primary results")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    python = str(root / ".venv/bin/python")
    env = dict(os.environ, PYTHONPATH=str(root / "src"), OMP_NUM_THREADS="8", MKL_NUM_THREADS="8")
    log_dir = root / "results/study_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    state = {"stage": "starting", "started_at": utc_now(), "completed_commands": [], "failed_commands": [],
             "manuscript_status": "Independent qPCR experiments pending; computational completion is not paper completion"}
    if args.resume:
        state = json.loads((root / "results/study_status.json").read_text())
    state.setdefault("execution_history", []).append({"at": utc_now(), "workers": args.workers, "resumed": args.resume,
        "threads_per_compute_process": 8, "deterministic_span_cache": True,
        "cache_wrapper_sha256": sha256(root / "scripts/evaluate_cached_spans.py"),
        "indexed_screen_wrapper_sha256": sha256(root / "scripts/screen_specificity_accelerated.py"),
        "indexed_pairing_sha256": sha256(root / "scripts/amplicon_pairing_fast.py"),
        "timing_policy": "Mixed fresh/cached per-gene timing is excluded from runtime comparisons"})
    previous = {row["label"]: row for row in state["completed_commands"]}

    def already_completed(label, command):
        if label not in previous:
            return False
        if previous[label]["command"] != command:
            raise ValueError("A completed command changed during resume: " + label)
        return True

    def update(stage):
        state["stage"] = stage
        state["updated_at"] = utc_now()
        write_json(root / "results/study_status.json", state)
        print(json.dumps({"stage": stage, "at": state["updated_at"]}), flush=True)

    def run(label, command):
        path = log_dir / (label + ".log")
        actual = [python, *command]
        if command[0] == "scripts/evaluate.py":
            actual = [python, "scripts/evaluate_cached_spans.py", "--evaluation-script", command[0], *command[1:]]
        elif command[0] == "scripts/screen_specificity.py":
            actual = [python, "scripts/screen_specificity_accelerated.py", *command[1:]]
        with path.open("w") as handle:
            result = subprocess.run(actual, cwd=root, env=env, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{label} failed (exit {result.returncode}); see {path}")
        return {"label": label, "command": command, "actual_command": actual, "log": str(path), "completed_at": utc_now()}

    def sequential(label, command):
        if already_completed(label, command):
            return
        update(label)
        state["completed_commands"].append(run(label, command))
        update(label + "_complete")

    def grouped(stage, jobs):
        update(stage)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run, label, command): label for label, command in jobs
                if not already_completed(label, command)}
            for future in as_completed(futures):
                label = futures[future]
                state["completed_commands"].append(future.result())
                update(stage + ":" + label + "_complete")

    plan = json.loads((root / "configs/analysis_plan.json").read_text())
    seeds = plan["seeds"]
    primary_dirs = []
    try:
        sequential("freeze_partial_dataset", ["scripts/prepare_data.py", "--freeze-partial", "--include-teacher"])
        sequential("verify_tests", ["-m", "pytest", "-q"])
        sequential("freeze_protocol", ["scripts/freeze_protocol.py"])
        train_jobs = [(f"train_{mode}_{seed}", ["scripts/train.py", "--mode", mode, "--seed", str(seed),
            "--dataset", "data/processed/dataset_partial_pilot.jsonl", "--allow-weak", "--partial-source"])
            for seed in seeds for mode in ["direct", "diffusion"]]
        grouped("training", train_jobs)

        def evaluation_jobs(split, bands, repairs=(None,), ablation=None):
            jobs, directories = [], []
            for amin, amax in bands:
                for repair in repairs:
                    extra = ["--repair", repair] if repair else []
                    suffix = "_fix_" + repair if repair else ""
                    for seed in seeds:
                        for mode in ["direct", "diffusion", "rule_sampler"]:
                            name = f"{mode}_{split}_{seed}_{amin}_{amax}{suffix}"
                            command = ["scripts/evaluate.py", "--split", split, "--seed", str(seed),
                                       "--amplicon-min", str(amin), "--amplicon-max", str(amax), *extra]
                            if mode == "rule_sampler":
                                command += ["--rule-sampler"]
                                stage = "pilot"
                            else:
                                command += ["--checkpoint", f"results/checkpoints/{mode}_weak_partial_pilot_seed{seed}.pt", "--oligo-guidance"]
                                stage = "weak_partial_pilot"
                            descriptor = f"{mode}_{stage}_{split}_seed{seed}_amp{amin}-{amax}{suffix}_oligo_guided"
                            directory = "results/evaluations/" + descriptor
                            jobs.append((name, command))
                            directories.append(directory)
                    seed = seeds[0]
                    name = f"primer3_{split}_{amin}_{amax}{suffix}"
                    command = ["scripts/evaluate.py", "--split", split, "--seed", str(seed),
                               "--amplicon-min", str(amin), "--amplicon-max", str(amax), *extra]
                    directory = f"results/evaluations/primer3_pilot_{split}_seed{seed}_amp{amin}-{amax}{suffix}"
                    jobs.append((name, command))
                    directories.append(directory)
            return jobs, directories

        validation_jobs, _ = evaluation_jobs("validation", [(80, 200)])
        grouped("validation_generation", validation_jobs)
        # No validation-based tuning is performed by this runner. Any changes to
        # frozen code/config must be documented before the final test is opened.
        sequential("verify_protocol_before_test", ["scripts/freeze_protocol.py"])
        test_jobs, primary_dirs = evaluation_jobs("test", [(80, 200)])
        grouped("primary_test_generation", test_jobs)
        screen_jobs = [("screen_" + Path(directory).name, ["scripts/screen_specificity.py", directory, "--genome"]) for directory in primary_dirs]
        grouped("primary_specificity", screen_jobs)
        sequential("aggregate_primary", ["scripts/aggregate.py", *primary_dirs, "--output", "results/primary_aggregate.json"])
        sequential("build_primary_paper", ["scripts/build_paper.py", "--aggregate", "results/primary_aggregate.json"])
        if args.extended:
            jobs, extra_dirs = evaluation_jobs("test", plan["condition_bands"])
            grouped("condition_bands", jobs)
            repair_jobs, repair_dirs = evaluation_jobs("test", [(80, 200)], repairs=("forward", "reverse"))
            grouped("fixed_primer_completion", repair_jobs)
            seed = seeds[0]
            ablations = [("geometry_only", []), ("raw", ["--raw"]), ("one_step", ["--oligo-guidance", "--steps", "1"])]
            ablation_jobs = []
            for label, extra in ablations:
                ablation_jobs.append(("diffusion_ablation_" + label, ["scripts/evaluate.py", "--split", "test", "--seed", str(seed),
                    "--checkpoint", f"results/checkpoints/diffusion_weak_partial_pilot_seed{seed}.pt", *extra]))
            grouped("decoder_ablations", ablation_jobs)
            state["extended_evaluation_directories"] = extra_dirs + repair_dirs
        state["primary_evaluation_directories"] = primary_dirs
        update("computational_pilot_complete_independent_qpcr_pending")
    except Exception as exc:
        state["failed_commands"].append({"error": str(exc), "at": utc_now()})
        update("pipeline_error")
        raise


if __name__ == "__main__":
    main()
