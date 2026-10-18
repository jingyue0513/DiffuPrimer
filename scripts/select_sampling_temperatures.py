#!/usr/bin/env python3
"""Select on complete validation screens, then evaluate all frozen test genes."""
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
directory = root / "results/temperature_tuning"
status_path = directory / "selection_status.json"
validation_path = root / "results/bowtie_audit/validation_RNA_summary.json"
python = str(root / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(root / "experiments/joint_tm/src"))


def state(stage, **extra):
    write_json(status_path, {"stage": stage, "at": utc_now(), **extra})
    print(json.dumps({"stage": stage, "at": utc_now(), **extra}), flush=True)


def temperature_group(temperature):
    return "joint_tm" if temperature == 1.0 else "joint_tm_temperature" + str(temperature)


try:
    state("waiting_for_complete_validation_screens")
    while True:
        sampling = json.loads((directory / "status.json").read_text())
        if sampling["stage"] == "validation_sampling_error":
            raise ValueError("Temperature validation generation failed")
        if validation_path.exists() and sampling["stage"] == "validation_sampling_complete_specificity_pending":
            summary = json.loads(validation_path.read_text())
            candidates = [r for r in summary["summaries"] if r["method"] in ["direct", "diffusion"] and
                r["study"] in [temperature_group(t) for t in [0.5, 1.0, 2.0, 4.0]] and
                r["budget"] == 20 and r["context"] == "amp80-200" and r["variant"] == "primary"]
            if len(candidates) == 8 and all(r["genes"] == 75 and r["seeds"] == [20261013, 20261014, 20261015] for r in candidates):
                break
        audit_status = root / "results/bowtie_audit/driver_status.json"
        if audit_status.exists() and json.loads(audit_status.read_text())["stage"] == "audit_error":
            raise ValueError("Independent validation screen failed")
        time.sleep(5)
    selected = {}
    for method in ["direct", "diffusion"]:
        own = [r for r in candidates if r["method"] == method]
        best = max(r["unique_RNA_pass_pairs"] for r in own)
        tied = [r for r in own if abs(r["unique_RNA_pass_pairs"] - best) <= 1e-12]
        ranked = []
        for row in tied:
            temperature = 1.0 if row["study"] == "joint_tm" else float(row["study"].split("temperature")[1])
            quality = row["best_RNA_pass_temperature_penalty"]
            ranked.append((quality if quality is not None else float("inf"), temperature, row))
        quality, temperature, row = min(ranked, key=lambda item: (item[0], item[1]))
        selected[method] = {"temperature": temperature, "validation_summary": row}
    selection_path = directory / "selection.json"
    planned = {"selected_at": utc_now(), "selected": selected, "validation_groups": candidates,
        "protocol_sha256": sha256(directory / "protocol.json"), "validation_summary_sha256": sha256(validation_path),
        "criterion": "Mean unique RNA-pass pairs across complete validation genes/seeds; ties by lower mean best RNA-pass temperature penalty, then lower sampling temperature",
        "test_outcomes_used_for_selection": False}
    if selection_path.exists():
        if json.loads(selection_path.read_text())["selected"] != selected:
            raise ValueError("Previously frozen temperature selection differs")
    else:
        write_json(directory / "selection_validation_summary.json", summary)
        write_json(selection_path, planned)
    state("selected_test_sampling", selected_temperatures={m: r["temperature"] for m,r in selected.items()})
    jobs, directories = [], []
    for method in ["direct", "diffusion"]:
        temperature = selected[method]["temperature"]
        label = {0.5: "05", 2.0: "2", 4.0: "4"}.get(temperature)
        experiment = root / ("experiments/temperature_" + label) if label else root / "experiments/joint_tm"
        for seed in [20261013, 20261014, 20261015]:
            output = experiment / f"results/evaluations/{method}_weak_partial_pilot_test_seed{seed}_amp80-200_oligo_guided_joint_tm"
            command = [python, str(root / "scripts/evaluate_cached_spans.py"), "--evaluation-script", str(experiment / "scripts/evaluate.py"),
                "--split", "test", "--seed", str(seed), "--joint-tm-guidance", "--oligo-guidance", "--temperature", str(temperature),
                "--checkpoint", str(root / f"results/checkpoints/{method}_weak_partial_pilot_seed{seed}.pt")]
            jobs.append((method, seed, temperature, output, command))
            directories.append(str(output))
    write_json(directory / "selected_test_protocol.json", {"frozen_at": utc_now(), "selection_sha256": sha256(selection_path),
        "commands": [job[-1] for job in jobs], "directories": directories, "temperature_1_policy": "Reuse owned joint-Tm test runs after their completion"})

    def run(job):
        method, seed, temperature, output, command = job
        if temperature == 1.0:
            while not (output / "completion.json").exists():
                if json.loads((root / "experiments/joint_tm/results/test_status.json").read_text())["stage"] == "test_sampling_error":
                    raise ValueError("Joint default test run failed")
                time.sleep(5)
        else:
            with (directory / f"selected_{method}_{seed}.log").open("a") as handle:
                subprocess.run(command, cwd=root, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT)
        return str(output)

    with ThreadPoolExecutor(max_workers=3) as executor:
        list(executor.map(run, jobs))
    state("selected_test_sampling_complete_specificity_pending", directories=directories)
    while not all((Path(path) / "specificity_bowtie_status.json").exists() for path in directories):
        if json.loads((root / "results/bowtie_audit/driver_status.json").read_text())["stage"] == "audit_error":
            raise ValueError("Selected test specificity failed")
        time.sleep(5)
    baselines = [str(root / f"experiments/joint_tm/results/evaluations/rule_sampler_pilot_test_seed{seed}_amp80-200_oligo_guided_joint_tm") for seed in [20261013, 20261014, 20261015]]
    baselines += [str(root / f"experiments/uniform_pairs/results/evaluations/uniform_pair_sampler_pilot_test_seed{seed}_amp80-200_oligo_guided_joint_tm") for seed in [20261013, 20261014, 20261015]]
    baselines.append(str(root / "results/evaluations/primer3_pilot_test_seed20261013_amp80-200"))
    while not all((Path(path) / "specificity_bowtie_status.json").exists() for path in baselines):
        time.sleep(5)
    with (directory / "selected_test_statistics.log").open("w") as handle:
        subprocess.run([python, "scripts/summarize_specificity_audit.py", *directories, *baselines,
            "--selection", str(selection_path), "--output", "results/temperature_tuning/selected_test_summary.json"],
            cwd=root, env=dict(env, PYTHONPATH=str(root / "src")), check=True, stdout=handle, stderr=subprocess.STDOUT)
    state("selected_test_analyzed", directories=directories)
except Exception as exc:
    state("selection_error", error=str(exc))
    raise
