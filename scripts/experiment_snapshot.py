#!/usr/bin/env python3
"""Snapshot actual completion ledgers for the running study."""
import json
from collections import Counter, defaultdict
from pathlib import Path

from diffuprimer.io import read_jsonl, utc_now, write_json, write_jsonl


root = Path(__file__).resolve().parents[1]
studies = {"original": (root / "results/evaluations", 81),
    "joint_tm": (root / "experiments/joint_tm/results/evaluations", 72),
    "epoch100": (root / "experiments/epoch100/results/evaluations", 12 if (root / "experiments/epoch100/three_seed_evaluation_protocol.json").exists() else 6 if (root / "experiments/epoch100/three_seed_protocol.json").exists() else 2)}
studies["uniform_pairs"] = (root / "experiments/uniform_pairs/results/evaluations", 18)
for experiment in sorted((root / "experiments").glob("temperature_*")):
    planned = 6
    selected_path = root / "results/temperature_tuning/selected_test_protocol.json"
    if selected_path.exists():
        planned += sum(str(experiment) + "/" in path for path in json.loads(selected_path.read_text())["directories"])
    condition_path = root / "results/temperature_tuning/selected_conditions_protocol.json"
    if condition_path.exists():
        planned += sum(str(experiment) + "/" in path for path in json.loads(condition_path.read_text())["directories"])
    studies[experiment.name] = (experiment / "results/evaluations", planned)
rows, summaries = [], []
for study, (base, planned) in studies.items():
    counts = Counter()
    for path in sorted(base.glob("*/manifest.json")):
        directory = path.parent
        manifest = json.loads(path.read_text())
        if manifest["stage"] == "calibration":
            continue
        arguments = manifest["arguments"]
        records = read_jsonl(directory / "per_gene.jsonl") if (directory / "per_gene.jsonl").exists() else []
        complete = (directory / "completion.json").exists()
        original_screen = (directory / "specificity_status.json").exists()
        audit_screen = (directory / "specificity_bowtie_status.json").exists()
        rows.append({"study": study, "method": manifest["method"], "split": arguments["split"],
            "seed": arguments["seed"], "amplicon_min": arguments["amplicon_min"], "amplicon_max": arguments["amplicon_max"],
            "fixed_side": arguments["repair"], "steps_override": arguments["steps"], "raw": arguments["raw"],
            "budget": arguments["budget"], "evaluated_genes": len(records), "expected_genes": len(manifest["expected_targets"]),
            "sampling_complete": complete, "original_specificity_complete": original_screen,
            "Bowtie_specificity_complete": audit_screen, "directory": str(directory)})
        counts["started_runs"] += 1
        counts["sampling_complete_runs"] += int(complete)
        counts["original_specificity_complete_runs"] += int(original_screen)
        counts["Bowtie_specificity_complete_runs"] += int(audit_screen)
    summaries.append({"study": study, "planned_runs": planned, **dict(counts)})
states = {}
for name, relative in {"core": "results/study_status.json", "joint_tm_test": "experiments/joint_tm/results/test_status.json",
    "step_sweep": "results/step_sweep/status.json", "Bowtie_audit": "results/bowtie_audit/driver_status.json",
    "result_analysis": "results/result_checks_status.json", "temperature_sampling": "results/temperature_tuning/status.json",
    "temperature_selection": "results/temperature_tuning/selection_status.json",
    "uniform_pairs": "experiments/uniform_pairs/results/status.json",
    "selected_conditions": "results/temperature_tuning/selected_conditions_status.json",
    "duration_training": "experiments/epoch100/results/three_seed_status.json",
    "duration_validation": "experiments/epoch100/results/three_seed_validation_status.json",
    "duration_reference_followup": "experiments/epoch100/results/followup_status.json"}.items():
    path = root / relative
    if path.exists():
        states[name] = json.loads(path.read_text())["stage"]
write_jsonl(root / "results/experiment_snapshot.runs.jsonl", rows)
snapshot = {"at": utc_now(), "summaries": summaries, "states": states,
    "manuscript_work_paused": (root / "results/paper_writing_paused.flag").exists(),
    "scope": "Original and joint-Tm studies, training-duration checks, uniform compatible-pair baseline, validation temperature selection and selected-temperature condition studies; reused evaluations counted once"}
write_json(root / "results/experiment_snapshot.json", snapshot)
print(json.dumps(snapshot, indent=2))
