#!/usr/bin/env python3
"""Preserve live alignment searches, then resume with equivalent fast pairing."""
import json
import os
import shlex
import signal
import subprocess
import time
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
check = json.loads((root / "results/checks/fast_pairing/check.json").read_text())
assert check["exact_ordered_output_equality"] and check["real_reference_candidate_reference_comparisons"] > 0
state_path = root / "results/study_status.json"
state = json.loads(state_path.read_text())


def processes():
    rows = []
    for line in subprocess.check_output(["ps", "-eo", "pid,ppid,args"], text=True).splitlines()[1:]:
        pid, parent, *command = line.split()
        rows.append((int(pid), int(parent), " ".join(command)))
    return rows


rows = processes()
parents = [pid for pid, _, command in rows if command.startswith(".venv/bin/python -u scripts/run_study.py")
           and Path(f"/proc/{pid}/cwd").resolve() == root]
if len(parents) != 1:
    raise ValueError("Cannot identify one original driver")
owner = parents[0]
screens = [(pid, shlex.split(command)[2]) for pid, parent, command in rows
    if parent == owner and " scripts/screen_specificity.py " in command]
if not screens:
    raise ValueError("No running original screening processes")
template_directory = root / "results/evaluations/diffusion_weak_partial_pilot_test_seed20261013_amp80-200_oligo_guided"
template_status = json.loads((template_directory / "specificity_status.json").read_text())
write_json(root / "results/study_status_before_indexed_pairing.json", state)
os.kill(owner, signal.SIGTERM)
write_json(root / "results/indexed_pairing_transition.json", {"at": utc_now(), "old_driver": owner,
    "screen_owners": screens, "policy": "Let original BLAST processes finish; reuse their byte-verified raw alignments. Replace only mathematically equivalent product enumeration."})
pending = dict(screens)
while pending:
    rows = processes()
    alive = {pid for pid, _, _ in rows}
    for pid, raw in list(pending.items()):
        if any(parent == pid and "blastn " in command for _, parent, command in rows):
            continue
        directory = root / raw
        files = [directory / f"blast_{name}.tsv" for name in ["transcriptome", "genome"]]
        if not all(path.exists() for path in files):
            raise ValueError("Inactive original search has missing raw files")
        if pid not in alive and not (directory / "specificity_status.json").exists():
            raise ValueError("Original screen exited without completed status")
        if pid in alive:
            # The original script reaches parsing/pairing only after check=True succeeds.
            # Both raw files exist and no blast child remains, so genome search succeeded.
            os.kill(pid, signal.SIGTERM)
        query = directory / "specificity_queries.fasta"
        for original_command in template_status["commands"]:
            command = [part.replace(str(template_directory), str(directory)) for part in original_command]
            output = Path(command[command.index("-out")+1])
            database = Path(command[command.index("-db")+1])
            write_json(output.with_suffix(".completed.json"), {"command": command, "query_sha256": sha256(query),
                "database_manifest_sha256": sha256(database.with_suffix(".manifest.json")), "output_sha256": sha256(output),
                "completed_at": utc_now(), "completion_evidence": "Original frozen script was alive beyond checked genome subprocess with no alignment child, or had completed final status",
                "original_screen_process": pid})
        print(json.dumps({"verified_original_alignments": raw, "at": utc_now()}), flush=True)
        pending.pop(pid)
    if pending:
        time.sleep(3)
with (root / "results/study.log").open("a") as handle:
    subprocess.run([str(root / ".venv/bin/python"), "-u", "scripts/run_study.py", "--extended", "--resume", "--workers", "6"],
        cwd=root, env=dict(os.environ, PYTHONPATH=str(root / "src")), check=True, stdout=handle, stderr=subprocess.STDOUT)
