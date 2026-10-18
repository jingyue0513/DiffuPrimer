#!/usr/bin/env python3
"""Reload orchestration after the current immutable screen finishes."""
import os
import signal
import subprocess
import time
from pathlib import Path

from diffuprimer.io import utc_now, write_json


root = Path(__file__).resolve().parents[1]
parent, screen = 919316, 919318
parent_cmd = Path(f"/proc/{parent}/cmdline").read_bytes()
screen_cmd = Path(f"/proc/{screen}/cmdline").read_bytes()
if b"scripts/run_bowtie_audit.py" not in parent_cmd or b"scripts/screen_bowtie.py" not in screen_cmd:
    raise ValueError("Expected audit parent/processor identity differs")
os.kill(parent,signal.SIGTERM)
write_json(root/"results/bowtie_audit/orchestration_reload.json",{
    "at":utc_now(),"stage":"waiting_for_preserved_screen","preserved_screen_pid":screen,
    "reason":"Include newly frozen uniform compatible-pair baseline; preserve all active alignments and ingestion"})
while Path(f"/proc/{screen}/cmdline").exists():
    command=Path(f"/proc/{screen}/cmdline").read_bytes()
    if command != screen_cmd:
        break
    time.sleep(5)
write_json(root/"results/bowtie_audit/orchestration_reload.json",{"at":utc_now(),"stage":"reloaded","preserved_screen_pid":screen})
with (root/"results/bowtie_audit_driver.log").open("a") as handle:
    subprocess.run([str(root/".venv/bin/python"),"-u","scripts/run_bowtie_audit.py"],cwd=root,
        env=dict(os.environ,PYTHONPATH=str(root/"src")),stdout=handle,stderr=subprocess.STDOUT,check=True)
