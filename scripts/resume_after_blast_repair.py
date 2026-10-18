#!/usr/bin/env python3
"""Recompute a truncated alignment and resume after preserved screens finish."""
import json
import os
import signal
import subprocess
import time
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root=Path(__file__).resolve().parents[1]
python=str(root/'.venv/bin/python')
env=dict(os.environ,PYTHONPATH=str(root/'src'))
state_path=root/'results/study_status.json'
prior_state=json.loads(state_path.read_text())
directory=root/'results/evaluations/direct_weak_partial_pilot_test_seed20261015_amp80-200_oligo_guided'
old_driver=916694
command=Path(f'/proc/{old_driver}/cmdline').read_bytes()
if b'scripts/run_study.py' not in command:
    raise ValueError('Expected original core driver differs')
os.kill(old_driver,signal.SIGTERM)
archive=root/'results/checks/truncated_blast_repair'
archive.mkdir(parents=True,exist_ok=True)
raw=directory/'blast_genome.tsv'
with raw.open('rb') as handle:
    handle.seek(max(0,raw.stat().st_size-1024));tail=handle.read()
if tail.endswith(b'\n'):
    raise ValueError('Expected truncated final line no longer present')
evidence={'at':utc_now(),'source':str(raw),'bytes':raw.stat().st_size,'raw_sha256':sha256(raw),
    'incomplete_last_line':tail.splitlines()[-1].decode('ascii'),'policy':'Preserve original truncated file; recompute complete genome search with unchanged query, database and BLAST flags. No records are deleted or completed by inference.'}
raw.rename(archive/'blast_genome.truncated.tsv')
(directory/'blast_genome.completed.json').rename(archive/'revoked_completion_marker.json')
write_json(archive/'repair.json',evidence)
write_json(state_path,{**prior_state,'stage':'primary_specificity_alignment_repair','updated_at':utc_now(),
    'reason':'Truncated genome file, rejected by frozen parser; unchanged candidates and screening parameters'})
with (archive/'rerun.log').open('a') as handle:
    subprocess.run([python,'scripts/screen_specificity_accelerated.py',str(directory),'--genome'],cwd=root,env=env,
        stdout=handle,stderr=subprocess.STDOUT,check=True)
remaining=[root/'results/evaluations'/name for name in ['rule_sampler_pilot_test_seed20261015_amp80-200_oligo_guided','primer3_pilot_test_seed20261013_amp80-200']]
while not all((p/'specificity_status.json').exists() for p in remaining):
    time.sleep(5)
write_json(archive/'completed.json',{'at':utc_now(),'new_genome_sha256':sha256(directory/'blast_genome.tsv'),
    'new_completion_marker_sha256':sha256(directory/'blast_genome.completed.json'),
    'specificity_status_sha256':sha256(directory/'specificity_status.json'),'original_candidates_sha256':sha256(directory/'candidates.jsonl')})
with (root/'results/study.log').open('a') as handle:
    subprocess.run([python,'-u','scripts/run_study.py','--extended','--resume','--workers','6'],cwd=root,env=env,
        stdout=handle,stderr=subprocess.STDOUT,check=True)
