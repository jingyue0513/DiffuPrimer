#!/usr/bin/env python3
"""Complete the memory-bounded replacement and resume the original study."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import sha256,utc_now,write_json


root=Path(__file__).resolve().parents[1]
python=str(root/'.venv/bin/python');env=dict(os.environ,PYTHONPATH=str(root/'src'))
directory=root/'results/evaluations/direct_weak_partial_pilot_test_seed20261015_amp80-200_oligo_guided'
archive=root/'results/checks/truncated_blast_repair'
guard=root/'results/checks/blast_batching/check.json'
while not guard.exists():time.sleep(5)
check=json.loads(guard.read_text())
if not check['exact_qualifying_hit_set_equality'] or check['script_sha256']!=sha256(root/'scripts/blast_query_batches.py'):
    raise ValueError('Batching equivalence check missing or stale')
raw=directory/'blast_genome.tsv'
if not (directory/'blast_genome.completed.json').exists() and raw.exists():
    raw.rename(archive/'blast_genome.monolithic_rerun_interrupted.tsv')
write_json(archive/'batch_execution.json',{'at':utc_now(),'equivalence_check_sha256':sha256(guard),
    'policy':'Same independent full-database queries and scoring; one query per subprocess bounds multiquery HSP accumulation',
    'script_sha256':sha256(root/'scripts/blast_query_batches.py')})
with (archive/'query_batched_rerun.log').open('a') as handle:
    subprocess.run([python,'scripts/screen_specificity_accelerated.py',str(directory),'--genome'],cwd=root,env=env,
        stdout=handle,stderr=subprocess.STDOUT,check=True)
remaining=[root/'results/evaluations'/name for name in ['rule_sampler_pilot_test_seed20261015_amp80-200_oligo_guided','primer3_pilot_test_seed20261013_amp80-200']]
while not all((p/'specificity_status.json').exists() for p in remaining):time.sleep(5)
write_json(archive/'completed.json',{'at':utc_now(),'new_genome_sha256':sha256(directory/'blast_genome.tsv'),
    'new_completion_marker_sha256':sha256(directory/'blast_genome.completed.json'),'query_batching_check_sha256':sha256(guard),
    'specificity_status_sha256':sha256(directory/'specificity_status.json'),'original_candidates_sha256':sha256(directory/'candidates.jsonl')})
with (root/'results/study.log').open('a') as handle:
    subprocess.run([python,'-u','scripts/run_study.py','--extended','--resume','--workers','6'],cwd=root,env=env,
        stdout=handle,stderr=subprocess.STDOUT,check=True)
