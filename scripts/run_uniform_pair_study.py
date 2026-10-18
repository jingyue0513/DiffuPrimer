#!/usr/bin/env python3
"""Run the frozen compatible-pair baseline under all primary conditions."""
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
exp = root/'experiments/uniform_pairs'
python = str(root/'.venv/bin/python')
env = dict(os.environ,PYTHONPATH=str(exp/'src'))
protocol = json.loads((exp/'planned_protocol.json').read_text())
for relative,checksum in protocol['sha256'].items():
    if sha256(exp/relative) != checksum:
        raise ValueError('Uniform-pair implementation changed')
check = json.loads((exp/'results/check.json').read_text())
if not check['exact_pool_equality'] or check['source_sha256'] != sha256(exp/'src/diffuprimer/guidance.py'):
    raise ValueError('Uniform-pair sampler check missing or stale')
jobs=[]
for seed in [20261013,20261014,20261015]:
    for lower,upper,role in [(80,200,None),(80,120,None),(121,160,None),(161,200,None),(80,200,'forward'),(80,200,'reverse')]:
        command=[python,str(root/'scripts/evaluate_cached_spans.py'),'--evaluation-script',str(exp/'scripts/evaluate.py'),
            '--split','test','--seed',str(seed),'--joint-tm-guidance','--rule-sampler','--amplicon-min',str(lower),'--amplicon-max',str(upper)]
        suffix=f'uniform_pair_sampler_pilot_test_seed{seed}_amp{lower}-{upper}'
        if role:
            command += ['--repair',role]
            suffix += '_fix_'+role
        suffix += '_oligo_guided_joint_tm'
        jobs.append((suffix,command,exp/'results/evaluations'/suffix))
state={'stage':'test_sampling','at':utc_now(),'completed':[],'directories':[str(j[2]) for j in jobs]}
write_json(exp/'results/execution_protocol.json',{'at':utc_now(),'planned_protocol_sha256':sha256(exp/'planned_protocol.json'),
    'commands':[j[1] for j in jobs],'wrapper_sha256':sha256(root/'scripts/evaluate_cached_spans.py')})

def run(job):
    label,command,directory=job
    with (exp/'results'/(label+'.log')).open('a') as handle:
        subprocess.run(command,cwd=root,env=env,check=True,stdout=handle,stderr=subprocess.STDOUT)
    if not (directory/'completion.json').exists():
        raise ValueError('Missing complete output '+str(directory))
    return label

try:
    with ThreadPoolExecutor(max_workers=2) as executor:
        for future in as_completed([executor.submit(run,j) for j in jobs]):
            state['completed'].append(future.result());state['at']=utc_now()
            write_json(exp/'results/status.json',state)
            print(json.dumps({'completed':state['completed'][-1],'at':state['at']}),flush=True)
    state['stage']='test_sampling_complete_specificity_pending'
    write_json(exp/'results/status.json',state)
except Exception as exc:
    state.update(stage='test_sampling_error',error=str(exc))
    write_json(exp/'results/status.json',state)
    raise
