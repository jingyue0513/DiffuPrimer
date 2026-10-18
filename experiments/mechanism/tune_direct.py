#!/usr/bin/env python3
"""Validation-only selection of full-support non-iterative score decoders.

Supplement to the frozen main study. Both trained networks receive the same
configuration search. No test screening is read until configurations are locked.
"""
import argparse
import dataclasses
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from run import (SEED,CHECKPOINT_SEED,Conditions,ConditionalDenoiser,torch,np,
                 guidance,install,read_jsonl,write_json,write_jsonl,sha256,
                 atomic_json,make_candidates,stable_seed,assess_pair,oligo_properties,summarize_candidates)
from controls import first_scores,all_support,weighted_order,pairs_from_order

CONFIGS={'topk':None,'temperature_0.5':0.5,'temperature_1':1.,'temperature_2':2.,'temperature_4':4.}


def screen(directories,out):
    with (ROOT/'results/bowtie_audit/mechanism_writer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        with (out/'screen.log').open('a') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/screen_bowtie.py'),*[str(p) for p in directories],
                            '--threads','8'],env=dict(os.environ,PYTHONPATH=str(ROOT/'src')),
                           cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)


def counts(directory,expected):
    completion=json.loads((directory/'completion.json').read_text())
    status=json.loads((directory/'specificity_bowtie_status.json').read_text())
    assert completion['all_expected_targets_evaluated'] and status['status']=='completed'
    assert completion['candidates_sha256']==status['candidate_sha256']==sha256(directory/'candidates.jsonl')
    assert status['output_sha256']==sha256(directory/'specificity_bowtie.jsonl')
    scores={r['candidate_id']:r for r in read_jsonl(directory/'specificity_bowtie.jsonl')}
    kept={str(g):set() for g in expected}
    rows=read_jsonl(directory/'candidates.jsonl')
    assert {r['candidate_id'] for r in rows}==set(scores)
    for r in rows:
        s=scores[r['candidate_id']]
        if r['assessment']['computational_feasible'] and s['gene_specificity_screen_pass'] and not s['genomic_DNA_product_detected']:
            kept[str(r['gene_id'])].add((r['forward_primer'],r['reverse_primer']))
    return {g:len(v) for g,v in kept.items()}


def generate(cohort,targets,configurations,models,conditions,out,protocol_path,selection_hash=None):
    gene_dir=out/'genes'/cohort;gene_dir.mkdir(parents=True,exist_ok=True)
    for number,target in enumerate(targets):
        path=gene_dir/f"{target['gene_id']}.json"
        if path.exists():
            old=json.loads(path.read_text())
            assert old['protocol_sha256']==sha256(protocol_path) and old.get('selection_sha256')==selection_hash
            continue
        template=target['template_sequence']
        span=guidance.feasible_oligo_spans(template,conditions)
        temperatures=guidance.compatible_span_temperatures(template,conditions,span)
        record={'protocol_sha256':sha256(protocol_path),'selection_sha256':selection_hash,'methods':{}}
        for mode in ['direct','diffusion']:
            begin=time.process_time()
            logits=first_scores(models[mode],template,conditions)
            coordinates,scores=all_support(logits,conditions,span[1],temperatures)
            common_cost=time.process_time()-begin
            for key in configurations[mode]:
                name=f'{mode}_{key}'
                oligo_properties.cache_clear()
                begin=time.process_time()
                order=np.argsort(-scores,kind='stable') if CONFIGS[key] is None else weighted_order(scores,stable_seed(cohort,str(target['gene_id'])),len(scores),CONFIGS[key])
                pairs=pairs_from_order(coordinates,order,template,20)
                assessments=[assess_pair(template,p,conditions) for p in pairs]
                elapsed=time.process_time()-begin+common_cost
                candidates=make_candidates(f'{cohort}_{name}',target,pairs,assessments)
                record['methods'][name]={'candidates':candidates,'ledger':{'gene_id':target['gene_id'],
                    **summarize_candidates(candidates,20),'generation_cpu_seconds':elapsed,
                    'support_size':len(scores),'common_score_support_cpu_seconds':common_cost}}
        atomic_json(path,record)
        atomic_json(out/'status.json',{'stage':'generation','cohort':cohort,'completed':number+1,'total':len(targets)})
        print(json.dumps({'cohort':cohort,'completed':number+1,'total':len(targets)}),flush=True)
    rows=[json.loads((gene_dir/f"{t['gene_id']}.json").read_text()) for t in targets]
    directories={}
    for mode in configurations:
        for key in configurations[mode]:
            name=f'{mode}_{key}'
            dest=out/'evaluations'/f'{cohort}_{name}';dest.mkdir(parents=True,exist_ok=True)
            write_json(dest/'manifest.json',{'method':name,'seed':SEED,'cohort':cohort,'config':key,
                       'protocol_sha256':sha256(protocol_path),'selection_sha256':selection_hash})
            write_jsonl(dest/'candidates.jsonl',[c for r in rows for c in r['methods'][name]['candidates']])
            write_jsonl(dest/'per_gene.jsonl',[r['methods'][name]['ledger'] for r in rows])
            write_json(dest/'completion.json',{'all_expected_targets_evaluated':True,'targets':len(targets),
                       'candidates_sha256':sha256(dest/'candidates.jsonl'),'per_gene_sha256':sha256(dest/'per_gene.jsonl')})
            directories[mode,key]=dest
    return directories


def main(out):
    out.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(4);torch.set_num_interop_threads(1)
    torch.manual_seed(SEED);torch.use_deterministic_algorithms(True)
    install()
    conditions=Conditions(**json.loads((ROOT/'configs/pilot.json').read_text())['conditions'])
    pilot=read_jsonl(ROOT/'data/processed/targets_pilot.jsonl')
    cohorts={'validation':[r for r in pilot if r['split']=='validation'],
             'test':[r for r in pilot if r['split']=='test'],
             'external_test':read_jsonl(ROOT/'data/processed/attention_external_targets.jsonl')}
    for k in cohorts:cohorts[k].sort(key=lambda r:str(r['gene_id']))
    models={};checkpoints={}
    for mode in ['direct','diffusion']:
        path=ROOT/f'results/checkpoints/{mode}_weak_partial_pilot_seed{CHECKPOINT_SEED}.pt'
        state=torch.load(path,map_location='cpu',weights_only=True)
        model=ConditionalDenoiser(**state['model_spec']);model.load_state_dict(state['state_dict']);models[mode]=model.eval()
        checkpoints[mode]=sha256(path)
    protocol={'seed':SEED,'checkpoint_hashes':checkpoints,'historical_training_seed':CHECKPOINT_SEED,
              'configurations':CONFIGS,'budget':20,'cohorts':{k:[r['gene_id'] for r in v] for k,v in cohorts.items()},
              'selection':'Per-network validation mean retained pairs; ties >=2 retained-pair coverage, then listed configuration order. No test metrics accessed.',
              'code_sha256':{p.name:sha256(p) for p in [Path(__file__),HERE/'controls.py',HERE/'run.py']},
              'target_sha256':{p.name:sha256(p) for p in [ROOT/'data/processed/targets_pilot.jsonl',ROOT/'data/processed/attention_external_targets.jsonl']}}
    protocol_path=out/'protocol.json'
    if protocol_path.exists():assert json.loads(protocol_path.read_text())==protocol
    else:atomic_json(protocol_path,protocol)
    selected_path=out/'selection.json'
    if not selected_path.exists():
        directories=generate('validation',cohorts['validation'],{m:list(CONFIGS) for m in models},models,conditions,out,protocol_path)
        atomic_json(out/'status.json',{'stage':'validation_screen'})
        screen(list(directories.values()),out)
        evidence=[];selection={}
        for mode in models:
            options=[]
            for index,key in enumerate(CONFIGS):
                vals=counts(directories[mode,key],protocol['cohorts']['validation'])
                mean=float(np.mean(list(vals.values())));two=float(np.mean([v>=2 for v in vals.values()]))
                evidence.append({'model':mode,'config':key,'retained_pairs':mean,'two_retained_coverage':two,
                                 'screen_sha256':sha256(directories[mode,key]/'specificity_bowtie.jsonl')})
                options.append((mean,two,-index,key))
            selection[mode]=max(options)[-1]
        atomic_json(selected_path,{'selected':selection,'validation_results':evidence,
                                  'protocol_sha256':sha256(protocol_path),'test_outcomes_consulted':False})
    selected=json.loads(selected_path.read_text())['selected']
    all_dirs={}
    for cohort in ['test','external_test']:
        ds=generate(cohort,cohorts[cohort],{m:[selected[m]] for m in models},models,conditions,out,protocol_path,sha256(selected_path))
        all_dirs.update({(cohort,m,k):d for (m,k),d in ds.items()})
    atomic_json(out/'status.json',{'stage':'test_screen','selected':selected})
    screen(list(all_dirs.values()),out)
    results=[]
    for (cohort,mode,key),directory in all_dirs.items():
        vals=counts(directory,protocol['cohorts'][cohort])
        results.append({'cohort':cohort,'model':mode,'config':key,'retained_pairs':float(np.mean(list(vals.values()))),
                        'one_retained':float(np.mean([v>=1 for v in vals.values()])),
                        'two_retained':float(np.mean([v>=2 for v in vals.values()])),
                        'per_gene_retained':vals,'directory':str(directory),
                        'screen_sha256':sha256(directory/'specificity_bowtie.jsonl')})
    atomic_json(out/'summary.json',{'selection_sha256':sha256(selected_path),'results':results})
    atomic_json(out/'status.json',{'stage':'complete','selected':selected})
    print(json.dumps({'selected':selected,'results':[{k:v for k,v in r.items() if k!='per_gene_retained'} for r in results]}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'tuned_seed42');args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:main(args.output.resolve())
        except Exception as exc:
            atomic_json(args.output/'status.json',{'stage':'failed','error':repr(exc)})
            raise
