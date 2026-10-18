#!/usr/bin/env python3
"""Execute the attention study with two generation workers and one reference screen."""
import os,sys,subprocess,time,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from diffuprimer.io import write_json,utc_now,sha256
BASE=ROOT/'experiments/attention';OUT=BASE/'results';PY=str(ROOT/'.venv_attention/bin/python');CPU=str(ROOT/'.venv/bin/python');env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
def status(stage,**extra):
 d=dict(stage=stage,at=utc_now(),seed=20261013,**extra);write_json(OUT/'status.json',d);print(json.dumps(d),flush=True)
def run(label,command):
 with (OUT/(label+'.log')).open('a') as f:subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 return label
jobs=[]
def add(mode,task='design',band=(80,200),side=None,variant='full',steps=16,split='test'):
 name=f'attention_{mode}_{split}_{task}_{side or variant}_amp{band[0]}-{band[1]}_steps{steps}_seed20261013'
 args=[PY,'experiments/attention/evaluate.py','--method',mode,'--task',task,'--amplicon-min',str(band[0]),'--amplicon-max',str(band[1]),'--split',split,'--steps',str(steps),'--variant',variant]
 if side:args+=['--fixed-side',side]
 jobs.append((name,args))
for mode in ['diffusion','direct']:
 add(mode,split='validation')
 for band in [(80,200),(80,120),(121,160),(161,200)]:add(mode,band=band)
 for task in ['completion','redesign']:
  for side in ['forward','reverse']:add(mode,task=task,side=side)
for steps in [1,4,8]:add('diffusion',steps=steps)
add('diffusion',task='ablation',variant='single_oligo')
protocol=dict(declared_at=utc_now(),seed=20261013,architecture=dict(type='conditional_transformer',width=96,layers=4,heads=4,ff_mult=2,position_encoding='rotary',conditioning='DNA/noisy-label embeddings + per-block FiLM design/time conditioning'),training=dict(epochs=30,batch_size=16,learning_rate=.001,weight_decay=.01,window_length=512,best_checkpoint='minimum frozen validation masked-CE',fixed_side_augmentation=.25),sampling=dict(steps=16,temperature=4,budget=20,decoder='same joint-Tm and anchor-prefilter support as final CNN analysis'),expected_runs=[n for n,_ in jobs],additional_baseline='qPrimer 2025; same target/reference/assessment',dataset_sha256=sha256(ROOT/'data/processed/dataset_partial_pilot.jsonl'),target_sha256=sha256(ROOT/'data/processed/targets_pilot.jsonl'),original_results='results/final_fixed_seed/study_summary.json',comparison='All actual outputs retained; shared full-reference screens and gene-level endpoints; one training/generation seed as user requested')
if not (BASE/'protocol.json').exists():write_json(BASE/'protocol.json',protocol)
if '--plan-only' in sys.argv:print(json.dumps(protocol,indent=2));sys.exit()
try:
 status('training')
 with ThreadPoolExecutor(max_workers=2) as ex:
  fs=[]
  for mode in ['diffusion','direct']:
   history=OUT/f'attention_{mode}_training.json'
   if history.exists() and len(json.loads(history.read_text())['history'])==30:continue
   fs.append(ex.submit(run,'train_'+mode,[PY,'experiments/attention/train.py','--mode',mode]))
  for f in as_completed(fs):status('training',completed=f.result())
 status('generation',expected_runs=len(jobs));finished=[]
 with ThreadPoolExecutor(max_workers=2) as ex:
  fs=[]
  for name,command in jobs:
   if (OUT/'evaluations'/name/'completion.json').exists():finished.append(name);continue
   fs.append(ex.submit(run,name,command))
  for f in as_completed(fs):finished.append(f.result());status('generation',completed_runs=len(finished),expected_runs=len(jobs),last_completed=finished[-1])
 dirs=[OUT/'evaluations'/n for n,_ in jobs]
 dirs+=sorted((OUT/'evaluations').glob('qprimer*'))
 status('reference_screening',runs=len(dirs))
 run('screening',[CPU,'scripts/screen_bowtie.py',*map(str,dirs),'--threads','8'])
 status('analysis');run('analysis',[CPU,'experiments/attention/analyze.py'])
 status('attention_experiments_complete',runs=len(dirs))
except Exception as e:
 status('error',error=str(e));raise
