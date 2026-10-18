import os,sys,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];out=ROOT/'experiments/attention_external/results';out.mkdir(exist_ok=True)
methods=sys.argv[1:] or ['hybrid_diffusion','hybrid_direct','cnn_diffusion','cnn_direct','rope_diffusion']
def run(m):
 with (out/(m+'_generation.log')).open('a') as f:subprocess.run([str(ROOT/'.venv_attention/bin/python'),'experiments/attention_external/evaluate_learned.py','--method',m],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
 return m
with ThreadPoolExecutor(max_workers=2) as ex:
 for f in as_completed([ex.submit(run,m) for m in methods]):print('Completed external learned model:',f.result(),flush=True)
