"""Screen completed new runs using the unchanged common Bowtie evaluator."""
import argparse,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('roots',nargs='+');p.add_argument('--threads',type=int,default=8);a=p.parse_args()
dirs=[]
for raw in a.roots:
 for d in sorted((ROOT/raw).iterdir()):
  if not d.is_dir() or not (d/'completion.json').exists():continue
  complete=json.loads((d/'completion.json').read_text());status=d/'specificity_bowtie_status.json'
  if status.exists() and json.loads(status.read_text()).get('candidate_sha256')==complete['candidates_sha256']:continue
  dirs.append(d)
print('New completed runs requiring reference screen:',len(dirs),flush=True)
if dirs:subprocess.run([str(ROOT/'.venv/bin/python'),'scripts/screen_bowtie.py',*map(str,dirs),'--threads',str(a.threads)],cwd=ROOT,check=True)
