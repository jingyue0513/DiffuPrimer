#!/usr/bin/env python3
"""Stage an explicit public release allowlist; never add private correspondence."""
from pathlib import Path
import json,subprocess
ROOT=Path(__file__).resolve().parents[2]
paths=[ROOT/n for n in ['.gitignore','.zenodo.json','LICENSE','CITATION.cff','README.md',
 'THIRD_PARTY_NOTICES.md','pyproject.toml','requirements-runtime.txt']]
for folder in ['src','tests','configs','examples','scripts']:
 for p in (ROOT/folder).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix in {'.py','.json','.fasta','.md','.cjs'}:paths.append(p)
for base in (ROOT/'experiments').iterdir():
 if base.is_dir():
  paths += [p for p in base.glob('*') if p.is_file() and p.suffix in {'.py','.json','.md'}]
  for name in ['src','data','configs']:
   item=base/name
   if item.is_symlink():paths.append(item)
   elif name=='src' and item.exists():paths += list(item.rglob('*.py'))
for p in (ROOT/'experiments/attention/baselines').glob('*'):
 if p.is_file() and p.suffix in {'.py','.json','.md'}:paths.append(p)
for name in ['reproduction.md','data_dictionary.md','independent_validation.md','attention_implementation_record.md',
             'attention_related_work.md','bioinformatics_submission_requirements.md']:
 paths.append(ROOT/'docs'/name)
for p in (ROOT/'release').iterdir():
 if p.is_file() and p.suffix in {'.json','.zip','.md'}:paths.append(p)
for name in ['bowtie_manifest.json','manifest.json']:paths.append(ROOT/'tools'/name)
for name in ['bowtie_exhaustive','sqlite_specificity']:
 paths.append(ROOT/'results/checks'/name/'check.json')
paths=sorted({p for p in paths if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts})
for p in paths:
 assert p.stat().st_size<95*1024**2,p
 if p.is_symlink():assert p.resolve().is_relative_to(ROOT) and not str(p.readlink()).startswith('/'),p
rel=[str(p.relative_to(ROOT)) for p in paths]
subprocess.run(['git','add','--force','--',*rel],cwd=ROOT,check=True)
(ROOT/'submission/public_release_files.json').write_text(json.dumps(rel,indent=2)+'\n')
print(len(rel),'files staged;',sum(p.stat().st_size for p in paths),'bytes')
