#!/usr/bin/env python3
"""Build review/upload files without submitting them or creating archive IDs."""
import hashlib,json,re,shutil,subprocess,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'paper/bioinformatics_submission'
OUT=ROOT/'submission'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zip_files(dest,files,base):
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(files)):z.write(p,str(p.relative_to(base)))
def main():
    OUT.mkdir(exist_ok=True);upload=OUT/'upload';upload.mkdir(exist_ok=True)
    for name in ['manuscript.pdf','supplementary.pdf']:
        shutil.copy2(PAPER/name,upload/name)
    for name in ['cover_letter.txt','author_metadata.json']:
        shutil.copy2(OUT/name,upload/name)
    figures=upload/'figures';figures.mkdir(exist_ok=True)
    for i,name in enumerate(['01_method.pdf','main_02_quality.pdf','main_03_flexibility.pdf'],1):
        shutil.copy2(PAPER/'figures'/name,figures/f'figure_{i}.pdf')
    for stem in ['manuscript','supplementary']:
        text=(PAPER/f'{stem}.tex').read_text()
        files=[PAPER/f'{stem}.tex',PAPER/'oup-authoring-template.cls',PAPER/'submission_formatting.tex']
        for name in re.findall(r'\\input\{([^}]+)\}',text):
            if name=='glyphtounicode':continue
            files.append(PAPER/(name if Path(name).suffix else name+'.tex'))
        for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',text):files.append(PAPER/name)
        if stem=='manuscript':files += [PAPER/'manuscript.bbl',PAPER/'oup-abbrvnat.bst',PAPER/'submission_references.bib',PAPER/'references.bib']
        files += [PAPER/'source_manifest.json',PAPER/'verification.json',PAPER/'reference_verification.json']
        zip_files(upload/f'{stem}_source.zip',files,PAPER)
    zip_files(upload/'supplementary_data.zip',[PAPER/'mechanism_evidence.json',PAPER/'mechanism_results.csv'],PAPER)
    commit=subprocess.check_output(['git','rev-parse','v0.1.0^{}'],cwd=ROOT,text=True).strip()
    subprocess.run(['git','archive','--format=zip','--output',str(OUT/'software_archive.zip'),commit],cwd=ROOT,check=True)
    subprocess.run(['git','bundle','create',str(OUT/'DiffuPrimer.git.bundle'),'--all'],cwd=ROOT,check=True)
    abstract=(PAPER/'manuscript.tex').read_text().split('\\abstract{',1)[1].split('\\keywords',1)[0]
    macros=dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^}]+)\}',(PAPER/'result_values.tex').read_text()))
    for key,value in macros.items():abstract=abstract.replace('\\'+key+'{}',value)
    abstract=re.sub(r'\\href\{([^}]+)\}\{([^}]+)\}',lambda m:m[2] if m[1].startswith('mailto:') else m[2]+' ('+m[1]+')',abstract)
    abstract=re.sub(r'\\url\{([^}]+)\}',r'\1',abstract)
    abstract=re.sub(r'\\textbf\{([^}]+)\}',r'\1',abstract)
    abstract=abstract.replace('\\\\','\n').rstrip('}\n')
    (upload/'abstract.txt').write_text(abstract+'\n')
    files=[p for p in upload.rglob('*') if p.is_file()]
    files += [OUT/'software_archive.zip',OUT/'DiffuPrimer.git.bundle',OUT/'archive_metadata.json',OUT/'README.md']
    report={'software_commit':commit,'public_repository':'https://github.com/jingyue0513/DiffuPrimer',
      'state':'prepared_author_review_required','files':[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]}
    (OUT/'package_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'SHA256SUMS.txt').write_text(''.join(f"{r['sha256']}  {r['path']}\n" for r in report['files']))
    zip_files(OUT/'DiffuPrimer_submission_package.zip',files+[OUT/'package_manifest.json',OUT/'SHA256SUMS.txt'],OUT)
    print(json.dumps({'software_commit':commit,'files':len(files),'package_bytes':(OUT/'DiffuPrimer_submission_package.zip').stat().st_size},indent=2))
if __name__=='__main__':main()
