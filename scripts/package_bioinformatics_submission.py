#!/usr/bin/env python3
"""Build a clean journal package; retain internal records outside upload files."""
import csv,hashlib,io,json,re,shutil,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'paper/bioinformatics_submission'
OUT=ROOT/'submission'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zip_files(dest,files,base):
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(files)):z.write(p,str(p.relative_to(base)))
def scientific(value):
    """Remove machine-local provenance fields, preserving scientific values."""
    if isinstance(value,dict):
        return {k:scientific(v) for k,v in value.items() if
                not any(term in k.lower() for term in ['run','path','source','directory','sha256'])}
    if isinstance(value,list):return [scientific(v) for v in value]
    return value

def supplementary_data(dest):
    payload={}
    for cohort,relative in [('original','experiments/attention/results/study_summary.json'),
                            ('additional','experiments/attention_external/results/study_summary.json')]:
        d=json.loads((ROOT/relative).read_text())
        payload['historical_'+cohort+'.json']=scientific({k:d[k] for k in
                    ['seed','evaluation_split','total_test_genes','endpoints','summaries','paired_comparisons']})
    evidence=json.loads((PAPER/'mechanism_evidence.json').read_text())
    payload['mechanism_controls.json']={k:evidence[k] for k in
        ['seed','groups','main_contrasts','validation_selected_contrasts','uncapped_contrasts','uncapped_budget_diagnostics']}
    payload['mechanism_validation_selection.json']=scientific(evidence['validation_selection'])
    payload['figure_values.json']=scientific(json.loads((ROOT/'paper/figures/cnn_publication_results/figure_values.json').read_text()))
    payload['biological_case.json']=scientific(json.loads((ROOT/'paper/figures/cnn_biological_analysis/illustrative_case.json').read_text()))
    bio=json.loads((ROOT/'paper/figures/cnn_biological_analysis/verification.json').read_text())
    payload['biological_geometry.json']={'cohort_genes':bio['cohort_genes'],'aggregate_values':bio['aggregate_values'],'category_values':bio['category_values']}
    description='''Supplementary numerical data for DiffuPrimer

historical_original.json and historical_additional.json contain the separate
75-gene cohort summaries and gene-paired contrasts used in the main benchmark,
Tables S1-S4 and Figures 2-3/S1-S4. Seed: 20261013.
mechanism_controls.json and mechanism_results.csv contain the 14 seed-42
mechanism configurations in each cohort (Table 1, Table S5 and Figure S5).
mechanism_validation_selection.json documents validation-only decoder selection.
figure_values.json contains the numerical values used in Figures 2-3/S1-S4;
biological_geometry.json and biological_case.json support Figure S3.

Retained pairs pass the shared computational constraints, intended-gene RNA
specificity and absence of a predicted genomic product under the reported screen.
Counts average over all 75 genes. Conditional means exclude undefined genes;
_n fields give denominators and ci95 fields give 95% gene-bootstrap intervals.
Fractions and coverage are stored on the 0-1 scale, with figure scale factors
recorded explicitly. Jaccard distance is 1 minus interval overlap/union.
Null numerical values indicate undefined endpoints, not missing experiments
or author placeholders. CPU times are process-CPU seconds for generation;
shared span preparation and reference screening are excluded.

Cohort labels test/original and external_test/additional refer to the original
and additional cohorts respectively. The historical and seed-42 studies use
different recorded generation seeds and must not be pooled. The cap-4096
control and uncapped matched-CPU control are distinct configurations.
Complete candidates, screening records, test sequences, checkpoints and
reproduction code are provided with the versioned repository release:
https://github.com/jingyue0513/DiffuPrimer
https://archive.softwareheritage.org/swh:1:rev:275d1c0acf64edba234d17c1a83094af597f6202/
'''
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,value in payload.items():z.writestr(name,json.dumps(value,indent=2,ensure_ascii=False)+'\n')
        z.write(PAPER/'mechanism_results.csv','mechanism_results.csv')
        z.writestr('README.txt',description)

def main():
    OUT.mkdir(exist_ok=True);upload=OUT/'upload';upload.mkdir(exist_ok=True)
    forms=OUT/'form_fields';forms.mkdir(exist_ok=True)
    for name in ['manuscript.pdf','supplementary.pdf']:shutil.copy2(PAPER/name,upload/name)
    shutil.copy2(OUT/'cover_letter.txt',upload/'cover_letter.txt')
    shutil.copy2(OUT/'author_metadata.json',forms/'author_metadata.json')
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
        zip_files(upload/f'{stem}_source.zip',files,PAPER)
    supplementary_data(upload/'supplementary_data.zip')
    abstract=(PAPER/'manuscript.tex').read_text().split('\\abstract{',1)[1].split('\\keywords',1)[0]
    macros=dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^}]+)\}',(PAPER/'result_values.tex').read_text()))
    for key,value in macros.items():abstract=abstract.replace('\\'+key+'{}',value)
    abstract=re.sub(r'\\href\{([^}]+)\}\{([^}]+)\}',lambda m:m[2] if m[1].startswith('mailto:') else m[2]+' ('+m[1]+')',abstract)
    abstract=re.sub(r'\\url\{([^}]+)\}',r'\1',abstract)
    abstract=re.sub(r'\\textbf\{([^}]+)\}',r'\1',abstract)
    abstract=abstract.replace('\\\\','\n').rstrip('}\n')
    (forms/'abstract.txt').write_text(abstract+'\n')
    # These are form-entry references, never peer-review attachments.
    for name in ['abstract.txt','author_metadata.json']:
        if (upload/name).exists():(upload/name).unlink()
    allowed={'manuscript.pdf','supplementary.pdf','cover_letter.txt','manuscript_source.zip',
             'supplementary_source.zip','supplementary_data.zip','figures/figure_1.pdf',
             'figures/figure_2.pdf','figures/figure_3.pdf'}
    assert {str(p.relative_to(upload)) for p in upload.rglob('*') if p.is_file()}==allowed
    files=[p for p in upload.rglob('*') if p.is_file()]+[forms/'abstract.txt',forms/'author_metadata.json',forms/'figure_alt_text.txt',OUT/'UPLOAD_GUIDE.txt']
    report={'software_commit':'275d1c0acf64edba234d17c1a83094af597f6202',
      'public_repository':'https://github.com/jingyue0513/DiffuPrimer',
      'state':'technical_checks_passed_author_submission_declarations_required',
      'files':[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]}
    (OUT/'package_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'SHA256SUMS.txt').write_text(''.join(f"{r['sha256']}  {r['path']}\n" for r in report['files']))
    # The clean handoff contains journal attachments and clearly separated form
    # references. Internal audit reports and this checksum manifest stay local.
    zip_files(OUT/'DiffuPrimer_submission_package.zip',files,OUT)
    print(json.dumps({'files':len(files),'journal_attachments':len(allowed),
         'package_bytes':(OUT/'DiffuPrimer_submission_package.zip').stat().st_size},indent=2))
if __name__=='__main__':main()
