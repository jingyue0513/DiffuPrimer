#!/usr/bin/env python3
"""Verify compiled manuscript layout and frozen numerical provenance."""
import hashlib,json,re
from pathlib import Path
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'paper/bioinformatics_submission'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    manifest=json.loads((OUT/'source_manifest.json').read_text())
    for source in manifest['sources']:assert sha(ROOT/source['path'])==source['sha256'],source
    old=json.loads((ROOT/'paper/bioinformatics_cnn_draft/source_manifest.json').read_text())
    datasets=[]
    for source in old['source_summaries']:
        path=ROOT/source['path'];assert sha(path)==source['sha256'];datasets.append(json.loads(path.read_text()))
    macros=dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^}]+)\}',(OUT/'result_values.tex').read_text()))
    checked=0
    for record in old['result_values']:
        data=datasets[0] if record.get('cohort')=='original' else datasets[1]
        if 'contrast' in record:
            c=record['contrast'];assert c in data['paired_comparisons']
            for suffix,v in [('Delta',c['mean_difference']),('Lo',c['ci95'][0]),('Hi',c['ci95'][1])]:
                assert macros[record['macro_prefix']+suffix]==f"{v*record['scale']:.{record['digits']}f}";checked+=1
            assert macros[record['macro_prefix']+'N']==str(c['n']);checked+=1
        elif 'derivation' in record:
            data=datasets[0] if record['macro'].startswith('Main') else datasets[1]
            gs=[g for g in data['summaries'] if (g['context'],g['variant'],g['budget'])==('amp80-200','full',20)]
            means={g['method']:g['unique_RNA_pass_pairs_without_genomic_flag'] for g in gs}
            v=100*(means['cnn_diffusion']/means['primer3']-1)
            assert macros[record['macro']]==f'{v:.1f}';checked+=1
        else:
            g=next(g for g in data['summaries'] if (g['method'],g['context'],g['variant'],g['budget'])==
                (record['method'],record['context'],record['variant'],record['budget']))
            v=g[record['endpoint']];assert v==record['value']
            assert macros[record['macro']]==f"{v*record['scale']:.{record['digits']}f}";checked+=1
    assert checked==len(macros)
    evidence=json.loads((OUT/'mechanism_evidence.json').read_text())
    assert evidence==json.loads((ROOT/'experiments/mechanism/uncapped_seed42/summary.json').read_text())
    table=(OUT/'mechanism_table.tex').read_text()
    methods=['direct_validation_selected','diffusion1_validation_selected','frozen_4','diffusion_4','frozen_16','diffusion_16','direct_equal_cpu_uncapped']
    for g in evidence['groups']:
        if g['method'] in methods:
            assert f"{g['retained_pairs']:.2f} & {100*g['one_retained']:.1f} & {g['generation_cpu_seconds']:.2f}" in table
    text=(OUT/'manuscript.tex').read_text();supp=(OUT/'supplementary.tex').read_text()
    assert 'Author contributions' not in text and 'CRediT' not in text
    assert 'Shiyue Mei' in text and 'lijingyue513@163.com' in text
    assert 'The authors declare no conflict of interest.' in text
    assert 'no specific grant' in text and 'manuscript drafting and revision' in text
    citations={k for b in re.findall(r'\\cite\w*\{([^}]+)\}',text) for k in b.split(',')}
    compiled=set(re.findall(r'\\bibitem\[.*?\]\{([^}]+)\}',(OUT/'manuscript.bbl').read_text(),re.S))
    assert citations==compiled and len(citations)==51,(len(citations),len(compiled))
    abstract=text.split('\\abstract{',1)[1].split('\\keywords',1)[0]
    abstract=re.sub(r'\\href\{[^}]+\}\{([^}]+)\}',r'\1',abstract)
    abstract=re.sub(r'\\url\{[^}]+\}','URL',abstract)
    abstract=re.sub(r'\\textbf\{[^}]+\}','',abstract)
    abstract=re.sub(r'\\[A-Za-z]+','',abstract)
    words=len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*",abstract));assert words<=150
    pdf=PdfReader(OUT/'manuscript.pdf');spdf=PdfReader(OUT/'supplementary.pdf');assert len(pdf.pages)<=7
    assert len(re.findall(r'\\begin\{figure\*\}',text))==3
    assert len(re.findall(r'\\begin\{figure\*\}',supp))==5
    for name in ['manuscript','supplementary']:
        log=(OUT/(name+'.log')).read_text(errors='replace')
        assert not re.search(r'^!',log,re.M)
        assert not re.search(r'(Citation .+ undefined|Reference .+ undefined|There were undefined references)',log)
        assert not [line for line in log.splitlines() if 'Overfull' in line and 'output is active' not in line]
    report={'passed':True,'main_pages':len(pdf.pages),'supplementary_pages':len(spdf.pages),
        'abstract_words_excluding_heading_labels':words,'cited_references':len(citations),
        'historical_numerical_macros_verified':checked,'mechanism_methods_per_cohort':14,
        'main_figures':3,'supplementary_complete_figures':5,'author_contribution_paragraph_in_main':False,
        'undefined_citations_or_references':False,'body_or_table_overflow':False,'source_hashes_verified':True,
        'main_pdf_sha256':sha(OUT/'manuscript.pdf'),'supplementary_pdf_sha256':sha(OUT/'supplementary.pdf'),
        'submission_ready':False,'pending':['GitHub publication verification','independent archive persistent URL',
            'researcher revision and journal LLM-policy resolution','author final approval and submission-system metadata']}
    (OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
