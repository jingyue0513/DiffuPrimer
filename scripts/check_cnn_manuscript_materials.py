#!/usr/bin/env python3
"""Verify frozen-source numbers, bibliography coverage and compiled draft layout."""
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'paper/bioinformatics_cnn_draft'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    manifest=json.loads((OUT/'source_manifest.json').read_text())
    sources=[]
    for source in manifest['source_summaries']:
        p=ROOT/source['path'];assert sha(p)==source['sha256'];sources.append(json.loads(p.read_text()))
    original,extra=sources
    macros=dict(re.findall(r'\\newcommand\{\\(\w+)\}\{([^}]+)\}',(OUT/'result_values.tex').read_text()))
    checked=0
    for record in manifest['result_values']:
        data=original if record.get('cohort')=='original' else extra
        if 'contrast' in record:
            c=record['contrast'];assert c in data['paired_comparisons']
            for suffix,v in [('Delta',c['mean_difference']),('Lo',c['ci95'][0]),('Hi',c['ci95'][1])]:
                assert macros[record['macro_prefix']+suffix]==f"{v*record['scale']:.{record['digits']}f}";checked+=1
            assert macros[record['macro_prefix']+'N']==str(c['n']);checked+=1
        elif 'derivation' in record:
            data=original if record['macro'].startswith('Main') else extra
            groups=[g for g in data['summaries'] if (g['context'],g['variant'],g['budget'])==('amp80-200','full',20)]
            means={g['method']:g['unique_RNA_pass_pairs_without_genomic_flag'] for g in groups}
            v=100*(means['cnn_diffusion']/means['primer3']-1);assert record['value']==v
            assert macros[record['macro']]==f'{v:.1f}';checked+=1
        else:
            g=next(g for g in data['summaries'] if (g['method'],g['context'],g['variant'],g['budget'])==
                   (record['method'],record['context'],record['variant'],record['budget']))
            v=g[record['endpoint']];assert v==record['value'];assert g[record['endpoint']+'_n']==record['n']
            assert macros[record['macro']]==f"{v*record['scale']:.{record['digits']}f}";checked+=1
    assert checked==len(macros)
    for fig in manifest['figures']:
        assert sha(ROOT/fig['source'])==fig['sha256']==sha(OUT/fig['copy'])
    for fig in manifest['main_figure_subsets']:
        assert sha(OUT/fig['source'])==fig['source_sha256'];assert sha(OUT/fig['main_output'])==fig['main_sha256']
        assert abs(float(PdfReader(OUT/fig['source']).pages[0].mediabox.width)-float(PdfReader(OUT/fig['main_output']).pages[0].mediabox.width))<.001
    tex=(OUT/'manuscript.tex').read_text();supp=(OUT/'supplementary.tex').read_text()
    cite_keys={k.strip() for b in re.findall(r'\\cite\w*\{([^}]+)\}',tex) for k in b.split(',')}
    full_keys=set(re.findall(r'^@\w+\{([^,]+),',(OUT/'references.bib').read_text(),re.M))
    print_keys=set(re.findall(r'^@\w+\{([^,]+),',(OUT/'submission_references.bib').read_text(),re.M))
    assert cite_keys==full_keys==print_keys and len(cite_keys)==52
    assert len(re.findall(r'\\bibitem', (OUT/'manuscript.bbl').read_text()))==52
    reference_check=json.loads((OUT/'reference_verification.json').read_text())
    assert reference_check['identity_verified_count']==52 and reference_check['references_bib_sha256']==sha(OUT/'references.bib')
    assert '\\nocite' not in tex
    rendered=tex
    for name,value in sorted(macros.items(),key=lambda x:-len(x[0])):rendered=rendered.replace('\\'+name+'{}',value)
    abstract=rendered.split('\\abstract{',1)[1].split('\\keywords',1)[0]
    abstract=re.sub(r'\\href\{[^}]+\}\{([^}]+)\}',r'\1',abstract)
    abstract=re.sub(r'\\rule\{[^}]+\}\{[^}]+\}','',abstract)
    abstract=re.sub(r'\\textbf\{([^}]+)\}',r'\1',abstract)
    abstract=re.sub(r'\\[A-Za-z]+','',abstract)
    abstract_words=len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*",abstract))
    assert abstract_words<=150
    main_pdf=PdfReader(OUT/'manuscript.pdf');supp_pdf=PdfReader(OUT/'supplementary.pdf')
    assert len(main_pdf.pages)<=7 and len(supp_pdf.pages)>0
    assert len(re.findall(r'\\begin\{figure\*\}',tex))==4
    assert len(re.findall(r'\\begin\{figure\*\}',supp))==4
    implementation=tex.split('\\section{Implementation and results}',1)[1].split('\\section{Discussion}',1)[0]
    assert len(re.findall(r'\\subsection\{',implementation))==8
    for name in ['manuscript','supplementary']:
        log=(OUT/(name+'.log')).read_text(errors='replace')
        assert not re.search(r'^!',log,re.M)
        assert not re.search(r'(Citation .+ undefined|Reference .+ undefined|There were undefined references)',log)
        assert not [line for line in log.splitlines() if 'Overfull' in line and 'output is active' not in line]
    report={'at':datetime.now(timezone.utc).isoformat(),'passed':True,'primary_method':'cnn_diffusion',
            'numerical_macros_verified':checked,'source_summaries_verified':2,'references_verified_and_cited':52,
            'abstract_word_count':abstract_words,'main_pages':len(main_pdf.pages),'supplementary_pages':len(supp_pdf.pages),
            'main_figures':4,'supplementary_complete_figures':4,'experimental_questions':8,
            'reference_style':'OUP author-date; compact print bibliography; full metadata preserved',
            'source_manifest_sha256':sha(OUT/'source_manifest.json'),'manuscript_pdf_sha256':sha(OUT/'manuscript.pdf'),
            'supplementary_pdf_sha256':sha(OUT/'supplementary.pdf'),'undefined_citations':False,'body_or_table_overflows':False,
            'class_output_box_warnings':'OUP template emits output-routine box warnings; page previews reviewed separately',
            'journal_llm_drafting_policy_compliance_asserted':False,'public_code_url':'','public_data_url':''}
    (OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'layout_report.json').write_text(json.dumps({k:v for k,v in report.items() if k not in ['public_code_url','public_data_url']},indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':main()
