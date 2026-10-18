#!/usr/bin/env python3
"""Typeset author-editable Bioinformatics result material, not an authored paper.

Numerical tables are loaded only from completed declared-protocol analyses.
--layout-preview reserves result areas without inserting unfinished results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from diffuprimer.io import read_jsonl,sha256,utc_now,write_json

ORDER=['primer3','qprimer','primercast','uniform_pair_sampler','rule_sampler','cnn_direct','cnn_diffusion',
       'attention_rope_direct','attention_rope_diffusion','attention_hybrid_direct','attention_hybrid_diffusion']
NAMES={'primer3':'Primer3','qprimer':'qPrimer','primercast':'PrimerCast','uniform_pair_sampler':'Uniform pairs','rule_sampler':'Forward-first',
       'cnn_direct':'CNN direct','cnn_diffusion':'CNN diffusion','attention_rope_direct':'RoPE direct','attention_rope_diffusion':'RoPE diffusion',
       'attention_hybrid_direct':'Hybrid direct','attention_hybrid_diffusion':'Hybrid diffusion'}
RNA='unique_RNA_pass_pairs';RETAINED='unique_RNA_pass_pairs_without_genomic_flag'


def tex(value):
    value=str(value).replace('–','--').replace('−','-').replace('≥',' >= ').replace('≤',' <= ').replace('°',' degrees ')
    replacements={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
    return ''.join(replacements.get(c,c) for c in value)


def load_completed(summary_path,protocol_path):
    raw=summary_path.read_bytes();summary=json.loads(raw)
    if summary['status']!='complete' or not summary['additional_experiments_complete']:
        raise ValueError('Numerical author material requires completed analyses: '+str(summary_path))
    if summary.get('protocol_sha256')!=sha256(protocol_path):
        raise ValueError('Analysis does not match declared protocol: '+str(protocol_path))
    if summary.get('pending_planned_runs'):
        raise ValueError('Declared experiment sources are pending')
    return summary,hashlib.sha256(raw).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--protocol',default='experiments/attention/final_protocol.json')
    p.add_argument('--summary',default='experiments/attention/results/study_summary.json')
    p.add_argument('--external-protocol',default='experiments/attention_external/protocol.json')
    p.add_argument('--external-summary',default='experiments/attention_external/results/study_summary.json')
    p.add_argument('--output',default='paper/bioinformatics_attention_materials')
    p.add_argument('--layout-preview',action='store_true',help='Compile structural fields and empty result areas; no unfinished numerical results')
    p.add_argument('--fifth-figure',choices=['generalization','budget'],default='generalization')
    p.add_argument('--no-compile',action='store_true')
    args=p.parse_args()
    out=ROOT/args.output;out.mkdir(parents=True,exist_ok=True);(out/'figures').mkdir(exist_ok=True)
    protocol_path=ROOT/args.protocol;protocol=json.loads(protocol_path.read_text())
    architecture=protocol.get('primary_architecture',protocol['architecture'])
    config=json.loads((ROOT/'configs/pilot.json').read_text());conditions=config['conditions']
    dataset_path=ROOT/'data/processed/dataset_partial_pilot.manifest.json';dataset=json.loads(dataset_path.read_text())
    targets=read_jsonl(ROOT/'data/processed/targets_pilot.jsonl')
    counts={s:sum(t['split']==s for t in targets) for s in ['train','validation','test']}
    external_manifest_path=ROOT/'experiments/attention_external/manifest.json'
    external_manifest=json.loads(external_manifest_path.read_text())
    summary=external=None;main_sha=external_sha=None
    if not args.layout_preview:
        summary,main_sha=load_completed(ROOT/args.summary,protocol_path)
        external,external_sha=load_completed(ROOT/args.external_summary,ROOT/args.external_protocol)
        if external.get('evaluation_split')!='external_test' or external.get('old_sources_included'):
            raise ValueError('Independent cohort must not contain original study sources')
        if external['targets_sha256']!=external_manifest['targets_sha256']:
            raise ValueError('External analysis is not the frozen additional-gene cohort')
    primary=protocol['primary_methods']['diffusion'];direct=protocol['primary_methods']['direct']
    data_values=[];figure_records=[]
    official=ROOT/'paper/bioinformatics_template/oup-authoring-template'
    for name in ['oup-authoring-template.cls','oup-abbrvnat.bst']:
        shutil.copyfile(official/name,out/name)
    shutil.copyfile(ROOT/'paper/bioinformatics_working/references.bib',out/'references.bib')

    def fact_table(rows,caption,label):
        content='\n'.join(tex(a)+' & '+b+r' \\' for a,b in rows)
        return r'\begin{table}[t]\caption{'+tex(caption)+r'}\label{'+label+r'}\centering\small\begin{tabular}{@{}p{.45\columnwidth}p{.48\columnwidth}@{}}\toprule Field & Recorded value \\\midrule'+'\n'+content+r'\bottomrule\end{tabular}\end{table}'

    architecture_table=fact_table([
        ('Hidden dimensions',str(architecture['width'])),('Local residual CNN blocks',str(architecture['layers'])),
        ('Dilation factors',r'$1,4,16,64$'),('Conditional attention blocks',str(architecture.get('attention_layers',architecture['layers']))),
        ('Attention heads',str(architecture['heads'])),('Position features','Relative transcript; rotary Q/K'),
        ('Training window',str(protocol['training']['window_length'])+' bp'),('Epochs / batch',str(protocol['training']['epochs'])+' / '+str(protocol['training']['batch_size'])),
        ('AdamW learning rate',r'$10^{-3}$'),('Anchor augmentation',r'25\%'),
        ('Reverse steps / proposal K',str(protocol['sampling']['steps'])+' / '+str(protocol['sampling']['budget'])),('Fixed seed',str(protocol['seed']))
    ],'Recorded primary architecture and training configuration.','tab:architecture')
    dataset_table=fact_table([
        ('Reference assembly',tex(config['reference_assembly'])),('Species','Mus musculus'),('Original genes',str(len(targets))),
        ('Train / validation / test',f"{counts['train']} / {counts['validation']} / {counts['test']}"),
        ('Teacher position records',str(dataset['by_evidence']['computational_candidate'])),('Mapped source records',str(dataset['by_evidence']['unverified'])),
        ('Training / validation labels',f"{dataset['by_split']['train']} / {dataset['by_split']['validation']}"),
        ('Test-source labels excluded',str(dataset['by_split']['test'])),('Additional held-out genes',str(external_manifest['selected_genes'])),
        ('New/original ID overlap',str(external_manifest['source_check']['gene_ID_overlap_with_original_500'])),
        ('New/training homology edges',str(external_manifest['source_check']['selected_training_homology_edges']))
    ],'Frozen public-reference cohorts and position-label records.','tab:data')
    conditions_table=fact_table([
        ('Primer length',f"{conditions['primer_min_length']}--{conditions['primer_max_length']} nt"),('Product length',f"{conditions['amplicon_min']}--{conditions['amplicon_max']} bp"),
        ('Conditional length bands','80--120; 121--160; 161--200 bp'),('Tm interval / optimum',r'$57$--$63\,^{\circ}\mathrm C$ / $60\,^{\circ}\mathrm C$'),
        ('Pair Tm difference',r'$\leq 2\,^{\circ}\mathrm C$'),('GC interval',r'35--65\%'),('Homopolymer limit','5'),('Structure Tm ceiling',r'$47\,^{\circ}\mathrm C$'),
        ('Monovalent / divalent',r'$50$ / $1.5$ mM'),('dNTP / oligo',r'$0.6$ mM / $50$ nM'),
        ('Reference alignment',r'Ungapped; $\leq2$ mismatches'),("Exact primer 3-prime end",'5 nt'),('Screened product interval','36--1000 bp')
    ],'Shared design and complete-reference assessment conditions.','tab:conditions')

    def group_map(doc):return {(g['context'],g['method'],g['variant'],g['budget']):g for g in doc['summaries']}
    def value(g,key,ci=False):
        if g is None or g[key] is None:return '--'
        data_values.append({'context':g['context'],'method':g['method'],'variant':g['variant'],'budget':g['budget'],'endpoint':key,'value':g[key],'n':g[key+'_n'],'ci95':g[key+'_ci95']})
        if ci:
            lo,hi=g[key+'_ci95'];return f"{g[key]:.2f} [{lo:.2f}, {hi:.2f}]"
        if key=='RNA_no_genome_fraction':return f"{g[key]*100:.1f} ({g[key+'_n']})"
        return f"{g[key]:.3f}" if key in ('RNA_pass_amplicon_distance','best_RNA_pass_temperature_penalty') else f"{g[key]:.2f}"
    if summary:
        gs,es=group_map(summary),group_map(external)
        rows=[]
        for method in ORDER:
            own=gs.get(('amp80-200',method,'full',20));extra=es.get(('amp80-200',method,'full',20))
            if own is None and extra is None:continue
            rows.append(tex(NAMES.get(method,method))+' & '+value(own,RNA)+' & '+value(own,RETAINED)+' & '+value(own,'RNA_no_genome_fraction')+' & '+value(extra,RNA)+' & '+value(extra,RETAINED)+' & '+value(extra,'RNA_no_genome_fraction')+r' \\')
        quality_table=r'''\begin{table*}[t]\caption{Recorded candidate supply and composition at $K=20$, with cohorts reported separately. RNA and retained columns count unique pairs per gene; retained pairs have no predicted genomic product. Fraction is the mean per-gene retained/RNA percentage, with its defined-gene count in parentheses. Count endpoints include all cohort genes; Figure~\ref{fig:02_quality} gives gene-bootstrap intervals.}\label{tab:quality}\centering\small
\begin{tabular}{@{}lrrrrrr@{}}\toprule & \multicolumn{3}{c}{Original ($N=__MAIN_N__$)} & \multicolumn{3}{c}{Additional ($N=__EXTERNAL_N__$)} \\
Method & RNA & Retained & Fraction (\%, $n$) & RNA & Retained & Fraction (\%, $n$) \\\midrule
'''+ '\n'.join(rows)+r'\bottomrule\end{tabular}\end{table*}'
        quality_table=quality_table.replace('__MAIN_N__',str(summary['total_test_genes'])).replace('__EXTERNAL_N__',str(external['total_test_genes']))
        paired=[]
        for context,title in [('amp80-200','Complete design'),('amp80-120','80--120 bp'),('amp121-160','121--160 bp'),('amp161-200','161--200 bp'),
                             ('amp80-200_fix_forward','Retain F: completion'),('amp80-200_fix_reverse','Retain R: completion'),('redesign_forward','Retain F: redesign'),('redesign_reverse','Retain R: redesign')]:
            c=next((c for c in summary['paired_comparisons'] if c['context']==context and c['left_method']==primary and c['right_method']==direct and c['right_variant']=='full' and c['budget']==20 and c['endpoint']==RNA),None)
            if c:
                paired.append(tex(title)+f" & {c['mean_difference']:+.2f} & [{c['ci95'][0]:.2f}, {c['ci95'][1]:.2f}] & {c['n']}"+r' \\')
                data_values.append({'type':'contrast',**c})
        paired_table=r'''\begin{table}[t]\caption{Within-gene RNA-specific count differences: hybrid diffusion minus matched hybrid direct prediction, $K=20$.}\label{tab:paired}\centering\small\begin{tabular}{@{}lrrr@{}}\toprule Task & Difference & 95\% CI & $n$ \\\midrule
'''+ '\n'.join(paired)+r'\bottomrule\end{tabular}\end{table}'
        component=[]
        for variant in ['single_oligo','steps1','steps4','steps8','full']:
            g=gs.get(('amp80-200',primary,variant,20))
            if g:
                component.append(tex({'single_oligo':'Single-oligo','steps1':'1 step','steps4':'4 steps','steps8':'8 steps','full':'Full, 16 steps'}[variant])+' & '+value(g,RETAINED)+' & '+value(g,'RNA_pass_amplicon_distance')+f" ({g['RNA_pass_amplicon_distance_n']})"+' & '+value(g,'best_RNA_pass_temperature_penalty')+f" ({g['best_RNA_pass_temperature_penalty_n']})"+r' \\')
        component_table=r'''\begin{table}[t]\caption{Recorded guidance and denoising variants using the same hybrid diffusion checkpoint. Distance and Tm entries report their defined-gene counts in parentheses.}\label{tab:components}\centering\small\begin{tabular}{@{}lrrr@{}}\toprule Variant & Retained & Distance ($n$) & Tm penalty ($n$) \\\midrule
'''+ '\n'.join(component)+r'\bottomrule\end{tabular}\end{table}'
    else:
        quality_table=paired_table=component_table=r'\AuthorField{Machine-filled table reserved until declared analyses complete.}'

    main_figure_dir=(ROOT/args.summary).parent/'figures'
    external_figure_dir=(ROOT/args.external_summary).parent/'figures'
    figure_specs=[('01_method',ROOT/protocol['method_figure_directory']/(protocol['method_figure_basename']+'.pdf'),
        'Recorded hybrid model: local residual sequence encoding, conditional attention, absorbing-mask position diffusion and constrained primer-pair sampling.',None,106.8),
        ('02_quality',main_figure_dir/'01_candidate_quality.pdf',None,'01_candidate_quality',72),
        ('03_flexible',main_figure_dir/'03_flexible_design.pdf',None,'03_flexible_design',122),
        ('04_architecture',main_figure_dir/'02_architecture_comparison.pdf',None,'02_architecture_comparison',112)]
    if args.fifth_figure=='generalization':
        figure_specs.append(('05_generalization',external_figure_dir/'05_external_generalization.pdf',None,'05_external_generalization',120))
    else:
        figure_specs.append(('05_budget',main_figure_dir/'04_budget_and_diversity.pdf',None,'04_budget_and_diversity',61))
    main_prov=external_prov={}
    if summary:
        main_prov=json.loads((main_figure_dir/'provenance.json').read_text());external_prov=json.loads((external_figure_dir/'provenance.json').read_text())
        if main_prov['summary_sha256']!=main_sha or external_prov['summary_sha256']!=external_sha:
            raise ValueError('Figures do not match the completed current analyses')
    figure_tex=[]
    for name,source,caption,record_name,height in figure_specs:
        if record_name and not args.layout_preview:
            provenance=external_prov if name=='05_generalization' else main_prov
            record=next((r for r in provenance['generated'] if r['name']==record_name),None)
            if record is None:raise ValueError('Required result figure missing: '+record_name)
            expected=next(f['sha256'] for f in record['files'] if f['path'].endswith(source.name))
            if sha256(source)!=expected:raise ValueError('Figure changed: '+source.name)
            caption=record['caption'];alt=record['alt_text']
        else:
            alt='Recorded architecture and generation schematic.' if name=='01_method' else 'Result figure area reserved for completed analysis.'
        if name=='01_method' or not args.layout_preview:
            shutil.copyfile(source,out/'figures'/(name+'.pdf'))
            content=r'\includegraphics[width=\textwidth]{figures/'+name+'.pdf}'
            figure_records.append({'name':name,'source':str(source),'sha256':sha256(source)})
        else:
            content=r'\fbox{\parbox[c]['+str(height)+r'mm][c]{.97\textwidth}{\centering '+tex(name.replace('_',' '))+r'\\Result figure reserved; no unfinished numerical results.}}'
            caption='Reserved author result figure area.'
        figure_tex.append(r'\begin{figure*}[t]\centering'+content+r'\caption{'+tex(caption)+r'}\label{fig:'+name+r'}\figalttext{'+tex(alt)+r'}\end{figure*}')
    mode='layout preview without numerical results' if args.layout_preview else 'completed-source author result material'
    template=r'''% Author result worksheet; prose fields remain for the researcher.
\documentclass[numsec,webpdf,modern,large,namedate]{oup-authoring-template}
\usepackage{booktabs}
\input{glyphtounicode}\pdfgentounicode=1
\pdfminorversion=5
\emergencystretch=2em
\newcommand{\societylogo}{}
\newcommand{\AuthorField}[1]{\par\noindent\textit{#1}\par\vspace{4pt}\hrule\vspace{9pt}}
\hypersetup{pdftitle={DiffuPrimer author result worksheet},pdfauthor={Jingyue Li}}
\begin{document}
\journaltitle{Bioinformatics}\DOI{}\copyrightyear{2026}\pubyear{2026}\vol{}\issue{}\access{}
\appnotes{Author result worksheet}\firstpage{1}\lastpage{1}
\title[DiffuPrimer author result material]{DiffuPrimer: Conditional Discrete Diffusion for Diverse, Constraint-Aware qPCR Primer-Pair Design}
\author[1,$\ast$]{Jingyue Li}
\address[1]{\orgname{Sun Yat-sen University}}
\corresp[$\ast$]{Corresponding author. \href{mailto:18837614038@163.com}{18837614038@163.com}}
\abstract{\textbf{Motivation:} [Author text.]\\
\textbf{Results:} [Author text based on completed results.]\\
\textbf{Availability and Implementation:}\\
\textbf{Contact:} \href{mailto:18837614038@163.com}{18837614038@163.com}\\
\textbf{Supplementary Information:} [Author field.]}
\keywords{qPCR, primer-pair design, conditional discrete diffusion, attention}
\maketitle
\section{Introduction: author fields}
\AuthorField{Research problem, related work and motivation.}
\AuthorField{Three non-overlapping contributions: model formulation; unified flexible sampling; public-reference evidence.}
__FIGURE1__
\section{System and methods: recorded material}
\subsection{Architecture and training}
__ARCHITECTURE_TABLE__
\AuthorField{Author explanation of local sequence features and paired-position attention.}
\subsection{Public data and independent genes}
__DATASET_TABLE__
\AuthorField{Author account of data acquisition, positional supervision, homology-aware partitions and additional genes.}
\subsection{Position representation and diffusion equations}
\begin{equation}\begin{aligned}
x_0&\in\{0,F,R\}^{L},\\
0&\leq s_F<e_F\leq s_R<e_R\leq L,\\
\ell_{\min}&\leq e_F-s_F,e_R-s_R\leq\ell_{\max},\\
a_{\min}&\leq e_R-s_F\leq a_{\max}.
\end{aligned}\end{equation}
\begin{equation}\begin{aligned}
q(x_t^i\mid x_0^i)&=\bar\alpha_t\mathbf1[x_t^i=x_0^i]\\
&\quad+(1-\bar\alpha_t)\mathbf1[x_t^i=M],\\
\bar\alpha_t&=1-t/T,\qquad T=16.
\end{aligned}\end{equation}
\begin{equation}\mathcal L=\mathbb E_{t,x_t}\left[\frac{T}{tL}\sum_i\mathbf1[x_t^i=M]\big(-\log p_\theta(x_0^i\mid x_t,A,c,t)\big)\right].\end{equation}
\begin{equation}\begin{aligned}S_\theta(z)&=\sum_{i\in F(z)}[u_i(F)-u_i(0)]\\&\quad+\sum_{i\in R(z)}[u_i(R)-u_i(0)].\end{aligned}\end{equation}
\AuthorField{Author definitions of fixed-template sequence extraction, conditional input, anchors, compatible support and label revelation.}
\subsection{Conditions and assessment}
__CONDITIONS_TABLE__
\AuthorField{Author definitions of RNA specificity, genomic-product flags, candidate prefixes, diversity, Tm penalty and gene-level statistical units.}
\section{Results: source-filled evidence}
__QUALITY_TABLE__
__FIGURE2__
\AuthorField{Author interpretation of candidate quality and independent-gene comparison.}
__PAIRED_TABLE__
__FIGURE3__
\AuthorField{Author interpretation of length control, completion and changed-condition redesign.}
__FIGURE4__
__COMPONENT_TABLE__
\AuthorField{Author interpretation of architecture, guidance and denoising depth.}
__FIGURE5__
\section{Discussion: author fields}
\AuthorField{Research implications and practical use.}
\AuthorField{Short limitations paragraph.}
\section*{Availability}
Code URL: \rule{.60\columnwidth}{.3pt}\\
Data URL: \rule{.60\columnwidth}{.3pt}
\section*{Author declarations}
\AuthorField{Funding, conflicts, CRediT, ORCID and assistance disclosure.}
\section*{Reference material}
\nocite{primer3,primerbank,d3pm,mdlm,mmseqs,bowtie}
\bibliographystyle{oup-abbrvnat}\bibliography{references}
\end{document}
'''
    replacements={'__ARCHITECTURE_TABLE__':architecture_table,'__DATASET_TABLE__':dataset_table,'__CONDITIONS_TABLE__':conditions_table,
                  '__QUALITY_TABLE__':quality_table,'__PAIRED_TABLE__':paired_table,'__COMPONENT_TABLE__':component_table}
    replacements.update({'__FIGURE'+str(i+1)+'__':f for i,f in enumerate(figure_tex)})
    for k,v in replacements.items():template=template.replace(k,v)
    (out/'author_results.tex').write_text(template)
    write_json(out/'reported_values.json',{'at':utc_now(),'mode':mode,'main_summary_sha256':main_sha,'external_summary_sha256':external_sha,'values':data_values})
    write_json(out/'source_manifest.json',{'at':utc_now(),'mode':mode,'protocol_path':str(protocol_path),'protocol_sha256':sha256(protocol_path),
        'main_summary_sha256':main_sha,'external_summary_sha256':external_sha,'primary_method':primary,'direct_method':direct,
        'method_architecture':architecture,'datasets_sha256':sha256(dataset_path),'figures':figure_records,
        'author':{'name':'Jingyue Li','institution':'Sun Yat-sen University','contact':'18837614038@163.com'},
        'code_url':'','data_url':'','author_prose_generated':False,'old_working_copy_modified':False})
    pages=None
    if not args.no_compile:
        for i,command in enumerate([
            ['pdflatex','-interaction=nonstopmode','-halt-on-error','author_results.tex'],['bibtex','author_results'],
            ['pdflatex','-interaction=nonstopmode','-halt-on-error','author_results.tex'],['pdflatex','-interaction=nonstopmode','-halt-on-error','author_results.tex']]):
            proc=subprocess.run(command,cwd=out,capture_output=True,text=True)
            (out/f'build_{i}.log').write_text(proc.stdout+'\n'+proc.stderr)
            if proc.returncode:raise RuntimeError('Worksheet typesetting failed; see '+str(out/f'build_{i}.log'))
        proc=subprocess.run(['pdfinfo','author_results.pdf'],cwd=out,capture_output=True,text=True,check=True)
        pages=int(re.search(r'^Pages:\s+(\d+)',proc.stdout,re.M).group(1))
        # The unmodified OUP class counts the final text page before deferred wide
        # floats. Populate its public lastpage field with the measured PDF count.
        final_tex=(out/'author_results.tex').read_text().replace(r'\firstpage{1}\lastpage{1}',r'\firstpage{1}\lastpage{'+str(pages)+'}')
        (out/'author_results.tex').write_text(final_tex)
        proc=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','author_results.tex'],cwd=out,capture_output=True,text=True)
        (out/'build_measured_lastpage.log').write_text(proc.stdout+'\n'+proc.stderr)
        if proc.returncode:raise RuntimeError('Measured-page metadata compilation failed')
        measured=subprocess.run(['pdfinfo','author_results.pdf'],cwd=out,capture_output=True,text=True,check=True)
        final_pages=int(re.search(r'^Pages:\s+(\d+)',measured.stdout,re.M).group(1))
        if final_pages!=pages:raise ValueError('Pagination changed while recording measured page count')
        write_json(out/'layout_report.json',{'at':utc_now(),'mode':mode,'pages':pages,'target_max_pages':7,'within_target':pages<=7,
            'template':'numsec,webpdf,modern,large,namedate','figures':5,'numerical_result_instances':len(data_values),
            'pdf_sha256':sha256(out/'author_results.pdf'),'note':'Measured author worksheet; researcher-written prose may alter final manuscript pagination'})
        if pages>7:raise ValueError(f'Author worksheet is {pages} template pages; consolidate layout before delivery')
    print(json.dumps({'output':str(out),'mode':mode,'pages':pages,'figures':5,'numerical_result_instances':len(data_values)}))


if __name__=='__main__':main()
