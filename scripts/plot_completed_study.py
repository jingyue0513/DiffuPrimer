#!/usr/bin/env python3
"""Export scientific figures from completed, matched reference-audited runs."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from diffuprimer.io import read_jsonl,sha256,utc_now,write_json


root=Path(__file__).resolve().parents[1]
base=root/'results/temperature_tuning'
files=[base/'selected_test_summary.json',base/'selected_conditions_summary.json']
summaries=[json.loads(p.read_text()) for p in files]
guide_file=root/'results/guide_ablation_summary.json'
guide_summary=json.loads(guide_file.read_text())
selection_path=base/'selection.json'
selection=json.loads(selection_path.read_text())
for summary in summaries:
    if summary['validation_selection_sha256'] != sha256(selection_path) or summary.get('RNA_only'):
        raise ValueError('Figures require the selected complete RNA/genome audit')
for summary in summaries+[guide_summary]:
    for item in summary['inputs']:
        directory=Path(item['directory'])
        if sha256(directory/'completion.json') != item['completion_sha256'] or sha256(directory/'specificity_bowtie.jsonl') != item['specificity_sha256']:
            raise ValueError('Figure inputs changed after statistical analysis')
records=[read_jsonl(p.with_suffix('.per_gene.jsonl')) for p in files]
methods=['primer3','uniform_pair_sampler','rule_sampler','direct','diffusion']
labels={'primer3':'Primer3','uniform_pair_sampler':'Uniform pairs','rule_sampler':'Forward-first','direct':'Direct','diffusion':'Diffusion'}
colors={'primer3':'#db9c27','uniform_pair_sampler':'#37926c','rule_sampler':'#a4a7a9','direct':'#3d7aa8','diffusion':'#855ab0'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
output=root/'results/figures/selected_study';output.mkdir(parents=True,exist_ok=True)
values=[]


def endpoint(rows,method,context,variant,budget,key,study=None):
    own=[r for r in rows if r['method']==method and r['context']==context and r['variant']==variant and r['budget']==budget and (study is None or r['study']==study)]
    if len(own)!=75 or len({r['gene_id'] for r in own})!=75:
        raise ValueError(f'Incomplete figure cohort: {method}/{context}/{variant}/{budget}')
    expected=[20261013] if method=='primer3' else [20261013,20261014,20261015]
    if any(r['seeds']!=expected for r in own):
        raise ValueError('Figure seed denominator differs')
    data=np.asarray([r[key] for r in own if r[key] is not None],float)
    if not len(data):return None
    rng=np.random.default_rng(20261013)
    bootstrap=np.concatenate([data[rng.integers(0,len(data),(1000,len(data)))].mean(1) for _ in range(10)])
    mean=float(data.mean());lo,hi=np.quantile(bootstrap,[.025,.975])
    values.append({'study':own[0]['study'],'method':method,'context':context,'variant':variant,'budget':budget,'endpoint':key,
        'mean':mean,'ci95_lower':float(lo),'ci95_upper':float(hi),'genes':len(data),'total_target_genes':75,'seed_count':len(expected)})
    return mean,float(lo),float(hi),len(data)


def save(fig,name):
    fig.savefig(output/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(output/(name+'.png'),dpi=220,bbox_inches='tight')
    plt.close(fig)


def interval(ax,x,lower,upper,width=.07):
    ax.vlines(x,lower,upper,color='black',lw=1)
    ax.hlines([lower,upper],x-width,x+width,color='black',lw=1)


fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
fields=[('unique_RNA_pass_pairs','Unique RNA-specific pairs',None),('two_RNA_pass_pairs','Targets with at least two pairs',(0,1)),
    ('RNA_pass_amplicon_distance','Amplicon Jaccard distance',(0,1)),('best_RNA_pass_temperature_penalty','Best paired Tm penalty',None)]
for ax,(key,label,limits) in zip(axes.flat,fields):
    for index,method in enumerate(methods):
        estimate=endpoint(records[0],method,'amp80-200','primary',20,key)
        if estimate is None:continue
        mean,lo,hi,n=estimate
        ax.bar(index,mean,color=colors[method],width=.65)
        interval(ax,index,lo,hi)
        if n!=75:ax.annotate(f'n={n}',(index,hi),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
    ax.set_xticks(range(5),[labels[m] for m in methods],rotation=20,ha='right')
    ax.set_ylabel(label)
    if limits:ax.set_ylim(limits[0],limits[1]*1.12)
    else:ax.set_ylim(bottom=0)
save(fig,'primary_RNA_endpoints')
fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for ax,(key,label) in zip(axes,[('unique_RNA_pass_pairs_without_genomic_flag','RNA-pass pairs without genomic product flag'),
                              ('two_RNA_pass_pairs_without_genomic_flag','Targets with at least two such pairs')]):
    for index,method in enumerate(methods):
        mean,lo,hi,n=endpoint(records[0],method,'amp80-200','primary',20,key)
        ax.bar(index,mean,color=colors[method],width=.65)
        interval(ax,index,lo,hi)
    ax.set_xticks(range(5),[labels[m] for m in methods],rotation=20,ha='right')
    ax.set_ylabel(label);ax.set_ylim(bottom=0)
fig.suptitle('Secondary genomic-product screen',fontsize=11)
save(fig,'genomic_product_screen')
exon_file=root/'results/exon_specificity_association.json'
exon_summary=json.loads(exon_file.read_text())
if sha256(root/'results/selected_exon_context_audit.per_candidate.jsonl')!=exon_summary['exon_records_sha256']:
    raise ValueError('Exon association annotation inputs changed')
exon_rows=read_jsonl(exon_file.with_suffix('.per_gene.jsonl'))
fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
bottom=np.zeros(len(methods))
for category,color,label in [('primer_spans_junction','#3d7aa8','Primer spans junction'),
                              ('amplicon_crosses_junction_only','#37926c','Amplicon crosses junction; primers within exons'),
                              ('terminal_polyA','#db9c27','Terminal poly(A)')]:
    estimates=[]
    for method in methods:
        own=[r for r in exon_rows if r['method']==method]
        if len(own)!=75:raise ValueError('Exon figure gene denominator differs')
        data=np.array([r['categories'][category]['RNA_pass_without_genomic_flag'] for r in own])
        estimates.append(float(data.mean()))
    ax.bar(range(len(methods)),estimates,bottom=bottom,width=.65,color=color,label=label)
    bottom+=estimates
for index,method in enumerate(methods):
    _,lo,hi,n=endpoint(records[0],method,'amp80-200','primary',20,'unique_RNA_pass_pairs_without_genomic_flag')
    interval(ax,index,lo,hi)
ax.set_xticks(range(len(methods)),[labels[m] for m in methods],rotation=20,ha='right')
ax.set_ylabel('Unique RNA-pass pairs without genomic product flag')
ax.set_ylim(bottom=0);ax.legend(frameon=False,fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.02))
save(fig,'exon_geometry_partition')
fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
for ax,(key,label) in zip(axes,[('unique_RNA_pass_pairs','Unique RNA-specific pairs'),('two_RNA_pass_pairs','Targets with at least two pairs')]):
    for method in methods:
        estimates=[endpoint(records[0],method,'amp80-200','primary',k,key) for k in [1,2,5,10,20]]
        y=np.asarray([e[0] for e in estimates]);lo=[e[1] for e in estimates];hi=[e[2] for e in estimates]
        ax.plot([1,2,5,10,20],y,'o-',label=labels[method],color=colors[method],lw=1.4,ms=3)
        ax.fill_between([1,2,5,10,20],lo,hi,color=colors[method],alpha=.12)
    ax.set_xlabel('Proposal budget');ax.set_ylabel(label);ax.set_ylim(bottom=0);ax.set_xticks([1,2,5,10,20])
axes[0].legend(frameon=False,fontsize=8)
save(fig,'proposal_budget')
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
contexts=['amp80-200','amp80-120','amp121-160','amp161-200','amp80-200_fix_forward','amp80-200_fix_reverse']
for ax,context in zip(axes.flat,contexts):
    for index,method in enumerate(methods):
        mean,lo,hi,n=endpoint(records[1],method,context,'primary',20,'unique_RNA_pass_pairs')
        ax.bar(index,mean,color=colors[method],width=.65)
        interval(ax,index,lo,hi)
    ax.set_title(context.replace('amp','').replace('_fix_forward','; fixed forward').replace('_fix_reverse','; fixed reverse'))
    ax.set_xticks(range(5),[labels[m] for m in methods],rotation=25,ha='right');ax.set_ylim(0,21)
    ax.set_ylabel('Unique RNA-specific pairs')
save(fig,'condition_and_completion')
fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
for ax,(key,label) in zip(axes,[('unique_RNA_pass_pairs','Unique RNA-specific pairs'),('RNA_pass_amplicon_distance','Amplicon Jaccard distance'),('best_RNA_pass_temperature_penalty','Best paired Tm penalty')]):
    estimates=[endpoint(records[1],'diffusion','amp80-200',v,20,key) for v in ['steps1','steps4','steps8','primary']]
    mean=[e[0] for e in estimates];lo=[e[1] for e in estimates];hi=[e[2] for e in estimates]
    ax.plot([1,4,8,16],mean,'o-',color=colors['diffusion']);ax.fill_between([1,4,8,16],lo,hi,color=colors['diffusion'],alpha=.15)
    ax.set_xlabel('Denoising steps');ax.set_ylabel(label);ax.set_xticks([1,4,8,16]);ax.set_ylim(bottom=0)
save(fig,'denoising_steps')
guide_rows=read_jsonl(guide_file.with_suffix('.per_gene.jsonl'))
fig,axes=plt.subplots(1,3,figsize=(11,3.8),layout='constrained')
guide_methods=['rule_sampler','direct','diffusion']
for ax,(key,label) in zip(axes,[('unique_RNA_pass_pairs','Unique RNA-pass pairs'),
                              ('RNA_pass_amplicon_distance','Amplicon Jaccard distance'),
                              ('best_RNA_pass_temperature_penalty','Best paired Tm penalty')]):
    for offset,study,label_study,alpha in [(-.18,'original','Single-oligo guidance',.4),(.18,'joint_tm','Joint Tm guidance',1)]:
        for index,method in enumerate(guide_methods):
            mean,lo,hi,n=endpoint(guide_rows,method,'amp80-200','primary',20,key,study=study)
            ax.bar(index+offset,mean,width=.32,color=colors[method],alpha=alpha,label=label_study if index==0 else None)
            interval(ax,index+offset,lo,hi,width=.035)
    mean,lo,hi,n=endpoint(guide_rows,'primer3','amp80-200','primary',20,key,study='original')
    ax.axhline(mean,color=colors['primer3'],lw=1.2,linestyle='--',label='Primer3')
    ax.set_xticks(range(3),[labels[m] for m in guide_methods]);ax.set_ylabel(label);ax.set_ylim(bottom=0)
axes[0].legend(frameon=False,fontsize=7,loc='lower center',bbox_to_anchor=(.5,1.02))
fig.suptitle('Guidance mechanism at fixed temperature 1',fontsize=11)
save(fig,'guidance_mechanism')
duration_file=root/'experiments/epoch100/results/duration_reference_comparison.json'
duration=json.loads(duration_file.read_text())
duration_rows=[]
for item in duration['inputs']:
    path=Path(item['summary'])
    if sha256(path)!=item['summary_sha256'] or sha256(path.with_suffix('.per_gene.jsonl'))!=item['per_gene_sha256']:
        raise ValueError('Duration figure inputs changed')
    own=[r for r in read_jsonl(path.with_suffix('.per_gene.jsonl')) if r['study']=='original' and r['context']=='amp80-200' and r['variant']=='primary' and r['budget']==20 and r['method'] in ['direct','diffusion']]
    duration_rows.extend([{**r,'study':f"duration_{item['split']}_epochs{item['epochs']}"} for r in own])
fig,axes=plt.subplots(2,3,figsize=(11,7),layout='constrained')
for row,split in enumerate(['validation','test']):
    for ax,(key,label) in zip(axes[row],[('unique_RNA_pass_pairs','Unique RNA-pass pairs'),
                                      ('RNA_pass_amplicon_distance','Amplicon Jaccard distance'),
                                      ('best_RNA_pass_temperature_penalty','Best paired Tm penalty')]):
        for offset,epochs,alpha in [(-.18,30,.4),(.18,100,1)]:
            for index,method in enumerate(['direct','diffusion']):
                mean,lo,hi,n=endpoint(duration_rows,method,'amp80-200','primary',20,key,study=f'duration_{split}_epochs{epochs}')
                ax.bar(index+offset,mean,width=.32,color=colors[method],alpha=alpha,label=f'{epochs} epochs' if index==0 else None)
                interval(ax,index+offset,lo,hi,width=.035)
                if n!=75:ax.annotate(f'n={n}',(index+offset,hi),xytext=(0,4),textcoords='offset points',ha='center',fontsize=7)
        ax.set_xticks([0,1],['Direct','Diffusion']);ax.set_ylabel(label);ax.set_title(split.capitalize());ax.set_ylim(bottom=0)
axes[0,0].legend(frameon=False,fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.12))
fig.suptitle('Training duration at fixed original decoding settings',fontsize=11)
save(fig,'training_duration_endpoints')
runtime_file=root/'results/runtime_benchmark/summary.json'
runtime=json.loads(runtime_file.read_text())
runtime_records=read_jsonl(root/'results/runtime_benchmark/per_gene.jsonl')
if runtime['targets']!=75 or runtime['repetitions']!=3 or not runtime['all_candidate_coordinates_match']:
    raise ValueError('Runtime figure requires all completed measurements')
runtime_values=[]
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for method in methods:
    own=[]
    genes=sorted({r['gene_id'] for r in runtime_records if r['method']==method})
    if len(genes)!=75:raise ValueError('Runtime figure target denominator differs')
    for gene in genes:
        repeated=[r for r in runtime_records if r['method']==method and r['gene_id']==gene]
        if sorted(r['repetition'] for r in repeated)!=[0,1,2]:
            raise ValueError('Runtime figure repeated measurements differ')
        record={'method':method,'gene_id':gene,'template_length':repeated[0]['template_length'],
            'repetitions':3,**{field:float(np.mean([r[field] for r in repeated]))
            for field in ['generation_seconds','assessment_seconds','shared_preparation_seconds','total_seconds_with_preparation']}}
        runtime_values.append(record);own.append(record)
    for ax,field in zip(axes,['generation_seconds','total_seconds_with_preparation']):
        ax.scatter([r['template_length'] for r in own],[r[field] for r in own],s=12,alpha=.65,
            color=colors[method],label=labels[method])
for ax,label in zip(axes,['Generation time (s)','Generation + assessment + preparation (s)']):
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Full transcript length (nt)');ax.set_ylabel(label)
axes[0].legend(frameon=False,fontsize=8)
save(fig,'generation_runtime')
with (output/'runtime_values.csv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(runtime_values[0]));writer.writeheader();writer.writerows(runtime_values)
with (output/'figure_values.csv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)
write_json(output/'provenance.json',{'at':utc_now(),'source_summaries':{str(p):sha256(p) for p in files+[guide_file,duration_file,exon_file,runtime_file]},
    'source_records':{str(p):sha256(p) for p in [*[p.with_suffix('.per_gene.jsonl') for p in files+[guide_file,exon_file]],
        root/'results/runtime_benchmark/per_gene.jsonl',root/'results/selected_exon_context_audit.per_candidate.jsonl']},
    'selection_sha256':sha256(selection_path),'script_sha256':sha256(Path(__file__)),
    'statistical_unit':'Target gene after averaging model seeds; 10000 gene bootstrap replicates; feature-specific available-gene counts retained',
    'temperature':{m:v['temperature'] for m,v in selection['selected'].items()},'manuscript_work_paused':True})
print(output,flush=True)
