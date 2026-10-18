#!/usr/bin/env python3
"""Redraw completed CNN research results, preserving every source and endpoint."""
from __future__ import annotations
import hashlib, json
from collections import Counter, defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'paper/figures/cnn_publication_results'
TEAL='#248A75'; DIRECT='#7D98A7'; OCHRE='#C5A05B'; CORAL='#D36E51'; PURPLE='#9B88AC'; INK='#293F4A'; GRID='#E4EAEA'
COLORS={'cnn_diffusion':TEAL,'cnn_direct':DIRECT,'primer3':OCHRE,'qprimer':CORAL,'primercast':PURPLE,'attention_rope_diffusion':PURPLE,'attention_hybrid_diffusion':'#5B95B1'}
LABELS={'cnn_diffusion':'CNN\ndiffusion','cnn_direct':'CNN\ndirect','primer3':'Primer3','qprimer':'qPrimer\ndesign','primercast':'PrimerCast','attention_rope_diffusion':'Attention','attention_hybrid_diffusion':'CNN +\nattention'}
COUNT='unique_RNA_pass_pairs_without_genomic_flag'; FRACTION='RNA_no_genome_fraction'; DIVERSITY='RNA_pass_amplicon_distance'
VALUES=[]; SOURCES={}; DERIVED={}; FIGURES=[]
PLOT_RATIONALE=[
 dict(figure='fig2_candidate_quality',panel='A/B',question='How do candidate supply and retained quality compare across methods in each independent cohort?',endpoint='Gene-mean retained pair count; conditional gene-mean retained/RNA-specific fraction',form='Cohort-offset mean points with existing 95% gene-bootstrap intervals',reason='Directly encodes estimates and uncertainty without adding filled bar area; marker shape identifies the independent cohort.'),
 dict(figure='fig2_candidate_quality',panel='C/D',question='Is the mean advantage spread across genes, and which genes favor either method?',endpoint='Per-gene retained count, CNN diffusion versus Primer3',form='Integer paired scatter with equality line and multiplicity key',reason='Preserves all 75 paired gene outcomes; overlapping integer coordinates are shown as scaled circles rather than hidden or jittered observations.'),
 dict(figure='fig3_flexible_design',panel='A',question='What amplicon lengths does each requested length condition actually produce?',endpoint='Exact integer amplicon lengths of unique feasible RNA-specific candidate pairs',form='Empirical cumulative distribution with shaded requested ranges',reason='Shows the observed discrete distributions without kernel smoothing, arbitrary histogram bins or implying continuous density. Each condition is normalized separately and candidate-weighted.'),
 dict(figure='fig3_flexible_design',panel='B/C',question='How many screened alternatives remain when either primer is fixed or the requested length changes?',endpoint='Mean retained pairs across all 75 genes',form='Grouped bars with existing 95% gene-bootstrap intervals',reason='The four methods and two fixed sides are finite task categories; a common zero baseline and common vertical scale make their counts comparable.'),
 dict(figure='fig3_flexible_design',panel='D/E',question='How do retained candidate supply and positional diversity change with proposal budget?',endpoint='Retained pair count; mean RNA-specific amplicon Jaccard distance',form='Recorded-prefix budget curves with existing uncertainty bands',reason='Proposal budgets are ordered and each point is an actual prefix, so curves expose the response to increasing budget. Diversity compares the matched predictors, using defined endpoint denominators.'),
 dict(figure='fig3_flexible_design',panel='F',question='Which primer stays fixed and where is the alternative partner placed on one transcript?',endpoint='Actual recorded F/R intervals, length and genomic-screen flag',form='Transcript coordinate tracks with directed primer arrows and anchor marks',reason='Position tracks display the actual geometry and preserved side directly; candidate-level genomic flags distinguish an illustrative RNA-specific candidate from a fully retained pair.'),
 dict(figure='fig5_mechanism',panel='A',question='What is the mean gain from diffusion relative to the matched CNN direct predictor?',endpoint='Gene-paired mean difference in retained candidate count',form='Horizontal paired mean differences with existing 95% paired bootstrap intervals and zero reference',reason='The inferential target is the paired mean effect, not a smoothed distribution of discrete gene differences; all gene differences remain in the source export.'),
 dict(figure='fig5_mechanism',panel='B',question='How do sampling constraints affect passage through three common assessment stages?',endpoint='Mean unique feasible, RNA-specific and retained pair counts',form='Annotated variant-by-assessment-stage heatmap with a common 0–20 scale',reason='The row-by-column matrix has a meaningful factorial structure; equal color limits and cell values support comparisons across both dimensions.'),
 dict(figure='fig5_mechanism',panel='C',question='How does retained supply vary with the number of reverse-diffusion steps?',endpoint='Mean retained pairs per gene',form='Measured schedule curve with existing 95% gene-bootstrap band',reason='The steps are ordered settings and the curve shows the observed schedule response without asserting monotonic improvement.'),
 dict(figure='fig5_mechanism',panel='D',question='Which diffusion denoiser backbone yields more screened candidates in each cohort?',endpoint='Mean retained pairs per gene for CNN, hybrid and attention diffusion',form='Grouped architectural mean bars with existing 95% intervals',reason='Three discrete backbones are compared at the same candidate budget; paired cohort bars retain cohort separation and the count zero baseline.'),
]

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(p):
 p=Path(p); digest=sha(p); expected=globals().get('HASH_EXPECTATIONS',{}).get(str(p)); assert expected is None or digest==expected, f'Source hash changed: {p}'; SOURCES[str(p.relative_to(ROOT))]=digest; return p

def read_jsonl(p): return [json.loads(x) for x in source(p).read_text().splitlines() if x.strip()]
def load():
 summaries={}; pergenes={}
 for cohort,p in [('original',ROOT/'experiments/attention/results/study_summary.json'),('external',ROOT/'experiments/attention_external/results/study_summary.json')]:
  s=json.loads(source(p).read_text()); assert s['status']=='complete'; assert s['total_test_genes']==75
  summaries[cohort]=s
  pergenes[cohort]=read_jsonl(p.with_suffix('.per_gene.jsonl'))
 return summaries,pergenes
S,P=load()
G={cohort:{(g['context'],g['method'],g['variant'],g['budget']):g for g in s['summaries']} for cohort,s in S.items()}
HASH_EXPECTATIONS={str(Path(i['directory'])/name):i[name+'_sha256'] for summary in S.values() for i in summary['inputs'] for name in ['candidates.jsonl','specificity_bowtie.jsonl'] if name+'_sha256' in i}

def value(method,key=COUNT,cohort='original',context='amp80-200',variant='full',budget=20):
 g=G[cohort][context,method,variant,budget]; scale=100 if key==FRACTION else 1
 v=float(g[key])*scale; ci=np.array(g[key+'_ci95'])*scale
 VALUES.append(dict(cohort=cohort,context=context,method=method,variant=variant,budget=budget,endpoint=key,mean=v,ci95=ci.tolist(),n=g[key+'_n'],scale=scale,run=g['run']))
 return v,ci

def records(method,cohort='original',context='amp80-200',variant='full',budget=20):
 return {r['gene_id']:r for r in P[cohort] if (r['method'],r['context'],r['variant'],r['budget'])==(method,context,variant,budget)}

def panel(ax,letter,title):
 ax.set_title(title,loc='left',fontsize=8.1,pad=11,color=INK)
 ax.text(-.12,1.09,letter,transform=ax.transAxes,fontsize=10,fontweight='bold',color=INK)
 ax.set_axisbelow(True); ax.grid(axis='y',color=GRID,lw=.5); ax.tick_params(length=2.5,width=.6)

def errbars(ax,x,ys,cis,**kw):
 yerr=np.array([[y-c[0],c[1]-y] for y,c in zip(ys,cis)]).T
 ax.errorbar(x,ys,yerr=yerr,fmt='none',color=INK,capsize=2,lw=.65,**kw)

def finish(fig,name,caption,alt,height_mm):
 fig.patch.set_facecolor('white'); fig.canvas.draw()
 # Never crop; actual output canvas remains exactly 178 mm wide.
 files=[]
 for ext in ['pdf','svg','png','tiff']:
  path=OUT/f'{name}.{ext}'; fig.savefig(path,dpi=350,facecolor='white')
  if ext in ['png','tiff']:
   with Image.open(path) as raster: rgb=raster.convert('RGB')
   rgb.save(path,dpi=(350,350),**({'compression':'tiff_lzw'} if ext=='tiff' else {}))
  files.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
 fig.savefig(OUT/f'{name}_print_review.png',dpi=100,facecolor='white')
 start=sum(len(f.get('value_indices',[])) for f in FIGURES)
 FIGURES.append(dict(name=name,value_indices=list(range(start,len(VALUES))),derived_value_keys=[k for k in DERIVED if k.startswith(name.split('_')[0].replace('fig','figure'))],plot_rationale=[r for r in PLOT_RATIONALE if r['figure']==name],width_mm=178,height_mm=height_mm,caption=caption,alt_text=alt,files=files))
 plt.close(fig)

def cohort_comparison(ax,key):
 methods=['cnn_diffusion','cnn_direct','primer3','qprimer','primercast']; x=np.arange(len(methods)); w=.23
 for j,cohort in enumerate(['original','external']):
  ys=[];cis=[]
  for m in methods:
   y,ci=value(m,key,cohort); ys.append(y);cis.append(ci)
  for i,m in enumerate(methods):
   y,ci=ys[i],cis[i]
   ax.errorbar(x[i]+(j-.5)*w,y,yerr=[[y-ci[0]],[ci[1]-y]],fmt='o' if j==0 else 's',color=COLORS[m],markerfacecolor=COLORS[m] if j==0 else 'white',markeredgewidth=.9,ms=4.1,capsize=2,lw=1,zorder=3)
 ax.set_xticks(x,[LABELS[m] for m in methods]); ax.set_ylim(bottom=0)
 ax.set_ylabel('Retained pairs per gene' if key==COUNT else 'Retained / RNA-specific pairs (%)')
 ax.legend(handles=[Line2D([],[],linestyle='none',marker='o',color='#617D79',markersize=4,label='Original test (75 genes)'),Line2D([],[],linestyle='none',marker='s',color='#617D79',markerfacecolor='white',markersize=4,label='Independent test (75 genes)')],frameon=False,fontsize=7,loc='upper right')

def paired_scatter(ax,cohort):
 left=records('cnn_diffusion',cohort);right=records('primer3',cohort); assert set(left)==set(right)
 xy=[(right[g][COUNT],left[g][COUNT]) for g in sorted(left)]
 multiplicities=Counter(xy)
 for (x,y),n in sorted(multiplicities.items()):
  ax.scatter(x,y,s=16+14*n,c=TEAL if y>x else CORAL if y<x else '#A4B0B4',alpha=.73,edgecolors='white',lw=.5,zorder=3)
 ax.plot([0,20],[0,20],color='#8F9DA3',lw=.8,ls='--');ax.set(xlim=(-.8,20.8),ylim=(-.8,20.8),xticks=[0,5,10,15,20],yticks=[0,5,10,15,20],xlabel='Primer3 retained pairs',ylabel='CNN diffusion retained pairs')
 wins=sum(y>x for x,y in xy);ties=sum(y==x for x,y in xy);losses=sum(y<x for x,y in xy)
 ax.text(.03,.97,f'CNN > Primer3: {wins}/75 genes\nEqual: {ties}  ·  CNN < Primer3: {losses}',transform=ax.transAxes,va='top',fontsize=7,color=INK,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
 DERIVED[f'figure2_{cohort}_paired_scatter']=dict(gene_ids=sorted(left),xy=xy,wins=wins,ties=ties,losses=losses,bubble_size='16 + 14 * number of genes at that integer coordinate')


def figure2():
 fig,axes=plt.subplots(2,2,figsize=(178/25.4,149/25.4));fig.subplots_adjust(left=.095,right=.98,top=.93,bottom=.155,wspace=.34,hspace=.53)
 panel(axes[0,0],'A','Reference-screened candidate supply');cohort_comparison(axes[0,0],COUNT)
 panel(axes[0,1],'B','Candidate quality');cohort_comparison(axes[0,1],FRACTION)
 for ax,c,l in zip(axes[1],['original','external'],'CD'):
  panel(ax,l,'Original test · 75 genes' if c=='original' else 'Independent test · 75 genes');paired_scatter(ax,c)
 multiplicity=[axes[1,0].scatter([],[],s=16+14*n,c='#A4B0B4',edgecolors='white',lw=.5) for n in [1,5,10]]
 fig.legend(multiplicity,['1 gene','5 genes','10 genes'],title='Genes sharing one integer coordinate',title_fontsize=7,fontsize=7,loc='lower center',bbox_to_anchor=(.55,.004),ncol=3,frameon=False,handletextpad=.6,columnspacing=2)
 finish(fig,'fig2_candidate_quality','Candidate quantity and quality under a shared proposal budget K = 20 and shared full-reference assessment. (A) Unique RNA-specific pairs with no predicted genomic product, averaged over all 75 genes in each cohort. (B) Per-gene ratio of these retained pairs to RNA-specific pairs; genes with no RNA-specific pairs are excluded only from this ratio. Points and whiskers show means and the existing 95% gene-bootstrap intervals; filled circles and open squares identify separate original and independent cohorts. Cohorts are never pooled. (C,D) Every test gene is shown in a CNN diffusion versus Primer3 comparison. One circle represents all genes sharing each integer pair of counts; the size key indicates this multiplicity, and the dashed line indicates equality. Circle area is 16 + 14n points squared for n overlapping genes. qPrimer denotes the upstream design module with common assessment. PrimerCast uses native quick-design proposals assessed under the same exact-template pair task. Retained means computational reference-screening outcomes, not measured amplification.','Four panels compare retained candidate count and retained fraction using mean points and confidence intervals across five methods and two independent cohorts, followed by per-gene CNN versus Primer3 scatter plots with a multiplicity key.',149)


def native(method,context='amp80-200',variant='full'):
 directory=Path(G['original'][context,method,variant,20]['run']); candidates=read_jsonl(directory/'candidates.jsonl');screens={r['candidate_id']:r for r in read_jsonl(directory/'specificity_bowtie.jsonl')}
 by=defaultdict(list)
 for c in candidates: by[c['gene_id']].append(c)
 output={}
 for gene,cs in by.items():
  seen=set(); chosen=[]
  for c in sorted(cs,key=lambda c:c['proposal_index'])[:20]:
   pair=(c.get('forward_primer'),c.get('reverse_primer'))
   sc=screens[c['candidate_id']]
   if c['assessment']['computational_feasible'] and sc['gene_specificity_screen_pass'] and pair not in seen:
    chosen.append(dict(c,screen=sc)); seen.add(pair)
  output[gene]=chosen
 return output

def task_bars(ax,contexts):
 methods=['cnn_diffusion','cnn_direct','primer3','qprimer']; x=np.arange(len(contexts));w=.18
 for j,m in enumerate(methods):
  ys=[];cis=[]
  for ctx in contexts:
   y,ci=value(m,context=ctx);ys.append(y);cis.append(ci)
  ax.bar(x+(j-1.5)*w,ys,w,color=COLORS[m],zorder=3);errbars(ax,x+(j-1.5)*w,ys,cis)
 ax.set_xticks(x,['Retain F','Retain R']);ax.set_ylabel('Retained pairs per gene');ax.set_ylim(0,6.6)

def curve(ax,key,methods):
 for m in methods:
  budgets=[1,2,5,10,20];ys=[];cis=[]
  for b in budgets:
   g=G['original']['amp80-200',m,'full',b]
   if g[key] is None: continue
   y,ci=value(m,key,budget=b);ys.append(y);cis.append(ci)
  xs=[b for b in budgets if G['original']['amp80-200',m,'full',b][key] is not None]
  ax.plot(xs,ys,'o-',color=COLORS[m],lw=1.5,ms=3,label=LABELS[m].replace('\n',' '));ax.fill_between(xs,[c[0] for c in cis],[c[1] for c in cis],color=COLORS[m],alpha=.10,lw=0)
 ax.set_xlabel('Proposal budget K');ax.set_xticks([1,5,10,20]);ax.set_ylim(bottom=0)


def footprint(ax):
 contexts=['amp80-200_fix_forward','amp80-200_fix_reverse','redesign_forward','redesign_reverse']; tasks={c:native('cnn_diffusion',c) for c in contexts};taskledger={r['gene_id']:r for r in read_jsonl(ROOT/'experiments/final_followup/redesign_tasks.jsonl')}
 # Deterministic illustration criterion: first numerical gene ID with any
 # RNA-specific candidate for each supported mode; first proposal in each.
 # No comparative performance ranking or genomic-retention filtering is used.
 common=set(taskledger)
 for cs in tasks.values(): common &= {g for g,rows in cs.items() if rows}
 gene=sorted(common,key=int)[0];anchor=taskledger[gene]['anchor'];selected=[rows[gene][0] for rows in tasks.values()]
 target=records('cnn_diffusion')[gene]
 for i,p in enumerate(selected):
  side='forward' if i in [0,2] else 'reverse'
  assert all(p[side+'_'+suffix]==anchor[side+'_'+suffix] for suffix in ['start','end','primer']), f'Anchor changed for {gene}: {contexts[i]}'
 positions=[anchor]+selected;lo=min(p['forward_start'] for p in positions)-8;hi=max(p['reverse_end'] for p in positions)+8
 for i,p in enumerate(positions):
  y=4-i;ax.plot([lo,hi],[y,y],color='#D6DEDF',lw=1)
  fixedside=None if i==0 else ('forward' if i in [1,3] else 'reverse')
  for side,color in [('forward',TEAL),('reverse',CORAL)]:
   start,end=p[side+'_start'],p[side+'_end'];ax.plot([start,end],[y,y],lw=2.4,color=color,solid_capstyle='butt',zorder=3)
   ax.annotate('',xy=(end if side=='forward' else start,y),xytext=(start if side=='forward' else end,y),arrowprops=dict(arrowstyle='-|>',color=color,lw=1,mutation_scale=12),zorder=4)
   if side==fixedside: ax.plot((start+end)/2,y+.13,'s',color=INK,ms=2.8)
  status='' if i==0 else (' · gDNA+' if p['screen']['genomic_DNA_product_detected'] else ' · gDNA−')
  ax.text(hi+4,y,f"{p['amplicon_length']} bp{status}",va='center',fontsize=7)
 ax.set_yticks(range(5),['Redesign · retain R','Redesign · retain F','Complete · retain R','Complete · retain F','Initial pair']);ax.set_ylim(-.45,4.6);ax.set_xlim(lo,hi+70);ax.set_xlabel('Transcript position (bp; zero-based)');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0,labelsize=7)
 ax.text(.99,.99,f"{target['gene_symbol']} · {selected[0]['transcript_accession_version']}   ■ = preserved primer",transform=ax.transAxes,va='top',ha='right',fontsize=7,color=INK)
 DERIVED['figure3_placement_case']=dict(gene_id=gene,gene_symbol=target['gene_symbol'],selection='first numerical gene ID with an RNA-specific candidate for all four tasks; first proposal per task; not a cohort filter',initial_pair=anchor,changed_length_band=taskledger[gene]['new_length_band'],candidates=[{k:v for k,v in p.items() if k!='screen'} for p in selected],genomic_flag=[p['screen']['genomic_DNA_product_detected'] for p in selected])


def figure3():
 fig=plt.figure(figsize=(178/25.4,183/25.4));gs=fig.add_gridspec(3,3,left=.10,right=.97,bottom=.075,top=.94,wspace=.75,hspace=.72,height_ratios=[1,1,1.15])
 ax=fig.add_subplot(gs[0,0]);panel(ax,'A','Length-conditioned design')
 bands=[(80,120),(121,160),(161,200)];band_colors=[TEAL,'#587C9E','#AB8B45']
 for i,(lo,hi) in enumerate(bands):
  ctx=f'amp{lo}-{hi}';cs=native('cnn_diffusion',ctx);flat=[p['amplicon_length'] for rows in cs.values() for p in rows]
  assert len(flat)==sum(r['unique_RNA_pass_pairs'] for r in records('cnn_diffusion',context=ctx).values())
  assert all(lo<=x<=hi for x in flat)
  length_counts=Counter(flat);xs=np.arange(75,206);cdf=np.array([sum(n for x,n in length_counts.items() if x<=threshold)/len(flat) for threshold in xs])
  ax.axvspan(lo,hi,color=band_colors[i],alpha=.075,lw=0)
  ax.step(xs,cdf,where='post',color=band_colors[i],lw=1.25,label=f'{lo}–{hi}')
  DERIVED['figure3_length_'+ctx]=dict(lengths=flat,n_unique_RNA_specific_pairs=len(flat),genes=75,definition='RNA-specific feasible unique sequence pairs in each gene prefix K=20; candidate-weighted length distribution, no genomic filtering',ecdf={'length_threshold_bp':xs.tolist(),'fraction_at_or_below':cdf.tolist()},requested_range=[lo,hi])
 ax.set(xlabel='Amplicon length (bp)',ylabel='Cumulative candidate fraction',xlim=(75,205),ylim=(-.02,1.03),xticks=[80,120,160,200],yticks=[0,.5,1]);ax.legend(title='Requested (bp)',title_fontsize=7,fontsize=7,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.24),ncol=3,handlelength=.8,labelspacing=.2,borderpad=0,columnspacing=.5,handletextpad=.3)
 ax=fig.add_subplot(gs[0,1]);panel(ax,'B','Fixed-primer completion');task_bars(ax,['amp80-200_fix_forward','amp80-200_fix_reverse'])
 ax=fig.add_subplot(gs[0,2]);panel(ax,'C','Changed-length redesign');task_bars(ax,['redesign_forward','redesign_reverse'])
 fig.legend(handles=[Patch(facecolor=COLORS[m],label=LABELS[m].replace('\n',' ')) for m in ['cnn_diffusion','cnn_direct','primer3','qprimer']],loc='upper center',bbox_to_anchor=(.52,.997),ncol=4,frameon=False,fontsize=7,columnspacing=1.2,handlelength=1.1)
 ax=fig.add_subplot(gs[1,:2]);panel(ax,'D','More proposals yield more retained alternatives');curve(ax,COUNT,['cnn_diffusion','cnn_direct','primer3','qprimer']);ax.set_ylabel('Retained pairs per gene')
 ax=fig.add_subplot(gs[1,2]);panel(ax,'E','Positional diversity');curve(ax,DIVERSITY,['cnn_diffusion','cnn_direct']);ax.set_ylabel('Amplicon Jaccard distance');ax.set_ylim(0,1);ax.legend(frameon=False,fontsize=7,loc='lower right')
 ax=fig.add_subplot(gs[2,:]); pos=ax.get_position(); ax.set_position([.235,pos.y0,.735,pos.height]); panel(ax,'F','One reference sequence · preserve either primer');footprint(ax)
 finish(fig,'fig3_flexible_design','Flexible usage with the same trained CNN diffusion model. (A) Empirical cumulative distributions of exact integer amplicon lengths, separately normalized within the three requested ranges. These candidate-weighted distributions contain 1,316, 1,324 and 1,327 unique feasible RNA-specific pairs from the 75 genes and K = 20; they are not per-gene independent samples or kernel-smoothed densities. Shaded vertical bands identify requested ranges. (B,C) Mean retained pair counts over all 75 genes, including failures and shortages, for fixed-primer completion and changed-length redesign; whiskers show existing 95% gene-bootstrap intervals. (D) Retained pair supply from the recorded prefixes K = 1, 2, 5, 10 and 20. (E) Diversity of RNA-specific amplicon intervals, compared with the matched CNN direct predictor; this endpoint requires at least two distinct intervals and source denominators are supplied. Shading in D,E shows existing 95% gene-bootstrap intervals. (F) Actual recorded RNA-specific candidates on one transcript for fixed-F/fixed-R completion and redesign; arrows point in primer extension directions and squares identify unchanged primers. gDNA+ and gDNA− respectively indicate presence and absence of a predicted genomic product. The initial pair has no genomic status assigned in this panel. This case illustrates placement, includes a gDNA+ redesign candidate, and does not claim every illustrated RNA-specific candidate passes genomic screening. Its deterministic selection does not filter any aggregate results.','Six panels show empirical amplicon length cumulative distributions, fixed-primer completion and redesign counts, candidate-budget accumulation, positional diversity, and real directed primer tracks with genomic screening flags.',183)


def figure5():
 fig,axes=plt.subplots(2,2,figsize=(178/25.4,141/25.4));fig.subplots_adjust(left=.16,right=.975,top=.925,bottom=.105,wspace=.40,hspace=.65)
 ax=axes[0,0];panel(ax,'A','Diffusion versus matched direct prediction');contrast_intervals=[]
 for i,c in enumerate(['original','external']):
  l=records('cnn_diffusion',c);r=records('cnn_direct',c);ids=sorted(l);delta=[l[g][COUNT]-r[g][COUNT] for g in ids]
  contrast=next(x for x in S[c]['paired_comparisons'] if (x['context'],x['left_method'],x['right_method'],x['left_variant'],x['right_variant'],x['budget'],x['endpoint'])==('amp80-200','cnn_diffusion','cnn_direct','full','full',20,COUNT))
  contrast_intervals.append(contrast['ci95'])
  ax.errorbar(contrast['mean_difference'],i,xerr=[[contrast['mean_difference']-contrast['ci95'][0]],[contrast['ci95'][1]-contrast['mean_difference']]],fmt='o',color=TEAL,capsize=3,ms=5,lw=1.4,zorder=4)
  ax.text(contrast['mean_difference'],i-.17,f"+{contrast['mean_difference']:.2f} pairs/gene",ha='center',va='bottom',fontsize=7,color=INK)
  DERIVED['figure5_diffusion_direct_'+c]=dict(gene_ids=ids,paired_difference=delta,paired_mean_difference=contrast['mean_difference'],paired_ci95=contrast['ci95'],n=contrast['n'])
 ax.axvline(0,color='#8F9DA3',lw=.8,ls='--');ax.set_yticks([0,1],['Original test\n75 genes','Independent test\n75 genes']);ax.set_xlabel('Paired mean gain: diffusion − direct');ax.set_ylim(1.55,-.65);ax.set_xlim(-.15,max(ci[1] for ci in contrast_intervals)+.3);ax.grid(False);ax.grid(axis='x',color=GRID,lw=.5)
 ax=axes[0,1];panel(ax,'B','Constraint-guided pair sampling')
 variants=['raw','geometry_only','single_oligo','full'];keys=['unique_feasible_pairs','unique_RNA_pass_pairs',COUNT]
 matrix=np.array([[value('cnn_diffusion',key,variant=v)[0] for key in keys] for v in variants]);ax.imshow(matrix,vmin=0,vmax=20,cmap='GnBu',aspect='auto')
 for i in range(4):
  for j in range(3):ax.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if matrix[i,j]>13 else INK)
 ax.set_xticks([0,1,2],['Feasible','RNA-specific','Retained']);ax.set_yticks(range(4),['Raw labels','Geometry','Single oligo','Full pair']);ax.tick_params(length=0);ax.set_xlabel('Mean unique pairs per gene · K = 20');ax.grid(False)
 DERIVED['figure5_constraint_guidance']=dict(variants=variants,endpoints=keys,means=matrix.tolist())
 ax=axes[1,0];panel(ax,'C','Denoising schedule')
 ys=[];cis=[]
 for step in [1,4,8,16]:
  y,ci=value('cnn_diffusion',variant='full' if step==16 else f'steps{step}');ys.append(y);cis.append(ci)
 ax.plot([1,4,8,16],ys,'o-',color=TEAL,lw=1.7,ms=4);ax.fill_between([1,4,8,16],[c[0] for c in cis],[c[1] for c in cis],color=TEAL,alpha=.12,lw=0)
 ax.set(xlabel='Reverse-diffusion steps',ylabel='Retained pairs per gene',xticks=[1,4,8,16],ylim=(0,10))
 ax=axes[1,1];panel(ax,'D','Denoiser architecture')
 methods=['cnn_diffusion','attention_hybrid_diffusion','attention_rope_diffusion'];x=np.arange(3);w=.31
 for j,c in enumerate(['original','external']):
  vals=[value(m,cohort=c) for m in methods];y=[v[0] for v in vals];cis=[v[1] for v in vals]
  ax.bar(x+(j-.5)*w,y,w,color=[COLORS[m] for m in methods],alpha=1 if j==0 else .36,edgecolor=[COLORS[m] for m in methods],lw=.7,zorder=3);errbars(ax,x+(j-.5)*w,y,cis)
 ax.set_xticks(x,['CNN','CNN +\nattention','Attention']);ax.set_ylabel('Retained pairs per gene');ax.set_ylim(0,10)
 ax.legend(handles=[Patch(facecolor='#617D79',label='Original test'),Patch(facecolor='#C9D7D5',edgecolor='#617D79',label='Independent test')],frameon=False,fontsize=7,loc='upper right')
 finish(fig,'fig5_mechanism','Mechanism and architectural analyses at K = 20. (A) Mean gene-paired gains in retained count from CNN diffusion relative to the matched CNN direct predictor; dots and horizontal whiskers show the frozen paired mean difference and 95% gene-paired bootstrap interval in each cohort. The dashed vertical line is no mean difference. All 75 underlying paired gene differences per cohort remain in the source export. (B) Staged constraint guidance: mean feasible, RNA-specific and retained unique pairs per gene across all 75 original test genes, including zero-yield outcomes. All cells use the common color scale of 0–20 pairs. (C) Retained candidate counts for recorded reverse-diffusion schedules of 1, 4, 8 and 16 steps; line and shading show mean and existing 95% gene-bootstrap intervals. (D) Diffusion denoiser backbones, shown with original and independent cohorts kept separate. Bars and whiskers show means and existing 95% gene-bootstrap intervals. Architectures and step counts are measured comparisons, not an assumption of monotonic improvement.','Four panels show paired mean differences between diffusion and direct prediction with confidence intervals, a constraint-guidance heatmap, the denoising schedule, and CNN, hybrid and attention architecture comparisons.',141)


def main():
 OUT.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.titlesize':8.1,'axes.labelsize':7.2,'xtick.labelsize':7.1,'ytick.labelsize':7.1,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#A7B4B9','axes.linewidth':.6,'text.color':INK,'axes.labelcolor':INK,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'path'})
 figure2();figure3();figure5()
 (OUT/'figure_values.json').write_text(json.dumps({'primary':'cnn_diffusion','values':VALUES,'derived':DERIVED},indent=2)+'\n')
 manifest=dict(primary='cnn_diffusion',seed=20261013,source_summary_status={c:s['status'] for c,s in S.items()},source_cohorts={c:s['total_test_genes'] for c,s in S.items()},source_files=SOURCES,derivation_code={'path':str(Path(__file__).relative_to(ROOT)),'sha256':sha(__file__)},figure_values_sha256=sha(OUT/'figure_values.json'),physical_width_mm=178,minimum_font_pt=7,output_dpi=350,cohort_policy='Original and external 75-gene cohorts never pooled or paired across cohorts',intervals='Reuse frozen gene-bootstrap/paired gene-bootstrap intervals; no new seeds or inferential resampling',plot_form_policy='Choose each form from the scientific question and endpoint; do not diversify types for appearance',question_endpoint_plot_rationale=PLOT_RATIONALE,figures=FIGURES)
 (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 provenance=dict(primary_method='cnn_diffusion',summary_sha256=sha(ROOT/'experiments/attention/results/study_summary.json'),external_summary_sha256=sha(ROOT/'experiments/attention_external/results/study_summary.json'),generated=FIGURES,figure_values_sha256=sha(OUT/'figure_values.json'),manifest_sha256=sha(OUT/'manifest.json'))
 (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 readme=['# CNN publication result figures','','Build: `.venv/bin/python experiments/attention/plot_publication_results.py`','','All source results and existing figures remain unchanged. Primary method is CNN diffusion. Figures are 178 mm wide, white background, vector PDF/SVG and 350 dpi PNG/TIFF. All original and external test cohorts have 75 genes and remain separate. `figure_values.json` preserves exact plotted source means, intervals, denominators, paired gene values, native amplicon lengths, and illustration candidate IDs. `manifest.json` records source and derivation hashes. Print reviews are exported at approximately 700 px width.','','The quantity endpoint counts unique RNA-specific pairs without a predicted genomic product over every gene. The quality endpoint is a per-gene ratio conditional on at least one RNA-specific pair; it is not a pooled ratio or measured amplification efficiency. qPrimer denotes the upstream design module plus shared assessment; PrimerCast proposals use the exact-reference task shared by all methods.','','']
 readme.extend(['## Scientific questions determine plot form','','The audit changes only three forms where the former form obscured the scientific target: cohort mean bars become mean points with intervals; length KDEs become empirical discrete-length CDFs; paired-difference KDEs become paired mean estimates with intervals. Task bars, budget curves, coordinate tracks, constraint heatmaps and architectural bars retain their existing roles.','','| Panel | Scientific question | Endpoint | Plot form and reason |','| --- | --- | --- | --- |'])
 for r in PLOT_RATIONALE:readme.append(f"| {r['figure']} {r['panel']} | {r['question']} | {r['endpoint']} | {r['form']}. {r['reason']} |")
 readme.append('')
 for f in FIGURES:readme.extend([f"## {f['name']}",'',f['caption'],'',f"Alt text: {f['alt_text']}",''])
 (OUT/'README.md').write_text('\n'.join(readme))
 # The figure builder verifies every plotted estimate against the frozen table.
 for r in VALUES:
  g=G[r['cohort']][r['context'],r['method'],r['variant'],r['budget']]
  assert r['mean']==g[r['endpoint']]*r['scale']
  assert r['ci95']==[x*r['scale'] for x in g[r['endpoint']+'_ci95']]
  assert r['n']==g[r['endpoint']+'_n']
 raster_checks=[]
 for f in FIGURES:
  for item in f['files']:
   path=ROOT/item['path']; assert sha(path)==item['sha256']
   if path.suffix in ['.png','.tiff']:
    with Image.open(path) as im:
     assert im.mode=='RGB' and im.width>=2450
     assert all(abs(dpi-350)<.1 for dpi in im.info['dpi'])
     raster_checks.append(dict(path=item['path'],size=list(im.size),mode=im.mode))
 verification=dict(status='passed',verified_source_values=len(VALUES),source_files_checked=len(SOURCES),source_summary_sha256=provenance['summary_sha256'],external_summary_sha256=provenance['external_summary_sha256'],provenance_sha256=sha(OUT/'provenance.json'),native_length_distribution_checks=[dict(context=k.replace('figure3_length_',''),pairs=r['n_unique_RNA_specific_pairs'],min=min(r['lengths']),max=max(r['lengths'])) for k,r in DERIVED.items() if k.startswith('figure3_length_')],preserved_anchor_geometry='All four illustrated fixed primers match source anchor sequence and interval',illustrated_genomic_flags=DERIVED['figure3_placement_case']['genomic_flag'],raster_checks=raster_checks,visual_review='pending for this build')
 (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
 print(json.dumps({'output':str(OUT),'figures':[f['name'] for f in FIGURES],'plotted_source_values':len(VALUES),'hashed_sources':len(SOURCES)},indent=2))
if __name__=='__main__':main()
