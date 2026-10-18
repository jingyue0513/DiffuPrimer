#!/usr/bin/env python3
"""Combine completed main controls with validation-selected strong baselines."""
import json
import sys
import time
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path[:0]=[str(ROOT/'src'),str(HERE)]
import numpy as np
from diffuprimer.io import read_jsonl,write_json,sha256
from diffuprimer.metrics import paired_bootstrap_difference


def main():
    base=HERE/'results_seed42'; tuned=HERE/'tuned_seed42'; out=HERE/'comparison_seed42'
    out.mkdir(exist_ok=True)
    write_json(out/'status.json',{'stage':'waiting_for_both_studies'})
    while True:
        states=[]
        for p in [base,tuned]:
            try:states.append(json.loads((p/'status.json').read_text()).get('stage'))
            except (FileNotFoundError,json.JSONDecodeError):states.append('starting')
        if 'failed' in states:
            raise RuntimeError(f'A prerequisite failed: {states}')
        if states==['complete','complete']:break
        time.sleep(5)
    summary=json.loads((base/'summary.json').read_text())
    check=json.loads((base/'verification.json').read_text())
    assert check['passed'] and check['summary_sha256']==sha256(base/'summary.json')
    selected=json.loads((tuned/'selection.json').read_text())
    supplement=json.loads((tuned/'summary.json').read_text())
    assert selected['test_outcomes_consulted'] is False
    assert supplement['selection_sha256']==sha256(tuned/'selection.json')
    genes={}
    groups=[]
    for r in read_jsonl(base/'per_gene.jsonl'):
        genes.setdefault((r['cohort'],r['method']),{})[r['gene_id']]=r['retained_pairs']
    for g in summary['groups']:
        groups.append({k:g[k] for k in ['cohort','method','retained_pairs','retained_pairs_ci95',
                       'one_retained','two_retained','generation_cpu_seconds']})
    for r in supplement['results']:
        method='direct_validation_selected' if r['model']=='direct' else 'diffusion1_validation_selected'
        own=r['per_gene_retained'];genes[r['cohort'],method]=own
        directory=Path(r['directory'])
        assert sha256(directory/'specificity_bowtie.jsonl')==r['screen_sha256']
        costs=[x['generation_cpu_seconds'] for x in read_jsonl(directory/'per_gene.jsonl')]
        stat=paired_bootstrap_difference(list(own.values()),np.zeros(len(own)),seed=42)
        groups.append({'cohort':r['cohort'],'method':method,'configuration':r['config'],
                       'retained_pairs':r['retained_pairs'],'retained_pairs_ci95':stat['ci95'],
                       'one_retained':r['one_retained'],'two_retained':r['two_retained'],
                       'generation_cpu_seconds':float(np.mean(costs))})
    contrasts=[]
    pairs=[('diffusion_16','direct_validation_selected'),('diffusion_4','direct_validation_selected'),
           ('diffusion_1','direct_validation_selected'),('diffusion1_validation_selected','direct_validation_selected'),
           ('diffusion_16','diffusion1_validation_selected')]
    for cohort in ['test','external_test']:
        for left,right in pairs:
            a,b=genes[cohort,left],genes[cohort,right]
            assert a.keys()==b.keys() and len(a)==75
            ids=sorted(a)
            stat=paired_bootstrap_difference([a[g] for g in ids],[b[g] for g in ids],seed=42)
            contrasts.append({'cohort':cohort,'left':left,'right':right,'endpoint':'retained_pairs',**stat})
    report={'seed':42,'groups':groups,'main_contrasts':summary['contrasts'],
            'validation_selected_contrasts':contrasts,'validation_selection':selected,
            'source_hashes':{'main_summary':sha256(base/'summary.json'),'main_verification':sha256(base/'verification.json'),
                             'tuned_summary':sha256(tuned/'summary.json'),'selection':sha256(tuned/'selection.json'),
                             'finalizer':sha256(Path(__file__))},
            'limitations':['One historical checkpoint pair; no training-seed uncertainty.',
                           'External cohort is same-species unseen genes, not cross-species.',
                           'Equal CPU allowance covers generation, not downstream full-reference search.',
                           'Deadline baseline has a 4096-candidate cap; unused CPU allowance is not proof of optimal direct search.',
                           'Gene-paired confidence intervals are unadjusted for multiple comparisons.']}
    write_json(out/'summary.json',report)
    lines=['Mechanism study and validation-selected full-support baselines; seed 42','']
    lines.append('Validation-selected configurations: '+json.dumps(selected['selected']))
    for cohort in ['test','external_test']:
        lines.extend(['',cohort])
        for g in groups:
            if g['cohort']==cohort:
                lines.append(f"{g['method']:32s} retained={g['retained_pairs']:.3f} >=1={g['one_retained']:.3f} >=2={g['two_retained']:.3f} CPU={g['generation_cpu_seconds']:.3f}s")
        for c in summary['contrasts']+contrasts:
            if c['cohort']==cohort and c['endpoint']=='retained_pairs':
                lines.append(f"{c['left']} - {c['right']}: {c['mean_difference']:+.3f}, 95% CI {c['ci95']}")
    (out/'summary.txt').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(16,7.5))
    for ax,cohort in zip(axes,['test','external_test']):
        own=[g for g in groups if g['cohort']==cohort]
        values=np.array([g['retained_pairs'] for g in own]);ci=np.array([g['retained_pairs_ci95'] for g in own])
        colors=['#9b6caa' if 'validation_selected' in g['method'] else '#df9c47' if 'frozen' in g['method'] else '#4c85a9' for g in own]
        ax.barh([g['method'] for g in own],values,xerr=np.maximum(0,np.array([values-ci[:,0],ci[:,1]-values])),color=colors,capsize=3)
        ax.invert_yaxis();ax.set_xlim(0,max(1.,ci[:,1].max())*1.15);ax.set_title(cohort+' (75 genes)');ax.set_xlabel('Unique retained pairs / gene (K=20)')
        ax.spines[['top','right']].set_visible(False)
        for i,v in enumerate(values):ax.text(ci[i,1]+.08,i,f'{v:.2f}',va='center',fontsize=9)
    fig.suptitle('Diffusion mechanism controls and validation-selected direct baselines | seed 42')
    fig.tight_layout();fig.savefig(out/'results.png',dpi=180);fig.savefig(out/'results.pdf');plt.close(fig)
    write_json(out/'status.json',{'stage':'complete','cohorts':2,'genes_per_cohort':75,'methods_per_cohort':13,
                                'summary_sha256':sha256(out/'summary.json')})
    print('\n'.join(lines),flush=True)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        write_json(HERE/'comparison_seed42/status.json',{'stage':'failed','error':repr(exc)})
        raise
