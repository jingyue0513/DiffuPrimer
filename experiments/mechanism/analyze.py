#!/usr/bin/env python3
"""Validate complete screens and report the predeclared mechanism contrasts."""
from __future__ import annotations
import argparse
import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from diffuprimer.io import read_jsonl,write_json,write_jsonl,sha256
from diffuprimer.metrics import paired_bootstrap_difference

CONTRASTS=[
    ('diffusion_1','direct_beam','denoising_training'),
    ('diffusion_4','diffusion_1','four_step_increment'),
    ('diffusion_16','diffusion_1','sixteen_step_increment'),
    ('frozen_4','diffusion_1','decoder_only_four_step_increment'),
    ('frozen_16','diffusion_1','decoder_only_sixteen_step_increment'),
    ('diffusion_4','frozen_4','learned_feedback_four_steps'),
    ('diffusion_16','frozen_16','learned_feedback_sixteen_steps'),
    ('direct_full_wor','direct_beam','stronger_direct_support'),
    ('diffusion1_full_wor','direct_full_wor','training_under_identical_full_support'),
    ('diffusion1_pool320','direct_pool320','training_under_identical_pool_selection'),
    ('diffusion_4','direct_full_wor','four_steps_vs_full_direct'),
    ('diffusion_16','direct_full_wor','sixteen_steps_vs_full_direct'),
    ('diffusion_16','direct_pool320','sixteen_steps_vs_direct_diverse_pool'),
    ('diffusion_16','direct_equal_cpu','equal_inference_cpu_budget'),
]
ENDPOINTS=['retained_pairs','rna_pairs','feasible_pairs','one_retained','two_retained',
           'retained_fraction','retained_diversity','rna_diversity','best_retained_tm_penalty',
           'generation_cpu_seconds','generation_wall_seconds']


def diversity(rows):
    spans=sorted({(r['forward_start'],r['reverse_end']) for r in rows})
    values=[]
    for (a,b),(c,d) in itertools.combinations(spans,2):
        intersection=max(0,min(b,d)-max(a,c))
        values.append(1-intersection/(b-a+d-c-intersection))
    return float(np.mean(values)) if values else None


def write_csv(path, rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=keys);writer.writeheader()
        writer.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)


def analyze(out):
    protocol=json.loads((out/'protocol.json').read_text())
    records,provenance=[],[]
    for cohort,cohort_spec in protocol['cohorts'].items():
        expected=set(map(str,cohort_spec['targets']))
        for method in protocol['methods']:
            directory=out/'evaluations'/f'{cohort}_{method}'
            manifest=json.loads((directory/'manifest.json').read_text())
            completion=json.loads((directory/'completion.json').read_text())
            status=json.loads((directory/'specificity_bowtie_status.json').read_text())
            assert manifest['protocol_sha256']==sha256(out/'protocol.json')
            assert completion['all_expected_targets_evaluated'] and status['status']=='completed'
            assert completion['candidates_sha256']==status['candidate_sha256']==sha256(directory/'candidates.jsonl')
            assert completion['per_gene_sha256']==sha256(directory/'per_gene.jsonl')
            assert status['output_sha256']==sha256(directory/'specificity_bowtie.jsonl')
            candidates=read_jsonl(directory/'candidates.jsonl')
            ledger=read_jsonl(directory/'per_gene.jsonl')
            screens=read_jsonl(directory/'specificity_bowtie.jsonl')
            assert len(ledger)==len(expected) and {str(r['gene_id']) for r in ledger}==expected
            assert len({r['candidate_id'] for r in candidates})==len(candidates)
            assert len({r['candidate_id'] for r in screens})==len(screens)
            assert {r['candidate_id'] for r in candidates}=={r['candidate_id'] for r in screens}
            screens={r['candidate_id']:r for r in screens}
            by_gene=defaultdict(list)
            for r in candidates:by_gene[str(r['gene_id'])].append(r)
            for l in ledger:
                gene=str(l['gene_id']); own=by_gene[gene]
                assert len(own)<=20
                feasible={}
                for r in own:
                    if r['assessment']['computational_feasible']:
                        assert 'gene_specificity_screen_pass' in screens[r['candidate_id']]
                        assert 'genomic_DNA_product_detected' in screens[r['candidate_id']]
                        feasible.setdefault((r['forward_primer'],r['reverse_primer']),r)
                rna=[r for r in feasible.values() if screens[r['candidate_id']]['gene_specificity_screen_pass']]
                retained=[r for r in rna if not screens[r['candidate_id']]['genomic_DNA_product_detected']]
                if method=='direct_equal_cpu':
                    assert l['accepted_checkpoint_cpu_seconds']<=l['cpu_budget_seconds']
                    assert all(r['cpu_seconds']<=l['cpu_budget_seconds'] for r in l['trace'] if r['accepted_before_deadline'])
                records.append({'cohort':cohort,'method':method,'gene_id':gene,
                    'retained_pairs':len(retained),'rna_pairs':len(rna),'feasible_pairs':len(feasible),
                    'one_retained':int(bool(retained)),'two_retained':int(len(retained)>=2),
                    'retained_fraction':len(retained)/len(rna) if rna else None,
                    'retained_diversity':diversity(retained),'rna_diversity':diversity(rna),
                    'best_retained_tm_penalty':min((r['assessment']['shared_temperature_penalty'] for r in retained),default=None),
                    'generation_cpu_seconds':l['generation_cpu_seconds'],
                    'generation_wall_seconds':l['generation_wall_seconds'],
                    'pool_size':l.get('pool_size',20),'support_size':l.get('support_size'),
                    'network_calls':l['network_calls'],
                    'cpu_budget_seconds':l.get('cpu_budget_seconds'),
                    'accepted_checkpoint_cpu_seconds':l.get('accepted_checkpoint_cpu_seconds'),
                    'overshoot_cpu_seconds':l.get('overshoot_cpu_seconds')})
            provenance.append({'directory':str(directory),'candidate_sha256':status['candidate_sha256'],
                'screen_sha256':status['output_sha256'],'screen_identity_sha256':status['identity_sha256']})
    assert len({p['screen_identity_sha256'] for p in provenance})==1
    grouped=defaultdict(list)
    for r in records:grouped[r['cohort'],r['method']].append(r)
    groups,contrasts=[],[]
    for (cohort,method),rows in grouped.items():
        g={'cohort':cohort,'method':method,'total_genes':len(rows)}
        for endpoint in ENDPOINTS:
            vals=[r[endpoint] for r in rows if r[endpoint] is not None]
            stat=paired_bootstrap_difference(vals,np.zeros(len(vals)),seed=42) if vals else None
            g[endpoint]=stat['mean_difference'] if stat else None
            g[endpoint+'_ci95']=stat['ci95'] if stat else None
            g[endpoint+'_n']=len(vals)
        g['generation_cpu_median']=float(np.median([r['generation_cpu_seconds'] for r in rows]))
        g['generation_wall_median']=float(np.median([r['generation_wall_seconds'] for r in rows]))
        g['mean_internal_pool_size']=float(np.mean([r['pool_size'] for r in rows]))
        if method=='direct_equal_cpu':
            g['mean_cpu_budget']=float(np.mean([r['cpu_budget_seconds'] for r in rows]))
            g['mean_cpu_overshoot']=float(np.mean([r['overshoot_cpu_seconds'] for r in rows]))
            g['max_cpu_overshoot']=max(r['overshoot_cpu_seconds'] for r in rows)
            g['genes_without_on_time_output']=sum(r['accepted_checkpoint_cpu_seconds']==0 for r in rows)
        groups.append(g)
    for cohort in protocol['cohorts']:
        for left,right,question in CONTRASTS:
            a={r['gene_id']:r for r in grouped[cohort,left]}
            b={r['gene_id']:r for r in grouped[cohort,right]}
            assert a.keys()==b.keys()
            for endpoint in ENDPOINTS:
                genes=sorted(g for g in a if a[g][endpoint] is not None and b[g][endpoint] is not None)
                if not genes:continue
                stat=paired_bootstrap_difference([a[g][endpoint] for g in genes],[b[g][endpoint] for g in genes],seed=42)
                contrasts.append({'cohort':cohort,'question':question,'left':left,'right':right,'endpoint':endpoint,**stat})
    result={'seed':42,'protocol_sha256':sha256(out/'protocol.json'),'analysis_code_sha256':sha256(Path(__file__)),'cohorts_separate':True,
            'bootstrap_unit':'gene; paired, 10000 repeats; CIs unadjusted for multiple contrasts',
            'training_seed_limit':'One reused checkpoint pair; confidence intervals do not estimate training-seed variation',
            'compute_limit':'CPU generation budgets include network, search and thermodynamics, exclude identical shared span preparation and downstream specificity. Actual overshoot reported; only on-time output accepted.',
            'provenance':provenance,'groups':groups,'contrasts':contrasts}
    write_json(out/'summary.json',result)
    write_jsonl(out/'per_gene.jsonl',records)
    write_csv(out/'groups.csv',groups);write_csv(out/'contrasts.csv',contrasts)
    write_json(out/'verification.json',{'passed':True,'runs':len(provenance),'gene_method_records':len(records),
              'screen_identity_sha256':provenance[0]['screen_identity_sha256'],
              'all_candidates_screened':True,'timed_output_deadlines_respected':True,
              'summary_sha256':sha256(out/'summary.json')})
    lines=['Single-seed mechanism study (seed=42; historical checkpoint seed retained as provenance)','']
    for cohort in protocol['cohorts']:
        lines.append(cohort)
        for g in groups:
            if g['cohort']==cohort:
                lines.append(f"{g['method']:24s} retained={g['retained_pairs']:.3f} >=1={g['one_retained']:.3f} >=2={g['two_retained']:.3f} CPU={g['generation_cpu_seconds']:.3f}s")
        lines.append('Predeclared retained-pair contrasts (left minus right):')
        for c in contrasts:
            if c['cohort']==cohort and c['endpoint']=='retained_pairs':
                lines.append(f"{c['left']} - {c['right']}: {c['mean_difference']:+.3f}, 95% CI {c['ci95']}")
        lines.append('')
    (out/'summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines),flush=True)
    plot(out,groups,contrasts,protocol)


def plot(out,groups,contrasts,protocol):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    methods=list(protocol['methods'])
    colors=['#607d8b','#4c78a8','#2878b5','#15517a','#f0a447','#ce772b',
            '#62a78f','#439478','#8f79b5','#7056a0','#bf5368']
    fig,axes=plt.subplots(1,len(protocol['cohorts']),figsize=(13,6.8),squeeze=False)
    for ax,cohort in zip(axes[0],protocol['cohorts']):
        lookup={g['method']:g for g in groups if g['cohort']==cohort}
        values=np.array([lookup[m]['retained_pairs'] for m in methods])
        ci=np.array([lookup[m]['retained_pairs_ci95'] for m in methods])
        ax.barh(methods,values,color=colors,xerr=np.maximum(0,np.array([values-ci[:,0],ci[:,1]-values])),capsize=3)
        ax.invert_yaxis();ax.set_xlabel('Unique retained pairs / gene (K=20)');ax.set_title(f"{cohort} (n={len(protocol['cohorts'][cohort]['targets'])})")
        ax.set_xlim(0,max(1.,ci[:,1].max())*1.15);ax.spines[['top','right']].set_visible(False)
        for i,v in enumerate(values):ax.text(ci[i,1]+.08,i,f'{v:.2f}',va='center',fontsize=9)
    fig.suptitle('Frozen-score and strong direct-decoding controls | sampling seed 42')
    fig.tight_layout();fig.savefig(out/'mechanism_results.png',dpi=180);fig.savefig(out/'mechanism_results.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,len(protocol['cohorts']),figsize=(12,5),squeeze=False)
    for ax,cohort in zip(axes[0],protocol['cohorts']):
        lookup={g['method']:g for g in groups if g['cohort']==cohort}
        for m,color in zip(methods,colors):
            g=lookup[m];ax.scatter(g['generation_cpu_seconds'],g['retained_pairs'],c=color,s=55,label=m)
        ax.set_xscale('log');ax.set_xlabel('Mean actual generation CPU seconds (log scale)')
        ax.set_ylabel('Retained pairs / gene');ax.set_title(cohort)
    axes[0,-1].legend(fontsize=8,loc='center left',bbox_to_anchor=(1,0.5))
    fig.tight_layout();fig.savefig(out/'quality_compute.png',dpi=180);fig.savefig(out/'quality_compute.pdf');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    analyze(p.parse_args().output.resolve())
