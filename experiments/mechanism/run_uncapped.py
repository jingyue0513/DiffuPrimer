#!/usr/bin/env python3
"""Deadline-driven direct search, removing the diagnostic 4096-pair cap.

Added because EVERY target hit that cap before its CPU deadline. Original
results remain intact. No quality-dependent hyperparameter is changed.
"""
import concurrent.futures
import fcntl
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from run import (torch,np,SEED,CHECKPOINT_SEED,ConditionalDenoiser,Conditions,guidance,
                 install,read_jsonl,write_jsonl,write_json,sha256,atomic_json,stable_seed,
                 make_candidates,assess_pair,oligo_properties,summarize_candidates)
from controls import first_scores,all_support,weighted_order,diverse_indices
from diffuprimer.decode import grammar_candidates
from diffuprimer.sequence import PrimerPair
from diffuprimer.metrics import paired_bootstrap_difference

MODEL=None
CONDITIONS=None


def initialize():
    global MODEL,CONDITIONS
    torch.set_num_threads(4);torch.set_num_interop_threads(1)
    torch.manual_seed(SEED);torch.use_deterministic_algorithms(True)
    install()
    ck=ROOT/f'results/checkpoints/direct_weak_partial_pilot_seed{CHECKPOINT_SEED}.pt'
    state=torch.load(ck,map_location='cpu',weights_only=True)
    MODEL=ConditionalDenoiser(**state['model_spec']);MODEL.load_state_dict(state['state_dict']);MODEL.eval()
    CONDITIONS=Conditions(**json.loads((HERE/'results_seed42/protocol.json').read_text())['conditions'])
    for _ in range(2):first_scores(MODEL,'ACGT'*128,CONDITIONS)


def search(template,allowed,temperatures,seed,budget):
    start=time.process_time();wall=time.perf_counter()
    pool,assessments,seen=[],[],set()
    accepted,accepted_assessments,trace=[],[],[]
    last_selection_cost=0.005

    def commit(stage):
        nonlocal accepted,accepted_assessments,last_selection_cost
        before=time.process_time()
        chosen=diverse_indices(pool,assessments)
        last_selection_cost=time.process_time()-before
        elapsed=time.process_time()-start
        on_time=elapsed<=budget
        trace.append({'stage':stage,'pool_size':len(pool),'cpu_seconds':elapsed,
                      'selection_cpu_seconds':last_selection_cost,'accepted_before_deadline':on_time})
        if on_time:
            accepted=[pool[i] for i in chosen];accepted_assessments=[assessments[i] for i in chosen]
        return on_time

    logits=first_scores(MODEL,template,CONDITIONS)
    beam,beam_scores=grammar_candidates(logits,CONDITIONS,span_allowed=allowed,span_temperatures=temperatures)
    for i in weighted_order(beam_scores,seed,len(beam_scores)):
        pair=beam[i];key=pair.sequences(template)
        if key not in seen:
            seen.add(key);pool.append(pair);assessments.append(assess_pair(template,pair,CONDITIONS))
        if len(pool)>=20 or time.process_time()-start>=budget:break
    commit('beam20')
    support_size=None;reason='deadline'
    if time.process_time()-start<budget:
        coordinates,scores=all_support(logits,CONDITIONS,allowed,temperatures)
        support_size=len(scores)
        if time.process_time()-start<budget:
            order=weighted_order(scores,seed,len(scores))
            exhausted=True
            for i in order:
                # Reserve measured selection cost so recently evaluated pairs
                # can actually contribute to an on-time committed output.
                remaining=budget-(time.process_time()-start)
                if remaining<=max(0.01,2*last_selection_cost):
                    exhausted=False
                    if remaining>0:commit('deadline_final')
                    break
                pair=PrimerPair(*map(int,coordinates[i]));key=pair.sequences(template)
                if key in seen:continue
                seen.add(key);pool.append(pair);assessments.append(assess_pair(template,pair,CONDITIONS))
                if len(pool)%512==0 and not commit('full_support'):
                    exhausted=False;break
            if exhausted:
                reason='entire_support_exhausted';commit('support_exhausted')
    elapsed=time.process_time()-start
    times=[r['cpu_seconds'] for r in trace if r['accepted_before_deadline']]
    return accepted,accepted_assessments,{'pool_size':len(pool),'support_size':support_size,
        'candidate_count_cap':None,'stop_reason':reason,'trace':trace,
        'generation_cpu_seconds':elapsed,'generation_wall_seconds':time.perf_counter()-wall,
        'cpu_budget_seconds':budget,'accepted_checkpoint_cpu_seconds':max(times,default=0.),
        'overshoot_cpu_seconds':max(0.,elapsed-budget),'budget_utilization':elapsed/budget,
        'network_calls':1,'network_examples':1}


def worker(job):
    cohort,target,budget,destination,protocol_hash=job
    path=Path(destination);template=target['template_sequence'];gene=str(target['gene_id'])
    span=guidance.feasible_oligo_spans(template,CONDITIONS)
    temperatures=guidance.compatible_span_temperatures(template,CONDITIONS,span)
    oligo_properties.cache_clear()
    pairs,assessments,details=search(template,span[1],temperatures,stable_seed(cohort,gene),budget)
    candidates=make_candidates(f'{cohort}_direct_equal_cpu_uncapped',target,pairs,assessments)
    row={'cohort':cohort,'gene_id':gene,'protocol_sha256':protocol_hash,'candidates':candidates,
         'ledger':{'gene_id':gene,**summarize_candidates(candidates,20),**details}}
    assert details['accepted_checkpoint_cpu_seconds']<=budget
    atomic_json(path,row)
    return {'cohort':cohort,'gene':gene,'pool':details['pool_size'],
            'cpu':round(details['generation_cpu_seconds'],3),'budget':round(budget,3)}


def summarize(out,protocol):
    from tune_direct import counts
    prior=HERE/'comparison_seed42/summary.json'
    report=json.loads(prior.read_text())
    report['source_hashes']['previous_comparison']=sha256(prior)
    report['source_hashes']['uncapped_protocol']=sha256(out/'protocol.json')
    report['source_hashes']['uncapped_analysis']=sha256(Path(__file__))
    report['uncapped_contrasts']=[];report['uncapped_budget_diagnostics']=[]
    base_rows=read_jsonl(HERE/'results_seed42/per_gene.jsonl')
    for cohort in protocol['cohorts']:
        dest=out/'evaluations'/f'{cohort}_direct_equal_cpu_uncapped'
        kept=counts(dest,protocol['cohorts'][cohort])
        ledger=read_jsonl(dest/'per_gene.jsonl')
        stats=paired_bootstrap_difference(list(kept.values()),np.zeros(len(kept)),seed=42)
        report['groups'].append({'cohort':cohort,'method':'direct_equal_cpu_uncapped',
            'retained_pairs':stats['mean_difference'],'retained_pairs_ci95':stats['ci95'],
            'one_retained':float(np.mean([v>=1 for v in kept.values()])),
            'two_retained':float(np.mean([v>=2 for v in kept.values()])),
            'generation_cpu_seconds':float(np.mean([r['generation_cpu_seconds'] for r in ledger]))})
        direct={r['gene_id']:r for r in ledger}
        for method in ['diffusion_4','diffusion_16','frozen_16']:
            left={r['gene_id']:r for r in base_rows if r['cohort']==cohort and r['method']==method}
            assert left.keys()==kept.keys() and len(left)==75
            ids=sorted(left)
            for metric in ['retained_pairs','one_retained','two_retained']:
                vals=[kept[g] if metric=='retained_pairs' else int(kept[g]>=(1 if metric=='one_retained' else 2)) for g in ids]
                s=paired_bootstrap_difference([left[g][metric] for g in ids],vals,seed=42)
                report['uncapped_contrasts'].append({'cohort':cohort,'left':method,'right':'direct_equal_cpu_uncapped','endpoint':metric,**s})
        report['uncapped_budget_diagnostics'].append({'cohort':cohort,
            'mean_pool_size':float(np.mean([r['pool_size'] for r in ledger])),
            'mean_utilization':float(np.mean([r['budget_utilization'] for r in ledger])),
            'min_utilization':min(r['budget_utilization'] for r in ledger),
            'max_overshoot_seconds':max(r['overshoot_cpu_seconds'] for r in ledger),
            'support_exhausted_genes':sum(r['stop_reason']=='entire_support_exhausted' for r in ledger),
            'all_deadlines_respected':all(r['accepted_checkpoint_cpu_seconds']<=r['cpu_budget_seconds'] for r in ledger)})
    write_json(out/'summary.json',report)
    lines=['Complete mechanism comparison, including uncapped CPU-deadline direct search. Seed 42.','']
    for cohort in protocol['cohorts']:
        lines.append(cohort)
        for g in report['groups']:
            if g['cohort']==cohort:
                lines.append(f"{g['method']:32s} retained={g['retained_pairs']:.3f} coverage1={g['one_retained']:.3f} coverage2={g['two_retained']:.3f} CPU={g['generation_cpu_seconds']:.3f}")
    lines.extend(['',json.dumps(report['uncapped_budget_diagnostics'],indent=2)])
    (out/'summary.txt').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    methods=['direct_validation_selected','diffusion1_validation_selected','diffusion_1','frozen_4',
             'diffusion_4','frozen_16','diffusion_16','direct_equal_cpu_uncapped']
    fig,axes=plt.subplots(1,2,figsize=(15,6.5))
    for ax,cohort in zip(axes,protocol['cohorts']):
        lookup={g['method']:g for g in report['groups'] if g['cohort']==cohort}
        v=np.array([lookup[m]['retained_pairs'] for m in methods]);ci=np.array([lookup[m]['retained_pairs_ci95'] for m in methods])
        ax.barh(methods,v,xerr=np.maximum(0,np.array([v-ci[:,0],ci[:,1]-v])),capsize=3,
                color=['#986ab0','#986ab0','#4c85a9','#e0a34c','#4c85a9','#e0a34c','#4c85a9','#ba5667'])
        ax.invert_yaxis();ax.set_xlim(0,max(ci[:,1])*1.15);ax.set_title(cohort+' (75 genes)')
        ax.set_xlabel('Unique retained pairs / gene (K=20)');ax.spines[['top','right']].set_visible(False)
        for i,x in enumerate(v):ax.text(ci[i,1]+.08,i,f'{x:.2f}',va='center')
    fig.tight_layout();fig.savefig(out/'results.png',dpi=180);fig.savefig(out/'results.pdf');plt.close(fig)
    print('\n'.join(lines),flush=True)


def main():
    out=HERE/'uncapped_seed42';out.mkdir(exist_ok=True)
    base=HERE/'results_seed42';base_protocol=json.loads((base/'protocol.json').read_text())
    targets={'test':[r for r in read_jsonl(ROOT/'data/processed/targets_pilot.jsonl') if r['split']=='test'],
             'external_test':read_jsonl(ROOT/'data/processed/attention_external_targets.jsonl')}
    caps={c:{str(r['gene_id']):r['generation_cpu_seconds'] for r in read_jsonl(base/'evaluations'/f'{c}_diffusion_16'/'per_gene.jsonl')} for c in targets}
    protocol={'seed':42,'workers':8,'threads_per_worker':4,'candidate_count_cap':None,
              'cohorts':{c:base_protocol['cohorts'][c]['targets'] for c in targets},
              'reason':'All 150 capped direct runs hit 4096 well before deadline. Remove that implementation cap without changing score model, selection objective, references, or test candidate budget.',
              'deadline':'Same stored per-gene diffusion16 process-CPU budget; checkpoint every512 pool additions plus final reserved selection. Only on-time output counts.',
              'source_hashes':{str(p.relative_to(ROOT)):sha256(p) for p in [Path(__file__),HERE/'controls.py',HERE/'run.py',base/'protocol.json']},
              'budget_source_hashes':{c:sha256(base/'evaluations'/f'{c}_diffusion_16'/'per_gene.jsonl') for c in targets}}
    pp=out/'protocol.json'
    if pp.exists():assert json.loads(pp.read_text())==protocol
    else:atomic_json(pp,protocol)
    jobs=[]
    for cohort,rows in targets.items():
        assert {str(r['gene_id']) for r in rows}==set(map(str,protocol['cohorts'][cohort]))
        (out/'genes'/cohort).mkdir(parents=True,exist_ok=True)
        # Longest CPU-budget jobs first to avoid a straggler at the end.
        for target in rows:
            gene=str(target['gene_id']);path=out/'genes'/cohort/f'{gene}.json'
            if path.exists():assert json.loads(path.read_text())['protocol_sha256']==sha256(pp)
            else:jobs.append((cohort,target,caps[cohort][gene],str(path),sha256(pp)))
    jobs.sort(key=lambda x:-x[2])
    atomic_json(out/'status.json',{'stage':'generation','remaining_genes':len(jobs)})
    with concurrent.futures.ProcessPoolExecutor(max_workers=8,mp_context=multiprocessing.get_context('spawn'),initializer=initialize) as ex:
        for number,res in enumerate(ex.map(worker,jobs,chunksize=1)):
            print(json.dumps(res),flush=True)
            # Count actual committed files; map returns in submitted order.
            done=len(list((out/'genes').glob('*/*.json')))
            atomic_json(out/'status.json',{'stage':'generation','completed_genes':done,'total_genes':150})
    directories=[]
    for cohort,rows in targets.items():
        records=[json.loads((out/'genes'/cohort/f"{r['gene_id']}.json").read_text()) for r in rows]
        dest=out/'evaluations'/f'{cohort}_direct_equal_cpu_uncapped';dest.mkdir(parents=True,exist_ok=True)
        write_json(dest/'manifest.json',{'method':'direct_equal_cpu_uncapped','seed':42,'cohort':cohort,'protocol_sha256':sha256(pp)})
        write_jsonl(dest/'candidates.jsonl',[c for r in records for c in r['candidates']])
        write_jsonl(dest/'per_gene.jsonl',[r['ledger'] for r in records])
        write_json(dest/'completion.json',{'all_expected_targets_evaluated':True,'targets':75,
                   'candidates_sha256':sha256(dest/'candidates.jsonl'),'per_gene_sha256':sha256(dest/'per_gene.jsonl')})
        directories.append(dest)
    atomic_json(out/'status.json',{'stage':'full_reference_screen'})
    with (ROOT/'results/bowtie_audit/mechanism_writer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        with (out/'screen.log').open('a') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts/screen_bowtie.py'),*[str(d) for d in directories],'--threads','8'],
                           env=dict(os.environ,PYTHONPATH=str(ROOT/'src')),stdout=log,stderr=subprocess.STDOUT,check=True)
    summarize(out,protocol)
    atomic_json(out/'status.json',{'stage':'complete','genes':150,'summary_sha256':sha256(out/'summary.json')})


if __name__=='__main__':
    out=HERE/'uncapped_seed42';out.mkdir(exist_ok=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:main()
        except Exception as exc:
            atomic_json(out/'status.json',{'stage':'failed','error':repr(exc)})
            raise
