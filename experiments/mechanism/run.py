#!/usr/bin/env python3
"""Frozen single-seed mechanism study, resumable by complete target gene."""
from __future__ import annotations
import argparse
import dataclasses
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'experiments/hybrid_attention/src'), str(ROOT / 'scripts'), str(HERE)]
import numpy as np
import torch
from diffuprimer.cnn_model import ConditionalDenoiser
from diffuprimer import guidance
from diffuprimer.constraints import Conditions, assess_pair, oligo_properties
from diffuprimer.io import read_jsonl, write_json, write_jsonl, sha256, utc_now
from diffuprimer.metrics import summarize_candidates
from span_cache import install
from controls import iterative, strong, timed_direct, first_scores

SEED = 42
CHECKPOINT_SEED = 20261013  # Historical file identity, NEVER a newly chosen seed.
METHODS = {
    'direct_beam': {'model': 'direct', 'kind': 'iterative', 'steps': 1},
    'diffusion_1': {'model': 'diffusion', 'kind': 'iterative', 'steps': 1},
    'diffusion_4': {'model': 'diffusion', 'kind': 'iterative', 'steps': 4},
    'diffusion_16': {'model': 'diffusion', 'kind': 'iterative', 'steps': 16},
    'frozen_4': {'model': 'diffusion', 'kind': 'iterative', 'steps': 4, 'frozen': True},
    'frozen_16': {'model': 'diffusion', 'kind': 'iterative', 'steps': 16, 'frozen': True},
    'direct_full_wor': {'model': 'direct', 'kind': 'strong', 'pool_size': 20},
    'diffusion1_full_wor': {'model': 'diffusion', 'kind': 'strong', 'pool_size': 20},
    'direct_pool320': {'model': 'direct', 'kind': 'strong', 'pool_size': 320, 'diverse': True},
    'diffusion1_pool320': {'model': 'diffusion', 'kind': 'strong', 'pool_size': 320, 'diverse': True},
    'direct_equal_cpu': {'model': 'direct', 'kind': 'timed', 'reference': 'diffusion_16', 'max_pool': 4096},
}


def stable_seed(cohort, gene):
    return int.from_bytes(hashlib.sha256(f'{SEED}:{cohort}:{gene}'.encode()).digest()[:8], 'little')


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    write_json(temp, value)
    temp.replace(path)


def source_hashes():
    paths = [HERE/'controls.py', HERE/'run.py', ROOT/'scripts/span_cache.py']
    paths += [ROOT/'experiments/hybrid_attention/src/diffuprimer'/f'{n}.py'
              for n in ['cnn_model','decode','guidance','constraints','sequence','metrics','io']]
    return {str(p.relative_to(ROOT)): sha256(p) for p in paths}


def make_candidates(name, target, pairs, assessments):
    output = []
    for i, (pair, assessment) in enumerate(zip(pairs, assessments)):
        row = {'gene_id': target['gene_id'], 'method': name, 'proposal_index': i,
               'candidate_id': f"{name}:{target['gene_id']}:{i}",
               'transcript_accession_version': target['transcript_accession_version'],
               'structural_valid': pair is not None, 'assessment': assessment,
               'requested_amplicon_min': 80, 'requested_amplicon_max': 200}
        if pair is not None:
            row.update(pair.as_dict(target['template_sequence']))
        output.append(row)
    return output


def run(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(args.threads)
    torch.set_num_interop_threads(1)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.use_deterministic_algorithms(True)
    install()
    conditions = Conditions(**json.loads((ROOT/'configs/pilot.json').read_text())['conditions'])
    target_paths = {'test': ROOT/'data/processed/targets_pilot.jsonl',
                    'external_test': ROOT/'data/processed/attention_external_targets.jsonl'}
    cohorts = {k: sorted([r for r in read_jsonl(p) if k != 'test' or r['split'] == 'test'],
                         key=lambda r: str(r['gene_id'])) for k, p in target_paths.items()}
    if args.smoke:
        # Smoke only uses validation genes and never contributes to final cohorts.
        cohorts = {'validation_smoke': [r for r in read_jsonl(target_paths['test']) if r['split'] == 'validation'][:2]}
        target_paths = {'validation_smoke': target_paths['test']}
    else:
        assert all(len(rows) == 75 for rows in cohorts.values())
        assert not ({r['gene_id'] for r in cohorts['test']} & {r['gene_id'] for r in cohorts['external_test']})
    checkpoints, models = {}, {}
    for mode in ['direct', 'diffusion']:
        path = ROOT/f'results/checkpoints/{mode}_weak_partial_pilot_seed{CHECKPOINT_SEED}.pt'
        state = torch.load(path, map_location='cpu', weights_only=True)
        model = ConditionalDenoiser(**state['model_spec'])
        model.load_state_dict(state['state_dict'])
        models[mode] = model.eval()
        checkpoints[mode] = {'path': str(path.relative_to(ROOT)), 'sha256': sha256(path),
                             'historical_training_seed': state['seed'], 'epoch': state.get('epoch')}
    protocol = {
        'version': 1, 'seed': SEED, 'seed_policy': 'Constant integer 42; per-gene seeds are SHA256 of seed/cohort/gene; never clock-derived',
        'checkpoint_policy': 'Reuse ONE existing direct/diffusion CNN checkpoint pair; historical training seed is provenance only; no retraining',
        'checkpoints': checkpoints, 'methods': METHODS, 'candidate_budget': 20,
        'temperature': 4.0, 'conditions': dataclasses.asdict(conditions), 'cpu_threads': args.threads,
        'device': 'cpu', 'precision': 'float32', 'code_sha256': source_hashes(),
        'cohorts': {k: {'targets': [r['gene_id'] for r in rs], 'target_file_sha256': sha256(target_paths[k])} for k, rs in cohorts.items()},
        'selection_policy': 'Predeclared defaults; no tuning on either test cohort; no specificity or exon annotation used for proposal selection',
        'equal_compute': 'Per-gene process CPU seconds for diffusion_16 network+decoder+pair assessment; direct output committed by that deadline. Shared span setup excluded from BOTH; cold pair-property caches per method; actual overrun reported separately. This is an inference CPU-budget comparison, not equal end-to-end screen cost or equal FLOPs.',
        'pool_selection': 'All compatible pairs; weighted sampling without replacement; sequence deduplication; pool320 uses feasible-only farthest-first Jaccard, first/minimum-Tm and deterministic ties',
        'statistics': 'Gene-paired bootstrap, seed42, 10000 resamples; two cohorts separate; one training-checkpoint pair, no training-seed uncertainty estimate',
        'primary_endpoint': 'Unique RNA-specific pairs without predicted genomic products per gene',
        'secondary_endpoints': ['one/two retained-pair gene coverage', 'retained-pair positional diversity', 'temperature penalty', 'CPU/wall costs'],
        'runtime_versions': {'python': sys.version, 'torch': torch.__version__, 'numpy': np.__version__, 'primer3': __import__('primer3').__version__},
    }
    protocol_path = out/'protocol.json'
    if protocol_path.exists():
        if json.loads(protocol_path.read_text()) != protocol:
            raise ValueError('Frozen protocol changed; use a new output directory, never overwrite prior runs')
    else:
        atomic_json(protocol_path, protocol)
    # Warm both models outside measured regions; no test-target-dependent tuning.
    for m in models.values():
        for _ in range(2):
            first_scores(m, 'ACGT'*128, conditions)
    started = time.monotonic()
    status_path = out/'status.json'
    all_dirs = []
    for cohort, targets in cohorts.items():
        gene_dir = out/'genes'/cohort
        gene_dir.mkdir(parents=True, exist_ok=True)
        for number, target in enumerate(targets):
            gene = str(target['gene_id'])
            record_path = gene_dir/f'{gene}.json'
            if record_path.exists():
                saved = json.loads(record_path.read_text())
                if saved['protocol_sha256'] != sha256(protocol_path) or set(saved['methods']) != set(METHODS):
                    raise ValueError('Incomplete or stale gene record')
                continue
            template = target['template_sequence']
            prep_wall, prep_cpu = time.perf_counter(), time.process_time()
            span = guidance.feasible_oligo_spans(template, conditions)
            temperatures = guidance.compatible_span_temperatures(template, conditions, span)
            prep = {'cached_span_setup_wall_seconds': time.perf_counter()-prep_wall,
                    'cached_span_setup_cpu_seconds': time.process_time()-prep_cpu}
            record = {'cohort': cohort, 'gene_id': gene, 'template_length': len(template),
                      'protocol_sha256': sha256(protocol_path), 'methods': {}, **prep}
            # Rotate ordinary method order by gene. Timed direct depends on the
            # completed reference cost and therefore runs after the rotation.
            ordinary = [k for k,v in METHODS.items() if v['kind'] != 'timed']
            shift = number % len(ordinary)
            order = ordinary[shift:] + ordinary[:shift] + ['direct_equal_cpu']
            for name in order:
                cfg = METHODS[name]
                oligo_properties.cache_clear()
                wall, cpu = time.perf_counter(), time.process_time()
                if cfg['kind'] == 'iterative':
                    pairs, details = iterative(models[cfg['model']], template, conditions, span[1], temperatures,
                                               stable_seed(cohort,gene), cfg['steps'], cfg.get('frozen',False))
                    assessments = [assess_pair(template,p,conditions) if p is not None else
                                   {'computational_feasible':False,'failures':['invalid_structure']} for p in pairs]
                elif cfg['kind'] == 'strong':
                    pairs, assessments, details = strong(models[cfg['model']], template, conditions, span[1], temperatures,
                                                          stable_seed(cohort,gene), cfg['pool_size'], cfg.get('diverse',False))
                else:
                    cap = record['methods']['diffusion_16']['ledger']['generation_cpu_seconds']
                    pairs, assessments, details = timed_direct(models['direct'], template, conditions, span[1], temperatures,
                                                               stable_seed(cohort,gene), cap, cfg['max_pool'])
                cpu_elapsed, wall_elapsed = time.process_time()-cpu, time.perf_counter()-wall
                candidates = make_candidates(f'{cohort}_{name}',target,pairs,assessments)
                ledger = {'gene_id':gene,'method':name,'conditions':dataclasses.asdict(conditions),
                          **summarize_candidates(candidates,20), 'generation_cpu_seconds':cpu_elapsed,
                          'generation_wall_seconds':wall_elapsed, **details}
                record['methods'][name] = {'candidates':candidates,'ledger':ledger}
            atomic_json(record_path,record)
            atomic_json(status_path, {'stage':'generation','cohort':cohort,'completed_in_cohort':number+1,
                                     'total_in_cohort':len(targets),'last_gene':gene,'elapsed_seconds':time.monotonic()-started,'updated_at':utc_now()})
            print(json.dumps({'cohort':cohort,'gene':gene,'completed':number+1,'total':len(targets),
                              'seconds':round(sum(x['ledger']['generation_wall_seconds'] for x in record['methods'].values()),2)}),flush=True)
        records = [json.loads((gene_dir/f"{t['gene_id']}.json").read_text()) for t in targets]
        for name,cfg in METHODS.items():
            dest = out/'evaluations'/f'{cohort}_{name}'
            dest.mkdir(parents=True,exist_ok=True)
            write_json(dest/'manifest.json',{'method':name,'arguments':{'seed':SEED,'budget':20,'split':cohort,'steps':cfg.get('steps',1),
                       'amplicon_min':80,'amplicon_max':200,'temperature':4.0},'conditions':dataclasses.asdict(conditions),
                       'expected_targets':[r['gene_id'] for r in records], 'protocol_sha256':sha256(protocol_path),
                       'checkpoint_sha256':checkpoints[cfg['model']]['sha256'],'architecture':'conditional_dilated_cnn'})
            write_jsonl(dest/'candidates.jsonl',[c for r in records for c in r['methods'][name]['candidates']])
            write_jsonl(dest/'per_gene.jsonl',[r['methods'][name]['ledger'] for r in records])
            write_json(dest/'completion.json',{'all_expected_targets_evaluated':True,'targets':len(records),
                       'candidates_sha256':sha256(dest/'candidates.jsonl'),'per_gene_sha256':sha256(dest/'per_gene.jsonl')})
            all_dirs.append(dest)
    if not args.skip_screen:
        atomic_json(status_path,{'stage':'full_reference_screen','runs':len(all_dirs),'updated_at':utc_now()})
        env = dict(os.environ,PYTHONPATH=str(ROOT/'src'))
        with (ROOT/'results/bowtie_audit/mechanism_writer.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            with (out/'screen.log').open('a') as log:
                subprocess.run([sys.executable,str(ROOT/'scripts/screen_bowtie.py'),*[str(p) for p in all_dirs],
                                '--threads',str(args.screen_threads)],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        atomic_json(status_path,{'stage':'analysis','updated_at':utc_now()})
        subprocess.run([sys.executable,str(HERE/'analyze.py'),'--output',str(out)],check=True,cwd=ROOT)
    atomic_json(status_path,{'stage':'generation_complete' if args.skip_screen else 'complete',
                            'runs':len(all_dirs),'elapsed_seconds':time.monotonic()-started,'updated_at':utc_now()})


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=HERE/'results_seed42')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--screen-threads',type=int,default=8)
    p.add_argument('--skip-screen',action='store_true')
    p.add_argument('--smoke',action='store_true')
    a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    with (a.output/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            run(a)
        except Exception as exc:
            atomic_json(a.output/'status.json',{'stage':'failed','error':repr(exc),'updated_at':utc_now()})
            raise
