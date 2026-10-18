#!/usr/bin/env python3
"""Complete the predeclared RNA-only selection endpoint while DNA audit continues."""
import hashlib
import json
import os
import sqlite3
import subprocess
import time
from pathlib import Path

from bowtie_audit import alignment_command
from bowtie_sqlite_cache import AlignmentCache
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.specificity import classify_transcript_products


root=Path(__file__).resolve().parents[1]
output=root/'results/bowtie_audit/rna_validation'
output.mkdir(parents=True,exist_ok=True)
status_path=output/'status.json'


def state(stage,**extra):
    write_json(status_path,{'stage':stage,'at':utc_now(),**extra})
    print(json.dumps({'stage':stage,'at':utc_now(),**extra}),flush=True)


try:
    state('waiting_for_complete_preserved_RNA_alignments')
    original=root/'results/bowtie_audit/alignments.sqlite'
    while True:
        connection=sqlite3.connect(f'file:{original}?mode=ro',uri=True)
        complete=connection.execute('SELECT count(*) FROM imports WHERE reference=0 AND complete=1').fetchone()[0]
        if complete:
            break
        connection.close();time.sleep(5)
    database=output/'alignments.sqlite'
    snapshot_marker=output/'snapshot.json'
    if not snapshot_marker.exists():
        state('snapshotting_complete_RNA_binding_sites')
        with sqlite3.connect(database) as destination:
            connection.backup(destination,pages=1024,sleep=.01)
        write_json(snapshot_marker,{'at':utc_now(),'source':str(original),'identity_sha256':sha256(root/'results/bowtie_audit/identity.json'),
            'scope':'Independent consistent SQLite snapshot; only complete RNA entries used; partial genome entries ignored'})
    connection.close()
    cache=AlignmentCache(database)
    check=json.loads((root/'results/checks/sqlite_specificity/check.json').read_text())
    if not check['exact_gene_specificity_decisions'] or check['cache_sha256'] != sha256(root/'scripts/bowtie_sqlite_cache.py'):
        raise ValueError('Exact specificity check stale')
    inputs=[];sequences=set()
    for method in ['direct','diffusion']:
        for seed in [20261013,20261014,20261015]:
            for label in [None,'05','2','4']:
                base=root/('experiments/temperature_'+label if label else 'experiments/joint_tm')
                directory=base/f'results/evaluations/{method}_weak_partial_pilot_validation_seed{seed}_amp80-200_oligo_guided_joint_tm'
                completion=json.loads((directory/'completion.json').read_text())
                if not completion['all_expected_targets_evaluated'] or completion['candidates_sha256'] != sha256(directory/'candidates.jsonl'):
                    raise ValueError('Validation designs incomplete or changed')
                candidates=read_jsonl(directory/'candidates.jsonl')
                sequences.update(s for r in candidates if r['assessment']['computational_feasible'] for s in [r['forward_primer'],r['reverse_primer']])
                inputs.append((directory,candidates))
    registered=cache.register(sequences)
    missing=sorted(s for s in sequences if not registered[s][1])
    if missing:
        qmap={'Q'+hashlib.sha256(s.encode()).hexdigest():s for s in missing}
        queries=output/'queries.fasta';queries.write_text(''.join(f'>{qid}\n{s}\n' for qid,s in qmap.items()))
        sam=output/'transcriptome.sam'
        command=alignment_command(root/'tools/bowtie-1.3.1-linux-x86_64/bowtie',root/'data/processed/bowtiedb/transcriptome',queries,sam,2)
        marker=output/'alignment_completed.json'
        valid=False
        if marker.exists():
            prior=json.loads(marker.read_text())
            valid=prior['command']==command and prior['queries_sha256']==sha256(queries) and prior['sam_sha256']==sha256(sam)
        if not valid:
            state('aligning_missing_validation_RNA_oligos',new_oligos=len(missing))
            with (output/'alignment.log').open('w') as handle:
                subprocess.run(command,check=True,stdout=handle,stderr=subprocess.STDOUT)
            write_json(marker,{'at':utc_now(),'command':command,'queries_sha256':sha256(queries),'sam_sha256':sha256(sam),
                'database_manifest_sha256':sha256(root/'data/processed/bowtiedb/transcriptome.manifest.json')})
        cache.ingest(sam,qmap,0)
        registered=cache.register(sequences)
    mapping=json.loads((root/'data/processed/blastdb/transcript_gene_map.json').read_text())
    state('classifying_complete_validation_RNA_pairs',runs=len(inputs))
    for directory,candidates in inputs:
        records=[]
        for row in candidates:
            if not row['assessment']['computational_feasible']:
                records.append({'candidate_id':row['candidate_id'],'status':'not_screened_computational_failure'})
                continue
            f,r=registered[row['forward_primer']],registered[row['reverse_primer']]
            if not f[1] or not r[1]:
                raise ValueError('Incomplete RNA binding sites')
            result=classify_transcript_products(cache.products(f[0],r[0],0),row['gene_id'],mapping)
            records.append({'candidate_id':row['candidate_id'],'status':'RNA_full_length_ungapped_mismatch_rule_assessed',**result})
        path=directory/'specificity_bowtie_rna.jsonl';write_jsonl(path,records)
        write_json(directory/'specificity_bowtie_rna_status.json',{'status':'completed','completed_at':utc_now(),
            'candidate_sha256':sha256(directory/'candidates.jsonl'),'output_sha256':sha256(path),
            'identity_sha256':sha256(root/'results/bowtie_audit/identity.json'),
            'product_evidence_format':'one_coordinate_witness_per_subject_and_primer_role','scope':'RNA-only endpoint; genome audit separate'})
    command=[str(root/'.venv/bin/python'),'scripts/summarize_specificity_audit.py',*[str(p) for p,_ in inputs],
        '--split','validation','--rna-only','--output','results/bowtie_audit/validation_RNA_summary.json']
    with (output/'statistics.log').open('w') as handle:
        subprocess.run(command,cwd=root,env=dict(os.environ,PYTHONPATH=str(root/'src')),check=True,stdout=handle,stderr=subprocess.STDOUT)
    state('validation_RNA_endpoints_complete',runs=len(inputs))
except Exception as exc:
    state('RNA_audit_error',error=str(exc))
    raise
