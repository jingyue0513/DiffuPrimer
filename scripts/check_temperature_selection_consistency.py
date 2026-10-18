#!/usr/bin/env python3
"""Verify independent RNA-only selection against the complete reference audit."""
import json
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json


root=Path(__file__).resolve().parents[1]
selection_path=root/'results/temperature_tuning/selection.json'
selection=json.loads(selection_path.read_text())
summary=json.loads((root/'results/bowtie_audit/validation_summary.json').read_text())
candidate_decisions=0
for source in json.loads((root/'results/temperature_tuning/selection_validation_summary.json').read_text())['inputs']:
    directory=Path(source['directory'])
    rna_path=directory/'specificity_bowtie_rna.jsonl'
    full_path=directory/'specificity_bowtie.jsonl'
    for path in [rna_path,full_path]:
        status=json.loads(path.with_name(path.stem+'_status.json').read_text())
        assert status['output_sha256']==sha256(path)
        assert status['candidate_sha256']==sha256(directory/'candidates.jsonl')
    rna={r['candidate_id']:r for r in read_jsonl(rna_path)}
    full={r['candidate_id']:r for r in read_jsonl(full_path)}
    assert rna.keys()==full.keys()
    for candidate in rna:
        a,b=rna[candidate],full[candidate]
        assert a.get('gene_specificity_screen_pass')==b.get('gene_specificity_screen_pass')
        for field in ['on_target_products','off_target_products','unresolved_products']:
            roles=lambda row:{(p['subject'],p['left_primer'],p['right_primer']) for p in row.get(field,[])}
            assert roles(a)==roles(b)
        candidate_decisions+=1
for group in selection['validation_groups']:
    matched=[r for r in summary['summaries'] if all(r[k]==group[k] for k in ['study','method','context','variant','budget','genes','seeds'])]
    assert len(matched)==1
    for key in ['unique_RNA_pass_pairs','one_RNA_pass_pair','two_RNA_pass_pairs','best_RNA_pass_temperature_penalty']:
        assert matched[0][key]==group[key]
write_json(root/'results/temperature_tuning/full_audit_consistency.json',{'at':utc_now(),
    'candidate_decisions_compared':candidate_decisions,'validation_groups_compared':len(selection['validation_groups']),
    'exact_RNA_decision_and_subject_role_agreement':True,'exact_selection_endpoint_agreement':True,
    'selection_sha256':sha256(selection_path),'complete_validation_summary_sha256':sha256(root/'results/bowtie_audit/validation_summary.json')})
print('RNA-only selection agrees exactly with complete RNA/genome audit',flush=True)
