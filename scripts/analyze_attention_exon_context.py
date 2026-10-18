"""Partition recorded RNA/genome-screened candidates by verified exon geometry."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src')]
from diffuprimer.io import read_jsonl, sha256, write_json, write_jsonl, utc_now
from experiments.attention.analyze import load_source, PendingRun

METHODS = [
    ('primer3', 'results/evaluations/primer3_pilot_test_seed20261013_amp80-200'),
    ('qprimer', 'experiments/attention/results/evaluations/qprimer_design_80-200_test_seed20261013'),
    ('primercast', 'experiments/attention/results/evaluations/primercast_design_80-200_test_seed20261013'),
    ('cnn_direct', 'experiments/temperature_4/results/evaluations/direct_weak_partial_pilot_test_seed20261013_amp80-200_oligo_guided_joint_tm'),
    ('cnn_diffusion', 'experiments/temperature_4/results/evaluations/diffusion_weak_partial_pilot_test_seed20261013_amp80-200_oligo_guided_joint_tm'),
    ('attention_rope_direct', 'experiments/attention/results/evaluations/attention_direct_test_design_full_amp80-200_steps16_seed20261013'),
    ('attention_rope_diffusion', 'experiments/attention/results/evaluations/attention_diffusion_test_design_full_amp80-200_steps16_seed20261013'),
    ('attention_hybrid_direct', 'experiments/hybrid_attention/results/evaluations/hybrid_attention_direct_test_design_full_amp80-200_steps16_seed20261013'),
    ('attention_hybrid_diffusion', 'experiments/hybrid_attention/results/evaluations/hybrid_attention_diffusion_test_design_full_amp80-200_steps16_seed20261013'),
]
CATEGORIES = ['primer_spans_junction', 'amplicon_crosses_junction_only', 'same_exon', 'terminal_polyA', 'unmapped']


def geometry(row, context):
    if 'exons' not in context:
        return {'category': 'unmapped', 'mapping_status': context['status']}
    fs, fe, rs, re = [row[k] for k in ('forward_start', 'forward_end', 'reverse_start', 'reverse_end')]
    boundaries = context['boundaries']
    fj = [j for j in boundaries if fs < j < fe]
    rj = [j for j in boundaries if rs < j < re]
    aj = [j for j in boundaries if fs < j < re]
    tail = fe > context['genomically_annotated_RNA_length'] or re > context['genomically_annotated_RNA_length']
    category = 'terminal_polyA' if tail else 'primer_spans_junction' if fj or rj else 'amplicon_crosses_junction_only' if aj else 'same_exon'
    return {'category': category, 'mapping_status': context['status'],
            'forward_junctions': fj, 'reverse_junctions': rj, 'amplicon_junctions': aj,
            'forward_terminal_five_bases_span_junction': any(fe - 5 < j < fe for j in fj),
            'reverse_terminal_five_bases_span_junction': any(rs < j < rs + 5 for j in rj)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='experiments/attention/results/exon_interpretation')
    parser.add_argument('--allow-pending', action='store_true')
    args = parser.parse_args()
    out = ROOT / args.output
    out.mkdir(parents=True, exist_ok=True)
    target_path = ROOT / 'data/processed/targets_pilot.jsonl'
    context_path = ROOT / 'results/frozen_target_exon_context.json'
    verification_path = ROOT / 'results/spliced_reference_verification.json'
    targets = {r['gene_id']: r for r in read_jsonl(target_path) if r['split'] == 'test'}
    contexts = json.loads(context_path.read_text())
    verification = json.loads(verification_path.read_text())
    assert len(targets) == 75
    assert verification['target_sha256'] == sha256(target_path)
    assert verification['exon_context_sha256'] == sha256(context_path)
    assert verification['reference_manifest_sha256'] == sha256(ROOT / 'data/processed/blastdb/genome.manifest.json')
    assert verification['total_differing_bases'] == 0
    sources, pending, per_gene, per_candidate, groups = [], [], [], [], []
    for method, raw in METHODS:
        directory = ROOT / raw
        try:
            manifest, candidates, ledger, screens, hashes = load_source(directory, targets, 20261013)
        except PendingRun as exc:
            pending.append({'method': method, 'directory': raw, 'reason': str(exc)})
            continue
        counts = {gene: {cat: {'RNA_pass': 0, 'retained': 0} for cat in CATEGORIES} for gene in targets}
        seen = {gene: set() for gene in targets}
        for row in sorted(candidates, key=lambda r: (int(r['gene_id']), r.get('proposal_index', 0), r['candidate_id'])):
            if not row['assessment']['computational_feasible']:
                continue
            screen = screens[row['candidate_id']]
            context = contexts[row['transcript_accession_version']]
            geo = geometry(row, context)
            pair = (row['forward_primer'], row['reverse_primer'])
            unique = pair not in seen[row['gene_id']]
            if screen['gene_specificity_screen_pass']:
                seen[row['gene_id']].add(pair)
            record = {'method': method, 'gene_id': row['gene_id'], 'gene_symbol': targets[row['gene_id']]['gene_symbol'],
                      'candidate_id': row['candidate_id'], 'proposal_index': row.get('proposal_index', 0),
                      'transcript_accession_version': row['transcript_accession_version'],
                      **{k: row[k] for k in ('forward_start', 'forward_end', 'reverse_start', 'reverse_end', 'forward_primer', 'reverse_primer', 'amplicon_length')},
                      **geo, 'RNA_pass': bool(screen['gene_specificity_screen_pass']),
                      'genomic_product_flag': bool(screen['genomic_DNA_product_detected']),
                      'retained': bool(screen['gene_specificity_screen_pass'] and not screen['genomic_DNA_product_detected']),
                      'unique_RNA_passing_synthesis_pair': bool(unique and screen['gene_specificity_screen_pass']),
                      'source_directory': raw}
            per_candidate.append(record)
            if screen['gene_specificity_screen_pass'] and unique:
                counts[row['gene_id']][geo['category']]['RNA_pass'] += 1
                counts[row['gene_id']][geo['category']]['retained'] += not screen['genomic_DNA_product_detected']
        for gene in sorted(targets, key=int):
            per_gene.append({'method': method, 'gene_id': gene, 'gene_symbol': targets[gene]['gene_symbol'], 'categories': counts[gene]})
        for category in CATEGORIES:
            total = sum(counts[gene][category]['RNA_pass'] for gene in targets)
            retained = sum(counts[gene][category]['retained'] for gene in targets)
            groups.append({'method': method, 'category': category, 'genes': 75,
                           'unique_RNA_pass_pairs': total, 'unique_retained_pairs': retained,
                           'mean_RNA_pass_pairs_per_gene': total / 75, 'mean_retained_pairs_per_gene': retained / 75,
                           'pooled_retained_fraction': retained / total if total else None})
        sources.append({'method': method, 'directory': raw, 'manifest_method': manifest['method'], **hashes})
    if pending and not args.allow_pending:
        raise PendingRun('Required screens still pending: ' + json.dumps(pending))
    write_jsonl(out / 'per_candidate.jsonl', per_candidate)
    write_jsonl(out / 'per_gene.jsonl', per_gene)
    write_json(out / 'summary.json', {
        'at': utc_now(), 'status': 'complete' if not pending else 'pending_sources', 'seed': 20261013,
        'cohort': 'original frozen 75 test genes', 'total_genes': 75, 'budget': 20,
        'groups': groups, 'inputs': sources, 'pending_sources': pending,
        'target_sha256': sha256(target_path), 'exon_context_sha256': sha256(context_path),
        'spliced_reference_verification_sha256': sha256(verification_path),
        'genome_manifest_sha256': verification['reference_manifest_sha256'],
        'reference_assembly': 'GCF_000001635.27 / GRCm39',
        'GFF_source': 'data/raw/refseq/GCF_000001635.27_GRCm39_genomic.gff.gz',
        'GFF_sha256': sha256(ROOT / 'data/raw/refseq/GCF_000001635.27_GRCm39_genomic.gff.gz'),
        'deduplication': 'Unique synthesis-sequence pairs within each gene; first RNA-passing occurrence in native proposal order determines descriptive geometry.',
        'screen_definition': 'Actual recorded RNA pass and genomic product flags. Exon geometry never substitutes for either flag.',
        'geometry_source': 'Frozen RefSeq GRCm39 exon mapping, with 500 genomic transcript portions previously verified base for base.',
        'analysis_role': 'Optional biological explanation within the existing reference-specificity experiment; all 75 genes retained.'
    })
    with (out / 'groups.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(groups[0]))
        writer.writeheader(); writer.writerows(groups)
    print(json.dumps({'completed_methods': len(sources), 'pending_methods': [r['method'] for r in pending], 'groups': len(groups)}))


if __name__ == '__main__':
    main()
