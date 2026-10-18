"""Select a traceable exon-layout illustration without altering cohort statistics."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from diffuprimer.io import read_jsonl, sha256, write_json
from diffuprimer.sequence import reverse_complement

OUT = ROOT / 'experiments/attention/results/exon_interpretation'
DEST = ROOT / 'paper/figures/cnn_biological_analysis'


def genomic_segments(start, end, exons):
    segments = []
    for exon in exons:
        lo, hi = max(start, exon['RNA_start']), min(end, exon['RNA_end'])
        if lo >= hi:
            continue
        if exon['strand'] == '+':
            gs = exon['genome_start'] + lo - exon['RNA_start']
            ge = exon['genome_start'] + hi - exon['RNA_start']
        else:
            gs = exon['genome_end'] - (hi - exon['RNA_start'])
            ge = exon['genome_end'] - (lo - exon['RNA_start'])
        segments.append({'RNA_start': lo, 'RNA_end': hi, 'genome_start': gs, 'genome_end': ge,
                         'chromosome': exon['chromosome'], 'strand': exon['strand']})
    return segments


def main():
    summary = json.loads((OUT / 'summary.json').read_text())
    if summary['status'] != 'complete':
        raise ValueError('Exon interpretation still awaits source screens')
    contexts = json.loads((ROOT / 'results/frozen_target_exon_context.json').read_text())
    targets = {r['gene_id']: r for r in read_jsonl(ROOT / 'data/processed/targets_pilot.jsonl')}
    records = read_jsonl(OUT / 'per_candidate.jsonl')
    screened, raw_candidates, sources = {}, {}, {}
    for source in summary['inputs']:
        if source['method'] not in {'primer3', 'qprimer', 'cnn_diffusion'}:
            continue
        directory = ROOT / source['directory']
        screened[source['method']] = {r['candidate_id']: r for r in read_jsonl(directory / 'specificity_bowtie.jsonl')}
        raw_candidates[source['method']] = {r['candidate_id']: r for r in read_jsonl(directory / 'candidates.jsonl')}
        sources[source['method']] = source
    pool = []
    for baseline in [r for r in records if r['method'] in {'primer3', 'qprimer'} and r['proposal_index'] == 0
                     and r['RNA_pass'] and r['genomic_product_flag'] and r['category'] == 'same_exon']:
        context = contexts[baseline['transcript_accession_version']]
        mapped = genomic_segments(baseline['forward_start'], baseline['reverse_end'], context['exons'])
        assert len(mapped) == 1
        segment = mapped[0]
        witnesses = [w for w in screened[baseline['method']][baseline['candidate_id']]['genomic_potential_amplicons']
                     if w['subject'] == segment['chromosome'] and w['start'] == segment['genome_start']
                     and w['end'] == segment['genome_end']]
        if not witnesses:
            continue
        for hybrid in [r for r in records if r['method'] == 'cnn_diffusion' and r['gene_id'] == baseline['gene_id']
                       and r['retained'] and r['unique_RNA_passing_synthesis_pair'] and r['category'] == 'primer_spans_junction']:
            if hybrid['transcript_accession_version'] != baseline['transcript_accession_version']:
                continue
            width = max(baseline['reverse_end'], hybrid['reverse_end']) - min(baseline['forward_start'], hybrid['forward_start'])
            pool.append({'baseline': baseline, 'hybrid': hybrid, 'RNA_union_width_bp': width,
                         'baseline_genomic_product_witness': witnesses[0]})
    if not pool:
        raise ValueError('No actual same-gene top-ranked baseline / retained-junction hybrid case available')
    # Prefer readable input-sized windows and nearly matched product lengths; these are display criteria, not effect-size criteria.
    pool.sort(key=lambda r: (0 if r['baseline']['method'] == 'primer3' else 1, 0 if r['RNA_union_width_bp'] <= 512 else 1, abs(r['baseline']['amplicon_length'] - r['hybrid']['amplicon_length']), r['RNA_union_width_bp'],
                             int(r['baseline']['gene_id']), r['hybrid']['proposal_index'], r['hybrid']['candidate_id']))
    chosen = pool[0]
    chosen['cnn_diffusion'] = chosen.pop('hybrid')
    gene = chosen['baseline']['gene_id']
    target = targets[gene]
    context = contexts[target['transcript_accession_version']]
    for role in ['baseline', 'cnn_diffusion']:
        row = chosen[role]
        raw = raw_candidates[row['method']][row['candidate_id']]
        assert target['template_sequence'][row['forward_start']:row['forward_end']] == row['forward_primer']
        assert reverse_complement(target['template_sequence'][row['reverse_start']:row['reverse_end']]) == row['reverse_primer']
        row['forward_genomic_segments'] = genomic_segments(row['forward_start'], row['forward_end'], context['exons'])
        row['reverse_genomic_segments'] = genomic_segments(row['reverse_start'], row['reverse_end'], context['exons'])
        row['thermodynamics'] = {k: raw['assessment'][k] for k in ['tm_difference', 'shared_temperature_penalty']}
        row['forward_tm'] = raw['assessment']['forward_properties']['tm']
        row['reverse_tm'] = raw['assessment']['reverse_properties']['tm']
        row['screen_record'] = screened[row['method']][row['candidate_id']]
    start = max(0, min(chosen['baseline']['forward_start'], chosen['cnn_diffusion']['forward_start']) - 25)
    end = min(len(target['template_sequence']), max(chosen['baseline']['reverse_end'], chosen['cnn_diffusion']['reverse_end']) + 25)
    write_json(DEST / 'illustrative_case.json', {
        'gene_id': gene, 'gene_symbol': target['gene_symbol'], 'transcript_accession_version': target['transcript_accession_version'],
        'strand': context['exons'][0]['strand'], 'chromosome': context['exons'][0]['chromosome'],
        'selection_rule': 'Top-ranked native Primer3 or qPrimer proposal (index0), RNA pass, same annotated exon, and recorded genomic witness exactly matching its annotated product. Pair with a retained CNN diffusion proposal spanning an exon junction on the same transcript. Prefer Primer3; prefer a combined RNA display window no wider than the 512-bp model input; choose smallest RNA product-length difference, then shortest combined RNA plot window, then numeric gene ID and CNN proposal order. No cohort records or aggregate outcomes are excluded.',
        'eligible_case_pairs': len(pool), 'eligible_genes': len({r['baseline']['gene_id'] for r in pool}),
        'comparison_role': 'Illustrative positional mechanism case, not an estimate of method-level superiority.',
        'RNA_plot_window': [start, end], 'RNA_window_sequence': target['template_sequence'][start:end],
        'all_annotated_exons': context['exons'], 'annotated_junctions': context['boundaries'],
        **chosen, 'source_inputs': [sources[chosen['baseline']['method']], sources['cnn_diffusion']],
        'exon_summary_sha256': sha256(OUT / 'summary.json'), 'exon_context_sha256': summary['exon_context_sha256'],
        'coordinate_convention': '0-based half-open RNA and genomic intervals; displayed genomic endpoints explicitly noted.'
    })
    print(json.dumps({'gene': target['gene_symbol'], 'accession': target['transcript_accession_version'],
                      'eligible_cases': len(pool), 'RNA_union_bp': chosen['RNA_union_width_bp']}))


if __name__ == '__main__':
    main()
