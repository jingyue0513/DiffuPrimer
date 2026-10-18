"""Publication-size exon partitions and one actual RNA/genome layout example."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Patch
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from diffuprimer.io import sha256, write_json
from select_attention_exon_case import genomic_segments

SOURCE = ROOT / 'experiments/attention/results/exon_interpretation'
OUT = ROOT / 'paper/figures/exon_interpretation'
LABELS = {'primer3': 'Primer3', 'qprimer': 'qPrimer', 'primercast': 'PrimerCast',
          'cnn_direct': 'CNN direct', 'cnn_diffusion': 'CNN diffusion',
          'attention_rope_direct': 'RoPE direct', 'attention_rope_diffusion': 'RoPE diffusion',
          'attention_hybrid_direct': 'Hybrid direct', 'attention_hybrid_diffusion': 'DiffuPrimer (hybrid)'}
CAT_NAMES = {'primer_spans_junction': 'Primer spans junction',
             'amplicon_crosses_junction_only': 'Product crosses junction only',
             'same_exon': 'Same exon', 'terminal_polyA': 'Terminal poly(A)', 'unmapped': 'Unmapped'}
CAT_COLORS = {'primer_spans_junction': '#2A9380', 'amplicon_crosses_junction_only': '#B6A16A',
              'same_exon': '#AEBEC6', 'terminal_polyA': '#9D8AB2', 'unmapped': '#D88D7E'}
F, R, INK = '#248A75', '#D36E51', '#193D50'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads((SOURCE / 'summary.json').read_text())
    case = json.loads((SOURCE / 'illustrative_case.json').read_text())
    if summary['status'] != 'complete':
        raise ValueError('Do not plot incomplete genotype screens as outcomes')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.titlesize': 8,
                         'axes.labelsize': 7.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .55, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'svg.fonttype': 'path', 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    generated = []

    def save(fig, name, caption, alt):
        files = []
        for ext in ('pdf', 'svg', 'png'):
            path = OUT / f'{name}.{ext}'
            fig.savefig(path, dpi=350)
            files.append({'path': str(path.relative_to(ROOT)), 'sha256': sha256(path)})
        # Render the exact178mm PDF for the nearest-pixel350dpi publication asset.
        subprocess.run(['pdftoppm', '-png', '-r', '350', '-singlefile',
                        str(OUT / f'{name}.pdf'), str(OUT / name)], check=True)
        with Image.open(OUT / f'{name}.png') as rendered:
            rendered.convert('RGB').save(OUT / f'{name}.png', dpi=(350, 350))
        files[-1]['sha256'] = sha256(OUT / f'{name}.png')
        generated.append({'name': name, 'files': files, 'caption': caption, 'alt_text': alt,
                          'physical_width_mm': 178})
        plt.close(fig)

    methods = [m for m in LABELS if any(g['method'] == m for g in summary['groups'])]
    categories = [c for c in CAT_NAMES if any(g['category'] == c and g['unique_RNA_pass_pairs'] > 0 for g in summary['groups'])]
    values = {(g['method'], g['category']): g for g in summary['groups']}
    fig, axes = plt.subplots(1, 2, figsize=(178 / 25.4, 88 / 25.4))
    fig.subplots_adjust(left=.195, right=.985, top=.825, bottom=.155, wspace=.18)
    for ax, endpoint, title in zip(axes, ['mean_RNA_pass_pairs_per_gene', 'mean_retained_pairs_per_gene'],
                                    ['a  RNA-specific candidates', 'b  Retained after genome screen']):
        left = np.zeros(len(methods))
        for cat in categories:
            means = [values[(m, cat)][endpoint] for m in methods]
            ax.barh(np.arange(len(methods)), means, left=left, height=.66,
                    color=CAT_COLORS[cat], edgecolor='white', linewidth=.45)
            left += means
        for i, total in enumerate(left):
            ax.text(total + .14, i, f'{total:.2f}', va='center', fontsize=7.5, color=INK)
        ax.set_yticks(range(len(methods)), [LABELS[m] for m in methods] if ax is axes[0] else [])
        ax.invert_yaxis(); ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
        ax.set_title(title, loc='left', fontsize=8, color=INK, pad=8)
        ax.set_xlabel('Unique pairs per gene (all 75 genes)')
        ax.set_xlim(0, max(left) * 1.16)
        ax.grid(axis='x', color='#E5E9EB', linewidth=.45); ax.set_axisbelow(True)
    fig.legend(handles=[Patch(facecolor=CAT_COLORS[c], label=CAT_NAMES[c]) for c in categories],
               loc='upper center', bbox_to_anchor=(.55, .99), ncol=len(categories), frameon=False,
               fontsize=7.2, handlelength=1.3, columnspacing=1.4)
    save(fig, 'exon_geometry_partition',
         'Exon geometry of unique RNA-specific primer pairs and the subset without recorded genomic-product flags. Counts are averaged over all 75 frozen test genes, with no-design genes included. Geometry categories are mutually exclusive and use verified transcript exon coordinates. Genomic flags come from the actual full-reference screen; categories do not replace the screen. Zero-count geometry categories are retained in the underlying data and omitted from the legend.',
         'Two horizontal stacked-bar plots partition RNA-specific and retained candidate counts into primer-junction, product-junction-only and same-exon geometry for native baselines and CNN, RoPE and hybrid neural methods.')

    # The case uses transcript-order genomic layout. Intron display widths are compressed.
    lo, hi = case['RNA_plot_window']
    visible = [(i + 1, e, max(lo, e['RNA_start']), min(hi, e['RNA_end']))
               for i, e in enumerate(case['all_annotated_exons']) if e['RNA_start'] < hi and lo < e['RNA_end']]
    gap = (hi - lo) * .13
    offsets = {}
    cursor = 0
    for index, exon, start, end in visible:
        offsets[index] = {'start': cursor, 'end': cursor + end - start, 'RNA_start': start, 'RNA_end': end}
        cursor += end - start + gap
    gxmax = cursor - gap
    exon_colors = ['#BFD9D3', '#CAD4E7', '#DDD1E5', '#DBDFBF']
    colors = {index: exon_colors[k % len(exon_colors)] for k, (index, *_rest) in enumerate(visible)}
    fig = plt.figure(figsize=(178 / 25.4, 100 / 25.4))
    axg = fig.add_axes([.095, .585, .875, .245])
    axr = fig.add_axes([.095, .17, .875, .30])
    fig.text(.04, .947, f"{case['gene_symbol']}  |  {case['transcript_accession_version']}  |  same-gene illustration",
             fontsize=10, fontweight='bold', color=INK)
    r_junctions = case['hybrid']['reverse_junctions']
    if r_junctions:
        split = case['hybrid']['reverse_end'] - r_junctions[0]
        rseq = case['hybrid']['reverse_primer']
        fig.text(.97, .947, f"R 5′ {rseq[:split]} | {rseq[split:]} 3′", ha='right',
                 fontsize=7.5, fontfamily='DejaVu Sans Mono', color=R)
    fig.text(.095, .872, f"a  Genomic locus ({case['chromosome']}, {case['strand']} strand; introns compressed)", fontsize=8, color=INK)
    fig.text(.095, .505, 'b  Spliced transcript and actual primer-binding intervals', fontsize=8, color=INK)
    axg.set_xlim(-gxmax * .02, gxmax * 1.02); axg.set_ylim(-.48, 1.30); axg.axis('off')
    for index, exon, start, end in visible:
        segment = genomic_segments(start, end, [exon])[0]
        display = offsets[index]
        axg.add_patch(Rectangle((display['start'], .57), end - start, .20,
                                facecolor=colors[index], edgecolor='#647F86', linewidth=.55))
        center = (display['start'] + display['end']) / 2
        axg.text(center, .96, f'Exon {index}', ha='center', fontsize=7.5, color=INK)
        axg.text(center, .40, f"{segment['genome_start']:,}–{segment['genome_end']:,}",
                 ha='center', va='top', fontsize=7.2, color='#526B79')
    for (_, e1, *_), (index2, e2, *_) in zip(visible, visible[1:]):
        index1 = index2 - 1
        xx = (offsets[index1]['end'] + offsets[index2]['start']) / 2
        intron = abs(e2['genome_start'] - e1['genome_end']) if case['strand'] == '+' else abs(e1['genome_start'] - e2['genome_end'])
        axg.plot([offsets[index1]['end'], offsets[index2]['start']], [.67, .67], color='#98A9AE', lw=.7)
        axg.text(xx, .665, '//', ha='center', va='center', color='#526B79', fontsize=8,
                 bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .15})
        axg.text(xx, 1.18, f'{intron:,} bp intron', ha='center', fontsize=7.2, color='#526B79')
    baseline = case['baseline']; hybrid = case['hybrid']
    baseline_exon = next(index for index, _e, start, end in visible if start <= baseline['forward_start'] < baseline['reverse_end'] <= end)
    mapper = lambda rna: offsets[baseline_exon]['start'] + rna - offsets[baseline_exon]['RNA_start']
    gs, ge = mapper(baseline['forward_start']), mapper(baseline['reverse_end'])
    axg.plot([gs, ge], [-.16, -.16], color='#C18670', linewidth=2.5)
    axg.text(ge, -.31,
             f"{LABELS[baseline['method']]}: predicted {case['baseline_genomic_product_witness']['amplicon_length']} bp genomic product",
             ha='right', va='top', fontsize=7.5, color='#A26049')

    axr.set_xlim(lo, hi); axr.set_ylim(-.25, 1.65)
    axr.spines[['left', 'right', 'top']].set_visible(False); axr.set_yticks([])
    axr.set_xlabel('RNA coordinates (0-based, half-open)', fontsize=7.5)
    for index, exon, start, end in visible:
        axr.add_patch(Rectangle((start, 1.26), end - start, .17, facecolor=colors[index], edgecolor='#647F86', linewidth=.55))
        axr.text((start + end) / 2, 1.53, f'Exon {index}', ha='center', fontsize=7.5, color=INK)
    for junction in case['annotated_junctions']:
        if lo < junction < hi:
            axr.axvline(junction, color='#A8B6BC', linewidth=.75, linestyle=(0, (2, 2)), zorder=0)
            axr.text(junction, -.06, f'junction {junction}', ha='center', va='top', fontsize=7.2, color='#526B79')
    for row, y, label in [(baseline, .92, LABELS[baseline['method']]), (hybrid, .42, 'DiffuPrimer (hybrid)')]:
        axr.plot([row['forward_start'], row['reverse_end']], [y, y], color='#D1D8D9', linewidth=1)
        for start, end, col, direction in [(row['forward_start'], row['forward_end'], F, 1),
                                          (row['reverse_start'], row['reverse_end'], R, -1)]:
            a, b = (start, end) if direction == 1 else (end, start)
            axr.add_patch(FancyArrowPatch((a, y), (b, y), arrowstyle='-|>', mutation_scale=6.5,
                                         linewidth=2.2, color=col, shrinkA=0, shrinkB=0))
        axr.text(lo, y + .14, f"{label}, proposal {row['proposal_index']}  ·  {row['amplicon_length']} bp RNA product",
                 fontsize=7.5, color=INK)
    junction_sides = ' / '.join(side for side, key in [('F', 'forward_junctions'), ('R', 'reverse_junctions')] if hybrid[key])
    junction_number = (hybrid['forward_junctions'] + hybrid['reverse_junctions'])[0]
    fig.text(.095, .042, f"Hybrid {junction_sides} crosses junction {junction_number}: RNA pass; no genomic product flagged.",
             fontsize=7.5, color=F)
    save(fig, 'same_gene_exon_case',
         f"Actual primer-pair illustration for {case['gene_symbol']} ({case['transcript_accession_version']}). The top-ranked native {LABELS[baseline['method']]} pair lies within one exon and has a coordinate-matched genomic-product witness. A retained hybrid alternative spans a transcript splice junction. Both RNA pass and genome flags are read from recorded full-reference screens. The case is selected by the documented display-width rule and does not alter aggregate method comparisons. Genomic display compresses introns; numeric intervals are actual zero-based half-open coordinates.",
         'A genomic exon layout with a recorded native-baseline product is aligned conceptually with a spliced transcript showing the baseline and an actual junction-spanning hybrid alternative. Native sequences, coordinates and screen witnesses are provided in the source case JSON.')

    write_json(OUT / 'provenance.json', {'source_summary': str((SOURCE / 'summary.json').relative_to(ROOT)),
              'source_summary_sha256': sha256(SOURCE / 'summary.json'),
              'case_source': str((SOURCE / 'illustrative_case.json').relative_to(ROOT)),
              'case_sha256': sha256(SOURCE / 'illustrative_case.json'), 'cohort_genes': 75,
              'raster_dpi': 350, 'figures': generated, 'minimum_text_font_pt': 7.2,
              'interpretation': 'Descriptive biological explanation using actual screen outcomes; no physical qPCR data.'})
    print(json.dumps({'figures': [f['name'] for f in generated], 'case_gene': case['gene_symbol']}))


if __name__ == '__main__':
    main()
