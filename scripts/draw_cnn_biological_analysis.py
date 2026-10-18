"""Draw a source-traceable exon interpretation composite for the CNN main model.

Uses completed screens and frozen exon annotations only. No inference or source
records are changed. Run select_cnn_biological_case.py first to select the case.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch, Patch
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from diffuprimer.io import sha256, write_json

SOURCE = ROOT / 'experiments/attention/results/exon_interpretation'
OUT = ROOT / 'paper/figures/cnn_biological_analysis'
INK, MUTED = '#193D50', '#59717C'
F, R, PURPLE = '#248A75', '#D36E51', '#9B88AC'
CAT = ['primer_spans_junction', 'amplicon_crosses_junction_only', 'same_exon', 'terminal_polyA', 'unmapped']
COL = ['#248A75', '#B9D7DC', '#E5D4B3', '#9B88AC', '#CCD4D8']
METHODS = ['cnn_diffusion', 'cnn_direct', 'primer3', 'qprimer']
LABELS = ['DiffuPrimer (CNN)', 'CNN direct', 'Primer3', 'qPrimer module']
EXON_COLORS = ['#C1DED8', '#C9DBEC', '#DCD4E7']


def label(fig, x, y, letter, title):
    fig.text(x, y, letter, fontsize=10, fontweight='bold', color=INK, va='baseline')
    fig.text(x + .025, y, title, fontsize=8.5, fontweight='bold', color=INK, va='baseline')


def primer_arrow(ax, start, end, y, forward=True, lw=2.5):
    a, b = (start, end) if forward else (end, start)
    ax.add_patch(FancyArrowPatch((a, y), (b, y), arrowstyle='-|>', mutation_scale=8,
                                linewidth=lw, color=F if forward else R, shrinkA=0, shrinkB=0))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads((SOURCE / 'summary.json').read_text())
    case = json.loads((OUT / 'illustrative_case.json').read_text())
    if summary['status'] != 'complete' or case['cnn_diffusion']['method'] != 'cnn_diffusion':
        raise ValueError('Completed CNN-specific evidence required')
    baseline, cnn = case['baseline'], case['cnn_diffusion']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.2,
                         'axes.labelsize': 7.2, 'xtick.labelsize': 7.2, 'ytick.labelsize': 7.2,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .55, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'svg.fonttype': 'path', 'savefig.facecolor': 'white', 'figure.facecolor': 'white'})
    fig = plt.figure(figsize=(178 / 25.4, 175 / 25.4))
    values = {(g['method'], g['category']): g for g in summary['groups']}
    label(fig, .035, .973, 'A', 'Exon geometry across the complete test cohort')
    fig.text(.975, .974, '75 genes · 20 proposals per gene', ha='right', color=MUTED, fontsize=7.2)
    # Small exon glyphs serve as the category legend, preserving a biological theme.
    glyph_labels = ['Primer spans a junction', 'Product crosses a junction', 'Same exon']
    for k, text in enumerate(glyph_labels):
        ax = fig.add_axes([.06 + k * .31, .921, .275, .035])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
        ax.add_patch(Rectangle((.01, .27), .105, .42, fc='#DDE9E7', ec='#829A9F', lw=.4))
        ax.add_patch(Rectangle((.117, .27), .105, .42, fc='#D5E1EC', ec='#829A9F', lw=.4))
        if k == 0:
            primer_arrow(ax, .075, .165, .8, True, 1.4)
        elif k == 1:
            primer_arrow(ax, .018, .072, .8, True, 1.4)
            primer_arrow(ax, .168, .217, .8, False, 1.4)
        else:
            primer_arrow(ax, .008, .05, .8, True, 1.4)
            primer_arrow(ax, .07, .112, .8, False, 1.4)
        ax.text(.26, .48, text, va='center', fontsize=7.2, color=INK)
        ax.plot([.26, .34], [.1, .1], color=COL[k], lw=3)
    fig.patches.append(Rectangle((.868, .933), .018, .009,
                                 transform=fig.transFigure, fc=PURPLE, ec='none'))
    fig.text(.894, .9375, 'Other', va='center', fontsize=7.2, color=INK)
    aa = [fig.add_axes([.23, .742, .335, .147]), fig.add_axes([.645, .742, .31, .147])]
    derived = []
    for ax, endpoint, title, xmax in zip(aa,
            ['mean_RNA_pass_pairs_per_gene', 'mean_retained_pairs_per_gene'],
            ['RNA-specific pairs', 'RNA-specific; no genomic product flag'], [20, 9]):
        for i, (method, name) in enumerate(zip(METHODS, LABELS)):
            total = 0
            for category, color in zip(CAT, COL):
                value = values[(method, category)][endpoint]
                ax.barh(i, value, left=total, height=.62, color=color, ec='white', lw=.6)
                derived.append({'panel': 'A', 'method': method, 'category': category,
                                'endpoint': endpoint, 'value': value})
                total += value
            ax.text(total + .13, i, f'{total:.2f}', va='center', color=INK, fontsize=7.2,
                    fontweight='bold' if method == 'cnn_diffusion' else 'normal')
        ax.set_xlim(0, xmax); ax.set_ylim(3.6, -.6)
        ax.set_yticks(range(4), LABELS if ax is aa[0] else [])
        ax.tick_params(axis='y', length=0, pad=7)
        for spine in ['left', 'top', 'right']: ax.spines[spine].set_visible(False)
        ax.set_title(title, fontsize=7.2, color=INK, loc='left', pad=7)
        ax.set_xlabel('Unique pairs per gene', color=MUTED)
        ax.grid(axis='x', lw=.45, color='#E7ECEE'); ax.set_axisbelow(True)
    fig.lines.append(plt.Line2D([.035, .965], [.652, .652], transform=fig.transFigure, color='#DDE5E7', lw=.6))

    label(fig, .035, .626, 'B', f"{case['gene_symbol']} · {case['transcript_accession_version']} · genomic locus")
    fig.text(.975, .628, f"{case['chromosome']} · {case['strand']} strand", ha='right', color=MUTED, fontsize=7.2)
    # The genomic panel is intentionally schematic; exact exons and gaps are labelled.
    lo, hi = case['RNA_plot_window']
    visible = [(i + 1, e) for i, e in enumerate(case['all_annotated_exons'])
               if e['RNA_start'] < hi and lo < e['RNA_end']]
    assert len(visible) == 3, 'This composed drawing expects a three-exon illustration'
    ag = fig.add_axes([.08, .497, .86, .105])
    ag.set_xlim(0, 1); ag.set_ylim(-.8, 1.8); ag.axis('off')
    starts = [.02, .365, .71]; width = .27
    for k, ((index, exon), start) in enumerate(zip(visible, starts)):
        ag.add_patch(Rectangle((start, .42), width, .42, fc=EXON_COLORS[k], ec='#718E99', lw=.7))
        ag.text(start + width / 2, 1.10, f'Exon {index}', ha='center', fontsize=7.2, color=INK)
        ag.text(start + width / 2, .12, f"{exon['genome_start']:,}\n{exon['genome_end']:,}",
                ha='center', va='top', fontsize=7.2, color=MUTED, linespacing=1.35)
        if k < 2:
            nxt = visible[k + 1][1]
            intron = nxt['genome_start'] - exon['genome_end'] if case['strand'] == '+' else exon['genome_start'] - nxt['genome_end']
            ag.plot([start + width, starts[k + 1]], [.63, .63], color='#9BAFB5', lw=1)
            mid = (start + width + starts[k + 1]) / 2
            ag.text(mid, .63, '//', ha='center', va='center', color=MUTED,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .0})
            ag.text(mid, 1.48, f'{intron:,} bp', ha='center', fontsize=7.2, color=MUTED)
    fig.text(.08, .478, 'Genomic layout is schematic; introns compressed. Endpoints: 0-based, half-open.',
             fontsize=7.2, color=MUTED)

    label(fig, .035, .449, 'C', 'Same transcript; similar RNA product lengths, different binding geometry')
    ar = fig.add_axes([.08, .292, .86, .13])
    ar.set_xlim(lo, hi); ar.set_ylim(-.18, 1.6); ar.set_yticks([])
    ar.spines[['left', 'right', 'top']].set_visible(False)
    ar.set_xlabel('RNA coordinate (bp; 0-based)', color=MUTED, labelpad=4)
    for k, (index, exon) in enumerate(visible):
        start, end = max(lo, exon['RNA_start']), min(hi, exon['RNA_end'])
        ar.add_patch(Rectangle((start, 1.23), end - start, .16, fc=EXON_COLORS[k], ec='#718E99', lw=.55))
        ar.text((start + end) / 2, 1.45, f'Exon {index}', ha='center', fontsize=7.2, color=INK)
    for junction in case['annotated_junctions']:
        if lo < junction < hi:
            ar.plot([junction, junction], [.19, 1.3], color=PURPLE, lw=.75, ls=(0, (2, 2)))
            ar.text(junction, .08, str(junction), ha='center', va='center', color=PURPLE, fontsize=7.2,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .15})
    for row, y, title in [(cnn, .91, 'DiffuPrimer'), (baseline, .36, 'Primer3')]:
        ar.plot([row['forward_start'], row['reverse_end']], [y, y], color='#AABBC2', lw=1.1)
        primer_arrow(ar, row['forward_start'], row['forward_end'], y, True)
        primer_arrow(ar, row['reverse_start'], row['reverse_end'], y, False)
        center = (row['forward_start'] + row['reverse_end']) / 2
        ar.text(center, y + .15, f"{title} · {row['amplicon_length']} bp", ha='center', fontsize=7.2,
                color=INK, fontweight='bold' if title == 'DiffuPrimer' else 'normal')
    ar.set_xticks([600, 700, 800, 900, 1000])

    label(fig, .035, .230, 'D', 'Actual primer sequences and screen outcomes')
    for x, row, method, color in [(.055, cnn, 'DiffuPrimer (CNN)', F), (.525, baseline, 'Primer3', '#A07B43')]:
        card = FancyBboxPatch((x, .038), .42, .166, boxstyle='round,pad=.009,rounding_size=.008',
                              transform=fig.transFigure, fc='#F4F9F7' if method.startswith('Diff') else '#FBF8F2',
                              ec='#D5E3DF' if method.startswith('Diff') else '#E8DDC9', lw=.6, zorder=-1)
        fig.patches.append(card)
        fig.text(x + .012, .183, f"{method} · proposal {row['proposal_index']}", color=color,
                 fontsize=7.5, fontweight='bold')
        fseq = row['forward_primer']; rseq = row['reverse_primer']
        if row['forward_junctions']:
            cut = row['forward_junctions'][0] - row['forward_start']; fseq = fseq[:cut] + '|' + fseq[cut:]
        if row['reverse_junctions']:
            cut = row['reverse_end'] - row['reverse_junctions'][0]; rseq = rseq[:cut] + '|' + rseq[cut:]
        fig.text(x + .012, .156, f"F 5′ {fseq} 3′", fontfamily='DejaVu Sans Mono', fontsize=7.2, color=F)
        fig.text(x + .012, .133, f"R 5′ {rseq} 3′", fontfamily='DejaVu Sans Mono', fontsize=7.2, color=R)
        fig.text(x + .012, .103, 'RNA specificity: pass', fontsize=7.2, color=INK)
        if row['genomic_product_flag']:
            witness = case['baseline_genomic_product_witness']
            fig.text(x + .012, .078, f"Genomic product: {witness['amplicon_length']} bp", fontsize=7.2, color=R)
            fig.text(x + .012, .053, f"[{witness['start']:,}, {witness['end']:,})", fontsize=7.2, color=MUTED)
        else:
            fig.text(x + .012, .078, 'Genomic product: none flagged', fontsize=7.2, color=F)
            side = 'F' if row['forward_junctions'] else 'R'
            junction = (row['forward_junctions'] + row['reverse_junctions'])[0]
            fig.text(x + .012, .053, f'{side} spans RNA junction {junction} (|)', fontsize=7.2, color=PURPLE)

    name = 'fig4_biological_specificity'
    for ext in ['pdf', 'svg']:
        fig.savefig(OUT / f'{name}.{ext}')
    plt.close(fig)
    subprocess.run(['pdftoppm', '-png', '-r', '350', '-singlefile',
                    str(OUT / f'{name}.pdf'), str(OUT / name)], check=True)
    with Image.open(OUT / f'{name}.png') as image:
        image.convert('RGB').save(OUT / f'{name}.png', dpi=(350, 350))
        image.convert('RGB').save(OUT / f'{name}.tiff', dpi=(350, 350), compression='tiff_lzw')
        image.thumbnail((700, 1000)); image.convert('RGB').save(OUT / f'{name}_review_700px.png')
    caption = (
        'Exon geometry and genomic-product screening of the CNN main model. (A) Unique RNA-specific candidates '
        'and their retained subset, partitioned by verified transcript exon geometry. Counts are means over all '
        '75 frozen test genes at K=20, including no-design genes; the geometry classes are mutually '
        'exclusive. Other consists of two terminal-poly(A) qPrimer candidates; there are no unmapped candidates in these plotted counts. '
        'Stacked counts display the number and exon composition together; the adjacent retained subset connects this geometry to the actual screen outcome. '
        '(B–D) Genomic and scaled RNA tracks locate exon boundaries and primer-binding intervals for an actual same-transcript illustration '
        f"for {case['gene_symbol']} ({case['transcript_accession_version']}). The top-ranked native Primer3 pair "
        f"has a {baseline['amplicon_length']}-bp RNA product and an exactly coordinate-matched genomic-product witness. "
        f"A CNN diffusion alternative yields a {cnn['amplicon_length']}-bp RNA product, with a forward primer spanning "
        'a splice junction, and has no genomic-product flag in the recorded screen. Primers are written 5′–3′; '
        '| marks the splice junction inside a primer; purple dashed lines and interior numbers mark RNA splice coordinates. RNA positions are drawn to scale, and the genomic layout '
        'is schematic with compressed introns. Candidate selection prioritizes Primer3 cases, display windows '
        'within 512 bp, minimal RNA product-length difference, then minimal combined window and stable IDs; '
        'this illustrative selection does not filter the cohort or alter aggregate outcomes. Screen: ≤2 ungapped '
        'mismatches, exact terminal five bases, potential products 36–1000 bp.'
    )
    alt = ('Four-panel biological interpretation figure. Stacked counts show that CNN diffusion produces 7.36 '
           'retained pairs per gene, versus 5.64 for CNN direct, 4.55 for Primer3 and 5.20 for qPrimer. A Susd5 '
           'case juxtaposes three annotated genomic exons with a spliced RNA coordinate track. The 103-bp CNN '
           'product has a forward primer across junction 588 and no genomic flag; the 104-bp same-exon Primer3 '
           'product has a 104-bp genomic witness. Actual primer sequences and genomic coordinates are shown.')
    write_json(OUT / 'provenance.json', {
        'figure': name, 'name': name, 'selected_case': {'gene_symbol': case['gene_symbol'], 'transcript_accession_version': case['transcript_accession_version'], 'CNN_candidate_id': cnn['candidate_id'], 'baseline_candidate_id': baseline['candidate_id'], 'eligible_case_pairs': case['eligible_case_pairs'], 'eligible_genes': case['eligible_genes']}, 'primary_method': 'cnn_diffusion', 'physical_width_mm': 178,
        'physical_height_mm': 175, 'minimum_font_pt': 7.2, 'raster_dpi': 350,
        'cohort_genes': 75, 'budget': 20, 'source_summary': str((SOURCE / 'summary.json').relative_to(ROOT)),
        'source_summary_sha256': sha256(SOURCE / 'summary.json'),
        'summary_sha256': sha256(ROOT / 'experiments/attention/results/study_summary.json'),
        'case_sha256': sha256(OUT / 'illustrative_case.json'), 'selection_rule': case['selection_rule'],
        'derived_values': derived, 'caption': caption, 'alt_text': alt,
        'files': [{'path': str((OUT / f'{name}.{ext}').relative_to(ROOT)), 'sha256': sha256(OUT / f'{name}.{ext}')}
                  for ext in ['pdf', 'svg', 'png', 'tiff']],
        'source_inputs': case['source_inputs'], 'geometry_source': summary['geometry_source'],
        'zero_count_categories': ['unmapped'],
        'plot_rationale': {'panel_A': 'Absolute stacked counts jointly answer how many candidates each geometry contributes and how many remain after the recorded screen, preserving the full gene denominator.',
                           'panels_B_C': 'Genomic and RNA tracks show actual exon boundaries and binding positions; a close RNA-product-length match isolates positional geometry in the illustration.',
                           'panel_D': 'Exact oligo sequences and screen outcomes make the coordinate interpretation traceable.'},
        'other_legend_categories': ['terminal_polyA', 'unmapped'],
        'interpretation_scope': 'Recorded computational screen; no physical assay measurements.'})
    record = json.loads((OUT / 'provenance.json').read_text())
    record['figures'] = [{'name': name, 'files': record['files'], 'caption': caption, 'alt_text': alt}]
    write_json(OUT / 'provenance.json', record)
    (OUT / 'caption.txt').write_text(caption + '\n\nAlt text: ' + alt + '\n')
    print(json.dumps({'figure': name, 'gene': case['gene_symbol'], 'RNA_products_bp':
                      [cnn['amplicon_length'], baseline['amplicon_length']], 'eligible_cases': case['eligible_case_pairs']}))


if __name__ == '__main__':
    main()
