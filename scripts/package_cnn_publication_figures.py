#!/usr/bin/env python3
"""Package completed publication graphics and the offline research review."""
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    page = ROOT / 'reports/cnn_publication'
    check = json.loads((page / 'browser_check.json').read_text())
    if not check['passed'] or check['page_sha256'] != sha(page / 'index.html'):
        raise ValueError('Current figure review page must pass browser verification')
    manifest = json.loads((page / 'build_manifest.json').read_text())
    for asset in manifest['asset_checks']:
        if sha(ROOT / asset['source']) != asset['source_sha256'] or sha(page / asset['download']) != asset['download_sha256']:
            raise ValueError('Review asset differs from its recorded source: ' + asset['source'])
    files = set()
    for directory in ('paper/figures/architecture_method', 'paper/figures/cnn_publication_results',
                      'paper/figures/cnn_biological_analysis', 'reports/cnn_publication'):
        for p in (ROOT / directory).rglob('*'):
            if p.is_file() and p.suffix in {'.pptx', '.pdf', '.svg', '.png', '.tiff', '.json', '.csv', '.md', '.html', '.txt', '.css', '.otf'}:
                if p.name in {'first_preview.png', 'refined_preview.png'}:
                    continue
                if directory.startswith('reports/') and p.name.startswith('preview_'):
                    continue
                files.add(p)
    for name in ('scripts/draw_diffuprimer_architecture_method.py', 'scripts/export_diffuprimer_architecture_method.py',
                 'scripts/select_cnn_biological_case.py', 'scripts/draw_cnn_biological_analysis.py',
                 'experiments/attention/plot_publication_results.py', 'scripts/build_cnn_publication_page.py',
                 'scripts/check_cnn_publication_page.cjs', 'scripts/package_cnn_publication_figures.py',
                 'docs/cnn_publication_figures.md', 'paper/author_metadata.json'):
        files.add(ROOT / name)
    archive = ROOT / 'paper/DiffuPrimer_publication_figures.zip'
    records = [{'path': str(p.relative_to(ROOT)), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(files)]
    record = {'at': datetime.now(timezone.utc).isoformat(), 'primary_method': 'cnn_diffusion',
              'method_drawing_agent': 'gpt-6-astra (explicit user request)', 'integrated_figures': 5,
              'main_summary_sha256': manifest['summary_sha256'], 'external_summary_sha256': manifest['external_summary_sha256'],
              'page_sha256': manifest['page_sha256'], 'files': records, 'formal_manuscript_completion_asserted': False}
    readme = ('DiffuPrimer publication figures and offline review\n\n'
              'Open reports/cnn_publication/index.html for the Chinese review.\n'
              'The editable English method figure is paper/figures/architecture_method/DiffuPrimer_method.pptx.\n'
              'All five figures have PDF/SVG/350 dpi PNG/TIFF and source/provenance records.\n'
              'This supplemental bundle contains graphics, derivation scripts and an offline page.\n'
              'It does not include model weights or all native candidate/reference files; full regeneration uses the complete project.\n'
              'See docs/cnn_publication_figures.md and the earlier full computational research archive.\n'
              'Original experimental outcomes, protocols and checkpoints are unchanged.\n')
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(files):
            z.write(p, 'DiffuPrimer_figures/' + str(p.relative_to(ROOT)))
        z.writestr('DiffuPrimer_figures/FIGURE_BUNDLE_README.txt', readme)
        z.writestr('DiffuPrimer_figures/figure_bundle_manifest.json', json.dumps(record, indent=2) + '\n')
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ValueError('Figure archive CRC failure')
        for r in records:
            payload = z.read('DiffuPrimer_figures/' + r['path'])
            if len(payload) != r['bytes'] or hashlib.sha256(payload).hexdigest() != r['sha256']:
                raise ValueError('Figure archive entry differs from source: ' + r['path'])
    receipt = {'archive': str(archive.relative_to(ROOT)), 'bytes': archive.stat().st_size, 'sha256': sha(archive),
               'source_files_checked': len(records), 'crc_passed': True, 'all_entry_hashes_passed': True,
               'figures': 5, 'primary_method': 'cnn_diffusion', 'page_sha256': manifest['page_sha256']}
    (ROOT / 'paper/cnn_publication_figure_bundle_check.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
