#!/usr/bin/env python3
"""Build the Chinese paper website from final, recorded scientific results."""
import base64
import json
import shutil
from pathlib import Path
from fontTools import subset
from fontTools.ttLib import TTFont
from diffuprimer.io import sha256, utc_now, write_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reports/chinese_research'
OUT.mkdir(parents=True, exist_ok=True)
summary_path = ROOT / 'results/final_fixed_seed/study_summary.json'
summary = json.loads(summary_path.read_text())
provenance = json.loads((ROOT / 'results/final_fixed_seed/figures/provenance.json').read_text())
assert provenance['summary_sha256'] == sha256(summary_path)
values = json.loads((ROOT / 'paper/reported_values.json').read_text())
assert values['summary_sha256'] == sha256(summary_path)
check = json.loads((ROOT / 'results/final_fixed_seed/completion_check.json').read_text())
assert check['passed']
data = {
    'seed': summary['seed'], 'groups': summary['summaries'],
    'contrasts': summary['paired_comparisons'], 'estimates': provenance['values'],
    'values': values['values'],
    'exon': json.loads((ROOT / 'results/final_fixed_seed/exon_partition.json').read_text())['rows'],
}
template = ROOT / 'reports/templates/chinese_research.html'
html = template.read_text().replace('__STUDY_DATA__', json.dumps(data,ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('<','\\u003c'))
method_dir = ROOT / 'paper/figures/method_ppt'
html = html.replace('__METHOD_PNG__', base64.b64encode((method_dir / 'DiffuPrimer_method.png').read_bytes()).decode())
for extension in ['pptx','pdf','svg','png']:
    shutil.copyfile(method_dir / f'DiffuPrimer_method.{extension}', OUT / f'DiffuPrimer_method.{extension}')
font_path = ROOT / 'reports/assets/report_chinese_subset.otf'
original_font = Path('/tmp/dfm-html-check/fonts/NotoSansCJKsc-Regular.otf')
if original_font.exists():
    font = TTFont(original_font)
    options = subset.Options()
    sub = subset.Subsetter(options=options)
    sub.populate(text=''.join(sorted(set(html))))
    sub.subset(font)
    font.save(font_path)
font_data = base64.b64encode(font_path.read_bytes()).decode()
html = html.replace('__CHINESE_FONT__', font_data)
(OUT / 'index.html').write_text(html)
for name in ['manuscript.pdf','supplementary.pdf']:
    shutil.copyfile(ROOT / 'paper' / name, OUT / name)
for name in ['study_summary.groups.csv','study_summary.contrasts.csv','study_summary.json']:
    shutil.copyfile(ROOT / 'results/final_fixed_seed' / name, OUT / name)
shutil.copyfile(ROOT / 'reports/assets/FONT_LICENSE.txt',OUT / 'FONT_LICENSE.txt')
write_json(OUT / 'build_manifest.json', {
    'built_at': utc_now(), 'presentation': 'paper_website',
    'page_sha256': sha256(OUT / 'index.html'), 'builder_sha256': sha256(Path(__file__)),
    'template_sha256': sha256(template), 'summary_sha256': sha256(summary_path),
    'reference_design': 'https://tactile-reactive-dexterous.github.io/',
    'sections': ['Abstract and contributions','Method and public data','Experiments','Limitations'],
    'contribution_count': 3, 'experiment_sections': 8,
    'method_figure': 'English editable PowerPoint with PDF, SVG and PNG exports', 'video_count': 0,
    'manuscript_complete': True, 'external_runtime_dependencies': [],
})
print(OUT / 'index.html')
