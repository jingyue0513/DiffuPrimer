#!/usr/bin/env python3
"""Build a Chinese figure-led review of the unchanged completed CNN results."""
from __future__ import annotations

import base64
import hashlib
import html
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reports/cnn_publication'
PRIMARY = 'cnn_diffusion'
DIRECT = 'cnn_direct'
RETAINED = 'unique_RNA_pass_pairs_without_genomic_flag'
FRACTION = 'RNA_no_genome_fraction'
RNA = 'unique_RNA_pass_pairs'
DIVERSITY = 'RNA_pass_amplicon_distance'
TM = 'best_RNA_pass_temperature_penalty'
LABELS = {RETAINED: '保留候选 / 基因', FRACTION: 'RNA 特异候选中的保留比例', RNA: 'RNA 特异候选 / 基因',
          DIVERSITY: '扩增子位置距离', TM: '最佳温度偏差 / °C', 'anchor_preservation_rate': '返回候选的固定区间保持率',
          'length_compliant_feasible_rate': '长度合规且可行候选 / K'}
NAMES = {'cnn_diffusion': 'DiffuPrimer · CNN 扩散', 'cnn_direct': 'CNN · 直接预测', 'primer3': 'Primer3',
         'qprimer': 'qPrimer 设计模块', 'primercast': 'PrimerCast', 'attention_rope_diffusion': '纯 Attention · 扩散',
         'attention_hybrid_diffusion': 'CNN + Attention · 扩散'}
VARIANTS = {'full': '完整联合约束 / 16 步', 'raw': '逐位置输出', 'geometry_only': '区间结构约束',
            'single_oligo': '结构 + 单引物约束', 'steps1': '1 步', 'steps4': '4 步', 'steps8': '8 步'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def esc(value):
    return html.escape(str(value), quote=True)


def fmt(value, key=None):
    if value is None:
        return '—'
    if key in (FRACTION, 'anchor_preservation_rate', 'length_compliant_feasible_rate'):
        return f'{100 * value:.1f}%'
    return f'{value:.3f}' if key in (DIVERSITY, TM) else f'{value:.2f}'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    main_path = ROOT / 'experiments/attention/results/study_summary.json'
    extra_path = ROOT / 'experiments/attention_external/results/study_summary.json'
    documents = {'original': read(main_path), 'external': read(extra_path)}
    for doc in documents.values():
        if doc['status'] != 'complete':
            raise ValueError('The figure review requires completed source analyses')
    groups = {cohort: {(g['context'], g['method'], g['variant'], g['budget']): g for g in doc['summaries']}
              for cohort, doc in documents.items()}
    if documents['external']['evaluation_split'] != 'external_test' or documents['external']['old_sources_included']:
        raise ValueError('Additional independent genes must remain separate')
    numeric, resources, figures = [], [], []
    asset_checks = []

    def value(method, key, cohort='original', context='amp80-200', variant='full', budget=20):
        return groups[cohort][context, method, variant, budget][key]

    def cell(g, key, cohort):
        v = g[key]
        record = {'cohort': cohort, 'method': g['method'], 'context': g['context'], 'variant': g['variant'],
                  'budget': g['budget'], 'endpoint': key, 'value': v}
        numeric.append(record)
        attrs = ' '.join(f'data-{k}="{esc(json.dumps(v) if k == "value" else v)}"' for k, v in record.items())
        ci = g.get(key + '_ci95')
        note = f'95% CI {fmt(ci[0], key)} – {fmt(ci[1], key)}' if ci is not None else ''
        if key in (FRACTION, DIVERSITY, TM):
            note += ' · n=' + str(g.get(key + '_n', 0))
        return f'<td {attrs}><span>{fmt(v, key)}</span><small>{esc(note)}</small></td>'

    def table(methods, keys=(RETAINED, FRACTION), cohort='original', context='amp80-200',
              variants=('full',), budget=20, table_id=''):
        rows = []
        for variant in variants:
            for method in methods:
                g = groups[cohort].get((context, method, variant, budget))
                if g is None:
                    continue
                label = NAMES[method] if variants == ('full',) else VARIANTS[variant]
                rows.append(f'<tr class="{"ours" if method == PRIMARY else ""}"><th scope="row">{esc(label)}</th>'
                            + ''.join(cell(g, key, cohort) for key in keys) + '</tr>')
        return f'<div class="table-wrap table-scroll"><table class="result-table" id="{table_id}"><thead><tr><th>方法 / 设置</th>' + ''.join(
            f'<th>{esc(LABELS[key])}</th>' for key in keys) + '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table></div>'

    def copy_asset(path, name=None):
        path = Path(path)
        dest = OUT / (name or path.name)
        if not path.is_file():
            raise FileNotFoundError(path)
        shutil.copyfile(path, dest)
        resources.append(dest.name)
        asset_checks.append({'source': str(path.relative_to(ROOT)), 'source_sha256': sha(path),
                             'download': dest.name, 'download_sha256': sha(dest)})
        return dest.name

    def figure(directory, basename, caption, alt, number):
        directory = ROOT / directory
        for ext in ('png', 'pdf', 'svg', 'tiff'):
            copy_asset(directory / f'{basename}.{ext}')
        source = directory / f'{basename}.png'
        encoded = base64.b64encode(source.read_bytes()).decode()
        figures.append({'number': number, 'basename': basename, 'source': str(directory.relative_to(ROOT)),
                        'caption': caption, 'alt_text': alt})
        return f'<figure class="result-figure" id="figure-{number}"><a href="{basename}.png" target="_blank">' \
            f'<img src="data:image/png;base64,{encoded}" alt="{esc(alt)}" loading="lazy"></a>' \
            f'<figcaption><strong>Figure {number}.</strong> {esc(caption)} ' \
            f'<a href="{basename}.pdf" target="_blank">PDF</a> · <a href="{basename}.svg" download>SVG</a> · ' \
            f'<a href="{basename}.tiff" download>TIFF</a></figcaption></figure>'

    quant_dir = ROOT / 'paper/figures/cnn_publication_results'
    biology_dir = ROOT / 'paper/figures/cnn_biological_analysis'
    for directory in (quant_dir, biology_dir):
        provenance = read(directory / 'provenance.json')
        if provenance.get('primary_method') != PRIMARY or provenance.get('summary_sha256') != sha(main_path):
            raise ValueError('Figure source/primary model mismatch: ' + str(directory))
        if directory == quant_dir and provenance.get('external_summary_sha256') != sha(extra_path):
            raise ValueError('Additional-cohort figure source mismatch')
        copy_asset(directory / 'provenance.json', directory.name + '_provenance.json')
    method_dir = 'paper/figures/architecture_method'
    method_base = 'DiffuPrimer_method'
    method_figure = figure(method_dir, method_base,
        '整体 CNN 计算结构与展开的残差块；固定 DNA 上的干净、带噪结合位置标签及恢复过程；下方展示四种设计用途。',
        'Macro and expanded CNN computations, clean and noisy binding maps on fixed DNA, and four primer-design modes.', 1)
    copy_asset(ROOT / method_dir / (method_base + '.pptx'))
    copy_asset(ROOT / method_dir / 'provenance.json', 'method_provenance.json')
    quality_figure = figure('paper/figures/cnn_publication_results', 'fig2_candidate_quality',
        '原测试集与新增独立队列分别展示保留候选数量、RNA 特异候选中的保留比例和逐基因比较。',
        'Candidate quality and paired gene-level comparisons in two separate 75-gene cohorts.', 2)
    flexibility_figure = figure('paper/figures/cnn_publication_results', 'fig3_flexible_design',
        '实际扩增子长度的经验累积分布、候选预算、单侧补全与局部重设计，使用同一 CNN 扩散模型。',
        'Recorded length distributions, candidate-budget curves and fixed-primer design tasks.', 3)
    biology_figure = figure('paper/figures/cnn_biological_analysis', 'fig4_biological_specificity',
        '全队列外显子结构统计与真实 CNN 引物案例，将结合位置、剪接结构和完整参考筛查结果联系起来。',
        'Exon-junction geometry and a real CNN versus Primer3 example on the same gene.', 4)
    mechanism_figure = figure('paper/figures/cnn_publication_results', 'fig5_mechanism',
        '匹配的直接预测与扩散、逐级采样约束、去噪步数，以及 CNN 与 attention 架构对照。',
        'Diffusion, constraint, iteration-depth and backbone comparisons from recorded sources.', 5)

    author = read(ROOT / 'paper/author_metadata.json')
    author_line = ' · '.join(esc(a['name']) + ' · ' + esc(a['institution']) +
        (f' · <a href="mailto:{esc(a["email"])}">{esc(a["email"])}</a>' if a.get('email') else '')
        for a in author['authors'])
    methods = [PRIMARY, DIRECT, 'primer3', 'qprimer', 'primercast']
    original = value(PRIMARY, RETAINED)
    external = value(PRIMARY, RETAINED, 'external')
    p3 = value('primer3', RETAINED)
    ep3 = value('primer3', RETAINED, 'external')
    intro_results = f'在相同 K=20 下，两个测试队列平均每基因保留 <strong>{original:.2f} / {external:.2f} 对</strong>' \
        f' RNA 特异且无基因组产物标记的候选；Primer3 为 <strong>{p3:.2f} / {ep3:.2f} 对</strong>。'

    def article(number, anchor, title, purpose, contents):
        return f'<article class="experiment" id="{anchor}"><h3><span class="index">3.{number}</span>{title}</h3>' \
               f'<p class="question"><strong>研究问题</strong>{purpose}</p>{contents}</article>'

    quality = quality_figure + '<h4>原测试集 · 75 基因</h4>' + table(methods) + \
        '<h4>新增独立队列 · 75 基因</h4>' + table(methods, cohort='external') + \
        '<p>保留数量衡量可用备选供给，保留比例衡量 RNA 特异候选经过基因组筛查后的组成质量。' \
        '图中每个基因均保留在计数分析中，配对散点展示同一基因上的方法差异。</p>'
    budget = flexibility_figure + f'<p>CNN 扩散在 K=20 时的 RNA 特异候选扩增子位置距离为 ' \
        f'<strong>{value(PRIMARY, DIVERSITY):.3f}</strong>，相同架构直接预测为 ' \
        f'<strong>{value(DIRECT, DIVERSITY):.3f}</strong>。下方可切换预算，查看备选积累和候选位置变化。</p>' \
        + '<div class="budget-switch"><label>候选预算 K</label><div class="budget-buttons">' \
        + ''.join(f'<button data-budget="{k}" aria-pressed="{str(k == 20).lower()}">{k}</button>' for k in (1, 2, 5, 10, 20)) \
        + '</div></div>' + table(methods, keys=(RETAINED, DIVERSITY, TM), table_id='budget-results')
    length = '<p>指定长度条件后，权重保持不变。Figure 3 的实际长度分布与目标区间对齐，' \
        '下表同时报告各区间的筛后候选供给。</p>'
    for context in ('amp80-120', 'amp121-160', 'amp161-200'):
        length += f'<h4>{context[3:].replace("-", "–")} bp</h4>' + table(
            methods[:4], keys=(RETAINED, 'length_compliant_feasible_rate'), context=context)
    completion = '<p>将已有 F 或 R 作为位置锚点输入，模型在同一套采样流程中寻找兼容伙伴。' \
        'Figure 3 同时展示固定侧与生成侧的实际位置。</p>'
    redesign = '<p>保留一侧引物，改变扩增子长度要求，再生成满足新条件的另一侧引物。' \
        '重设计沿用既有权重与锚点表示。</p>'
    for side, label in (('forward', '保留 F，设计 R'), ('reverse', '保留 R，设计 F')):
        completion += f'<h4>{label}</h4>' + table(methods[:4], keys=(RETAINED, RNA, 'anchor_preservation_rate'),
                                                context='amp80-200_fix_' + side)
        redesign += f'<h4>{label}</h4>' + table(methods[:4], keys=(RETAINED, RNA, 'anchor_preservation_rate'),
                                              context='redesign_' + side)
    biology = biology_figure + '<p>外显子分区统计覆盖完整原测试队列；案例使用真实 CNN 候选与同基因 Primer3 候选。' \
        '转录本上的结合位置、剪接连接和基因组产物记录共同解释筛后候选的生物学差异。</p>'
    architecture = mechanism_figure + table([PRIMARY, DIRECT, 'attention_rope_diffusion', 'attention_hybrid_diffusion'],
                                            keys=(RETAINED, DIVERSITY, TM)) + \
        '<p>CNN 扩散在两个独立测试队列获得最多的保留候选，因此本次图表以它作为主模型；' \
        '纯 attention 与混合 attention 保留为架构对照。直接预测使用匹配的 CNN，分离迭代生成与架构的作用。</p>'
    guidance = '<p>依次考察逐位置输出、成对区间结构、单引物资格与联合温度约束；' \
        '再比较 1、4、8、16 步恢复的候选数量和多样性。Figure 5 汇总这些组件的作用。</p>' \
        + '<h4>采样约束</h4>' + table([PRIMARY], keys=(RETAINED, RNA), variants=('raw', 'geometry_only', 'single_oligo', 'full')) \
        + '<h4>去噪步数</h4>' + table([PRIMARY], keys=(RETAINED, DIVERSITY, TM), variants=('steps1', 'steps4', 'steps8', 'full'))
    experiments = '\n'.join([
        article(1, 'exp-quality', '候选质量与独立泛化', '共同预算与筛查条件下，能获得多少通过筛查的备选？候选组成质量如何？', quality),
        article(2, 'exp-budget', '多解生成与候选预算', '随着采样预算增加，保留候选与位置多样性如何变化？', budget),
        article(3, 'exp-length', '扩增子长度控制', '同一个模型能否响应不同长度区间？', length),
        article(4, 'exp-completion', '固定单侧引物的补全', '已有一条引物时，能否生成兼容伙伴并保持固定位置？', completion),
        article(5, 'exp-redesign', '条件改变后的局部重设计', '保留已有引物，能否适应新的扩增子要求？', redesign),
        article(6, 'exp-biology', '基因结构与检测特异性', '结合位置和剪接结构如何解释不同设计方案的筛查结果？', biology),
        article(7, 'exp-architecture', '模型结构与扩散机制', 'CNN、attention 和直接预测分别表现如何？', architecture),
        article(8, 'exp-mechanism', '联合约束与恢复深度', '兼容性采样与迭代恢复各自贡献什么？', guidance)])

    draft_links = ''
    draft_dir = ROOT / 'paper/bioinformatics_cnn_draft'
    if (draft_dir / 'verification.json').is_file():
        draft_check = read(draft_dir / 'verification.json')
        if draft_check.get('passed') and draft_check['manuscript_pdf_sha256'] == sha(draft_dir / 'manuscript.pdf') \
                and draft_check['supplementary_pdf_sha256'] == sha(draft_dir / 'supplementary.pdf'):
            copy_asset(draft_dir / 'manuscript.pdf', 'manuscript_draft.pdf')
            copy_asset(draft_dir / 'supplementary.pdf', 'supplementary_draft.pdf')
            draft_links = '<a href="manuscript_draft.pdf" target="_blank">英文研究初稿 · PDF</a>' \
                '<a href="supplementary_draft.pdf" target="_blank">补充材料 · PDF</a>'

    style = (OUT / 'base_style.css').read_text()
    style += '\n.result-figure img{display:block;width:100%;height:auto}.result-figure{margin:30px 0 42px}' \
        '.hero .author-line{margin:16px 0 26px;font-size:14px;color:#526778}' \
        '.hero .author-line a{color:inherit}.experiment h4{margin-top:26px}.question{margin-top:14px}' \
        '.table-wrap{margin:20px 0 30px}.result-figure figcaption{line-height:1.7}'
    navigation = ''.join(f'<a href="#{a}">{n:02d}　{t}</a>' for n, a, t in (
        (1, 'abstract', '摘要与贡献'), (2, 'method', '方法与数据'), (3, 'experiments', '实验'), (4, 'limitations', '局限')))
    budgets = {'summaries': [g for g in documents['original']['summaries'] if g['context'] == 'amp80-200'
                            and g['method'] in methods and g['variant'] == 'full']}
    document = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DiffuPrimer · 条件离散扩散引物对设计</title><style>{style}</style></head><body>
<header class="masthead"><a class="masthead-brand" href="#title">DiffuPrimer</a><span class="masthead-meta">CONDITIONAL DISCRETE DIFFUSION</span></header>
<div class="layout"><aside class="outline"><p class="outline-title">CONTENTS</p>{navigation}</aside><nav class="mobile-nav">{navigation}</nav><main class="article">
<div class="title-block hero" id="title"><p class="kicker">LEARNING PRIMER-PAIR BINDING POSITIONS</p><h1 class="wordmark">Diffu<em>Primer</em></h1>
<p class="paper-title">Conditional Discrete Diffusion for Diverse,<br>Constraint-Aware qPCR Primer-Pair Design</p>
<p class="chinese-title">基于条件离散扩散的多样化约束感知 qPCR 引物对设计</p><p class="author-line">{author_line}</p>
<div class="resource-links"><a href="{method_base}.pptx" download>英文方法图 · 可编辑 PPT</a><a href="study_summary.json" download>完整结果数据</a>
<a href="study_summary.contrasts.csv" download>基因配对比较</a>{draft_links}</div></div>
<section id="abstract"><h2><span class="section-no">01</span>摘要<span class="en">Abstract</span></h2>
<p class="abstract-copy">qPCR 引物设计需要同时协调序列特异性、热力学兼容性和扩增子长度，并支持同一靶标的多个备选方案及已有引物的复用。
DiffuPrimer 将正、反向引物表示为固定转录本上的成对结合区间，通过条件离散扩散恢复背景、F 与 R 标签。
条件 CNN 预测位置标签，联合区间采样将引物资格、成对温度匹配与固定锚点融入恢复过程，统一支持完整设计、长度控制、单侧补全和局部重设计。</p>
<p class="abstract-copy">{intro_results}</p><p class="contribution-heading">OUR CONTRIBUTIONS</p><ul class="contributions">
<li><strong>引物对位置分布建模。</strong>在固定生物序列上构建 0/F/R 标签的条件离散扩散，学习同一靶标的多种兼容结合位置方案。</li>
<li><strong>支持灵活使用的统一约束采样。</strong>将成对结构、单引物资格、联合温度约束与固定位置锚点融入同一恢复过程，支持四种设计模式。</li>
<li><strong>公开参考基准与生物结构评估。</strong>通过独立基因队列、完整 RNA／基因组筛查和外显子结构解释，系统评价候选质量、设计能力与作用机制。</li></ul></section>
<section id="method"><h2><span class="section-no">02</span>方法与数据<span class="en">Method &amp; data</span></h2>
<h3>2.1　条件 CNN 与成对位置扩散</h3><div id="method-figure">{method_figure}</div>
<p>输入为 DNA 碱基、带噪 0/F/R/M 标签、相对位置、设计条件及扩散时间。96 通道的四个残差卷积块使用
1、4、16、64 的膨胀率，逐位置输出背景／F／R 分数。训练对标签施加吸收态掩码，以干净位置标签监督恢复；推理时，
每一步从分数中采样兼容的 F/R 区间，再恢复部分掩码位置。固定引物和已恢复标签持续保留。</p>
<p>九个归一化条件编码引物及产物长度范围、Tm 范围、最大成对 Tm 差与 GC 范围。成对采样满足区间顺序、长度、
单引物资格、联合温度和位置锚点要求。最终从模板提取 F，并将 R 区间反向互补，随后评估异源二聚体及完整参考特异性。</p>
<h3>2.2　公开数据与独立评估</h3><p>使用小鼠 GRCm39／RefSeq 参考序列与注释。原队列包含 500 个基因，按相似簇划分为
350 个训练、75 个验证和 75 个测试基因。监督结合已缓存的 285 条映射公开引物位置记录与 7,410 条 Primer3 计算位置记录；
模型拟合使用训练基因标签，检查点按验证表现选择。</p><p>额外 75 个基因排除原队列，按固定随机顺序选择；所有 NM_ 异构体与训练基因异构体进行了
双向 MMseqs2 相似性筛查，使用 80% 一致性与双侧 80% 覆盖规则。两个测试队列分别分析，训练与生成 seed 固定为 20261013。
所有方法使用共同候选预算和参考评估规则。qPrimer 使用实际设计模块；PrimerCast 使用实际开源实现和模型权重。</p></section>
<section id="experiments"><h2><span class="section-no">03</span>实验<span class="en">Experiments</span></h2>{experiments}</section>
<section id="limitations"><h2><span class="section-no">04</span>局限<span class="en">Limitations</span></h2><p class="limitations">
当前基准围绕小鼠参考转录本及既定比对、热力学规则构建。后续可以扩展至更多物种与等位变异背景，并进一步提高灵活设计任务的候选覆盖。</p></section>
</main></div><script id="study-data" type="application/json">__STUDY_DATA__</script><script>
const D=JSON.parse(document.getElementById('study-data').textContent),NAME=__NAMES__,METHODS=__METHODS__,KEYS=__KEYS__;
function fmt(v,k){{if(v===null)return '—';return Number(v).toFixed(k==='unique_RNA_pass_pairs_without_genomic_flag'?2:3)}}
function updateBudget(k){{const b=document.querySelector('#budget-results tbody');b.replaceChildren();for(const method of METHODS){{
const g=D.summaries.find(g=>g.method===method&&g.budget===k);if(!g)continue;const tr=document.createElement('tr');if(method==='cnn_diffusion')tr.className='ours';
const th=document.createElement('th');th.scope='row';th.textContent=NAME[method];tr.append(th);for(const key of KEYS){{const td=document.createElement('td');
Object.assign(td.dataset,{{cohort:'original',method,context:g.context,variant:g.variant,budget:String(k),endpoint:key,value:JSON.stringify(g[key])}});
const span=document.createElement('span');span.textContent=fmt(g[key],key);const small=document.createElement('small');const ci=g[key+'_ci95'];
small.textContent=ci?'95% CI '+fmt(ci[0],key)+' – '+fmt(ci[1],key):'';if(key!=='unique_RNA_pass_pairs_without_genomic_flag')small.textContent+=' · n='+g[key+'_n'];
td.append(span,small);tr.append(td)}}b.append(tr)}}document.querySelectorAll('[data-budget].budget-button,.budget-buttons button').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.budget)===k)))}}
document.querySelectorAll('.budget-buttons button').forEach(b=>b.addEventListener('click',()=>updateBudget(Number(b.dataset.budget))));window.researchReportReady=true;
</script></body></html>'''
    document = document.replace('__STUDY_DATA__', json.dumps(budgets, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')) \
        .replace('__NAMES__', json.dumps(NAMES, ensure_ascii=False)).replace('__METHODS__', json.dumps(methods)) \
        .replace('__KEYS__', json.dumps([RETAINED, DIVERSITY, TM]))
    font_path = Path('/tmp/dfm-html-check/fonts/NotoSansCJKsc-Regular.otf')
    if not font_path.is_file():
        font_path = Path('/root/.local/share/fonts/NotoSansCJKsc-Regular.otf')
    from fontTools import subset
    from fontTools.ttLib import TTFont
    font = TTFont(font_path)
    sub = subset.Subsetter(options=subset.Options())
    sub.populate(text=re.sub(r'data:image/[^;]+;base64,[A-Za-z0-9+/=]+', '', document))
    sub.subset(font)
    font_output = OUT / 'report_chinese_subset.otf'
    font.save(font_output)
    font_data = 'data:font/otf;base64,' + base64.b64encode(font_output.read_bytes()).decode()
    document = document.replace('__CHINESE_FONT_DATA__', font_data)
    (OUT / 'index.html').write_text(document)
    for source, dest in ((main_path, 'study_summary.json'), (extra_path, 'external_study_summary.json'),
        (main_path.with_suffix('.groups.csv'), 'study_summary.groups.csv'),
        (main_path.with_suffix('.contrasts.csv'), 'study_summary.contrasts.csv'),
        (ROOT / 'reports/assets/FONT_LICENSE.txt', 'FONT_LICENSE.txt')):
        copy_asset(source, dest)
    manifest = {'built_at': datetime.now(timezone.utc).isoformat(), 'status': 'complete', 'primary_method': PRIMARY,
        'primary_selection': 'Post-evaluation figure presentation; original experiment protocols and outcomes unchanged',
        'summary_sha256': sha(main_path), 'external_summary_sha256': sha(extra_path), 'page_sha256': sha(OUT / 'index.html'),
        'builder_sha256': sha(__file__), 'sections': 4, 'contributions': 3, 'experimental_questions': 8,
        'figures': figures, 'numeric_instances': numeric, 'resources': sorted(set(resources)), 'asset_checks': asset_checks,
        'cohorts_pooled': False, 'author_metadata_sha256': sha(ROOT / 'paper/author_metadata.json'),
        'manuscript_completion_asserted': False}
    (OUT / 'build_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({'output': str(OUT), 'page_bytes': (OUT / 'index.html').stat().st_size,
                      'primary_method': PRIMARY, 'figures': len(figures), 'numeric_instances': len(numeric)}))


if __name__ == '__main__':
    main()
