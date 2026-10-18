#!/usr/bin/env python3
"""Build a separate offline technical-review website from real attention results.

This script does not change a journal manuscript or the completed CNN website.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import html
import json
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json

METHODS = ["primer3", "qprimer", "primercast", "uniform_pair_sampler", "rule_sampler", "cnn_direct",
           "cnn_diffusion", "attention_rope_direct", "attention_rope_diffusion", "attention_direct", "attention_diffusion",
           "attention_hybrid_direct", "attention_hybrid_diffusion"]
NAMES = {"primer3": "Primer3", "qprimer": "qPrimer", "primercast": "PrimerCast", "uniform_pair_sampler": "均匀兼容对",
         "rule_sampler": "正向优先", "cnn_direct": "CNN · 直接预测", "cnn_diffusion": "CNN · 扩散",
         "attention_direct": "Attention · 直接预测", "attention_diffusion": "DiffuPrimer · Attention"}
NAMES.update({"attention_rope_direct":"RoPE Attention · 直接预测", "attention_rope_diffusion":"RoPE Attention · 扩散",
              "attention_hybrid_direct":"CNN + Attention · 直接预测", "attention_hybrid_diffusion":"DiffuPrimer · CNN + Attention"})
RETAINED = "unique_RNA_pass_pairs_without_genomic_flag"
RNA = "unique_RNA_pass_pairs"
FRACTION = "RNA_no_genome_fraction"
DIVERSITY = "RNA_pass_amplicon_distance"
TM = "best_RNA_pass_temperature_penalty"
LABELS = {RETAINED: "保留候选 / 基因", RNA: "RNA 特异候选 / 基因", FRACTION: "保留比例", DIVERSITY: "扩增子 Jaccard 距离", TM: "最佳温度偏差 / °C",
          "geometry_valid_rate": "结构有效率", "anchor_preservation_rate": "固定区间保持率", "length_compliant_feasible_rate": "长度合规候选 / K"}


def esc(value):
    return html.escape(str(value), quote=True)


def fmt(value, key=None):
    if value is None:
        return "—"
    if key in (FRACTION, "geometry_valid_rate", "anchor_preservation_rate", "length_compliant_feasible_rate"):
        return f"{value * 100:.1f}%"
    return f"{value:.3f}" if key in (DIVERSITY, TM) else f"{value:.2f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-partial", action="store_true", help="Build an explicitly labeled technical working preview")
    parser.add_argument("--summary", default="experiments/attention/results/study_summary.json")
    parser.add_argument("--protocol", default="experiments/attention/protocol.json", help="Declared experiment/primary architecture protocol")
    parser.add_argument("--output", default="reports/attention_research")
    parser.add_argument("--primary-method", help="Primary diffusion label (default: protocol/summary declaration)")
    parser.add_argument("--direct-method", help="Matched primary direct predictor")
    parser.add_argument("--method-figure-directory", help="Directory containing editable English method assets")
    parser.add_argument("--method-figure-basename", help="Override the protocol method_figure_basename")
    parser.add_argument("--external-summary", help="Optional separately analyzed new-gene cohort; displayed inside experimental question 3.1 without pooling")
    args = parser.parse_args()
    summary_path = ROOT / args.summary
    summary_bytes = summary_path.read_bytes()
    summary = json.loads(summary_bytes)
    summary_sha = hashlib.sha256(summary_bytes).hexdigest()
    partial = summary["status"] != "complete"
    if partial and not args.allow_partial:
        parser.error("Complete analysis is required; --allow-partial is available for a technical draft preview")
    out = ROOT / args.output
    if out.resolve() == (ROOT / "reports/chinese_research").resolve():
        parser.error("This builder must preserve the original completed CNN research website")
    out.mkdir(parents=True, exist_ok=True)
    groups = {(g["context"], g["method"], g["variant"], g["budget"]): g for g in summary["summaries"]}
    contrasts = summary["paired_comparisons"]
    protocol_path = ROOT / args.protocol
    protocol = json.loads(protocol_path.read_text())
    if not args.allow_partial and summary.get("protocol_sha256") != sha256(protocol_path):
        raise ValueError("Summary was analyzed under a different protocol; analyze the declared final protocol before final review build")
    observed = {g["method"] for g in summary["summaries"]}
    primary = args.primary_method or protocol.get("primary_methods", {}).get("diffusion") or summary.get("primary_methods", {}).get("diffusion") or next(
        (m for m in ("attention_hybrid_diffusion", "attention_rope_diffusion", "attention_diffusion") if m in observed), "attention_diffusion")
    direct = args.direct_method or protocol.get("primary_methods", {}).get("direct") or summary.get("primary_methods", {}).get("direct") or primary.removesuffix("diffusion") + "direct"
    method_order = METHODS + sorted(observed - set(METHODS))
    for m in observed:
        NAMES.setdefault(m,m.replace("_"," "))
    data_path = ROOT / "data/processed/dataset_partial_pilot.manifest.json"
    dataset = json.loads(data_path.read_text())
    targets = read_jsonl(ROOT / "data/processed/targets_pilot.jsonl")
    split_counts = {s: sum(t["split"] == s for t in targets) for s in ("train", "validation", "test")}
    method_dir = ROOT / (args.method_figure_directory or protocol.get("method_figure_directory", "paper/figures/attention_style_draft"))
    method_basename = args.method_figure_basename or protocol.get("method_figure_basename", "DiffuPrimer_attention_method")
    method_png = method_dir / (method_basename + ".png")
    figure_dir = ROOT / "experiments/attention/results/figures"
    provenance_path = figure_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text()) if provenance_path.is_file() else {}
    current_figures = provenance.get("summary_sha256") == summary_sha
    if not args.allow_partial and not current_figures:
        raise ValueError("Run plot_results.py against the complete current summary before building the final review website")
    numeric_instances, used_figures, resource_files = [], [], []
    author_path = ROOT / "paper/author_metadata.json"
    author_metadata = json.loads(author_path.read_text())
    author_line = ' · '.join(esc(a['name']) + ' · ' + esc(a['institution'])
        + (f' · <a href="mailto:{esc(a["email"])}">{esc(a["email"])}</a>' if a.get('corresponding') and a.get('email') else '')
        for a in author_metadata['authors'])
    external, external_path, external_sha = None, None, None
    if args.external_summary:
        external_path = ROOT / args.external_summary
        external_bytes = external_path.read_bytes()
        external = json.loads(external_bytes)
        external_sha = hashlib.sha256(external_bytes).hexdigest()
        if external.get("evaluation_split") != "external_test" or external.get("old_sources_included"):
            raise ValueError("External block requires a separate external_test analysis without original sources")
        if not args.allow_partial and external["status"] != "complete":
            raise ValueError("The requested external cohort analysis has not completed")
        (out / "external_study_summary.json").write_bytes(external_bytes)
        resource_files.append("external_study_summary.json")

    def cell(g, key, cohort="primary"):
        value = g[key]
        numeric_instances.append({"cohort":cohort,"method": g["method"], "context": g["context"], "variant": g["variant"], "budget": g["budget"], "endpoint": key, "value": value})
        ci = g.get(key + "_ci95")
        annotation = "" if ci is None else f"95% CI {fmt(ci[0], key)}–{fmt(ci[1], key)}"
        n = g.get(key + "_n")
        if key in (FRACTION, DIVERSITY, TM) and n is not None:
            annotation += f" · n={n}"
        return f'<td data-cohort="{cohort}" data-method="{g["method"]}" data-context="{g["context"]}" data-variant="{g["variant"]}" data-budget="{g["budget"]}" data-endpoint="{key}" data-value="{esc(json.dumps(value))}"><span>{fmt(value, key)}</span><small>{esc(annotation)}</small></td>'

    def result_table(context="amp80-200", keys=(RETAINED, RNA, FRACTION), methods=None, variants=("full",), table_id="", group_lookup=None, cohort="primary"):
        rows = []
        methods = method_order if methods is None else methods
        lookup = groups if group_lookup is None else group_lookup
        for variant in variants:
            for method in methods:
                g = lookup.get((context, method, variant, 20))
                if g is None:
                    continue
                label = NAMES.get(method, method)
                if len(variants) > 1:
                    label += " · " + {"full": "完整联合引导", "single_oligo": "单条引物引导", "steps1": "1 步", "steps4": "4 步", "steps8": "8 步"}.get(variant, variant)
                rows.append(f'<tr class="{"ours" if method == primary else ""}"><th scope="row">{esc(label)}</th>{"".join(cell(g,k,cohort) for k in keys)}</tr>')
        if not rows:
            return '<p class="pending">该实验将在生成与完整参考评估完成后显示真实结果。</p>'
        return f'<div class="table-scroll"><table class="result-table" id="{table_id}"><thead><tr><th>方法</th>{"".join("<th>"+LABELS[k]+"</th>" for k in keys)}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'

    def contrast_text(right, key=RETAINED, context="amp80-200", left=None, variant="full"):
        left = primary if left is None else left
        found = next((c for c in contrasts if c["context"] == context and c["budget"] == 20 and c["left_method"] == left and c["right_method"] == right and c["right_variant"] == variant and c["endpoint"] == key), None)
        if found is None:
            return ""
        delta, ci = found["mean_difference"], found["ci95"]
        factor = 100 if key == FRACTION else 1
        unit = " 个百分点" if key == FRACTION else ""
        numeric_instances.append({"type": "contrast", "left": left, "right": right, "context": context, "endpoint": key, "difference": delta, "ci95": ci, "n": found["n"]})
        return f'<p class="finding">{esc(NAMES[left])} 相对 {esc(NAMES.get(right,right))} 的<strong>{esc(LABELS[key])}</strong>基因配对差为 <strong>{delta * factor:+.3f}{unit}</strong>（95% CI {ci[0] * factor:.3f}–{ci[1] * factor:.3f}，n={found["n"]}）。</p>'

    def figure(name, caption):
        record = next((g for g in provenance.get("generated", []) if g["name"] == name), None)
        png = figure_dir / (name + ".png")
        if not current_figures or record is None or not png.is_file():
            return ""
        for ext in ("png", "pdf", "svg"):
            source = figure_dir / (name + "." + ext)
            expected = next((f["sha256"] for f in record["files"] if f["path"].endswith(source.name)), None)
            if expected != sha256(source):
                raise ValueError("Figure changed after recorded provenance: " + source.name)
            shutil.copyfile(source, out / source.name)
            resource_files.append(source.name)
        used_figures.append(name)
        encoded = base64.b64encode(png.read_bytes()).decode()
        return f'<figure class="result-figure"><a href="{name}.png" target="_blank"><img src="data:image/png;base64,{encoded}" alt="{esc(record["alt_text"])}"></a><figcaption>{esc(caption)} <a href="{name}.pdf" target="_blank">PDF</a> · <a href="{name}.svg" download>SVG</a></figcaption></figure>'

    def pending_attention(context="amp80-200", variant="full"):
        if (context, primary, variant, 20) in groups:
            return ""
        return '<p class="pending">Attention 版本结果待完整参考评估；下表如有数值，均为已完成来源，CNN 参照单独标注。</p>'

    def experiment(index, anchor, title, purpose, contents):
        return f'<article class="experiment" id="{anchor}"><h3><span class="index">3.{index}</span>{title}</h3><p class="question"><strong>目的</strong>{purpose}</p>{contents}</article>'

    abstract_result = ""
    attn = groups.get(("amp80-200", primary, "full", 20))
    p3 = groups.get(("amp80-200", "primer3", "full", 20))
    if attn and p3:
        abstract_result = f'<p class="abstract-copy">在相同候选预算与完整参考评估下，Attention 版本每个测试基因获得 <strong>{fmt(attn[RETAINED])} 对</strong> RNA 特异且无预测基因组产物的候选，Primer3 为 <strong>{fmt(p3[RETAINED])} 对</strong>。正文同时展示架构、长度条件、单侧补全、局部重设计、联合引导与去噪深度的实际结果。</p>'
        numeric_instances.extend([{"type":"abstract","method":primary,"endpoint":RETAINED,"value":attn[RETAINED]}, {"type":"abstract","method":"primer3","endpoint":RETAINED,"value":p3[RETAINED]}])
    else:
        abstract_result = '<p class="abstract-copy">我们使用公开参考与来源引物构建同源感知的基因级基准，并从候选质量、神经架构、设计条件和灵活使用方式展开系统比较。实验数值在相应来源完成全参考评估后进入本页。</p>'
    architecture = protocol.get("primary_architecture",protocol["architecture"])
    attention_layers = architecture.get("attention_layers",architecture["layers"])
    local_encoder = architecture.get("local_encoder", architecture.get("cnn_encoder"))
    local_encoder_note = f"前端采用 {architecture['layers']} 层一维残差 CNN 提取局部碱基特征，并保留相对转录本位置，再将这些特征输入 Attention；" if local_encoder else ""
    train = protocol["training"]
    method_encoded = base64.b64encode(method_png.read_bytes()).decode()
    for ext in ("png", "pdf", "svg", "pptx"):
        name = "DiffuPrimer_attention_method." + ext
        shutil.copyfile(method_dir / (method_basename + "." + ext), out / name)
        resource_files.append(name)
    quality = pending_attention() + figure("01_candidate_quality", "图 2 · 相同 K=20 下的候选质量、RNA 特异候选供给与保留比例。点为基因均值，线为 95% bootstrap 区间。") + result_table() + contrast_text("primer3") + contrast_text("qprimer") + contrast_text("primercast")
    exon_provenance_path = ROOT / "paper/figures/exon_interpretation/provenance.json"
    exon_case_sha = None
    if exon_provenance_path.is_file():
        exon_provenance = json.loads(exon_provenance_path.read_text())
        case_path = ROOT / exon_provenance["case_source"]
        source_path = ROOT / exon_provenance["source_summary"]
        if sha256(case_path) != exon_provenance["case_sha256"] or sha256(source_path) != exon_provenance["source_summary_sha256"]:
            raise ValueError("Exon illustration sources changed after provenance")
        record = next(r for r in exon_provenance["figures"] if r["name"] == "same_gene_exon_case")
        for asset in record["files"]:
            source = ROOT / asset["path"]
            if sha256(source) != asset["sha256"]:
                raise ValueError("Exon illustration asset changed after provenance")
            shutil.copyfile(source, out / source.name)
            resource_files.append(source.name)
        exon_case_sha = sha256(case_path)
        shutil.copyfile(case_path, out / "illustrative_case.json")
        resource_files.append("illustrative_case.json")
        encoded = base64.b64encode((out / "same_gene_exon_case.png").read_bytes()).decode()
        quality += '<h4>基因检测实例：Ndp 的外显子与引物位置</h4><p>在 Ndp（NM_010883.3）中，Primer3 的首位引物对位于同一外显子，完整参考搜索记录了对应基因组产物；DiffuPrimer 的一组保留候选跨越转录本剪接连接。下图展示实际序列坐标与搜索记录，帮助解释多位置候选如何提供不同的检测方案。</p>'
        quality += f'<figure class="result-figure"><a href="same_gene_exon_case.png" target="_blank"><img src="data:image/png;base64,{encoded}" alt="{esc(record["alt_text"])}"></a><figcaption>Ndp 实例 · 基因组示意压缩了内含子，标注坐标来自实际参考。该实例按记录的显示宽度规则选择，整体比较使用全部测试基因。<a href="same_gene_exon_case.pdf" target="_blank">PDF</a> · <a href="same_gene_exon_case.svg" download>SVG</a> · <a href="illustrative_case.json" download>序列与搜索记录</a></figcaption></figure>'
        used_figures.append("same_gene_exon_case")
    external_data_note = ""
    if external:
        ext_groups = {(g["context"],g["method"],g["variant"],g["budget"]):g for g in external["summaries"]}
        quality += f'<h4>额外独立基因队列：N={external["total_test_genes"]}</h4><p>在同一公开 RefSeq 来源中预先冻结额外基因队列，独立计算候选质量与基因配对差；该队列与原测试结果分开报告。</p>'
        quality += result_table(group_lookup=ext_groups,cohort="external")
        quality += '<p class="subtle"><a href="external_study_summary.json" download>下载独立队列的完整统计来源</a></p>'
        ext_figures = external_path.parent / "figures"
        ext_provenance_path = ext_figures / "provenance.json"
        if ext_provenance_path.is_file():
            ext_provenance = json.loads(ext_provenance_path.read_text())
            ext_record = next((r for r in ext_provenance.get("generated",[]) if r["name"] == "05_external_generalization"),None)
            if ext_record and ext_provenance.get("summary_sha256") == external_sha:
                for extension in ("png","pdf","svg"):
                    source = ext_figures / ("05_external_generalization." + extension)
                    expected = next(f["sha256"] for f in ext_record["files"] if f["path"].endswith(source.name))
                    if sha256(source) != expected:
                        raise ValueError("External figure changed after provenance")
                    shutil.copyfile(source,out / source.name);resource_files.append(source.name)
                encoded = base64.b64encode((ext_figures / "05_external_generalization.png").read_bytes()).decode()
                quality += f'<figure class="result-figure"><a href="05_external_generalization.png" target="_blank"><img src="data:image/png;base64,{encoded}" alt="{esc(ext_record["alt_text"])}"></a><figcaption>额外基因泛化评估 · 两个队列分别计算均值与 95% CI，不跨队列配对或合并基因。<a href="05_external_generalization.pdf" target="_blank">PDF</a> · <a href="05_external_generalization.svg" download>SVG</a></figcaption></figure>'
                used_figures.append("05_external_generalization")
        external_data_note = f'<p class="subtle">另设 {external["total_test_genes"]} 个额外基因的独立评估队列，结果在 3.1 内分开报告；目标清单与哈希记录随独立统计来源提供。</p>'
    neural = [m for m in method_order if m.endswith(("_direct","_diffusion"))]
    architecture_result = pending_attention() + figure("02_architecture_comparison", "图 3 · CNN、纯 RoPE Attention 以及 CNN + Attention 的架构及生成机制比较；下排使用同一基因的配对差，只有完成来源进入图中。") + result_table(keys=(RETAINED, DIVERSITY, TM), methods=neural) + contrast_text("cnn_diffusion") + contrast_text(direct, DIVERSITY) + contrast_text(direct, TM)
    if primary == "attention_hybrid_diffusion":
        architecture_result += contrast_text("attention_rope_diffusion")
    budget = pending_attention() + figure("04_budget_and_diversity", "图 4 · 同一批候选的前 K 个方案。多样性和温度指标使用相应有效基因，分母在表中列明。") + '<div class="budget-switch"><label>查看候选预算 K</label><div class="budget-buttons" role="group" aria-label="候选预算">' + ''.join(f'<button data-budget="{k}" aria-pressed="{str(k==20).lower()}">{k}</button>' for k in (1,2,5,10,20)) + '</div></div>' + result_table(keys=(RETAINED,DIVERSITY,TM), table_id="budget-results")
    length = figure("03_flexible_design", "图 5 · 同一模型的长度控制、固定一侧补全与改变长度要求后的局部重设计。")
    for context in ("amp80-120", "amp121-160", "amp161-200"):
        length += '<h4>' + context.replace("amp", "").replace("-", "–") + ' bp</h4>' + pending_attention(context) + result_table(context, keys=(RETAINED,RNA)) + contrast_text("primer3", context=context)
    completion, redesign = "", ""
    for side, label in (("forward", "固定 F：补全 R"), ("reverse", "固定 R：补全 F")):
        context = "amp80-200_fix_" + side
        completion += '<h4>' + label + '</h4>' + pending_attention(context) + result_table(context, keys=(RETAINED,RNA,"anchor_preservation_rate")) + contrast_text(direct, RNA, context)
        context = "redesign_" + side
        redesign += '<h4>' + label.replace("补全", "重设计") + '</h4>' + pending_attention(context) + result_table(context, keys=(RETAINED,RNA,"anchor_preservation_rate")) + contrast_text(direct, RNA, context)
    guidance = pending_attention(variant="single_oligo") + result_table(keys=(RETAINED,RNA,TM), methods=[primary], variants=("single_oligo","full")) + contrast_text(primary, RNA, variant="single_oligo")
    steps = pending_attention(variant="steps1") + result_table(keys=(RETAINED,DIVERSITY,TM), methods=[primary], variants=("steps1","steps4","steps8","full")) + contrast_text(primary, DIVERSITY, variant="steps1")
    experiments = [
        experiment(1,"exp-quality","候选质量与同期方法比较","在共同的模板、设计条件、候选预算和 RNA／基因组参考规则下，比较可用备选的数量及组成质量。",quality),
        experiment(2,"exp-architecture","Attention 架构与迭代生成的作用","分别对照 CNN 与 Attention、直接预测与扩散，考察联合位置建模对候选供给、多样性和温度匹配的影响。",architecture_result),
        experiment(3,"exp-budget","多解生成：候选预算、位置多样性与温度匹配","随 K 增长观察有效方案的积累，并测量备选覆盖不同结合区域的程度和最佳候选的成对温度偏差。",budget),
        experiment(4,"exp-length","同一模型的产物长度控制","保持模型权重不变，分别指定三个产物长度区间，检查条件响应与各长度要求下的候选质量。",length),
        experiment(5,"exp-completion","固定任一侧引物的补全","保留已有 F 或 R 区间，保持原始长度要求，生成兼容伙伴并检查锚点是否原样保留。",completion),
        experiment(6,"exp-redesign","改变条件后的局部重设计","将原始引物对移到不兼容的新产物长度要求中，固定一侧并重新生成另一侧，测试同一模型的修改能力。",redesign),
        experiment(7,"exp-guidance","单条与联合热力学引导","采用相同检查点与采样预算，比较单条引物过滤和进一步加入成对 Tm 差约束后的候选质量。",guidance),
        experiment(8,"exp-steps","去噪深度与候选集形成","使用同一权重比较 1、4、8、16 步恢复过程，分析迭代生成对保留候选、位置分散度与温度偏差的贡献。",steps),
    ]
    template_path = ROOT / "reports/templates/chinese_research.html"
    css = re.search(r'<style>(.*?)</style>',template_path.read_text(), re.S).group(1)
    css += '''
:root{--gold:#21877f;--gold-bg:#e7f0ec;--paper:#f7f7f1}.masthead-meta{letter-spacing:.06em}.draft-label{font-size:12px;color:#836347;border-left:2px solid #c5985f;padding-left:14px;margin-top:20px}.result-figure{margin:28px 0}.result-figure img,#method-figure img{width:100%;height:auto;display:block;border-radius:4px;background:#fff;border:1px solid var(--line)}figcaption{font-size:12px;color:#77746d;line-height:1.85;margin-top:12px}.table-scroll{overflow-x:auto}.result-table{border-collapse:collapse;width:100%;font-size:12px;line-height:1.75;margin:18px 0 25px;min-width:530px}.result-table th{text-align:left;font-weight:500;border-bottom:1px solid #bdb8ac;padding:10px 12px 10px 0;white-space:nowrap}.result-table td{border-bottom:1px solid var(--line);padding:12px 12px 12px 0;vertical-align:top;font-variant-numeric:tabular-nums}.result-table tbody th{border-bottom:1px solid var(--line);vertical-align:top;padding-top:12px}.result-table td small{display:block;color:#8b8780;font-size:9px;margin-top:3px;white-space:nowrap}.result-table .ours{color:#176d68;background:#eef3ef}.pending{color:#8c7a64;font-size:13px;line-height:1.9}.experiment h4{font-weight:550;font-size:15px;margin:25px 0 12px}.method-detail p{font-size:15px}.article{padding-top:70px}.main-note{font-size:12px;color:#77746d;border-left:2px solid #d0c6b4;padding-left:15px}.experiment{scroll-margin-top:75px}@media(max-width:650px){.result-table{font-size:11px;min-width:500px}.result-table td small{font-size:8px}.article{padding-top:42px}.method-detail p{font-size:14px}}
'''
    css += '.author-line{font-size:14px;color:#77746d;line-height:1.8;margin:15px 0 0}.author-line a{color:#21877f;text-decoration:none;overflow-wrap:anywhere}'
    draft = '<p class="draft-label">技术审阅草稿 · 仅展示已完成评估来源；CNN、纯 RoPE Attention 与混合架构分别标注。</p>' if partial or args.allow_partial else ''
    headings = [("abstract","摘要与贡献"),("method","方法与数据"),("experiments","实验"),("limitations","局限")]
    sidebar = ''.join(f'<a href="#{a}">0{i} &nbsp; {b}</a>' for i,(a,b) in enumerate(headings,1))
    mobile = ''.join(f'<a href="#{a}">{b}</a>' for a,b in headings)
    split = ''.join(f'<div class="split-part {"val" if s=="validation" else s}"><b>{split_counts[s]}</b>{label}</div>' for s,label in (("train","训练基因"),("validation","验证基因"),("test","独立测试基因")))
    document = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DiffuPrimer · Attention 研究审阅</title><style>{css}</style></head><body>
<header class="masthead"><a class="masthead-brand" href="#title">DiffuPrimer</a><span class="masthead-meta">ATTENTION · CONDITIONAL DIFFUSION</span></header>
<nav class="outline" aria-label="研究目录"><div class="outline-title">CONTENTS</div>{sidebar}</nav><nav class="mobile-nav" aria-label="分节导航">{mobile}</nav>
<div class="article"><div class="title-block" id="title"><div class="kicker">Learning primer-pair binding positions</div><h1 class="wordmark">Diffu<em>Primer</em></h1><p class="paper-title">Conditional Discrete Diffusion for Diverse,<br>Constraint-Aware qPCR Primer-Pair Design</p><p class="chinese-title">基于条件离散扩散的多样化约束感知 qPCR 引物对设计</p><p class="author-line">{author_line}</p>{draft}<div class="resource-links"><a href="DiffuPrimer_attention_method.pptx" download>英文方法图 · 可编辑 PPT</a><a href="study_summary.groups.csv" download>结果数据</a><a href="study_summary.contrasts.csv" download>配对比较</a></div></div>
<main><section id="abstract"><h2><span class="section-no">01</span>摘要<span class="en">Abstract</span></h2><p class="abstract-copy">qPCR 引物对设计需要同时协调序列特异性、热力学兼容性和扩增子长度，并支持同一靶标上的多个备选方案及已有引物的复用。DiffuPrimer 在固定转录本上表示正、反向引物的结合位置，通过条件离散扩散恢复背景、F 与 R 标签；Attention 去噪器连接不同序列位置，结合结构化采样、联合热力学引导与位置锚定，将完整设计、长度控制、单侧补全和局部重设计统一到同一模型中。</p>{abstract_result}
<div class="contribution-heading" id="contributions">OUR CONTRIBUTIONS</div><ul class="contributions"><li><strong>引物对的位置分布建模。</strong>在固定生物序列上构建 0／F／R 标签的条件离散扩散，采用 Attention 建模结合区域之间的关系，学习同一靶标的多种可行位置方案。</li><li><strong>支持灵活使用的统一约束采样。</strong>将成对区间结构、热力学兼容性与固定引物位置融入恢复过程，以同一模型处理完整设计、长度控制、固定任一端补全以及条件变化下的局部重设计。</li><li><strong>公开参考基准与系统质量评价。</strong>通过同源感知的独立基因划分和完整 RNA／基因组评估，将候选质量、架构机制与任务能力联系起来，并比较成熟方法及适用的同期工作。</li></ul></section>
<section id="method"><h2><span class="section-no">02</span>方法<span class="en">Method</span></h2><p>输入为目标转录本、九维设计条件和可选的固定引物。模型生成结合位置，从已知模板提取 F 序列及 R 区间的反向互补序列，形成可合成的引物对。</p><figure id="method-figure"><a href="DiffuPrimer_attention_method.png" target="_blank"><img src="data:image/png;base64,{method_encoded}" alt="English biological method schematic: DNA helix, transcript nucleotide tokens, clean F/R position labels, absorbing-mask corruption, conditional attention denoising, constrained paired-primer sampling and four flexible design modes."></a><figcaption>图 1 · 英文方法总览：生物模板、位置标签加噪与恢复、条件 Attention 和四种设计场景。点击图片可查看原图；<a href="DiffuPrimer_attention_method.pptx" download>可编辑 PPT</a> · <a href="DiffuPrimer_attention_method.pdf" target="_blank">矢量 PDF</a> · <a href="DiffuPrimer_attention_method.svg" download>SVG</a>。</figcaption></figure>
<div class="method-detail" id="method-detail"><h3><span class="index">2.1</span>Attention 条件去噪与成对约束采样</h3><h4>序列与位置标签共同输入</h4><p>每个核苷酸与其当前标签分别嵌入，标签包括背景 0、正向结合区间 F、反向结合区间 R 和吸收态掩码 M。训练对可生成位置随机遮盖，固定引物的位置保持可见；模型学习从带噪标签预测干净的位置类别。</p><div class="formula">x<sub>0</sub> ∈ {{0, F, R}}<sup>L</sup> &nbsp; → &nbsp; x<sub>t</sub> ∈ {{0, F, R, M}}<sup>L</sup> &nbsp; → &nbsp; (F, R)</div>
<h4>Attention 连接正、反向结合区域</h4><p>{local_encoder_note}去噪器使用 {attention_layers} 个 Transformer block、{architecture['heads']} 个注意力头与 {architecture['width']} 维隐藏表示。Q／K 采用旋转位置编码；九维设计条件与扩散时间共同嵌入，并以逐 Attention block 的缩放与偏移调制注入。多头 Attention 允许一个位置结合远处的序列、带噪标签与已有引物信息，再输出 0／F／R 概率。</p>
<h4>在兼容的引物对空间中恢复</h4><p>每个反向步骤根据网络得分采样兼容的连续 F／R 区间，同时满足先后顺序、引物长度、扩增子长度以及已经恢复的标签。单条引物候选通过 GC、Tm、同聚物与二级结构筛选，成对候选进一步满足联合 Tm 差要求。随后以 1/t 概率揭示剩余掩码位置，逐步恢复完整引物对。固定 R 时先过滤与锚点兼容的 F 支持，再进行 beam 排序。</p>
<h4>训练与使用方式</h4><p>按冻结训练基因构建每轮四个 {train['window_length']} bp 窗口，{train['fixed_side_augmentation']*100:.0f}% 样本保留一个引物区间；AdamW、学习率 {train['learning_rate']:.0e}、batch size {train['batch_size']}，训练 {train['epochs']} 轮并按验证损失选择检查点。设计过程使用 {protocol['sampling']['steps']} 步恢复、候选预算 K={protocol['sampling']['budget']}，所有模式复用同一模型权重。</p><div class="task-modes"><div class="mode"><div class="mode-name">完整设计与多解</div><p>在同一模板与条件下生成多个不同结合位置的备选。</p></div><div class="mode"><div class="mode-name">长度条件控制</div><p>修改扩增子范围，使配对位置响应新的检测要求。</p></div><div class="mode"><div class="mode-name">固定任一侧补全</div><p>保留已有 F 或 R，为其生成兼容伙伴。</p></div><div class="mode"><div class="mode-name">局部重设计</div><p>保留可复用引物，在改变长度条件后重设计另一端。</p></div></div></div>
<div id="data"><h3><span class="index">2.2</span>公开数据收集与独立基因划分</h3><p>参考来自 NCBI 小鼠 GRCm39 RefSeq 的完整 RNA、基因组与基因注释。基准包含 {len(targets)} 个蛋白编码基因；来源引物经模板坐标映射，Primer3 在训练及验证模板上提供计算位置标签。</p><table class="data-table"><thead><tr><th>数据来源</th><th>记录</th><th>用途</th></tr></thead><tbody><tr><td>NCBI RefSeq GRCm39</td><td>{len(targets)} 个基因模板</td><td>序列、转录本注释和完整参考评估</td></tr><tr><td>PrimerBank 来源记录</td><td>{dataset['by_evidence']['unverified']} 个映射位置记录</td><td>公开引物结合位置监督</td></tr><tr><td>Primer3 计算设计</td><td>{dataset['by_evidence']['computational_candidate']} 个位置标签</td><td>训练及验证基因的计算监督</td></tr></tbody></table><div class="split-strip">{split}</div><p class="subtle">按所有 NM_ isoform 的 MMseqs2 核苷酸相似关系划分基因簇，阈值为 80% 一致性及双向 80% 覆盖。训练、验证、测试基因分别用于参数拟合、检查点选择和独立结果评价；测试来源位置记录不进入训练。</p>{external_data_note}</div></section>
<section id="experiments"><h2><span class="section-no">03</span>实验<span class="en">Experiments</span></h2><p class="setup">所有方法使用同一冻结测试基因、设计条件与候选预算。RNA 特异性和基因组潜在产物采用共同完整参考规则；均值与 95% CI 以基因为统计单位。数字只来自完成生成及完整筛查的记录。保留候选指 RNA 特异且无预测基因组产物的唯一合成序列对。</p>{''.join(experiments)}</section>
<section id="limitations"><h2><span class="section-no">04</span>局限<span class="en">Limitations</span></h2><p class="limitations">当前基准围绕小鼠参考转录本和既定比对、热力学规则构建。后续可在更多物种与等位变异背景中扩展目标覆盖，并优化长转录本上的 Attention 推理成本。</p></section></main></div>
<script id="study-data" type="application/json">__STUDY_DATA__</script><script>
'use strict';const D=JSON.parse(document.getElementById('study-data').textContent);const NAME=__METHOD_NAMES__;const METHODS=__METHOD_ORDER__;const PRIMARY=__PRIMARY_METHOD__;
const KEYS=['{RETAINED}','{DIVERSITY}','{TM}'];function fmt(v,k){{return v===null?'—':Number(v).toFixed(k==='{RETAINED}'?2:3)}}
function updateBudget(k){{const tbody=document.querySelector('#budget-results tbody');if(!tbody)return;tbody.replaceChildren();for(const method of METHODS){{const g=D.summaries.find(x=>x.context==='amp80-200'&&x.method===method&&x.variant==='full'&&x.budget===k);if(!g)continue;const tr=document.createElement('tr');if(method===PRIMARY)tr.className='ours';const th=document.createElement('th');th.scope='row';th.textContent=NAME[method];tr.append(th);for(const key of KEYS){{const td=document.createElement('td');Object.assign(td.dataset,{{method,context:g.context,variant:g.variant,budget:String(k),endpoint:key,value:JSON.stringify(g[key])}});const span=document.createElement('span');span.textContent=fmt(g[key],key);const small=document.createElement('small');const ci=g[key+'_ci95'];small.textContent=ci?'95% CI '+fmt(ci[0],key)+'–'+fmt(ci[1],key):'';if(key!=='{RETAINED}'&&g[key+'_n']!==undefined)small.textContent+=' · n='+g[key+'_n'];td.append(span,small);tr.append(td)}}tbody.append(tr)}}document.querySelectorAll('.budget-buttons button[data-budget]').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.budget)===k)))}}
document.querySelectorAll('.budget-buttons button[data-budget]').forEach(b=>b.addEventListener('click',()=>updateBudget(Number(b.dataset.budget))));window.researchReportReady=true;
</script></body></html>'''
    # Only the candidate-budget switch needs client-side statistics. Full source
    # records and paired comparisons remain in the unchanged downloadable JSON.
    budget_keys = (RETAINED, DIVERSITY, TM)
    budget_fields = {"context", "method", "variant", "budget"} | {
        key + suffix for key in budget_keys for suffix in ("", "_ci95", "_n")}
    interactive_data = {"summaries": [
        {key:value for key,value in g.items() if key in budget_fields}
        for g in summary["summaries"]
        if g["context"] == "amp80-200" and g["variant"] == "full"
        and g["budget"] in (1,2,5,10,20)]}
    data = json.dumps(interactive_data,ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('<','\\u003c')
    document = document.replace('__STUDY_DATA__',data).replace('__METHOD_NAMES__',json.dumps(NAMES,ensure_ascii=False)).replace('__METHOD_ORDER__',json.dumps(method_order)).replace('__PRIMARY_METHOD__',json.dumps(primary))
    font_candidates = [Path('/tmp/dfm-html-check/fonts/NotoSansCJKsc-Regular.otf'), Path('/root/.local/share/fonts/NotoSansCJKsc-Regular.otf')]
    original_font = next((p for p in font_candidates if p.is_file()),None)
    font_output = out / 'report_chinese_subset.otf'
    if original_font:
        from fontTools import subset
        from fontTools.ttLib import TTFont
        font = TTFont(original_font)
        sub = subset.Subsetter(options=subset.Options())
        sub.populate(text=''.join(sorted(set(document))))
        sub.subset(font)
        font.save(font_output)
    else:
        shutil.copyfile(ROOT / 'reports/assets/report_chinese_subset.otf',font_output)
    document = document.replace('__CHINESE_FONT__',base64.b64encode(font_output.read_bytes()).decode())
    (out / 'index.html').write_text(document)
    for suffix in ('.json','.groups.csv','.contrasts.csv'):
        source = summary_path.with_suffix(suffix)
        if suffix == '.json':
            (out / 'study_summary.json').write_bytes(summary_bytes)
        else:
            shutil.copyfile(source,out / ('study_summary' + suffix))
        resource_files.append('study_summary' + suffix)
    shutil.copyfile(ROOT / 'reports/assets/FONT_LICENSE.txt',out / 'FONT_LICENSE.txt')
    write_json(out / 'build_manifest.json', {'built_at':utc_now(),'presentation':'attention_technical_review','status':summary['status'],'presentation_status':'working_preview' if args.allow_partial else 'completed_results_review',
        'page_sha256':sha256(out / 'index.html'),'builder_sha256':sha256(Path(__file__)), 'summary_sha256':summary_sha,
        'protocol_path':str(protocol_path),'protocol_sha256':sha256(protocol_path),'dataset_manifest_sha256':sha256(data_path),'method_figure_sha256':sha256(method_png),
        'method_figure_directory':str(method_dir),'method_figure_basename':method_basename,
        'sections':['Abstract and contributions','Method and public data','Eight experimental questions','Limitations'],
        'contribution_count':3,'experiment_sections':8,'used_result_figures':used_figures,'resources':sorted(set(resource_files)),
        'numeric_instances':numeric_instances,'external_runtime_dependencies':[],'original_website_preserved':True,
        'journal_manuscript_written':False,'method_architecture':architecture,'primary_method':primary,'direct_method':direct,
        'external_summary_sha256':external_sha,'external_total_genes':external.get('total_test_genes') if external else None,'cohorts_pooled':False,
        'exon_case_sha256':exon_case_sha,
        'author_metadata_sha256':sha256(author_path),'author_metadata':author_metadata,
        'interactive_data_bytes':len(data.encode()),'complete_statistics_download':'study_summary.json',
        'source_policy':'Completed source runs from the unified summary only; old CNN and rotary-attention/hybrid architectures explicitly labeled'})
    (out / 'README.txt').write_text('DiffuPrimer Attention 技术审阅页\n直接打开 index.html 可离线查看。方法图、字体和结果图已嵌入。\n本页独立于原 CNN 研究网页，不改写投稿论文。\n复现：.venv/bin/python scripts/build_attention_research_page.py --protocol experiments/attention/final_protocol.json --external-summary experiments/attention_external/results/study_summary.json\n计算尚未完成时加 --allow-partial，页面会明确标为技术审阅草稿。\n')
    print(json.dumps({'page':str(out / 'index.html'),'status':summary['status'],'sections':4,'contributions':3,'experiment_questions':8,'result_figures':used_figures},ensure_ascii=False))


if __name__ == '__main__':
    main()
