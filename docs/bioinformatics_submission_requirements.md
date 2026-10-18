# Bioinformatics 投稿要求与本项目实施清单

当前英文研究初稿另存于 `paper/bioinformatics_cnn_draft/`，正文实测七模板页、138 词结构摘要、三点贡献、八个实验问题、四幅核心图、52 条已核实且实际引用的参考文献；补充材料八页，保存四幅完整图和四张完整对照表。参考师姐论文的科学叙述与图组层次，使用目标期刊的作者年份引用。此项更新不把 AI 辅助研究初稿标为已经满足该刊的 LLM 起草政策，也不覆盖历史正文或结果工作表。代码与数据网址按作者要求继续留空。

核对日期：2026-10-13。本轮重新检索官方指南并读取本地 OUP 类与当前出版图表清单；官方页面直接打开超时，条款由 Oxford Academic 官方页面的搜索缓存核对。目标为 **Bioinformatics Original Paper**；建议投稿主题分类为 **Sequence analysis**。后者是结合引物结合位点建模内容作出的项目建议，并非期刊指定分类。

## 官方要求速查

以下是 [Bioinformatics author guidelines](https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines) 的压缩提要；具体条款位置列在括号中：

- Original Paper: 7 template pages, approximately 5,000 words excluding figures; actual biological data. Application Note: 4 pages, approximately 2,600 words. (Types of Manuscript)
- Abstract: recommended maximum 150 words; Motivation, Results, Availability and Implementation, Contact, Supplementary Information. (Structured Abstract)
- ML: dedicated training-data subsection; independent test; homology-aware sequence partitioning. (Machine learning)
- Software: working implementation, stable HTTPS source URL, documented installation, parameters and test data. (Software)
- Data: availability statement, persistent public deposits, full dataset citations. (Supporting Data; Data availability)
- Figures: embedded, white background, adjacent captions, 86/178 mm width; alt text under each legend; line art 1,200 dpi, colour/halftone 350 dpi. (Figures; How to prepare text and figures)
- Supplement: consolidated text/figures/tables, separately submitted and cited. Cover letter、真实作者／单位／联系信息、基金／利益声明和 CRediT 贡献需准备；提交作者需提供 ORCID。 (Supplementary Data; Cover letter; Authorship)
- LLM policy: code assistance/correction/translation allowed; prompt-based paper drafting prohibited. Exact excerpt: “It is not acceptable to use LLMs or related technologies to draft papers … from a prompt text.” Accepted use requires cover-letter/manuscript disclosure and supplementary documentation. (Policy for Acceptable Use of Large Language Models)

省略号仅省去该条款中关于文本、图、表和参考文献的范围列举；不是新增例外。原文禁止的是通过 prompt 草拟论文及其组成内容；并未禁止代码辅助、结果异常核查或作者已有文本的语言纠正。最终作者稿应由研究者撰写，不能把 AI 草拟稿标成已符合该政策。AI 辅助的实际用途、工具与使用过程应如实记录。此处记录期刊政策及正式投稿的要求，不限制已经授权的代码、计算实验和结果核查工作；当前研究初稿与完全投稿合规的状态分别记录，不以这条说明重复中断工作。

## 官方 LaTeX 模板已获取

已保存到 `paper/bioinformatics_template/`：

- `oup-authoring-template/`：OUP 原始 LaTeX 类、样例、参考文献样式、说明 PDF，保留原许可证与版权；版本 1.5，2026-07-14。
- `oup-authoring-template.zip`：原始包；下载 URL、实际镜像、SHA-256 与大小见 `source_manifest.json`。
- `OUP_LaTeX_supported_jnls.xlsx`：OUP 官方期刊选项表。

官方选项表中 Bioinformatics 行为 **author-date** 引文以及：

```tex
\documentclass[numsec,webpdf,modern,large]{oup-authoring-template}
```

样例文件默认是 `contemporary`，不能直接照用。PREDAC-FluB 的 Briefings in Bioinformatics 排版在此表中是 `contemporary,large,numbered`，与 Bioinformatics 不同；可借鉴其科学图设计语言，同时使用目标期刊的排版和引文方式。

获取依据与路径：

- [OUP 准备和提交稿件：LaTeX 模板与期刊选项](https://academic.oup.com/pages/for-authors/journals/preparing-and-submitting-your-manuscript)，其中链接官方期刊表和下载包。
- [OUP 官方期刊选项表](https://static.primary.prod.gcms.the-infra.com/static/site/journals/document/OUP_LaTeX_supported_jnls.xlsx?node=bc6b6c99a53a8b765b76)。本地已直接读取 XLSX 验证 Bioinformatics 行。
- [OUP General Template（Oxford University Press 提供）](https://www.overleaf.com/latex/templates/oup-general-template/ybpypwncdxyb)。
- [CTAN：OUP authoring template](https://ctan.org/pkg/oup-authoring-template)。

## 当前完成状态与作者待办

最终统计来源仍为 `experiments/attention/results/study_summary.json` 和额外独立基因的 `experiments/attention_external/results/study_summary.json`。目录名中的 attention 不表示当前论文主模型。当前出版图表以 **条件 CNN + 离散扩散（`cnn_diffusion`）** 为主模型；CNN 直接预测、纯 RoPE Attention、CNN + Attention 保留为机制和架构对照。这是**完整评估后作出的展示选择**，不是预先指定 CNN 获胜，也不改写冻结协议、权重、候选或统计。两个 75 基因测试队列各自分析，训练集为原划分的 350 个基因、验证集为 75 个基因；额外队列不跨队列配对或混合。

主模型的实际结构见 `src/diffuprimer/model.py` 与 `paper/figures/architecture_method/README.md`：DNA 与带噪位置标签各自嵌入，加入两维相对位置特征和九维设计条件／四维时间特征；**条件只在输入处相加**。四层宽度 96 的残差膨胀 Conv1d 使用 1、4、16、64 膨胀率，三分类头输出背景／F／R logits。主模型没有 Attention、FiLM 或 adaLN；不能沿用旧混合 Attention 图或说明。

| 项目 | 当前状态 | 作者下一步 |
| --- | --- | --- |
| 计算研究 | 模型训练、两队列完整 RNA／基因组参考筛查、统计及来源检查均已完成 | 用最终统一来源撰写方法与结果 |
| 比较方法 | Primer3、qPrimer、PrimerCast、CNN／纯 Attention／混合 Attention、direct／diffusion 及采样对照分别保留真实标签 | 交代适配范围、共同条件和候选预算 |
| 灵活使用 | 产物长度、保留 F 或 R 的补全、改变条件后的重设计、候选预算和多样性已评估 | 将灵活使用提炼为独立贡献，配合对应实验 |
| 作者元数据 | `paper/author_metadata.json`：Jingyue Li；Sun Yat-sen University；18837614038@163.com | 提供实际完整单位地址、ORCID、基金、利益与贡献信息 |
| 当前研究审阅页 | `reports/cnn_publication/index.html`，四节、三点贡献、八个实验；`analysis.html` 可单独查看分析图 | 作者核对最终论证和图注 |
| 当前出版图 | 下表五张整合图及其 PDF／SVG／350 dpi PNG／TIFF；主方法图另有可编辑 PPTX | 按期刊最终物理尺寸检查和导出 |
| 历史排版测试 | `paper/bioinformatics_attention_materials/author_results.pdf` 是旧混合 Attention 结果工作表，实测六页；`paper/bioinformatics_working/` 是更早旧 CNN 八页布局测试 | 这两个测量均不能证明当前 CNN 五图作者稿已满足七页；须重排并重新测量 |
| 正式作者稿与发布 | 科学作者稿正在准备；代码／数据 URL 按用户要求留空 | 写完正文及声明后测量；提交前填真实稳定公开链接和存档 |

## 当前五图与八个实验问题

路径及内容以 `docs/cnn_publication_figures.md`、当前图表 manifest／provenance 为准；旧 Attention 素材留档。

| 主图 | 当前 PDF 路径 | 所承载的实验问题 |
| --- | --- | --- |
| Figure 1 | `paper/figures/architecture_method/DiffuPrimer_method.pdf` | 条件 CNN 计算图及展开残差块；位置标签加噪／恢复；四种设计用途 |
| Figure 2 | `paper/figures/cnn_publication_results/fig2_candidate_quality.pdf` | ①候选质量与独立基因泛化；两队列分别报告；保留全部基因的配对散点 |
| Figure 3 | `paper/figures/cnn_publication_results/fig3_flexible_design.pdf` | ②候选预算／多样性；③长度控制；④单侧补全；⑤条件改变后的局部重设计 |
| Figure 4 | `paper/figures/cnn_biological_analysis/fig4_biological_specificity.pdf` | ⑥外显子结构与检测特异性；完整队列组成及实际 Susd5 引物／基因组产物记录 |
| Figure 5 | `paper/figures/cnn_publication_results/fig5_mechanism.pdf` | ⑦架构与扩散机制；⑧联合约束与恢复深度 |

实体 qPCR 不属于用户确认的研究范围；论文使用公开真实生物数据与计算评估。固定训练 seed 保持 20261013。主文应通过这八个问题串联科学结论，而不是计数实验运行次数。Fig. 4 的 Susd5 只是与完整队列统计配套的真实位置机制示例，不能替代总体方法评估。

## 章节、引用和提交材料

官方章节清单的顺序是：标题页 → 结构摘要 → 引言 → 系统与方法 → 算法 → 实现 → 讨论 → 参考文献。它与网页四节不是同一个结构。建议把数据、监督来源、独立基因划分和评价定义写入系统与方法；位置扩散、条件化、联合区间采样写入算法；实现包括训练／采样参数、基线适配、复现与计算评估。若独立设置 Results 节，应明确这是作者对评估内容的组织安排，而不是声称官方清单原样列出了 Results。 [官方章节与参考文献说明](https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines)

正文引用采用本地官方期刊选项表的 author-date；当前官方指南另指向 Oxford SCIMED。使用 `\citet{key}` 作句内作者引用、`\citep{key}` 作括号引用；让 BibTeX 与 `oup-abbrvnat` 生成参考文献，勿手工编号。同期比较方法、软件版本、参考组装／数据存档都需实际引用；摘要不放引用、公式或对其他章节的依赖。后者也见已缓存 OUP 样例 `oup-authoring-template.tex`。

额外正文、图和表合并成独立补充文件，并在正文恰当位置引用；主要模型定义和关键结果留正文。作者需准备基金（无基金亦据实说明）、利益冲突、CRediT、提交作者 ORCID、cover letter、代码／数据可获得性及实际辅助用途记录。补充资料经期刊永久保存，不能只给一个随时可能消失的个人网页。图注邻近图，Alt text 放图注之后。 [官方投稿、补充资料、作者声明和图说明](https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines)

## 可直接采用的 LaTeX 骨架

以下仅是已缓存官方模板的语法及填写区，不包含代拟摘要或正文。OUP 期刊选项表指定 `numsec,webpdf,modern,large`；本地类默认数值引用，补上它公开支持的 `namedate` 以稳定启用 author-date。它内部加载 natbib，勿再以冲突选项重复加载。

```tex
\documentclass[numsec,webpdf,modern,large,namedate]{oup-authoring-template}
\usepackage{booktabs}
\pdfminorversion=5
\newcommand{\societylogo}{}
\graphicspath{{figures/}}
\begin{document}
\journaltitle{Bioinformatics}
\appnotes{Original Paper}
\firstpage{1}
\title[DiffuPrimer]{DiffuPrimer: Conditional Discrete Diffusion for Diverse,
Constraint-Aware qPCR Primer-Pair Design}
\author[1,$\ast$]{Jingyue Li}
\address[1]{\orgname{Sun Yat-sen University}}
% 补充作者提供的真实 department / city / postcode / country。
\corresp[$\ast$]{Corresponding author.
\href{mailto:18837614038@163.com}{18837614038@163.com}}
\abstract{\textbf{Motivation:} [Author text.]\\
\textbf{Results:} [Author text.]\\
\textbf{Availability and Implementation:} [URL pending.]\\
\textbf{Contact:} \href{mailto:18837614038@163.com}{18837614038@163.com}\\
\textbf{Supplementary Information:} [Author text.]}
\keywords{qPCR, primer-pair design, conditional discrete diffusion}
\maketitle
\section{Introduction}
\section{System and methods}
\section{Algorithm}
\section{Implementation}
\section{Discussion}
% 在正文末尾补作者确认的 availability / acknowledgements /
% funding / conflicts / contribution 声明及补充材料指引。
\bibliographystyle{oup-abbrvnat}
\bibliography{references}
\end{document}
```

类文件与 `oup-abbrvnat.bst` 从 `paper/bioinformatics_template/oup-authoring-template/` 获取，保留许可证。标题、作者元数据来自本项目；摘要所有实质文本等待作者填写。`DOI`、期号、卷号、正式发布日期等由出版流程确定，不写虚构值。官方样例提供 `\figalttext{...}`，但本地 v1.5 类将它定义为 `\relax`（不显示任何内容）。因此仅调用它并不能满足“图注下方可见 Alt text”的稿件要求；还要在图注之后排出 `\par\small\noindent\textbf{Alt text:} ...` 等实际文本，并核对上传 PDF。

## 七页版式注意点与项目建议

**已核对的本地类几何：**`modern,large` 的成品页是 **210 × 276 mm**，不是普通 A4 的 210 × 297 mm；版心宽 178 mm、栏距 6 mm，因此单栏为 86 mm。保持默认双栏，不加 `geometry` 改页边距、改 `\textwidth`／`\textheight`、强行缩正文或启用 `\onecolumn` 来凑页数。官方原类不修改。 [OUP 模板与期刊选项入口](https://academic.oup.com/pages/for-authors/journals/preparing-and-submitting-your-manuscript)

五图是项目的整合选择，不是期刊强制图数。多面板图使用 `figure*` 配 `width=\textwidth`；真正的单栏图用 `figure` 配 `width=\columnwidth`。当前方法图成品 178 × 116 mm，文字至少 7 pt；缩放后须重新检查，尤其不要把七面板或长坐标轨迹压到单栏。图保持白背景，色彩图／半色调投稿资产为 350 dpi、纯线图按 1200 dpi；正文保留矢量 PDF。 [官方图资产要求](https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines)

建议摘要控制在约 140–150 词，正文先以约 3,000–3,500 词规划，给五图、必要表、图注和引用留出空间。约 5,000 词是页面限制的近似说明，不能理解为带五张大图仍保证七页。论文论证建议用一个紧凑主表，避免再复制五幅图中每个面板的详表；完整基因值、所有条件、输入参数、设计候选及扩展组件表放补充材料，关键不利结果仍按科学问题如实解释。

编译建议 `pdflatex → bibtex → pdflatex → pdflatex`，再用 `pdfinfo` 数实际页面，逐页检查 PDF。双栏跨栏浮动可能延后，避免在正文密集加 `FloatBarrier`、`[H]` 或手工 `\clearpage` 导致空白页；依据首次引用与图高调节顺序，保持图注邻近。参考文献、浮动图、声明全部加入后再认定七页。官方 `\lastpage{N}` 仅记录真实最后页，不能替代物理页数测量；旧工作表的六页不包含当前论文作者正文，也不能直接作为新稿页数。

## 最后核对清单

- [x] 条件 CNN 主模型及 Attention 架构对照、同期方法的真实结果进入统一统计来源，展示选择明确标为 post-evaluation。
- [x] 原测试及额外独立基因队列完成，来源与序列／坐标一致性核查通过。
- [x] 官方模板、author-date 技术样式与真实作者元数据已获取并核对；当前作者稿可采用本文骨架。
- [x] 当前 CNN 主线五图具备实际尺寸、图注、Alt text 和矢量／光栅资产。
- [ ] 作者完成约 150 词结构摘要、科学正文、结果解释和三点贡献。
- [ ] 当前 CNN 五图正式作者稿在全部正文填写后按期刊模板测量不超过七页。
- [ ] 代码、权重、数据／参考清单、split 与结果的真实稳定 HTTPS 地址可用（当前按用户要求留空）。
- [ ] 发布版本的软件许可、安装和小规模示例实际可运行。
- [ ] 统一补充 PDF、Funding、conflicts、CRediT、cover letter 和作者 ORCID 完成。
- [ ] ScholarOne 转换后的 PDF 内容、字体和图正确。

本清单记录已完成计算研究与投稿交付的实际边界；不把技术工作表等同于正式作者稿或投稿完成。
