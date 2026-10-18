#!/usr/bin/env python3
"""Prepare an auditable Bioinformatics manuscript candidate from frozen evidence.

This does not certify journal LLM-policy compliance or submit anything.
"""
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / 'paper/bioinformatics_cnn_draft'
OUT = ROOT / 'paper/bioinformatics_submission'
MECHANISM = ROOT / 'experiments/mechanism/uncapped_seed42/summary.json'
URL = 'https://github.com/jingyue0513/DiffuPrimer'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_json(p, data):
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def replace_once(text, old, new):
    assert text.count(old) == 1, old[:100]
    return text.replace(old, new, 1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ['oup-authoring-template.cls', 'oup-abbrvnat.bst', 'draft_formatting.tex',
                 'references.bib', 'submission_references.bib', 'reference_verification.json',
                 'result_values.tex', 'supplement_tables.tex']:
        shutil.copy2(DRAFT / name, OUT / name)
    (OUT/'submission_formatting.tex').write_text(r'''% Submission-only headers; official OUP class is unmodified.
\makeatletter
\def\ps@opening{
  \def\@oddhead{\hbox to \textwidth{\sffamily\small Bioinformatics\hfill\@appnotes}}
  \let\@evenhead\@oddhead
  \def\@oddfoot{\hfil\thepage\hfil}
  \let\@evenfoot\@oddfoot
}
\makeatother
''')
    shutil.copytree(DRAFT / 'figures', OUT / 'figures', dirs_exist_ok=True)
    shutil.copy2(MECHANISM.parent / 'results.pdf', OUT / 'figures/mechanism_controls.pdf')
    s = json.loads(MECHANISM.read_text())
    groups = {(g['cohort'], g['method']): g for g in s['groups']}
    methods = ['direct_validation_selected', 'diffusion1_validation_selected', 'frozen_4',
               'diffusion_4', 'frozen_16', 'diffusion_16', 'direct_equal_cpu_uncapped']
    names = ['Direct, validation-selected', 'Denoising, one step, validation-selected',
             'Frozen scores, 4 steps', 'Full diffusion, 4 steps', 'Frozen scores, 16 steps',
             'Full diffusion, 16 steps', 'Direct, matched CPU deadline']
    table = [r'\begin{table*}[t]\centering\small',
             r'\caption{Mechanism controls with seed 42. Counts are unique retained pairs per gene; coverage is the percentage of genes with at least one retained pair. Each cohort contains 75 genes; $K=20$. CPU seconds cover generation only.}\label{tab:mechanism}',
             r'\begin{tabular}{lrrrrrr}\toprule',
             r'&\multicolumn{3}{c}{Original cohort}&\multicolumn{3}{c}{Additional cohort}\\',
             r'Method & Pairs & Coverage & CPU (s) & Pairs & Coverage & CPU (s)\\\midrule']
    for method, name in zip(methods, names):
        a, b = groups['test', method], groups['external_test', method]
        vals = [f"{g['retained_pairs']:.2f} & {100*g['one_retained']:.1f} & {g['generation_cpu_seconds']:.2f}" for g in [a,b]]
        table.append(name + ' & ' + ' & '.join(vals) + r'\\')
    table.extend([r'\bottomrule\end{tabular}', r'\end{table*}'])
    (OUT / 'mechanism_table.tex').write_text('\n'.join(table) + '\n')
    with (OUT / 'mechanism_results.csv').open('w', newline='') as f:
        fields=['cohort','method','retained_pairs','one_retained','two_retained','generation_cpu_seconds']
        w=csv.DictWriter(f,fields);w.writeheader()
        for g in s['groups']:w.writerow({k:g[k] for k in fields})
    write_json(OUT / 'mechanism_evidence.json', s)
    authors = r'''\author[1,2,3]{Jingyue Li}
\author[2,3,$\ast$]{Shiyue Mei}
\address[1]{Sun Yat-sen University, Guangzhou, China}
\address[2]{Zhengzhou University, Zhengzhou, China}
\address[3]{Henan Children's Hospital, Zhengzhou, China}
\corresp[$\ast$]{Corresponding author: Shiyue Mei. \href{mailto:18837614038@163.com}{18837614038@163.com}. Jingyue Li: \href{mailto:lijingyue513@163.com}{lijingyue513@163.com}.}'''
    abstract = r'''\abstract{\textbf{Motivation:} Primer-pair design requires compatible alternatives and adaptation to changing constraints.\\
\textbf{Results:} DiffuPrimer models paired binding positions using conditional discrete diffusion. Constrained recovery supports length control, fixed-primer completion and redesign. Two separate 75-gene mouse cohorts show increased reference-screened retention relative to Primer3. Strong direct decoding and frozen-score controls attribute most count gains to denoising-trained scores and constrained sampling. Sixteen-step network updates provide no consistent count gain over frozen scores, while improving retained-pair positional diversity. Equal-CPU direct search provides greater gene coverage but fewer retained pairs.\\
\textbf{Availability and Implementation:} Python source, checkpoints and test data: \url{https://github.com/jingyue0513/DiffuPrimer}.\\
\textbf{Contact:} \href{mailto:18837614038@163.com}{18837614038@163.com}.\\
\textbf{Supplementary Information:} Supplementary methods, tables and figures accompany this article.}'''
    text=(DRAFT/'manuscript.tex').read_text()
    text=re.sub(r'^%.*\n','',text,flags=re.M)
    text=text.replace('pdfauthor={Jingyue Li}', 'pdfauthor={Jingyue Li and Shiyue Mei}')
    text=text.replace(r'\appnotes{Research draft for author review}',r'\appnotes{Original Paper}')
    start=text.index(r'\author[');end=text.index(r'\abstract{')
    text=text[:start]+authors+'\n'+text[end:]
    start=text.index(r'\abstract{');end=text.index(r'\keywords')
    text=text[:start]+abstract+'\n'+text[end:]
    start=text.index('This study makes three contributions.')
    end=text.index('\\begin{figure*}',start)
    text=text[:start]+('DiffuPrimer represents each pair as forward and reverse binding-position labels on a fixed transcript. Conditional recovery supports alternative generation, product-length control, fixed-primer completion and local redesign. We evaluate these uses on homology-aware mouse-gene partitions and a separate additional-gene cohort, then use frozen-score and strong-decoding controls to distinguish denoising training from iterative network updates.\n\n')+text[end:]
    text=text.replace(' Panel grouping follows the biological/computational hierarchy exemplified by \\citet{predac}.', '')
    text=text.replace('Generative primer work, including Primer C-VAE, motivates learned design \\citep{primercvae}.', 'Learned sequence-based design has been explored in Primer C-VAE \\citep{primercvae}.')
    text=text.replace('The fixed training/generation seed is 20261013. Main comparisons use 16 recovery steps.',
                      'Historical training and main-benchmark generation use seed 20261013 and 16 recovery steps. The added mechanism study reuses one historical Direct/diffusion checkpoint pair with constant generation and bootstrap seed 42. No new models were trained.')
    text=text.replace('Iterative generation thus improves the matched learned proposals beyond merely returning more items.',
                      'These matched-model differences combine training and sampling effects; the controls below test the role of iterative updates.')
    # Keep the complete biological panels in the unified supplement to allow
    # the new causal controls in the seven-page main article.
    text=re.sub(r'\\begin\{figure\*\}\[t\]\\centering\\includegraphics\[width=\\textwidth\]\{figures/main_04_biology.pdf\}.*?\\end\{figure\*\}', '',text,flags=re.S)
    text=text.replace(r'Fig.~\ref{fig:biology}A--C',r'Supplementary Fig.~S3B--D')
    start=text.index(r'\subsection{Diffusion and constraint-guidance mechanisms}')
    end=text.index(r'\section{Discussion}',start)
    text=text[:start]+r'''\subsection{Constraint guidance and architecture}
Raw independent labels yield no valid pairs. Geometry support yields \GeometryRet{} retained pairs, single-oligo guidance \OligoRet{}, and full pair compatibility \MainCNNRet{} (Supplementary Fig.~S4). CNN diffusion exceeds the tested pure-attention and hybrid models: \MainCNNRet{} versus \MainAttentionRet{} and \MainHybridRet{} in the original cohort, and \ExtraCNNRet{} versus \ExtraAttentionRet{} and \ExtraHybridRet{} in the additional cohort. This architecture comparison is specific to the available training data and configurations.

\subsection{Frozen scores, strong direct decoding and CPU budgets}
The added mechanism study separates learned scores from iterative updates (Table~\ref{tab:mechanism}; Supplementary Fig.~S5). Frozen-score recovery calls the diffusion network once and reuses its initial scores while retaining the full decoder, label-reveal process and step count. Strong static decoders enumerate compatible pairs, sample without replacement, remove synthesis-sequence duplicates and optionally select diverse candidates. Both networks receive the same validation-only search over top-$K$ and temperatures 0.5, 1, 2 and 4; selected settings are locked before test evaluation.

Validation-selected one-step denoising exceeds validation-selected Direct by 1.83 pairs (95\% paired CI $[-0.21,3.97]$) and 2.51 ($[0.59,4.52]$) in the original and additional cohorts. Sixteen-step updates versus frozen scores change counts by only 0.03 ($[-0.24,0.28]$) and $-0.24$ ($[-0.55,0.05]$), while increasing retained-pair Jaccard diversity by 0.061 ($[0.025,0.103]$, $n=43$ jointly defined genes) and 0.045 ($[0.026,0.065]$, $n=49$). Counts therefore do not establish a need for sixteen learned updates; diversity identifies an additional benefit.

The deadline-driven Direct baseline samples without a candidate-count cap and selects by temperature compatibility and positional dispersion. It uses 99.2\% of the per-gene sixteen-step process-CPU allowance and retains 4.28 pairs in each cohort, compared with 7.43 and 8.21 for full diffusion. Direct covers 80.0\% of genes in both cohorts, versus 70.7\% and 73.3\% for diffusion. This comparison demonstrates a count--coverage trade-off under the specified selection objective. It measures generation CPU time, excluding shared span preparation and downstream reference screening; it establishes neither equal GPU FLOPs nor end-to-end speed. Gene-paired intervals use 10,000 resamples and seed 42 and are unadjusted; a single checkpoint pair cannot measure training-seed uncertainty.
\input{mechanism_table.tex}

''' +text[end:]
    start=text.index(r'\section{Discussion}');end=text.index(r'\section*{Data and code availability}',start)
    text=text[:start]+r'''\section{Discussion}
DiffuPrimer unifies paired binding positions, interpretable oligos and partial observations. Compatible recovery supports alternative generation, length control, completion and redesign without task-specific retraining. Public-reference screening identifies candidate supply and composition gains across two disjoint mouse-gene cohorts.

The mechanism controls refine the attribution of these gains. Denoising-trained scores retain utility under strong static decoding. Repeated constrained sampling with frozen scores accounts for much of the retained-count performance; sixteen learned updates chiefly add positional diversity and some coverage benefits at higher CPU cost. Greater counts do not imply greater gene coverage, as the deadline-driven Direct control illustrates. The study therefore supports a useful conditional design formulation without establishing that long diffusion trajectories are necessary for yield.

The benchmark evaluates reference discrimination under specified alignment rules rather than amplification efficiency or experimental assay success. Weak supervision is largely computational. Primer3 and qPrimer retain advantages in temperature matching, dispersion and several anchored tasks. Same-species unseen genes assess transfer within the reference; cross-species and variant-background generalization remain untested. Scaling data or model size should retain Direct, one-step, frozen-score and full-update controls to determine which component benefits.

''' +text[end:]
    start=text.index(r'\section*{Data and code availability}')
    text=text[:start]+r'''\section*{Data and code availability}
Reference sequences and annotations are from NCBI RefSeq assembly GCF\_000001635.27, annotation RS\_2024\_02. Source code, checkpoints, test data, processed benchmark records and reproduction instructions are prepared at \url{https://github.com/jingyue0513/DiffuPrimer}. The release bundles test sequences, model weights and machine-readable evaluation records.

\section*{Acknowledgements}
\subsection*{Funding}
This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.
\subsection*{Conflict of interest}
The authors declare no conflict of interest.
\subsection*{AI assistance}
OpenAI Codex assisted with code, analyses, figures, documentation and manuscript drafting and revision. Details are provided in the Supplementary Material.

{\small\setlength{\bibsep}{0pt}\bibliographystyle{oup-abbrvnat}\bibliography{submission_references}}
\end{document}
'''
    text=text.replace(r'\input{draft_formatting.tex}',r'\input{submission_formatting.tex}').replace(r'\copyrightyear{2026}\pubyear{2026}',r'\copyrightyear{}\pubyear{}')
    (OUT/'manuscript.tex').write_text(text)
    supp=(DRAFT/'supplementary.tex').read_text()
    start=supp.index(r'\author[');end=supp.index(r'\abstract{')
    supp=supp[:start]+authors+'\n'+supp[end:]
    supp=supp.replace('research draft','manuscript candidate').replace(r'\appnotes{Supplementary research material}',r'\appnotes{Supplementary Material}')
    start=supp.index(r'\section{Assistance and author-completion record}')
    end=supp.index(r'\clearpage',start)
    supp=supp[:start]+r'''\section{Added mechanism controls and fixed-seed protocol}
The added study reuses the historical CNN Direct and diffusion checkpoints trained with seed 20261013. Generation and gene bootstrap use the constant 42, never a timestamp. Each cohort contains all 75 frozen genes. Eleven fixed controls, two validation-selected static decoders and one uncapped deadline baseline give 14 configurations per cohort. Original historical outcomes retain their own seed; mechanism results are not substituted into historical figures.

Frozen-score controls evaluate the network only on the fully masked initial state; subsequent compatible sampling and revelation use the same seed and step schedule as full recovery. Full-support static decoders enumerate all geometry-, single-oligo- and pair-Tm-compatible pairs and remove sequence duplicates. Candidate selection never uses RNA/genome screen outcomes or exon annotations. For the pool controls, the first pair minimizes the shared Tm penalty; subsequent pairs maximize minimum amplicon Jaccard distance, with Tm and pool order as ties. The five static configurations are deterministic top-$K$ and weighted sampling without replacement at temperatures 0.5, 1, 2 and 4. Each network is selected independently on all 75 validation genes by retained counts, then at-least-two coverage and configuration order. Direct selects temperature 2; diffusion selects top-$K$. Selection precedes test generation and screening.

Generation CPU accounting includes network prediction, compatible support, sampling, deduplication, pair assessment and final selection. Both sides exclude model loading and shared cached span preparation. Each sixteen-step diffusion gene supplies its measured process-CPU budget. The initial bounded Direct control reached its 4,096-pair cap on every target, so a separately recorded follow-up removes the cap without changing the score model, seed or selection objective. Eight worker processes each use four PyTorch threads; the process CPU clock isolates per-target consumption. Selection checkpoints occur every 512 pool additions with a final reserved selection. Only on-time committed outputs are evaluated. Mean CPU utilization is 99.20\% and 99.15\%; minimum utilization is 98.71\% and 98.62\% for the original and additional cohorts. Maximum actual non-preemptible overrun is 0.00765 s; late work never changes the evaluated output. All 150 outputs respect their commit deadlines. Mean explored pools contain 28,014 and 25,808 distinct sequence pairs. No target exhausts its compatible support.

The uncapped Direct baseline uses a diversity-oriented selection objective and is not an exhaustive optimal static search. Its mean retained counts and coverage describe that specified objective. Screening uses the same full references and identity as the main study. All 3,000 uncapped proposals were screened; hash and deadline checks pass. Gene-paired resampling uses seed 42 and 10,000 replicates. Confidence intervals are descriptive and unadjusted. Retained diversity compares only jointly defined genes; training-seed robustness is not measured. Complete machine-readable tables and contrasts accompany the code release.
\input{mechanism_table.tex}

\section{AI assistance and remaining author declarations}
OpenAI Codex assisted with code implementation, analysis, figure preparation, documentation and manuscript drafting and revision. Historical assistance records also document model-assisted method-figure preparation. This disclosure covers substantive drafting as well as editing. Jingyue Li is first author and Shiyue Mei is corresponding author. The authors report no specific funding and no conflict of interest.

''' +supp[end:]
    supp=supp.replace(r'\path{paper/bioinformatics_cnn_draft}',r'\path{paper/bioinformatics_submission}')
    supp=supp.replace(r'\path{scripts/build_cnn_manuscript_materials.py}',r'\path{scripts/prepare_bioinformatics_submission.py}')
    supp=supp.replace(r'.venv/bin/python',r'python')
    pos=supp.index(r'\end{document}')
    supp=supp[:pos]+r'''\begin{figure*}[p]\centering\includegraphics[width=\textwidth]{figures/mechanism_controls.pdf}
\caption{Added fixed-seed mechanism controls. Each cohort includes all 75 genes. Bars show retained-pair means with marginal 95\% gene-bootstrap intervals. Strong static decoding is selected on validation only. Frozen scores reuse one network call. The uncapped Direct search uses the stored per-gene sixteen-step CPU budget. Paired contrasts are reported in the main text and machine-readable evidence.}\label{fig:mechanismcontrols}
\par\noindent\textit{Alt text:} One-step and frozen-score models retain most count performance; full updates add diversity, and deadline Direct search returns fewer pairs while covering more genes.\end{figure*}
''' +supp[pos:]
    supp=supp.replace(r'\input{draft_formatting.tex}',r'\input{submission_formatting.tex}').replace(r'\copyrightyear{2026}\pubyear{2026}',r'\copyrightyear{}\pubyear{}')
    (OUT/'supplementary.tex').write_text(supp)
    sources=[DRAFT/'manuscript.tex',DRAFT/'supplementary.tex',DRAFT/'source_manifest.json',MECHANISM,
             MECHANISM.parent/'verification.json',ROOT/'experiments/mechanism/tuned_seed42/selection.json']
    write_json(OUT/'source_manifest.json',{'sources':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in sources],
        'authors':[{'name':'Jingyue Li','affiliations':[1,2,3],'email':'lijingyue513@163.com','first_author':True},
                   {'name':'Shiyue Mei','affiliations':[2,3],'email':'18837614038@163.com','corresponding_author':True}],
        'funding':'No specific funding','conflicts':'None declared by user','license':'MIT for original code',
        'generation_seed':42,'historical_checkpoint_training_seed':20261013,
        'status':'author_review_candidate_not_yet_submission_ready',
        'unresolved':['researcher writing and LLM-policy resolution','detailed affiliation addresses','submitting-author ORCID','independent code/test-data archival DOI','GitHub publication verification']} )
    print(OUT)


if __name__=='__main__':main()
