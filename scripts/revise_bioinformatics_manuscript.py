#!/usr/bin/env python3
"""Revise the submission prose around scientific methods and verified results."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'paper/bioinformatics_submission'

def main():
    main=P/'manuscript.tex';text=main.read_text()
    abstract=r'''\abstract{\textbf{Motivation:} qPCR primer design often requires alternative pairs or replacement of one primer as design requirements change. These tasks require compatible binding positions, suitable thermodynamic properties and target specificity.\\
\textbf{Results:} We present DiffuPrimer, a conditional discrete diffusion model of paired primer-binding positions on a fixed transcript. A constrained decoder supports alternative generation, product-length control, fixed-primer completion and local redesign. With a budget of 20 candidates per gene, DiffuPrimer yielded \MainCNNRet{} and \ExtraCNNRet{} RNA-specific pairs per gene without predicted genomic products on two disjoint 75-gene mouse test sets, compared with \MainPrimerRet{} and \ExtraPrimerRet{} for Primer3. Frozen-score controls preserved most candidate-count gains, while iterative network updates increased positional diversity.\\
\textbf{Availability and implementation:} Code and data are available at \url{https://github.com/jingyue0513/DiffuPrimer}, with an \href{https://archive.softwareheritage.org/swh:1:rev:275d1c0acf64edba234d17c1a83094af597f6202/}{archived release}.}'''
    a=text.index(r'\abstract{');b=text.index(r'\keywords',a);text=text[:a]+abstract+'\n'+text[b:]
    data_note='Supplementary data provide machine-readable cohort summaries, mechanism-control results and figure data.'
    if data_note not in text:
        text=text.replace(r'\section*{Conflict of interest}',data_note+'\n\n'+r'\section*{Conflict of interest}')
    disclosure='OpenAI Codex assisted with code and analyses.'
    a=text.index(r'\section*{Acknowledgements}') if r'\section*{Acknowledgements}' in text else text.index(r'\section*{Conflict of interest}');b=text.index(r'{\small\setlength{\bibsep}',a)
    text=text[:a]+r'''\section*{Conflict of interest}
The authors declare no conflict of interest.

'''+text[b:]
    anchor=r'\subsection{Reference-screened candidate quality}'
    text=text.replace(disclosure+'\n\n','')
    text=text.replace(anchor,disclosure+'\n\n'+anchor,1)
    text=re.sub(r'\\par\\noindent\\textit\{Alt text:\}[^\n]*?(?=\\end\{figure\*\})','',text)
    main.write_text(text)
    supp=(P/'supplementary.tex').read_text();head=supp[:supp.index(r'\section{')]
    head=head.replace(r'\lastpage{8}',r'\lastpage{10}')
    head=re.sub(r'\\abstract\{[^}]*\}',r'\\abstract{}',head)
    if r'\renewcommand{\thesection}' not in head:
        head=head.replace(r'\begin{document}',r'\renewcommand{\thesection}{S\arabic{section}}'+'\n'+r'\renewcommand{\thetable}{S\arabic{table}}'+'\n'+r'\begin{document}')
    figures=supp[supp.index(r'\clearpage\n'):] if r'\clearpage\n' in supp else supp[supp.index('\\clearpage\n'):]
    methods=r'''\section{Reference sequences, supervision and gene partitions}
\subsection{Reference and template selection}
All experiments use mouse GRCm39, RefSeq assembly GCF\_000001635.27 and annotation release RS\_2024\_02. The complete RNA reference, genome and GFF annotation provide target sequences, specificity references and exon coordinates. The initial cohort comprises 500 protein-coding genes with annotated NM transcripts, sampled using seed 20261013. Each gene is represented at inference by its lexicographically first annotated NM accession. Template selection precedes primer-source outcomes. Sequence reconstruction from genomic exons verifies the transcript sequences, with terminal poly(A) handled separately.

\subsection{Positional supervision}
The acquired PrimerBank snapshot contains 311 pairs from 339 queried genes. Exact remapping finds compatible binding positions for 289 pairs; 22 have no compatible remap. Removing four mapped duplicates leaves 285 public-source positional records. Acquisition stopped when the source denied further access. Consequently, these records represent a partial acquisition rather than complete coverage of PrimerBank.

Primer3 supplies an additional 7,410 weak-supervision pairs on training and validation genes. For each transcript, up to five evenly spaced 512-base windows are used; the complete transcript is also queried when longer than 512 bases. Each region requests up to four Primer3 pairs. Coordinates are translated back to the transcript, duplicate synthesis-sequence pairs are removed, and common thermodynamic criteria are applied. The combined 7,695 records comprise 6,260 training records, 1,378 validation records and 57 historical test-gene records. The 57 test-gene records are excluded from fitting. The training labels specify primer binding positions; computational designs do not provide amplification-efficiency measurements.

\subsection{Homology control and additional test genes}
All annotated NM isoforms are searched using MMseqs2 in nucleotide mode. Two genes are connected if any isoform alignment has at least 80\% identity and covers at least 80\% of both sequences. Gene-level connected components define a 350/75/75 training/validation/test partition. No edge crosses partitions at that threshold. Sensitivity searches at 70\% identity/70\% bidirectional coverage and 90\%/90\% find no training--test edge; the lower threshold identifies one training--validation edge.

A second test cohort comprises 75 genes and 169 NM isoforms from the same reference. These genes are disjoint from the original 500 and have no qualifying 80\%/80\% homology edge to the 781 training isoforms. Their selection precedes candidate generation. The original and additional test cohorts are analysed separately.

\section{Model architecture and training}
\subsection{Inputs and position labels}
DNA bases use five input codes for A, C, G, T and N. Position labels distinguish background, forward binding, reverse binding and an absorbing mask. A clean pair occupies two non-overlapping intervals ordered along the forward transcript. The synthesized reverse primer is the reverse complement of its interval. Padding is masked in the network and loss; oligos containing ambiguous bases are ineligible for generation.

The nine condition values are primer minimum/maximum length, product minimum/maximum length, Tm minimum/maximum, maximum pair Tm difference and GC minimum/maximum. Their normalization divisors are 25, 25, 350, 350, 65, 65, 5, 100 and 100, respectively. Four time features encode normalized time, squared time, and its sine and cosine. The condition/time vector is projected through two linear layers with SiLU activation.

\subsection{CNN and architecture controls}
The CNN has width 96 and four residual blocks. Each block uses eight-group normalization, SiLU, a width-five convolution with dilation 1, 4, 16 or 64, SiLU and a pointwise convolution. Relative position and its sinusoid are projected and added to the summed DNA and label embeddings. A normalized pointwise head predicts three clean-label logits per base.

The pure-attention comparator uses width 96, four attention blocks, four heads and a feed-forward expansion factor of two. Rotary embeddings encode relative position; block-wise affine modulation incorporates the condition/time embedding. The hybrid comparator adds two such attention blocks after the four CNN blocks. Direct and diffusion versions share the corresponding architecture. The CNN is reported following the architecture comparison; all completed architecture results are included in Table~\ref{tab:allmethods}.

\subsection{Windows, augmentation and optimization}
Each eligible gene contributes four conditional windows per epoch. A labelled pair is sampled within the gene and placed inside a random 512-base transcript window. Validation windows and pair choices are deterministic. Eligible training pairs have 18--25 nt oligos and products no longer than 400 bp. Training product bands are sampled from 60--120, 80--200, 150--300 and 200--400 bp when they contain the labelled pair; otherwise, a band extending up to 25 bp around its length is used. One labelled primer interval remains observed in 25\% of augmented examples, with forward or reverse chosen at random.

Training uses AdamW for 30 epochs, learning rate $10^{-3}$, weight decay 0.01 and batch size 16. Diffusion uses 16 absorbing-mask levels. At level $t$, each non-anchored valid label is independently masked with probability $t/16$. Cross-entropy on masked sites is normalized by the valid sequence length and weighted by $16/t$. Direct prediction reconstructs the clean map from fully masked unknown positions with the same eligible data and anchor augmentation. Checkpoints minimize deterministic validation reconstruction loss. Historical training and the main benchmark use seed 20261013. The mechanism controls reuse one existing Direct/diffusion CNN checkpoint pair and use generation seed 42.

OpenAI Codex assisted with code and analyses.

\section{Thermodynamic criteria and paired-interval sampling}
\subsection{Shared requirements}
Standard design requests 18--25 nt oligos and 80--200 bp products. GC content is 35--65\%, Tm is 57--63\,$^{\circ}$C with optimum 60\,$^{\circ}$C, and pair Tm difference is at most 2\,$^{\circ}$C. Homopolymers are at most five bases. Hairpin, homodimer and heterodimer Tm are each at most 47\,$^{\circ}$C. Primer3 thermodynamic functions use 50 mM monovalent salt, 1.5 mM divalent salt, 0.6 mM dNTP and 50 nM oligo. Primer-junction and exon annotations do not enter proposal scoring.

\subsection{Compatible recovery}
For a proposed pair, the decoder sums forward-minus-background logits over the forward interval and reverse-minus-background logits over the reverse interval. Support enforces interval order, oligo/product lengths, individual-oligo eligibility, pair Tm difference, observed labels and supplied anchors. The beam contains the 64 highest-scoring eligible forward spans and up to eight reverse partners per span; the preceding compatible pair remains available. Reverse-anchored tasks filter compatible forward support before beam ranking.

At reverse step $t$, a pair is sampled with probability proportional to $\exp(S/\tau)$ and each unknown label is revealed with probability $1/t$. Previously revealed labels and anchors persist. Pair extraction is followed by heterodimer assessment and common reference screening. Main-benchmark temperature is selected on validation from 0.5, 1, 2 and 4; the CNN methods use $\tau=4$. Recovery-depth experiments use 1, 4, 8 or 16 steps with the same checkpoint and normalized-time schedule.

\section{Comparator configurations and conditional tasks}
\subsection{Conventional and learned comparators}
Primer3 receives the common oligo/product bounds and thermodynamic settings. The qPrimer design module is pinned to commit \href{https://github.com/swu1019lab/qPrimer/tree/34f047e5d77ec0353846a9a18b0b94f74013db3b}{\texttt{34f047e5d77e}}. Its internal-probe option is disabled for pair-only design; common biological bounds and concentrations are matched, while additional native scoring and prefilters retain their defaults.

PrimerCast is pinned to commit \href{https://github.com/broadinstitute/PrimerCast/tree/50d357b43a727196b5ae2e1c2e7d041519243673}{\texttt{50d357b43a72}}. Its native quick-design workflow performs primer tiling, thermodynamic filtering, Bowtie2 mapping, feature preparation, pretrained prediction and pair filtering. The evaluated environment uses RNAduplex 2.4.17, Bowtie2 2.4.4 and sam2pairwise 1.0.0. Supported length, Tm and GC settings are passed through the native interface. Native duplex free-energy and activity/coverage thresholds retain their defaults. At most 20 highest-ranked returned pairs enter the common assessment. PrimerCast provides complete-design results; anchored tasks are evaluated only for methods with the corresponding interface.

The uniform-pair control samples from compatible pairs. The forward-first control samples a forward span and then a compatible reverse partner. They share the eligible-oligo support used by the learned samplers. All comparators are evaluated on identical target templates and downstream references. Invalid outputs, duplicate pairs and tool shortages consume the proposal budget.

\subsection{Length control, completion and redesign}
Length-controlled design requests 80--120, 121--160 or 161--200 bp. Completion fixes either oligo from the broad-range Primer3 top pair and requests a partner yielding an 80--200 bp product. Redesign changes the product band cyclically: 80--120 to 121--160, 121--160 to 161--200, and 161--200 to 80--120 bp. All 75 original partners fall outside their redesigned band. A valid anchored output preserves the supplied interval exactly and satisfies the new product bounds.

\section{Transcriptome and genome assessment}
Bowtie 1.3.1 searches the complete transcriptome and genome using all permitted hit strata and both orientations. Alignments are full-length and ungapped, allow at most two substitutions, and require the five terminal 3$'$ primer bases to match exactly. Inward-facing products of 36--1,000 bp are assessed for F/R, R/F, F/F and R/R primer-role combinations. Thus genome assessment is not restricted to the requested RNA-product band.

An RNA-specific pair produces an intended-gene RNA product and no other-gene or unresolved RNA product; annotated isoforms from the same gene are accepted. Any potential genomic product supplies a genomic flag. A retained pair passes the shared computational constraints, RNA specificity and absence of a genomic flag. Counts use unique synthesized forward/reverse sequence pairs. Exact-sequence hit caching shares alignments across methods without changing the matching or product rules.

\section{Frozen-score, static-decoding and compute controls}
\subsection{Matched frozen-score recovery}
Full and frozen-score recovery use identical checkpoints, candidate budgets, decoder support, initial states, random seeds and reveal schedules. Frozen recovery evaluates the network once on the fully masked state and reuses that score map at every step. Four- and sixteen-step frozen controls therefore isolate state-dependent network updates from repeated constrained sampling.

\subsection{Full-support and diverse-pool decoding}
Static decoding enumerates every pair satisfying geometry, individual-oligo eligibility and pair Tm compatibility, then removes synthesis-sequence duplicates. Weighted sampling without replacement uses Gumbel perturbations of the temperature-scaled pair scores. The direct and one-step diffusion networks use the same decoder. Pool controls draw 320 distinct sequence pairs and select at most 20 computationally feasible pairs. The first pair minimizes the shared Tm penalty; later pairs maximize minimum amplicon Jaccard distance to the selected set. Ties use Tm penalty and then pool order.

A separate validation search compares deterministic top-$K$ with weighted sampling without replacement at temperatures 0.5, 1, 2 and 4. Each network is selected on all 75 validation genes by mean retained count, then at-least-two-pair coverage and configuration order. Direct selects temperature 2; one-step diffusion selects top-$K$. Settings are fixed before either test cohort is evaluated. RNA/genome outcomes and exon annotations never enter test-time proposal selection.

\subsection{CPU-matched direct search}
Each gene's measured sixteen-step diffusion process-CPU time supplies the Direct search budget. Timing includes network prediction, support construction, sampling, sequence deduplication, thermodynamic assessment and selection. Both methods exclude model loading and shared span preparation. Pair-property caches are cleared before each method; ordinary method order rotates across genes. Screening follows proposal selection and is excluded from generation timing.

An initial 4,096-candidate cap is reported as a bounded-pool control. The uncapped comparison keeps the same Direct scores and diversity-oriented selection while continuing to the per-gene deadline. Selection checkpoints occur every 512 pool additions, reserving the measured selection cost for a final checkpoint. Only outputs committed before the deadline are evaluated. Eight worker processes each use four intra-operation threads, and budgets use each worker's process-CPU clock.

Mean budget utilization is 99.20\% in the original cohort and 99.15\% in the additional cohort; minima are 98.71\% and 98.62\%. The largest non-preemptible operation overrun is 0.00765 s and does not change the on-time committed output. Mean explored pools contain 28,014 and 25,808 distinct sequence pairs. No target exhausts its compatible support. The final output budget remains 20 pairs per gene. These measurements compare generation CPU cost under the stated selection objective.

\section{Endpoints and statistical analysis}
For $K=20$, all 75 genes in each cohort enter count and coverage analyses, including genes with no returned or retained pair. Budget curves use recorded prefixes $K\in\{1,2,5,10,20\}$. The retained fraction is the per-gene ratio of retained to RNA-specific unique pairs and is defined when its denominator is nonzero.

For product intervals $I$ and $J$, positional distance is $1-|I\cap J|/|I\cup J|$. Diversity averages this distance over distinct interval pairs and requires at least two intervals. Historical diversity summaries use RNA-specific pairs; the mechanism study additionally reports diversity among retained pairs. Temperature matching uses the minimum shared penalty $|T_F-60|+|T_R-60|+|T_F-T_R|$ among eligible RNA-specific or retained pairs, as labelled in the corresponding analysis.

Marginal confidence intervals resample genes within each cohort. Paired intervals resample within-gene differences over jointly defined genes, using 10,000 replicates and percentile endpoints. Historical analyses use seed 20261013; mechanism analyses use seed 42. Cohorts are not pooled. Undefined conditional endpoints are omitted from that endpoint's paired comparison. Intervals are descriptive and unadjusted for multiple comparisons.

\section{Exon annotation and illustrative cases}
Exon coordinates classify RNA-specific candidates as primer-junction, amplicon-junction-only, same-exon, terminal-poly(A), or unmapped. Genome flags are obtained independently. Two qPrimer terminal-poly(A) candidates remain in RNA-specific totals and appear as Other in the aggregate figure.

The Susd5 example uses NM\_001101510.1. Its 103-bp retained CNN product has a forward primer spanning splice coordinate 588. Primer3 proposal 0 produces a 104-bp RNA product and a genomic witness on NC\_000075.7 at the zero-based half-open interval [113924803,113924907). The example is selected from 225 same-transcript pairings with a same-exon comparator genomic witness and a junction-spanning retained CNN alternative. Selection prioritizes Primer3, a combined RNA window within 512 bp, minimum product-length difference, minimum combined window, and stable identifiers.

The C4bp example uses NM\_001419982.1. The initial RNA product is 193 bp; forward-fixed and reverse-fixed completion give 188 and 90 bp, while forward-fixed and reverse-fixed redesign give 113 and 90 bp. The 113-bp redesign candidate has a genomic flag. The figure displays these outcomes alongside the observed primer positions.

\section{Supplementary results}
Tables~\ref{tab:allmethods}--\ref{tab:contrasts} provide the complete-design, conditional-task, component and paired-endpoint results. Table~\ref{tab:allmechanism} reports all 14 mechanism configurations in both cohorts. Numerical results refer to the stated seeds and candidate budgets.
\input{supplement_tables.tex}
\input{mechanism_complete_table.tex}

'''
    figures=figures.replace('Added fixed-seed mechanism controls.','Fixed-seed mechanism controls.')
    figures=re.sub(r'\\par\\noindent\\textit\{Alt text:\}[^\n]*?(?=\\end\{figure\*\})','',figures)
    (P/'supplementary.tex').write_text(head+methods+figures)
    results=json.loads((P/'mechanism_evidence.json').read_text())
    groups={(g['cohort'],g['method']):g for g in results['groups']}
    names={'direct_beam':'Direct, original beam','diffusion_1':'Diffusion, 1 step','diffusion_4':'Diffusion, 4 steps','diffusion_16':'Diffusion, 16 steps','frozen_4':'Frozen scores, 4 steps','frozen_16':'Frozen scores, 16 steps','direct_full_wor':'Direct, full-support sampling','diffusion1_full_wor':'One-step diffusion, full support','direct_pool320':'Direct, pool 320','diffusion1_pool320':'One-step diffusion, pool 320','direct_equal_cpu':'Direct, CPU allowance, cap 4,096','direct_validation_selected':'Direct, validation-selected','diffusion1_validation_selected':'One-step diffusion, validation-selected','direct_equal_cpu_uncapped':'Direct, uncapped CPU deadline'}
    lines=[r'\begin{table*}[t]\centering\small\setlength{\tabcolsep}{4pt}',r'\caption{Complete mechanism controls with seed 42 and $K=20$. Values are mean unique retained pairs, gene coverage with at least one retained pair (\%), and generation CPU seconds. Each cohort includes all 75 genes. The capped control and uncapped deadline comparison are reported separately.}\label{tab:allmechanism}',r'\begin{tabular}{lrrrrrr}\toprule',r'&\multicolumn{3}{c}{Original cohort}&\multicolumn{3}{c}{Additional cohort}\\',r'Method & Pairs & Coverage & CPU (s) & Pairs & Coverage & CPU (s)\\\midrule']
    for method,name in names.items():
        row=[]
        for cohort in ['test','external_test']:
            g=groups[cohort,method];row.append(f"{g['retained_pairs']:.2f} & {100*g['one_retained']:.1f} & {g['generation_cpu_seconds']:.2f}")
        lines.append(name+' & '+' & '.join(row)+r'\\')
    lines.extend([r'\bottomrule\end{tabular}\end{table*}'])
    (P/'mechanism_complete_table.tex').write_text('\n'.join(lines)+'\n')
    cited={k for block in re.findall(r'\\cite\w*\{([^}]+)\}',text) for k in block.split(',')}
    for name in ['references.bib','submission_references.bib']:
        path=P/name;bib=path.read_text().replace('https://iclr.cc/Downloads/2019','https://openreview.net/forum?id=Bkg6RiCqY7').replace('https://www.iclr.cc/archive/www/2016.html','https://arxiv.org/abs/1511.07122').replace('https://openreview.net/pdf?id=PEQFHRUFca','https://openreview.net/forum?id=PEQFHRUFca')
        blocks=re.findall(r'(@\w+\{([^,]+),.*?)(?=\n@|\Z)',bib,re.S)
        path.write_text('\n'.join(block.rstrip() for block,key in blocks if key in cited)+'\n')
    audit_path=ROOT/'docs/reference_verification_20261018.json'
    if audit_path.exists():
        import hashlib
        audit=json.loads(audit_path.read_text())
        fields={key:{k:v for k,v in re.findall(r'^\s*(\w+)\s*=\s*\{(.*)\},?$',block,re.M)} for key,block in re.findall(r'@\w+\{([^,]+),(.*?)(?=\n@|\Z)',(P/'references.bib').read_text(),re.S)}
        for row in audit['records']:
            assert fields[row['key']]['title']==row['title']
            assert fields[row['key']].get('year')==row['year']
            assert fields[row['key']].get('doi')==row.get('doi')
        audit['main_bibliography_sha256']=hashlib.sha256((P/'references.bib').read_bytes()).hexdigest()
        (P/'reference_verification.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
    import subprocess,sys
    subprocess.run([sys.executable,str(ROOT/'scripts/plot_submission_mechanism.py')],check=True)
    print('Revised abstract, declarations and supplementary experimental methods')

if __name__=='__main__':main()
