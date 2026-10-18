#!/usr/bin/env python3
"""Insert recorded final values and compile the manuscript and supplement."""
import json
import re
import subprocess
from pathlib import Path
from diffuprimer.io import sha256, utc_now, write_json

root = Path(__file__).resolve().parents[1]
paper = root / "paper"
state = json.loads((root / "experiments/final_followup/results/status.json").read_text())
if not state.get("anchor_prefilter_tasks_complete"):
    raise RuntimeError("The final experiment barrier is not complete")
summary_path = root / "results/final_fixed_seed/study_summary.json"
data = json.loads(summary_path.read_text())


def group(method="diffusion", context="amp80-200", variant="full"):
    return next(r for r in data["summaries"] if r["method"] == method and r["context"] == context and r["variant"] == variant and r["budget"] == 20)


def contrast(endpoint, right="primer3_full", context="amp80-200"):
    return next(r for r in data["paired_comparisons"] if r["context"] == context and r["right"] == right and r["endpoint"] == endpoint and r["budget"] == 20)


d, p, direct, uniform = [group(m) for m in ["diffusion", "primer3", "direct", "uniform_pair_sampler"]]
count = "unique_RNA_pass_pairs_without_genomic_flag"
rna = "unique_RNA_pass_pairs"
quality = contrast("RNA_no_genome_fraction")
c = contrast(count)
values = {"DCount": d[count], "PCount": p[count], "CountGain": (d[count]/p[count]-1)*100,
    "CountDelta": c["mean_difference"], "CountLo": c["ci95"][0], "CountHi": c["ci95"][1],
    "QualityGap": quality["mean_difference"]*100, "QualityLo": quality["ci95"][0]*100, "QualityHi": quality["ci95"][1]*100,
    "DRNA": d[rna], "PRNA": p[rna], "DirectRNA": direct[rna], "DirectQualityCount": direct[count], "UniformQualityCount": uniform[count],
    "DirectCountGap": contrast(rna, "direct_full")["mean_difference"],
    "DiversityGap": contrast("RNA_pass_amplicon_distance", "direct_full")["mean_difference"],
    "ThermalGap": -contrast("best_RNA_pass_temperature_penalty", "direct_full")["mean_difference"]}
gains = []
for letter, context in zip("ABC", ["amp80-120", "amp121-160", "amp161-200"]):
    value = group("diffusion", context)[count]
    gain = (value/group("primer3", context)[count]-1)*100
    values["BandCount"+letter], values["BandGain"+letter] = value, gain
    gains.append(gain)
values.update(BandGainMin=min(gains), BandGainMax=max(gains))
for context, names in [("amp80-200_fix_forward", ("FCompletion", "FDirect")),
                       ("amp80-200_fix_reverse", ("RCompletion", "RDirect")),
                       ("redesign_forward", ("FRedesign", "FRedDirect")), ("redesign_reverse", ("RRedesign", "RRedDirect"))]:
    for method, name in zip(["diffusion", "direct"], names):
        values[name] = group(method, context)[rna]
decimals = {"CountGain": 1, "BandGainA": 1, "BandGainB": 1, "BandGainC": 1, "BandGainMin": 1, "BandGainMax": 1, "DiversityGap": 3, "ThermalGap": 3}
(paper / "final_values.tex").write_text("\n".join("\\newcommand{\\"+name+"}{"+f"{value:.{decimals.get(name,2)}f}"+"}" for name, value in values.items())+"\n")
write_json(paper / "reported_values.json", {"at": utc_now(), "seed": 20261013, "summary_sha256": sha256(summary_path), "values": values})
supplement = r"""\documentclass[11pt]{article}
\usepackage[margin=0.9in]{geometry}
\usepackage{graphicx,hyperref,amsmath,placeins}
\input{glyphtounicode}\pdfgentounicode=1
\emergencystretch=2em
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\title{Supplementary material: Conditional Discrete Diffusion for Diverse, Constraint-Aware qPCR Primer-Pair Design}
\author{}\date{13 October 2026}\begin{document}\maketitle
\section{Experimental records and statistical files}
The final fixed-seed analysis comprises 47 distinct target-level runs, each assessing 75 test genes under a 20-proposal budget. It includes five methods across broad-band design, three length conditions, two completion directions, two changed-condition redesign directions, component variants, and three reduced denoising depths. Existing candidate files are reused for completed matching tasks. Eighteen new evaluation runs were conducted, including the four runs implementing corrected reverse-anchor support. The earlier versions of affected outputs are retained as an algorithm audit and do not enter the final matched comparison.

The file \texttt{results/final\_fixed\_seed/study\_summary.groups.csv} contains 235 groups by method, task, variant and budget. The accompanying \texttt{study\_summary.contrasts.csv} contains all 2,599 computed paired endpoint comparisons, including their gene denominators and confidence intervals. JSON supplies identical group information and source hashes; per-gene JSONL records preserve individual targets and no-design outcomes. These files include all five methods for the anchored tasks and all defined endpoints, beyond the matched neural panels in the main flexibility figure. Candidate lists retain duplicates and invalid proposals rather than replacing them.

\section{Public data construction and source evidence}
Reference identifiers, raw-file URLs, compressed sizes, MD5 values, and SHA-256 values are recorded in \texttt{data/processed/reference\_manifest.json}. Frozen targets and gene partitions are \texttt{targets\_pilot.jsonl} and \texttt{data/splits/pilot.jsonl}. The positional dataset manifest gives 7,410 computational candidates and 285 unverified mapped source pairs, with split counts of 6,260 training, 1,378 validation, and 57 excluded test-source records. Fifteen previously published source-evidence examples are retained for evidence classification separately from model fitting and evaluation.

The 80\% identity / 80\% bidirectional-coverage partition has seven cross-gene edges and 493 connected components. The sensitivity audit retains the original partition. At 70\%/70\%, it finds one training--validation edge and no training--test edge; at 90\%/90\%, it finds no cross-partition edge. The search uses nucleotide alignments and the declared coverage rules. Spliced-reference verification reconstructs 500 templates with zero differing bases, including explicit terminal-poly(A) handling.

\section{Training histories and validation selection}
The main analysis uses the original master seed 20261013 and checkpoints selected within the 30-epoch optimization schedule. Existing records also contain seeds 20261014 and 20261015 for both neural models and continuations to 100 epochs. These continuations reproduce the original first 30 epochs exactly. Their objectives are compared within each model, because direct reconstruction and weighted denoising have different scales.
\begin{figure}[htbp]\centering\includegraphics[width=\linewidth]{../results/figures/training/training_curves.pdf}\caption{Recorded training and validation trajectories, including the additional existing seeds and 100-epoch continuations.}\end{figure}

Temperature selection used validation RNA-specific candidate counts among 0.5, 1, 2, and 4, with temperature penalty and then lower temperature as tie-breakers. Both neural models selected 4. The stored selection identifies validation sources and records that test outcomes were not used for temperature selection. The 100-epoch duration comparison holds the earlier single-oligo guidance and temperature 1 fixed; it is separate from the main joint-guidance temperature-4 study. Longer training improves best-temperature matching, gives no clear candidate-count increase, and reduces test-set positional diversity under that duration-comparison recipe.
\begin{figure}[htbp]\centering\includegraphics[width=\linewidth]{../results/figures/selected_study/training_duration_endpoints.pdf}\caption{Existing duration experiment under its fixed recipe, with validation and test endpoints shown separately.}\end{figure}

\section{Additional-seed primary evidence}
The original three-seed aggregate is preserved separately from the fixed-seed main results. It gives 6.60 retained RNA-specific pairs per gene for diffusion and 4.55 for Primer3, a 45.1\% increase. Its paired count difference is 2.05 (95\% CI 1.02--3.08). Diffusion and direct prediction give 17.53 and 16.33 RNA-specific pairs per gene, respectively. These historical unanchored comparisons use the same full-reference rules and are unaffected by the reverse-anchor sampling correction. Complete summaries remain in \texttt{results/temperature\_tuning/selected\_test\_summary.json}.

\section{Anchored-support correction and compliance}
With a reverse primer fixed at $[s_R,e_R)$, eligible forward spans satisfy $e_F\leq s_R$ and $a_{\min}\leq e_R-s_F\leq a_{\max}$; joint Tm guidance also constrains their temperature difference. These filters are applied before retaining the highest-scoring forward spans. Original reverse-anchored redesign truncated the forward beam before compatibility filtering and returned structural-validity rates of 49.33\% for diffusion and 56.00\% for direct prediction. The corrected support-first decoder gives 100\% for both across all 75 tasks. It uses the same beam width, temperature, seed, and weights. Valid anchored outputs preserve exact binding coordinates; candidate assessment additionally checks sequences and the changed product range.

The original evaluator is preserved in \texttt{experiments/final\_followup/script\_history/original\_evaluate.py}. Corrected decoder code and its protocol are \texttt{anchored\_decode.py} and \texttt{anchor\_prefilter\_protocol.json}. The correction affects only fixed-reverse neural sampling; the unanchored and fixed-forward paths retain their original behavior.

\section{Runtime and complete-reference assessment}
The measured runtime benchmark comprises 75 full transcripts, five methods, and three timing repetitions, for 1,125 measurements. Reconstructed coordinates match the stored master-seed outputs. Mean total times including shared constraint preparation and candidate assessment are 7.81 s for diffusion, 4.69 s for direct prediction, 5.02 s for uniform-pair sampling, 4.34 s for forward-first sampling, and 0.064 s for Primer3. Model loading, disk I/O, and database searches are excluded from these generation timings. Measurements use eight CPU threads; per-target timing and environment records accompany the summary.
\begin{figure}[htbp]\centering\includegraphics[width=\linewidth]{../results/figures/selected_study/generation_runtime.pdf}\caption{Measured unanchored generation times for the existing fixed-seed candidates. Timing repetitions are repeated measurements, not additional training seeds.}\end{figure}

The reference assessment uses Bowtie's all-alignment mismatch search and tests original-primer 3$'$ terminal matching after parsing alignments. Cached alignment hits are indexed by oligo sequence and reference; cached product witnesses retain the existence decisions for all four primer-role combinations. Computationally infeasible candidates receive an explicit unassessed status. Full candidate and assessment files accompany each evaluation directory. Existing aligner/parser and product-pairing checks compare actual outputs with brute-force reference constructions and coordinate cases.

\section{Reproduction commands}
The project uses its pinned Python environment and the supplied source tree. Existing target selection, labels, checkpoints, and completed evaluations are retained. The additional experiments and final analysis can be reproduced with:
\begin{verbatim}
PYTHONPATH=src .venv/bin/python experiments/final_followup/run.py
PYTHONPATH=src .venv/bin/python experiments/final_followup/run_anchor_prefilter.py
PYTHONPATH=src .venv/bin/python scripts/analyze_fixed_seed_study.py
PYTHONPATH=src .venv/bin/python scripts/plot_fixed_seed_study.py
PYTHONPATH=src .venv/bin/python scripts/finalize_paper.py
\end{verbatim}
Analysis inputs and figure provenance identify exact candidate, assessment, and summary files. All final target-count endpoints retain the complete frozen test cohort.\FloatBarrier\end{document}
"""
supplement = re.sub(r"\\texttt\{([^{}]+)\}", lambda m: "\\path{" + m.group(1).replace("\\_", "_") + "}", supplement)
(paper / "supplementary.tex").write_text(supplement)
for name in ["manuscript", "supplementary"]:
    for pass_number in range(2):
        with (paper / f"{name}.build.{pass_number}.log").open("w") as log:
            subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", name+".tex"], cwd=paper,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
write_json(paper / "build_manifest.json", {"built_at": utc_now(), "status": "formal_computational_manuscript_complete",
    "seed": 20261013, "physical_qpcr_required": False, "summary_sha256": sha256(summary_path),
    "source_sha256": sha256(paper / "manuscript.tex"), "pdf_sha256": sha256(paper / "manuscript.pdf"),
    "supplementary_source_sha256": sha256(paper / "supplementary.tex"), "supplementary_pdf_sha256": sha256(paper / "supplementary.pdf"),
    "reported_values_sha256": sha256(paper / "reported_values.json"), "figures": json.loads((root / "results/final_fixed_seed/figures/provenance.json").read_text())["figures"]})
print(paper / "manuscript.pdf")
print(paper / "supplementary.pdf")
