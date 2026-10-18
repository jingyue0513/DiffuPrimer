# Contemporary primer-design work and comparison plan

Search date: 13 October 2026. Sources below are original papers and author repositories. Publication dates distinguish online publication from issue year. This document records related work; benchmark results live under `experiments/attention/results/`.

## Recommended comparisons

1. **qPrimer / qPrimerDB 2.0** is the strongest readily runnable contemporary same-task tool. It designs gene-specific qPCR primer pairs from transcripts with thermodynamic constraints and optional reference/exon/SNP assessment. The paper appeared online **9 August 2024**, in the **2025 Nucleic Acids Research database issue**, DOI **10.1093/nar/gkae684**. Its Python CLI is open source and provides separate design, check, annotate and visualize modules. The native design engine is Primer3; consequently, label it as a **qPCR tool pipeline**, rather than as an independent learned generator. [Paper](https://academic.oup.com/nar/article/53/D1/D205/7730531), [author repository](https://github.com/swu1019lab/qPrimer).

2. **Primer PICKR**, *Nature Communications*, online **28 May 2026**, DOI **10.1038/s41467-026-73648-2**, retrieves literature-derived RT–qPCR primers and ranks recombined pairs using publication evidence, biophysics and pair synergy. Ranked mouse records can provide a **retrieval baseline** on the same frozen genes, and its released raw literature primers provide an independent source for evaluation. Coverage depends on literature availability; evaluate all frozen genes, plus a clearly identified common-covered subset for quality. No test-gene PICKR record should enter attention-model training. [Paper](https://www.nature.com/articles/s41467-026-73648-2), [code](https://github.com/tmolley2/PICKR), [database](https://github.com/tmolley2/PICKR_DB).

3. **PrimerCast** provides learned amplification evaluation and a current ML-guided qPCR design pipeline. A public manuscript titled *PrimerCast: Predictive Modeling of PCR Amplification with an AI-Ready Experimental Dataset* describes a NeurIPS **2025 submission**. Its final conference/journal acceptance is not established by the retrieved sources; do not label it an accepted 2026 paper without verification. The current Broad Institute repository exposes design/evaluate commands and packaged pretrained classifiers/regressors. It is useful as a **learned design/ranking comparison**, after confirming that its viral-trained predictors and transcript input interface can run on the frozen mouse references. A predictor score is an additional model-based endpoint, not a measured amplification outcome. [Public manuscript](https://openreview.net/pdf/946f5e08be95e3d6ab9d72c2d1e493f7bdc1a74a.pdf), [author repository](https://github.com/broadinstitute/PrimerCast).

## Additional related work

| Method | Verified date/status | What it does | Role in this paper |
|---|---|---|---|
| Primer C-VAE | arXiv v1 3 March 2025; v2 16 December 2025 | CNN/VAE classification and sequence reconstruction to identify discriminative forward primers, followed by reverse-primer design | Closest contemporary generative-primer citation; its variant-population input differs from single-transcript conditional pair generation |
| varVAMP | Nature Communications 2025, DOI 10.1038/s41467-025-60175-9 | Degenerate qPCR and tiled primers from diverse viral multiple-sequence alignments | Strong related tool; compare in a separate multi-variant extension if that task is added |
| PrimeGen | Online 30 July 2025; Nature Biomedical Engineering issue February 2026, DOI 10.1038/s41551-025-01455-z | GPT-4o-controlled multi-agent workflow for targeted sequencing panel design and laboratory automation | Workflow-level comparison in Related Work; not an interchangeable primer-pair generator on the present benchmark |
| PRISM | Bioinformatics 42(7), 2026, DOI 10.1093/bioinformatics/btag478 | Constrained submodular multiplex primer-set selection balancing coverage and dimer risk | Especially relevant to target-journal context; separate multiplex-set objective from per-target alternative-pair design |
| PrimerAST | Scientific Reports, online 25 March 2026, DOI 10.1038/s41598-026-38238-8 | Primer3 generation and supervised scoring from 16 engineered primer features | Related learned-quality assessment; evaluate only if its actual scoring model/interface can be reproduced |
| BioInnovate AI | Diagnostics, 6 June 2025, DOI 10.3390/diagnostics15121445 | Primer/probe assay development supported by random-forest and gradient-boosting amplification prediction | Contemporary learned assay-quality reference; its data availability statement is by author request |

Primary sources: [Primer C-VAE paper](https://arxiv.org/abs/2503.01459), [Primer C-VAE repository](https://github.com/harrywang9917/Primer_C-VAE); [varVAMP paper](https://www.nature.com/articles/s41467-025-60175-9), [varVAMP repository](https://github.com/jonas-fuchs/varVAMP); [PrimeGen paper](https://www.nature.com/articles/s41551-025-01455-z); [PRISM paper](https://academic.oup.com/bioinformatics/article/42/7/btag478/8722296), [PRISM repository](https://github.com/yhhan19/PRISM-new); [PrimerAST paper](https://www.nature.com/articles/s41598-026-38238-8).

Additional primary sources: [BioInnovate AI full article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12191490/). Earlier learned specificity work includes [eDNAssay, Molecular Ecology Resources 2022](https://pubmed.ncbi.nlm.nih.gov/35778862/), with [author-released scripts, trained models and empirical SYBR/TaqMan data](https://github.com/NationalGenomicsCenter/eDNAssay). eDNAssay predicts cross-amplification rather than generating paired binding positions.

The search found multiple contemporary learned, optimization and retrieval approaches. It did not establish a prior method performing the precise combination of absorbing-mask diffusion over paired binding-position labels, shared conditional/anchor inference, and full-reference transcript/genome assessment. Phrase the contribution around this concrete formulation; do not claim to be the first machine-learning or generative primer-design method.

## Runnable qPrimer adapter

`experiments/attention/baselines/run_qprimer.py` calls the **actual upstream** `qPrimer/design.py::design_primers`, pinned at commit `34f047e5d77ec0353846a9a18b0b94f74013db3b` (package 1.0.6). It uses the shared primer/product length, Tm, GC, structure-Tm, homopolymer and salt-concentration bounds and K=20. It retains qPrimer's additional alignment-score/end-stability rules, native optimal size and GC weighting. Internal probe generation is disabled to match the common forward/reverse-pair task. That adjustment is recorded in each manifest.

```bash
git clone https://github.com/swu1019lab/qPrimer.git experiments/attention/baselines/vendor/qPrimer
git -C experiments/attention/baselines/vendor/qPrimer checkout 34f047e5d77ec0353846a9a18b0b94f74013db3b
.venv/bin/python experiments/attention/baselines/run_qprimer.py
.venv/bin/python experiments/attention/baselines/run_qprimer.py --amplicon-min 80 --amplicon-max 120
.venv/bin/python experiments/attention/baselines/run_qprimer.py --task completion --fixed-side forward
.venv/bin/python experiments/attention/baselines/run_qprimer.py --task redesign --fixed-side reverse
```

The output has the same candidates/per-gene/completion schema as the existing common Bowtie assessment. It does **not** replace common specificity checks with upstream-specific checks. The broad band, three narrow bands, two-sided completion and two-sided redesign provide eight compatible settings. Training data and neural checkpoints are untouched.

## PrimerCast runnable interfaces

The author repository documents:

```bash
conda env create -f environment.yml
conda activate primercast
pip install .
adapt design --params params.txt
adapt evaluate --for ATCGATCGATCG --rev GCTAGCTAGCTA
adapt evaluate --pset my_primers.fa
```

Design inputs are local unaligned FASTA files under `target_seqs/original/`, with target and offtarget files; `workflows/params.txt.template` controls primer length, Tm/GC, product lengths and other settings. Generated output is `{target}_final.csv` with candidates/scores/coverage. Pretrained weights are bundled under `src/primercast/data/`; training scripts exist under `training/`, while the raw experimental dataset is stated as available on request. No fetch/monitor/email workflow is needed for local evaluation. [Repository usage](https://github.com/broadinstitute/PrimerCast), [parameter template](https://raw.githubusercontent.com/broadinstitute/PrimerCast/main/workflows/params.txt.template), [internal CLI documentation](https://raw.githubusercontent.com/broadinstitute/PrimerCast/main/docs/primercast_cli.md).

An initial shallow clone failed with a GitHub TLS transport interruption. The actual source and bundled model files have now been downloaded from author-hosted raw URLs at pinned commit `50d357b43a727196b5ae2e1c2e7d041519243673` to `experiments/attention/baselines/vendor/PrimerCast/`. `source_snapshot.json` records 42 files (875,085 bytes). Missing dependencies were installed in a separate `.venv_primercast`; the pretrained classifier and regressor load successfully. Official RNAduplex 2.4.17 and Bowtie2 2.4.4 were installed, and sam2pairwise v1.0.0 was built from the author's source. `run_primercast.py` now executes unchanged upstream quick-design on the frozen transcripts and preserves native logs, alignments, feature tables and neural scores.

The author environment recommends RNAduplex >=2.6 and Bowtie2 >=2.5, whereas the runnable distribution tools here are 2.4.17 and 2.4.4. Exact observed versions are recorded in manifests. No placeholder energies or artificial neural scores are used.

The public quick-design interface supports length/Tm/GC-upper-bound controls; unsupported common bounds receive shared final assessment. It has no fixed-primer argument, so anchored completion/redesign are not attributed to it. Native output permits imperfect on-target matches: five audited incompatible exact-template examples showed 4–7 reverse mismatches or a forward insertion, yielding native 129–188 bp predictions while exact binding positions formed 407–2938 bp products. These are genuine native proposals rather than adapter orientation errors. `adapter_failure_audit.json` preserves actual sequences and native evidence. RefSeq accession FASTA headers avoid pandas converting all-numeric gene identifiers to integers; the early numeric-header trial is archived and excluded.

## Acquired PICKR source data

The database README advertises per-species scored CSV directories, but the actual Git tree at commit `f72760013722c2d698c9c7f3b80ac1640ea8b4c7` currently contains only README, license, citation metadata and `raw_primers_3.csv.gz`. The raw file was successfully collected under `data/raw/pickr/` (59,476,480 bytes; SHA-256 `89ba90da80c92de47e42e5b6a49bff3b4ec70a9dbae8b1da9dcecc409e3f8c9d`). Its fields are `Sequence`, `PMCID`, `Orientation`, `Gene`, `Citation Count`; species is absent and some gene/orientation fields are unknown. This is a real independent literature sequence resource that requires exact template mapping. It is not the platform's pre-ranked mouse pair table. A PICKR retrieval comparison should use actual scored records from the platform or a verified release, rather than relabeling reconstruction from raw strings as the published ranker.

`map_pickr_literature.py` has processed 5,852,944 original rows. Of 3,042,411 ACGT rows with lengths 18–25, there are 2,147,155 unique sequences, of which 97,378 appear at least three times. Exactly 490 source sequences match the 500 frozen transcripts; all matched intervals retain PMCID provenance. Among the 75 test genes, 35 have matched source oligos, 23 have source oligos in both orientations, and 18 have at least one oligo satisfying the individual common thermodynamic constraints. Outputs are under `experiments/attention/results/pickr_literature/`. This source enters no training and supplies no fabricated pair-level outcome labels.

## Concrete experimental protocol

Use the existing similarity-clustered 350/75/75 split, seed 20261013, full reference transcript templates, original RNA/genome assemblies, common 3′-exact mismatch rules and K=20. Keep 75 targets in count/coverage endpoints. Report RNA-specific pairs, RNA-specific/no-genomic-product pairs, retained fraction, two-or-more coverage, interval diversity, best shared Tm penalty and runtime. Use the same frozen anchors and changed-length redesign tasks. Preserve both tool-specific proposal settings and common final assessment.

Main comparison: attention diffusion, matched direct attention, original CNN diffusion, Primer3, qPrimer, rule sampler and uniform compatible-pair sampling. PICKR can add a separate retrieval panel; PrimerCast can add a learned-pipeline panel once its local execution is verified. The attention experiment needs the direct attention comparator so architecture gain and iterative-diffusion gain can be distinguished. Prefer these concrete tests to adding task-incompatible multiplex tools to the main table.

## Suggested Introduction wording

“Recent primer-design systems address complementary aspects of assay development. qPrimer supports gene-specific qPCR design with thermodynamic and reference-aware analysis, whereas Primer PICKR prioritizes literature-derived primer pairs. Primer C-VAE explores generative learning for variant-discriminative primers, and PrimerCast learns primer–target amplification behavior. Multi-variant and multiplex systems, including varVAMP, PrimeGen and PRISM, optimize conservation, panel coverage or cross-primer interactions. We formulate transcript-level alternative-pair design as conditional generation of coupled binding intervals, and use attention-based discrete diffusion to connect de novo design, length control, anchored completion and local redesign within one inference procedure.”

Use the final sentence only after the attention implementation and corresponding experiments are completed.
