# Attention revision: implementation and recorded execution

Author: Jingyue Li, Sun Yat-sen University. Correspondence: 18837614038@163.com. Repository and data-archive URLs are deliberately left blank at the author's request.

All models use training/generation seed 20261013. Supervision and the 350/75/75 gene split are reused without adding test-gene teacher pairs. The model generates binding-position labels on fixed reference sequences. Masking and recovery apply to labels, not nucleotides.

| Backbone | Sequence encoding and coupling | Parameters | Training |
|---|---|---:|---|
| Original CNN | 4 residual dilated Conv1d blocks, 96 channels, dilation 1/4/16/64 | 235,011 | Existing 30-epoch checkpoints |
| Pure attention | 4 conditional Transformer blocks, width 96, 4 heads, rotary positions, feed-forward multiplier 2 | 385,635 | 30 epochs, diffusion and direct objectives |
| Local encoding + attention | 4 local residual Conv1d blocks, relative-position features, then 2 conditional Transformer blocks with rotary positions | 421,827 | 30 epochs, diffusion and direct objectives |

All runs use 512-bp training windows, batch size 16, AdamW learning rate 0.001 and weight decay 0.01. Known primer labels are retained in 25% of training examples. Checkpoints minimize the respective validation loss. The diffusion objective uses absorbing masks and masked cross-entropy weighted by T/t, with T=16. Direct prediction uses identical backbone and conditioning but predicts from masked positions in one pass.

The pure attention trial is retained as an architecture comparison. Its initial broad-band candidate audit did not improve the retained RNA-specific/no-genomic-product count over the existing CNN. The local-encoding attention model is therefore evaluated on the frozen validation genes before full task evaluation; its scientific outcome is not inferred from lower training loss.

Generation shares the corrected joint-span decoder, individual oligo constraints, joint Tm compatibility, optional exact anchors, temperature 4 and K=20. No failed or duplicate proposal is replaced. The temperature is inherited from the existing validation-selected protocol, rather than claimed to be separately optimized for each attention backbone.

The GPU environment is separate from the existing CPU analysis environment: Python 3.10, CUDA PyTorch 2.7.1+cu128, NumPy 2.2.6 and primer3-py 2.3.1. The shared thermodynamic evaluator uses the same primer3 version. Model inference uses bfloat16 autocasting; candidate assessment and full-reference Bowtie screening use the original numerical and reference rules. Hardware-dependent timings must identify the execution environment and must not be compared as if old CPU and new GPU runs used identical hardware.

Contemporary tool comparisons execute the actual pinned upstream qPrimer and PrimerCast implementations. Their native proposals and constraints are preserved in run provenance and assessed using the same reference-based endpoints. PrimerCast's permissive mismatch/indel target scoring is distinguished from exact compatible binding intervals in this task; native model scores are not experimental amplification measurements.

A second, independent cohort comprises 75 additional RefSeq mouse genes excluded from the original 500. Selection uses the public reference and a fixed random gene order, without model outputs or Primer3 feasibility. All NM isoforms were checked against training-gene isoforms at 80% identity and 80% coverage in both directions using MMseqs2. Cohort files and source checks are under experiments/attention_external/. External results are analyzed separately from the original test cohort.

Experimental results, gene-level paired estimates and paper figures are produced from completed candidates plus their complete shared-reference screens. The old completed study remains intact. The Bioinformatics technical working copy retains old CNN statements until the author revises it using the new results; it does not certify a new attention paper as submission-ready.

## Completed evaluation

All planned training, generation and full-reference screens are complete. The final primary protocol is `experiments/attention/final_protocol.json`. The main analysis uses 75 frozen test genes; the additional analysis uses 75 separate genes. Read-only verification passed for both analyses, including recomputed endpoint means, paired estimates, original sequence geometry and exact-anchor preservation. There were zero anchor violations among returned anchored candidates.

The endpoint below counts distinct candidates passing the common thermodynamic and RNA-specificity checks without a genomic-product flag, at a proposal budget of 20. Cohorts are analyzed separately.

| Method | Original 75 test genes: mean pairs/gene | Additional 75 genes: mean pairs/gene |
| --- | ---: | ---: |
| Local encoding + attention diffusion | 5.6933 | 6.4667 |
| Matched local encoding + attention direct | 4.8800 | 5.1733 |
| Pure attention diffusion | 4.0800 | 4.4667 |
| Original CNN diffusion | 7.3600 | 8.1067 |
| Primer3 | 4.5467 | 4.7200 |
| qPrimer design module | 5.2000 | 5.0667 |

On the additional cohort, the hybrid diffusion model's paired mean difference versus Primer3 is 1.7467 pairs/gene (95% gene-bootstrap interval 0.4267–3.0533); versus qPrimer it is 1.4000 (−0.1733–2.9467). These correspond to mean relative increases of 37.0% and 27.6%, respectively. On the original cohort the differences are 1.1467 (−0.0800–2.4403) and 0.4933 (−0.9467–1.9333). These intervals describe gene resampling of a fixed trained model and do not represent between-training variability.

Attention improved the pure-attention trial when combined with local sequence encoding; it did not surpass the original CNN diffusion in retained candidate count. The complete architecture comparison is therefore retained. Flexible tasks, joint constraints and diffusion-depth comparisons use the same frozen genes, proposal budgets and common assessment; values and paired intervals are in the final summary CSVs. The exon interpretation uses recorded screened candidates and reference exon geometry, including the Ndp example, rather than replacing whole-reference outcomes with a structural proxy.
