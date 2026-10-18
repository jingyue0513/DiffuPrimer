# Diffusion mechanism controls

This experiment reuses one existing CNN direct/diffusion checkpoint pair. New
generation and bootstrap seeds are the constant **42**, never the current time.
The checkpoint filenames retain the historical training seed `20261013`; this
is provenance, not an additional training run or a newly selected seed.

The original 75 test genes and the separate 75 external genes are evaluated
separately. No wet experiments or scaling runs are part of this round.

## Run and resume

From the repository root:

```bash
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  .venv_attention/bin/python experiments/mechanism/test_controls.py

OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  .venv_attention/bin/python -u experiments/mechanism/run.py \
  > experiments/mechanism/run_seed42.log 2>&1
```

The runner freezes the configuration, source hashes, checkpoint hashes, cohorts
and reference-independent selection rules before generation. It writes complete
per-gene records atomically, resumes completed genes, exports 22 method/cohort
directories, runs the existing full-reference Bowtie audit, and generates the
analysis. A changed frozen implementation or protocol requires a new output
directory; it cannot silently resume old records.

Run the validation-selected supplemental baselines and final comparison with:

```bash
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  .venv_attention/bin/python -u experiments/mechanism/tune_direct.py \
  > experiments/mechanism/tuned_seed42.log 2>&1

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv_attention/bin/python -u experiments/mechanism/finalize.py \
  > experiments/mechanism/finalize_seed42.log 2>&1
```

The supplemental run compares full-support deterministic top-k with weighted
without-replacement temperatures 0.5/1/2/4, separately for the direct and diffusion
checkpoints. It selects on all 75 original validation genes using the same full
RNA/genome screen and the predeclared retained-pair count endpoint (ties: >=2
coverage, then configuration order). Selection is saved before generating or
screening its test outputs. This retains the original fixed-temperature controls
and adds a matched, validation-selected strong baseline. It uses the same single
seed 42; configuration search is not a multiple-seed experiment.

`tuned_seed42/selection.json` records all validation configurations and their
screen hashes. `comparison_seed42/summary.json`, `summary.txt`, and `results.png`
combine the 11 main methods with the two validation-selected methods, separately
on both cohorts. `finalize.py` waits for both studies and stops on a failure.

### Uncapped CPU-deadline follow-up

The completed runtime audit found that all 150 `direct_equal_cpu` runs reached
the 4,096-candidate cap before their CPU deadlines. Preserve that result as a
bounded-pool control, but use the following follow-up for a genuinely
deadline-driven comparison:

```bash
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  .venv_attention/bin/python -u experiments/mechanism/run_uncapped.py \
  > experiments/mechanism/uncapped_seed42.log 2>&1
```

Run this after both original studies and `finalize.py` have completed. It reuses
the exact measured diffusion16 CPU budgets, weights, seed, full-support sampling
and diversity-selection objective. It removes the candidate-count cap and
commits selection checkpoints every 512 pool additions, reserving the measured
selection cost near the deadline. A target stops early only when its finite
compatible support has been exhausted. Eight independent worker processes
parallelize targets; each target's budget is measured by its own process CPU
clock. Only the output committed before the deadline is evaluated; actual
overrun is reported.

`uncapped_seed42/summary.json` and `results.pdf` are the final combined report,
including all original controls plus `direct_equal_cpu_uncapped` (14 methods per
cohort). The inherited 4,096-cap limitation in the original comparison applies
to `direct_equal_cpu`, not to `direct_equal_cpu_uncapped`; the latter has separate
budget-utilization and deadline diagnostics. The code change was prompted by
runtime-budget utilization, not selection of a favorable test quality outcome.

Additional deadline and finite-support checks:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv_attention/bin/python experiments/mechanism/test_uncapped.py
```

`results_seed42/status.json` records the active stage. The final products are
`summary.json`, `groups.csv`, `contrasts.csv`, `per_gene.jsonl`, `summary.txt`,
`verification.json`, and PNG/PDF result figures. All raw proposals and screens are
retained under `evaluations/`. Logs and timestamps document execution only and
never determine a random seed.

## Predeclared methods

| Method | Network | Decoder / selection |
| --- | --- | --- |
| direct_beam | Direct training, one step | Original beam 64 × 8 |
| diffusion_1 / 4 / 16 | Diffusion training, indicated steps | Original decoder and reveal process |
| frozen_4 / 16 | Diffusion checkpoint, first scores reused | Same decoder/reveal process, network called once |
| direct_full_wor | Direct, one cached score map | All compatible pairs, weighted without replacement, 20 outputs |
| diffusion1_full_wor | Diffusion, one cached score map | Identical full-support decoder |
| direct_pool320 | Direct, one cached score map | 320 sequence-distinct samples, feasible farthest-first selection of 20 |
| diffusion1_pool320 | Diffusion, one cached score map | Identical pool/selection control |
| direct_equal_cpu | Direct, one cached score map | Anytime beam then full-support pools, up to 4,096 candidates |

Full support enumerates every pair satisfying geometry, individual-oligo and
joint-Tm support, with the same pair-level assessment as the existing experiment.
Sequence duplicates are removed from the stronger direct samplers. No genome or
transcriptome screening result, exon annotation, or test outcome enters proposal
selection. Temperature is fixed at the existing value 4, without new test tuning.

Pool selection starts at the smallest shared Tm penalty, then greedily maximizes
minimum amplicon Jaccard distance. Ties use Tm penalty then original pool order.
The matched diffusion-one-step versions isolate training effects under the same
strengthened decoder. Pools and dropped candidates are counted as internal work;
the final output budget remains 20 and shortages remain failures.

## Compute accounting

The equal-budget endpoint uses **process CPU seconds**, including network,
support construction, random sampling, sequence deduplication, thermodynamic
assessment, and selection. Each target's full 16-step diffusion cost supplies the
budget. Both sides exclude the same shared, cached span preparation and model
loading. Pair-property caches are cleared before each method. Ordinary method
order rotates by gene; the dependent timed method follows its reference run.

The direct method commits completed output checkpoints at pool sizes 20, 40, 80,
160, 320, 640, 1,280, 2,560 and 4,096. Only checkpoints completed before the CPU
deadline are used. Non-preemptible operation overrun is recorded separately;
actual CPU and wall time are retained, not replaced by the requested budget.
The cap can leave unused compute, which is reported through actual costs and
pool size. Saved output/checkpoint traces support audit of time-dependent work.

This is an inference CPU-budget comparison, **not** equal GPU FLOPs or an
end-to-end runtime benchmark. Full-reference screening is applied uniformly
after selection, reusing exact-sequence alignment caches. Its batch runtime is
not attributed to individual methods. No inference efficiency claim should be
extended to total screening cost without a separate benchmark.

## Statistics and interpretation

The primary endpoint is the number of unique RNA-specific primer pairs without
a predicted genomic product. Report all 75 genes per cohort, including no-output
genes. Secondary endpoints include >=1 and >=2 retained-pair coverage, diversity
among retained pairs, temperature penalty, and CPU/wall cost. Bootstrap target
genes as paired units using seed 42 and 10,000 resamples. Confidence intervals
are unadjusted descriptive intervals; one reused training seed does not measure
training-seed robustness.

The key contrasts separate:

1. Diffusion-one-step versus direct: benefit of the trained score distribution.
2. Full versus frozen at the same step count: benefit of learned state updates.
3. Frozen versus one step: benefit of repeated constrained sampling alone.
4. Matched full-support/pool decoders: training benefits after strengthening the
   non-iterative baseline.
5. Full diffusion versus deadline-committed direct output: utility under the
   same inference CPU allowance.

## Scaling and generalization after this decision

Do not start scaling automatically. If learned iterative updates retain a
meaningful advantage, compare matched direct/diffusion models at CNN widths
96/192/384 (235,011 / 930,819 / 3,704,835 parameters), then nested training sets
of 350/1,000/3,000/10,000 genes. If the benefit is chiefly denoising training,
retain the one-step model as a primary efficient method rather than attributing
its benefit to sixteen-step inference.

Use one explicit constant seed, frozen gene-family-separated test sets, and
training-compute accounting for all future scaling runs. Generalization should
include a new held-out gene cohort and, after adding other reference species,
an entirely held-out species with orthology/homology controls. The present
external cohort measures same-species unseen-gene transfer, not cross-species
generalization. Keep data-source mixtures and per-gene supervision controlled
when varying dataset size.
