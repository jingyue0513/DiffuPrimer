# Contemporary baseline execution

## qPrimer (qPrimerDB 2.0, NAR 2025 issue)

The actual upstream `design.py` is pinned at `34f047e5d77ec0353846a9a18b0b94f74013db3b`.

```bash
git clone https://github.com/swu1019lab/qPrimer.git experiments/attention/baselines/vendor/qPrimer
git -C experiments/attention/baselines/vendor/qPrimer checkout 34f047e5d77ec0353846a9a18b0b94f74013db3b
.venv/bin/python experiments/attention/baselines/run_qprimer.py
```

The adapter matches common biological bounds/concentrations and disables internal probe generation for a pair-only comparison. Other native qPrimer rules are retained. Completion and redesign accept `--task completion|redesign --fixed-side forward|reverse`. Length-band runs accept `--amplicon-min` and `--amplicon-max`.

## PrimerCast (public manuscript and Broad Institute pretrained software)

Source and pretrained weights are pinned at `50d357b43a727196b5ae2e1c2e7d041519243673`. `run_primercast.py` calls its **unchanged native quick-design function**: primer tiling, native Tm/GC/duplex filters, Bowtie2 mapping, native feature preparation, pretrained classifier/regressor and pair filtering. Supply RefSeq accession FASTA headers so native pandas parsing preserves string identifiers.

Needed runtime: Python >=3.10; PyTorch, Biopython, pandas, NumPy, joblib, einops, scikit-learn, SciPy, openpyxl; official RNAduplex, Bowtie2, sam2pairwise. The execution here uses CPU and a separate `.venv_primercast`, plus official RNAduplex 2.4.17, Bowtie2 2.4.4 and source-built sam2pairwise v1.0.0. Actual tool versions and checkpoint/source hashes are in each manifest. The author's conda environment suggests newer RNAduplex/Bowtie2 versions; package constraints and the observed execution versions are distinct.

```bash
git clone https://github.com/broadinstitute/PrimerCast.git experiments/attention/baselines/vendor/PrimerCast
git -C experiments/attention/baselines/vendor/PrimerCast checkout 50d357b43a727196b5ae2e1c2e7d041519243673
python -m venv .venv_primercast
.venv_primercast/bin/pip install experiments/attention/baselines/vendor/PrimerCast
PATH="$PWD/.venv_primercast/bin:$PATH" CUDA_VISIBLE_DEVICES='' .venv_primercast/bin/python experiments/attention/baselines/run_primercast.py
```

If using a git clone rather than the downloaded source snapshot, add `source_snapshot.json` with the repository URL and pinned commit; the runner uses that acquisition metadata. The source-fetch script below creates the metadata automatically without downloading git history.

```bash
.venv/bin/python experiments/attention/baselines/fetch_primercast.py
```

Supported length/Tm/GC maximum controls use the native parameter interface. Native ΔG>=−6 kcal/mol and activity/coverage thresholds retain upstream defaults. K=20 selects the highest native-scoring returned pairs; shortages, incompatible exact positions and failures of common final constraints remain in the budget. It exposes no fixed-primer interface, so completion/redesign are not attributed to this tool.

## Independent PICKR literature-source mapping

The original CSV is pinned at database commit `f72760013722c2d698c9c7f3b80ac1640ea8b4c7`, with acquisition metadata under `data/raw/pickr/`. Run the mapper with `pyahocorasick` installed:

```bash
.venv_primercast/bin/python experiments/attention/baselines/map_pickr_literature.py
```

Outputs are exact oligo binding sites, gene coverage and PMCID evidence under `experiments/attention/results/pickr_literature/`. The original resource contains no per-record species field and many ambiguous gene/orientation hints, so mapping uses sequence identity against the frozen references. It supplies neither experimentally labelled pairs nor the published PICKR ranked-pair output, and enters no model training.

All candidate-generating baseline outputs share the `candidates.jsonl`, `per_gene.jsonl`, `manifest.json`, `completion.json` schema and receive the same downstream full-reference RNA/genome assessment as the learned methods.
