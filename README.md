# DiffuPrimer

Conditional discrete diffusion over paired binding-position labels for qPCR primer design. The target DNA stays fixed; a constrained interval decoder produces forward and reverse oligos. One model supports alternative pairs, product-length control, fixed-primer completion and local redesign.

Authors: **Jingyue Li** (Sun Yat-sen University; Zhengzhou University; Henan Children's Hospital) and **Shiyue Mei** (Zhengzhou University; Henan Children's Hospital; corresponding author, 18837614038@163.com). Jingyue Li: lijingyue513@163.com.

## Install and try

Python 3.10 or later, Linux. The provided checkpoints run on CPU. From a clean checkout:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
python scripts/verify_release.py
python scripts/restore_release_data.py
python scripts/design_cnn.py --fasta examples/attention_reference.fasta \
  --output example_design.json --seed 42 --steps 4 --budget 20
```

The example does not require downloading whole reference genomes. It reports oligos, zero-based half-open binding intervals and thermodynamic assessments. Full transcriptome/genome specificity screening is a separate step; the example output explicitly records that it has not been screened. See `python scripts/design_cnn.py --help` for every option and default. Set `--steps 1` or `--steps 16` for the tested recovery depths. Supply `--fixed-forward START:END` or `--fixed-reverse START:END` for anchoring. Narrow product bands are controlled by `--amplicon-min` and `--amplicon-max`.

The constant generation seed is **42**. Historical checkpoints and benchmark runs retain their training seed 20261013 in filenames and provenance. No clock-derived seed is used. The mechanism study reuses one Direct/diffusion checkpoint pair.

## What the experiments support

Both held-out cohorts contain 75 mouse genes. Retained pairs satisfy common computational constraints, intended-gene RNA specificity and no predicted genomic product under the recorded full-reference screen. These are reference-screen endpoints.

| Seed-42 mechanism control | Original genes: pairs/gene | Additional genes: pairs/gene |
|---|---:|---:|
| Direct, validation-selected strong decoding | 5.65 | 5.48 |
| Denoising-trained one-step, validation-selected | 7.48 | 7.99 |
| Frozen initial scores, 16 recovery steps | 7.40 | 8.45 |
| Full diffusion, 16 network updates | 7.43 | 8.21 |
| Direct, uncapped matched generation-CPU budget | 4.28 | 4.28 |

Denoising-trained scores remain useful with strong static decoding. Sixteen network updates do not consistently improve retained counts over frozen scores, but increase retained-pair positional diversity. Matched-CPU Direct search covers 80% of genes in each cohort versus 70.7%/73.3% for full diffusion, illustrating a count/coverage trade-off. Full details, paired intervals and CPU accounting are in [the mechanism study](experiments/mechanism/README.md), and after restoring the release data in `experiments/mechanism/uncapped_seed42/summary.json`.

## Reproduce and inspect

```bash
python scripts/verify_release.py
python -m pytest -q
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 python experiments/mechanism/test_controls.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 python experiments/mechanism/test_uncapped.py
```

The `release/` archives restore frozen data, checkpoints, complete numerical summaries, candidate outputs and screen records to their historical paths. `release/data_manifest.json` records every file and SHA-256. Offline verification checks all archive members and completed candidate/screen identities. The original absolute paths inside historical evidence are preserved for provenance; restore into this checkout and use the current repository-relative drivers for new work.

For full reference screening, obtain NCBI RefSeq GRCm39 (GCF_000001635.27; annotation RS_2024_02) with `scripts/download_reference.py`, verify the reference manifest, and build the Bowtie indexes. Whole reference files, hit caches and alignment databases are excluded from Git. See [reference/reproduction details](docs/reproduction.md) and [the data dictionary](docs/data_dictionary.md). The partial PrimerBank snapshot is included; replay does not require contacting a blocked service.

The CNN implementation used by the final controls is preserved in `experiments/hybrid_attention/src/diffuprimer/cnn_model.py`; `scripts/design_cnn.py` selects it explicitly. `src/diffuprimer/` retains the original package and training code. This separation preserves frozen source identities. Model/architecture comparison drivers, protocols and source snapshots are retained under `experiments/`.

## License and citation

Original code is [MIT licensed](LICENSE). Reference data and third-party tools retain their own terms; see [third-party notices](THIRD_PARTY_NOTICES.md). [CITATION.cff](CITATION.cff) records software authorship. The manuscript is being prepared for Bioinformatics; no acceptance, publication DOI or independent archival DOI is claimed. An independent software/test-data archive must accompany the journal submission.
