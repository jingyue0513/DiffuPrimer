# Additional mouse RefSeq gene cohort

`data/processed/attention_external_targets.jsonl` contains 75 new protein-coding mouse genes, selected from the existing official GRCm39 RNA/GFF snapshot and disjoint from every gene in the original 500-gene study. These are additional genes from the same species and source reference, for a new held-out gene evaluation.

Selection uses seed 20261013 to shuffle the full gene-ID-sorted catalog after excluding the original 500 IDs. Eligibility uses complete NM_ sequence availability and MMseqs2 similarity to original training genes; no primer-design, thermodynamic-feasibility or model-output criterion enters selection. The template is the lexicographically first annotated NM_ accession version for each gene.

All NM_ isoforms in the first 128 random candidates were compared with all 781 NM_ isoforms of the original 350 training genes in both search directions. There were no hits satisfying 80% nucleotide identity with at least 80% query **and** target coverage. Thus, the first 75 candidates in the frozen random order were retained. The selected cohort has 169 isoforms; template lengths range from 512 to 8200 bases, with median 2258. The `split` value is `external_test`.

```bash
.venv/bin/python experiments/attention_external/prepare_cohort.py
```

The cohort was prepared without running any model or Primer3. Evaluation should retain all 75 targets, including targets producing no valid candidates. This cohort has not yet been evaluated by this preparation script.

Provenance and selection files:

- `protocol.json`: selection rules, reference hashes and original partition identity.
- `candidate_order.jsonl`: complete randomized candidate order.
- `selection_ledger.jsonl`: availability, homology decisions and selected flag for searched genes.
- `gene_selection.jsonl`, `external_split.jsonl`: selected metadata and partition designation.
- `selected_all_isoforms.fasta`: all selected gene isoforms.
- `homology/`: actual MMseqs2 commands, logs and hit tables for both directions.
- `source_check.json`, `manifest.json`: completed checks and source/output hashes.

Target SHA-256: `c1ea14d95321fd57f3b73928a45592f1e0f90a1018a78516f2f946cb5e638b47`.

External tool generation uses all 75 frozen targets, K=20, seed 20261013 and the common 80–200 bp conditions:

```bash
.venv/bin/python experiments/attention_external/run_baselines.py --method primer3
.venv/bin/python experiments/attention_external/run_baselines.py --method qprimer
.venv/bin/python experiments/attention_external/run_baselines.py --method primercast
```

Outputs are stored separately under `experiments/attention_external/results/evaluations/`, with explicit `external_test` metadata. PrimerCast executes the same native pretrained quick-design procedure as the original-cohort comparison. Common constraint failures and shortages remain in the proposal budget. All tools still require the shared full-reference RNA/genome assessment after generation.
