# Frozen benchmark and evidence archives

Run `python scripts/verify_release.py` to check every archive member and completed screen record, then `python scripts/restore_release_data.py` to restore paths. Existing different files are never silently overwritten.

- `benchmark_data.zip`: training/validation labels and provenance, frozen targets/splits, small source records and the two historical seed-20261013 CNN checkpoints.
- `original_evidence.zip`: original-cohort published summaries and candidate/per-gene/full-screen evidence from all source runs used by those summaries.
- `additional_evidence.zip`: additional-cohort summary and candidate/full-screen evidence.
- `mechanism_seed42.zip`: fixed-seed mechanism protocols, results, validation selection, per-gene records and complete screens, including the uncapped deadline follow-up.
- `data_manifest.json`: SHA-256 and byte sizes of archives and all members. Historical files are preserved byte for byte.

Whole RefSeq references, aligner indexes, SQLite caches, raw SAM hit files, third-party vendor code and binaries are excluded. Public reference accessions and manifests identify those inputs. `THIRD_PARTY_NOTICES.md` explains provenance and license boundaries.

The source release includes reconstruction and evaluation code. Numerical results are computational reference-screen outcomes; the study has no wet-lab assay outcome dataset. A dedicated persistent archive identifier is not yet assigned.
