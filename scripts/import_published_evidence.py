#!/usr/bin/env python3
"""Extract explicit published assay outcomes; retain sequence-level verification."""
import re
from pathlib import Path

from diffuprimer.io import sha256, write_jsonl


root = Path(__file__).resolve().parents[1]
source = root / "data/raw/primerbank/2008_S4_springer.txt"
pdf = source.with_suffix(".pdf")
rows, evidence, failure = [], None, None
for line in source.read_text().splitlines():
    if "5 successful" in line:
        evidence, failure = "experimental_success", None
    if "5 failed" in line:
        evidence, failure = "experimental_failure", "agarose_gel" if "agarose" in line else "sequence_BLAST"
    match = re.match(r"\s*[a-e]\s+(\d+[ac]\d+)\s+(.+?)\s+([ACGT]+)/([ACGT]+)\s+(\d+)\s*$", line)
    if match:
        rows.append({"primerbank_id": match[1], "forward_primer": match[3], "reverse_primer": match[4],
                     "amplicon_length": int(match[5]), "evidence_type": evidence,
                     "experimental_result": "success" if evidence == "experimental_success" else "failure",
                     "failure_reason": failure, "source": "Spandidos et al. 2008 Additional file 4",
                     "source_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2F1471-2164-9-633/MediaObjects/12864_2008_1826_MOESM4_ESM.pdf",
                     "source_pdf_sha256": sha256(pdf),
                     "scope": "Source-assay outcome only; not independent validation of DiffuPrimer candidates"})
assert len(rows) == 15 and sum(row["evidence_type"] == "experimental_success" for row in rows) == 5
write_jsonl(root / "data/processed/published_evidence_examples.jsonl", rows)
print("Imported 15 explicit source-assay outcomes for evidence-parser calibration")
