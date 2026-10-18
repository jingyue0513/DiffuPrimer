#!/usr/bin/env python3
"""Compare every target RNA base to genomic exons, including strand direction."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.sequence import reverse_complement


root = Path(__file__).resolve().parents[1]
context_path = root / "results/frozen_target_exon_context.json"
contexts = json.loads(context_path.read_text())
targets = {r["transcript_accession_version"]: r for r in read_jsonl(root / "data/processed/targets_pilot.jsonl")}
by_chromosome = defaultdict(list)
for accession, context in contexts.items():
    if "exons" not in context:
        raise ValueError("Target lacks exon context: " + accession)
    by_chromosome[context["exons"][0]["chromosome"]].append(accession)
rows = []
fasta = root / "data/processed/blastdb/genome.fasta"
for record in SeqIO.parse(fasta, "fasta"):
    if record.id not in by_chromosome:
        continue
    # Lowercase genomic bases indicate repeat soft masking, not substitutions.
    genome = str(record.seq).upper()
    for accession in by_chromosome[record.id]:
        context = contexts[accession]
        pieces = []
        for exon in context["exons"]:
            piece = genome[exon["genome_start"]:exon["genome_end"]]
            if exon["strand"] == "-":
                piece = reverse_complement(piece)
            pieces.append(piece)
        spliced = "".join(pieces)
        target = targets[accession]
        rna = target["template_sequence"]
        length = context["genomically_annotated_RNA_length"]
        if len(spliced) != length:
            raise ValueError("Genomic exon extraction length differs")
        mismatches = [i for i, (a,b) in enumerate(zip(spliced, rna[:length])) if a != b]
        rows.append({"gene_id": target["gene_id"], "accession": accession, "split": target["split"],
            "status": "exact_spliced_reference_match" if not mismatches else "spliced_reference_sequence_difference",
            "genomic_RNA_length": length, "terminal_polyA_length": context["terminal_polyA_length"],
            "mismatch_count": len(mismatches), "mismatch_positions": mismatches,
            "RNA_genomic_portion_sha256": hashlib.sha256(rna[:length].encode()).hexdigest(),
            "spliced_genome_sha256": hashlib.sha256(spliced.encode()).hexdigest()})
    print(json.dumps({"chromosome": record.id, "verified_targets": len(by_chromosome[record.id])}), flush=True)
if {r["accession"] for r in rows} != set(targets):
    raise ValueError("Some target chromosomes were not found")
write_jsonl(root / "results/spliced_reference_verification.per_target.jsonl", rows)
write_json(root / "results/spliced_reference_verification.json", {"verified_at": utc_now(), "targets": len(rows),
    "status_counts": dict(Counter(r["status"] for r in rows)), "total_differing_bases": sum(r["mismatch_count"] for r in rows),
    "reference_manifest_sha256": sha256(root / "data/processed/blastdb/genome.manifest.json"),
    "exon_context_sha256": sha256(context_path), "target_sha256": sha256(root / "data/processed/targets_pilot.jsonl"),
    "policy": "All annotated exon bases compared in transcript strand order. Appended RNA polyA tails excluded from genomic comparison and retained explicitly."})
