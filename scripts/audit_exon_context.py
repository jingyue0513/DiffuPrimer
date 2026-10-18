#!/usr/bin/env python3
"""Attach annotated exon/junction context to primer intervals on frozen targets."""
import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="*")
    parser.add_argument("--output", default="results/exon_context_audit.json")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target_path = root / "data/processed/targets_pilot.jsonl"
    targets = {r["transcript_accession_version"]: r for r in read_jsonl(target_path)}
    annotation = root / "data/raw/refseq/GCF_000001635.27_GRCm39_genomic.gff.gz"
    exon_rows = defaultdict(set)
    with gzip.open(annotation, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "exon":
                continue
            attributes = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
            accession = unquote(attributes.get("transcript_id", ""))
            if accession in targets:
                exon_rows[accession].add((fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]))
    contexts = {}
    for accession, target in targets.items():
        exons = sorted(exon_rows[accession], key=lambda e: e[1])
        chromosomes = {e[0] for e in exons}
        strands = {e[3] for e in exons}
        if not exons or len(chromosomes) != 1 or len(strands) != 1:
            contexts[accession] = {"status": "missing_or_multiple_loci"}
            continue
        if exons[0][3] == "-":
            exons.reverse()
        length = sum(e[2] - e[1] for e in exons)
        tail_length = len(target["template_sequence"]) - length
        appended_poly_a = tail_length > 0 and set(target["template_sequence"][-tail_length:]) == {"A"}
        if length != len(target["template_sequence"]) and not appended_poly_a:
            contexts[accession] = {"status": "exon_sum_reference_RNA_length_mismatch", "exon_sum_length": length,
                "RNA_length": len(target["template_sequence"])}
            continue
        boundaries, offsets = [], []
        offset = 0
        for chrom, start, end, strand in exons:
            offsets.append({"RNA_start": offset, "RNA_end": offset + end - start, "chromosome": chrom,
                "genome_start": start, "genome_end": end, "strand": strand})
            offset += end - start
            boundaries.append(offset)
        contexts[accession] = {"status": "annotated_exons_plus_terminal_polyA" if appended_poly_a else "annotated_exon_lengths_match_RNA",
            "boundaries": boundaries[:-1], "exons": offsets, "genomically_annotated_RNA_length": length,
            "terminal_polyA_length": tail_length if appended_poly_a else 0}
    records, summaries = [], []
    for raw in args.directories:
        directory = root / raw
        completion = json.loads((directory / "completion.json").read_text())
        if not completion["all_expected_targets_evaluated"] or completion["candidates_sha256"] != sha256(directory / "candidates.jsonl"):
            raise ValueError("Incomplete or changed experiment")
        mapped = eligible = spanning = crossed = same_exon = 0
        for row in read_jsonl(directory / "candidates.jsonl"):
            if not row["assessment"]["computational_feasible"]:
                continue
            eligible += 1
            accession = row["transcript_accession_version"]
            context = contexts[accession]
            record = {"candidate_id": row["candidate_id"], "gene_id": row["gene_id"], "run_directory": str(directory),
                "transcript_accession_version": accession, "mapping_status": context["status"]}
            if context["status"] in {"annotated_exon_lengths_match_RNA", "annotated_exons_plus_terminal_polyA"}:
                mapped += 1
                fs, fe, rs, re = [row[key] for key in ["forward_start", "forward_end", "reverse_start", "reverse_end"]]
                junctions = context["boundaries"]
                f_junctions = [j for j in junctions if fs < j < fe]
                r_junctions = [j for j in junctions if rs < j < re]
                amp_junctions = [j for j in junctions if fs < j < re]
                genomic_length = context["genomically_annotated_RNA_length"]
                same_annotated_exon = not amp_junctions and re <= genomic_length
                record.update(forward_spans_junction=bool(f_junctions), reverse_spans_junction=bool(r_junctions),
                    forward_terminal_five_bases_span_junction=any(fe - 5 < j < fe for j in f_junctions),
                    reverse_terminal_five_bases_span_junction=any(rs < j < rs + 5 for j in r_junctions),
                    amplicon_crosses_annotated_junction=bool(amp_junctions), entire_amplicon_within_one_exon=same_annotated_exon,
                    forward_overlaps_terminal_polyA=fe > genomic_length, reverse_overlaps_terminal_polyA=re > genomic_length)
                spanning += bool(f_junctions or r_junctions)
                crossed += bool(amp_junctions)
                same_exon += same_annotated_exon
                if not f_junctions and not r_junctions and re <= genomic_length:
                    def genomic_base(position):
                        exon = next(e for e in context["exons"] if e["RNA_start"] <= position < e["RNA_end"])
                        shift = position - exon["RNA_start"]
                        return exon["genome_start"] + shift if exon["strand"] == "+" else exon["genome_end"] - 1 - shift
                    record["annotated_target_genomic_product_bp"] = abs(genomic_base(re - 1) - genomic_base(fs)) + 1
            records.append(record)
        summaries.append({"run_directory": str(directory), "feasible_proposals": eligible, "annotation_mapped_proposals": mapped,
            "primer_junction_spanning_proposals": spanning, "amplicon_junction_crossing_proposals": crossed,
            "same_exon_amplicon_proposals": same_exon, "mapping_denominator": mapped,
            "candidate_sha256": completion["candidates_sha256"]})
    output = root / args.output
    write_jsonl(output.with_suffix(".per_candidate.jsonl"), records)
    write_json(output, {"generated_at": utc_now(), "target_sha256": sha256(target_path), "GFF_sha256": sha256(annotation),
        "target_mapping_status_counts": {status: sum(c["status"] == status for c in contexts.values()) for status in {c["status"] for c in contexts.values()}},
        "runs": summaries, "interpretation": "Annotated transcript exon geometry; exact appended all-A terminal tails are recorded as non-genomic RNA sequence. Genomic binding is evaluated separately by the recorded full-reference alignment screen."})
    write_json(root / "results/frozen_target_exon_context.json", contexts)
    print(json.dumps({"targets": len(contexts), "mapping_status": json.loads(output.read_text())["target_mapping_status_counts"], "runs": len(summaries)}))


if __name__ == "__main__":
    main()
