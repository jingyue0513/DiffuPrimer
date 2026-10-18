#!/usr/bin/env python3
"""Frozen random gene pilot, resumable cached sequential acquisition."""
import argparse
import json
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.acquisition import CachedHTTP, RateLimiter, SourceAccessBlocked, acquire_gene, load_reference_sequences, reference_gene_catalog
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", type=int, default=500)
    parser.add_argument("--limit", type=int, help="Run only first N frozen genes for parser calibration")
    parser.add_argument("--workers", type=int, default=1, choices=[1, 2, 3], help="Shared global limit: at most one new HTTP request per second")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if (root / "results/source_access_block.json").exists():
        raise SourceAccessBlocked("A source-access block is recorded. Use cached data and continue other study work until access is restored; do not auto-retry.")
    config = json.loads((root / "configs/pilot.json").read_text())
    gff = root / "data/raw/refseq/GCF_000001635.27_GRCm39_genomic.gff.gz"
    fasta = root / "data/raw/refseq/GCF_000001635.27_GRCm39_rna.fna.gz"
    selection_path = root / "data/processed/pilot_gene_selection.jsonl"
    if selection_path.exists():
        selection = read_jsonl(selection_path)
        if len(selection) != args.genes:
            raise ValueError("Pilot selection already frozen with a different size")
    else:
        catalog = reference_gene_catalog(gff)
        rng = random.Random(config["seed"])
        selection = rng.sample(catalog, min(args.genes, len(catalog)))
        write_jsonl(selection_path, selection)
        write_json(selection_path.with_suffix(".manifest.json"), {
            "selected_at": utc_now(), "seed": config["seed"], "catalog_genes": len(catalog),
            "selected_genes": len(selection), "selection_method": "uniform random sample of protein-coding genes with NM_ reference transcripts, before querying PrimerBank outcomes",
            "gff_sha256": sha256(gff), "fasta_sha256": sha256(fasta), "selection_sha256": sha256(selection_path)})
    sequence_ids = {acc for gene in selection for acc in gene["reference_transcripts"]}
    references = load_reference_sequences(fasta, sequence_ids)
    limiter = RateLimiter(delay=1.0)
    destination = root / "data/processed/primerbank_pilot.jsonl"
    rows = read_jsonl(destination) if destination.exists() else []
    ledger_path = root / "data/processed/acquisition_ledger.jsonl"
    ledger = read_jsonl(ledger_path) if ledger_path.exists() else []
    latest = {x["gene_id"]: x["status"] for x in ledger}
    completed = {gene for gene, status in latest.items() if status in {"complete", "no_records_returned"}}
    scheduled = [(index, gene) for index, gene in enumerate(selection[:args.limit] if args.limit else selection)
                 if gene["gene_id"] not in completed]

    def worker(gene):
        # Separate sessions; shared limiter caps requests across workers.
        http = CachedHTTP(root / "data/raw/primerbank/cache", limiter=limiter)
        return acquire_gene(gene, http, references)

    executor = ThreadPoolExecutor(max_workers=args.workers)
    futures = {executor.submit(worker, gene): (index, gene) for index, gene in scheduled}
    for future in as_completed(futures):
        index, gene = futures[future]
        try:
            records = future.result()
            rows.extend(records)
            write_jsonl(destination, rows)
            ledger.append({"gene_id": gene["gene_id"], "status": "complete" if records else "no_records_returned", "records": len(records), "at": utc_now(),
                           "query_field": "NCBI Gene Symbol", "gene_symbol": gene["gene_symbol"]})
            completed.add(gene["gene_id"])
            print(json.dumps({"index": index + 1, "gene_id": gene["gene_id"], "pairs": len(records),
                              "mapped": sum(r["mapping_status"] == "unique_exact_pair_remap" for r in records)}, ensure_ascii=False), flush=True)
        except SourceAccessBlocked as exc:
            entry = {"gene_id": gene["gene_id"], "status": "source_access_blocked", "error": str(exc), "at": utc_now()}
            ledger.append(entry)
            write_jsonl(ledger_path, ledger)
            write_json(root / "results/source_access_block.json", entry)
            for pending in futures:
                pending.cancel()
            print(json.dumps(entry), flush=True)
            break
        except Exception as exc:
            ledger.append({"gene_id": gene["gene_id"], "status": "error", "error": str(exc), "at": utc_now()})
            print(json.dumps(ledger[-1]), flush=True)
        write_jsonl(ledger_path, ledger)
        write_json(root / "results/acquisition_status.json", {
            "updated_at": utc_now(), "planned_genes": len(selection), "completed_genes": len(completed),
            "pairs": len(rows), "mapping_status": dict(Counter(r["mapping_status"] for r in rows)),
            "evidence_type": dict(Counter(r["evidence_type"] for r in rows)),
            "experimental_record_available": sum(r.get("experimental_record_available", False) for r in rows)})
    executor.shutdown(wait=True, cancel_futures=True)


if __name__ == "__main__":
    main()
