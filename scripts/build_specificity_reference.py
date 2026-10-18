#!/usr/bin/env python3
import argparse
import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


parser = argparse.ArgumentParser()
parser.add_argument("--genome", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
bin_dir = root / "tools/ncbi-blast-2.17.0+/bin"
directory = root / "data/processed/blastdb"
directory.mkdir(parents=True, exist_ok=True)
names = [("transcriptome", "GCF_000001635.27_GRCm39_rna.fna.gz")]
if args.genome:
    names.append(("genome", "GCF_000001635.27_GRCm39_genomic.fna.gz"))
for label, name in names:
    source = root / "data/raw/refseq" / name
    fasta = directory / (label + ".fasta")
    manifest_path = directory / (label + ".manifest.json")
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest["compressed_reference_sha256"] != sha256(source):
            raise ValueError("Existing BLAST database uses different reference bytes")
        print(f"Existing database retained: {label}", flush=True)
        continue
    if not source.exists():
        raise ValueError(f"Reference missing: {source}; download and verify first")
    with gzip.open(source, "rb") as handle, fasta.open("wb") as output:
        shutil.copyfileobj(handle, output)
    database = directory / label
    command = [str(bin_dir / "makeblastdb"), "-in", str(fasta), "-dbtype", "nucl", "-parse_seqids", "-out", str(database)]
    subprocess.run(command, check=True, stdout=(directory / (label + ".build.log")).open("w"), stderr=subprocess.STDOUT)
    write_json(manifest_path, {"built_at": utc_now(), "command": command,
                              "compressed_reference_sha256": sha256(source), "fasta_sha256": sha256(fasta),
                              "blast_version": subprocess.check_output([str(bin_dir / "blastn"), "-version"], text=True).strip()})
transcript_map = {}
with gzip.open(root / "data/raw/refseq/GCF_000001635.27_GRCm39_genomic.gff.gz", "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) != 9:
            continue
        attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
        gene = re.search(r"(?:^|,)GeneID:(\d+)(?:,|$)", attrs.get("Dbxref", ""))
        accession = attrs.get("transcript_id")
        if accession and gene:
            transcript_map[accession] = gene.group(1)
write_json(directory / "transcript_gene_map.json", transcript_map)
print(f"Transcript gene mapping: {len(transcript_map)}", flush=True)
