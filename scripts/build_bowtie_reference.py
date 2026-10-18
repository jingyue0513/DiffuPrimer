#!/usr/bin/env python3
"""Build an independent short-read index from the frozen reference FASTA."""
import argparse
import json
import subprocess
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


parser = argparse.ArgumentParser()
parser.add_argument("reference", choices=["transcriptome", "genome"])
parser.add_argument("--threads", type=int, default=6)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
directory = root / "data/processed/bowtiedb"
directory.mkdir(parents=True, exist_ok=True)
fasta = root / f"data/processed/blastdb/{args.reference}.fasta"
manifest_path = directory / f"{args.reference}.manifest.json"
fasta_hash = sha256(fasta)
if manifest_path.exists():
    previous = json.loads(manifest_path.read_text())
    if previous["fasta_sha256"] != fasta_hash:
        raise ValueError("Existing index reference differs")
    for filename, checksum in previous["index_sha256"].items():
        if sha256(directory / filename) != checksum:
            raise ValueError("Existing index bytes differ")
    print("Existing verified index retained", flush=True)
else:
    binary = root / "tools/bowtie-1.3.1-linux-x86_64/bowtie-build"
    command = [str(binary), "--threads", str(args.threads), "--offrate", "3", str(fasta), str(directory / args.reference)]
    started = utc_now()
    with (directory / f"{args.reference}.build.log").open("w") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
    files = list(directory.glob(args.reference + ".*.ebwt*"))
    if len(files) != 6:
        raise ValueError("Incomplete Bowtie index")
    write_json(manifest_path, {"started_at": started, "completed_at": utc_now(), "command": command,
        "fasta_sha256": fasta_hash, "frozen_blast_reference_manifest_sha256": sha256(root / f"data/processed/blastdb/{args.reference}.manifest.json"),
        "tool_manifest_sha256": sha256(root / "tools/bowtie_manifest.json"),
        "index_sha256": {p.name: sha256(p) for p in files}})
    print(json.dumps({"index": args.reference, "completed_at": utc_now()}), flush=True)
