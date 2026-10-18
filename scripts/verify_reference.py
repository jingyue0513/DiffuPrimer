#!/usr/bin/env python3
import gzip
import hashlib
import json
import argparse
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json
from diffuprimer.acquisition import REFSEQ_BASE


root = Path(__file__).resolve().parents[1]
directory = root / "data/raw/refseq"
parser = argparse.ArgumentParser()
parser.add_argument("--genome", action="store_true")
args = parser.parse_args()
checksums = {}
for line in (directory / "md5checksums.txt").read_text().splitlines():
    checksum, name = line.split(maxsplit=1)
    checksums[name.removeprefix("./")] = checksum
manifest = {"verified_at": utc_now(), "reference_assembly": "GCF_000001635.27",
            "annotation": "GCF_000001635.27-RS_2024_02", "files": []}
names = ["GCF_000001635.27_GRCm39_rna.fna.gz", "GCF_000001635.27_GRCm39_genomic.gff.gz"]
if args.genome:
    names.append("GCF_000001635.27_GRCm39_genomic.fna.gz")
for name in names:
    path = directory / name
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            digest.update(block)
    if digest.hexdigest() != checksums[name]:
        raise ValueError(f"Reference MD5 mismatch for {name}")
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(2**20), b""):
            pass
    manifest["files"].append({"name": name, "url": REFSEQ_BASE + name,
                              "md5": digest.hexdigest(), "sha256": sha256(path), "bytes": path.stat().st_size})
write_json(root / "data/processed/reference_manifest.json", manifest)
print(json.dumps(manifest, indent=2))
