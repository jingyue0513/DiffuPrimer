#!/usr/bin/env python3
"""Atomic reference downloads with official compressed-file checksums."""
import hashlib
import argparse
import time
from pathlib import Path

import requests

from diffuprimer.acquisition import REFSEQ_BASE


directory = Path(__file__).resolve().parents[1] / "data/raw/refseq"
directory.mkdir(parents=True, exist_ok=True)
parser = argparse.ArgumentParser()
parser.add_argument("--genome", action="store_true")
args = parser.parse_args()
response = requests.get(REFSEQ_BASE + "md5checksums.txt", timeout=60)
response.raise_for_status()
(directory / "md5checksums.txt").write_bytes(response.content)
checksums = {name.removeprefix("./"): checksum for checksum, name in
             (line.split(maxsplit=1) for line in response.text.splitlines())}
names = ["GCF_000001635.27_GRCm39_rna.fna.gz", "GCF_000001635.27_GRCm39_genomic.gff.gz"]
if args.genome:
    names.append("GCF_000001635.27_GRCm39_genomic.fna.gz")
for name in names:
    target = directory / name
    if target.exists():
        with target.open("rb") as handle:
            existing_digest = hashlib.md5()
            for block in iter(lambda: handle.read(2**20), b""):
                existing_digest.update(block)
        if existing_digest.hexdigest() == checksums[name]:
            print(f"Already verified: {name}", flush=True)
            continue
    for attempt in range(5):
        try:
            digest = hashlib.md5()
            with requests.get(REFSEQ_BASE + name, stream=True, timeout=(30, 120)) as response:
                response.raise_for_status()
                with target.with_suffix(target.suffix + ".part").open("wb") as handle:
                    for block in response.iter_content(2**20):
                        digest.update(block)
                        handle.write(block)
            if digest.hexdigest() != checksums[name]:
                raise ValueError(f"MD5 mismatch for {name}")
            target.with_suffix(target.suffix + ".part").replace(target)
            print(f"Downloaded and verified: {name}", flush=True)
            break
        except (requests.RequestException, ValueError) as exc:
            print(f"Attempt {attempt + 1} failed: {exc}", flush=True)
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
