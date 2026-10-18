#!/usr/bin/env python3
"""Execute an unchanged evaluation script with deterministic span caching."""
import argparse
import runpy
import sys
from pathlib import Path

from span_cache import install
from diffuprimer.io import sha256


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--evaluation-script", required=True)
args, remaining = parser.parse_known_args()
root = Path(__file__).resolve().parents[1]
import json
check = json.loads((root / "results/checks/span_cache/check.json").read_text())
if not check["exact_output_equality"] or check["cache_sha256"] != sha256(root / "scripts/span_cache.py"):
    raise ValueError("Span cache equivalence check missing or stale")
install()
sys.argv = [args.evaluation_script, *remaining]
runpy.run_path(args.evaluation_script, run_name="__main__")
