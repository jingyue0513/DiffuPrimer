#!/usr/bin/env python3
"""Compatibility entry point for the final manuscript and supplement."""
import argparse
import runpy
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--aggregate", help="Historical argument; final fixed-seed results populate the paper")
parser.parse_args()
runpy.run_path(str(Path(__file__).with_name("finalize_paper.py")), run_name="__main__")
