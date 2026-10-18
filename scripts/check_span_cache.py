#!/usr/bin/env python3
"""Compare real model samples before/after caching deterministic calculations."""
import time
from pathlib import Path

import numpy as np
import torch

from diffuprimer import guidance
from diffuprimer.constraints import Conditions
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.sampling import generate
from span_cache import install


root = Path(__file__).resolve().parents[1]
torch.set_num_threads(8)
template = read_jsonl(root / "data/processed/targets_pilot.jsonl")[0]["template_sequence"][:512]
conditions = Conditions()
start = time.monotonic()
original_spans = guidance.feasible_oligo_spans(template, conditions)
original_temperatures = guidance.compatible_span_temperatures(template, conditions, original_spans)
assert all(int(flags.sum()) > 0 for flags in original_spans[1])
fresh_seconds = time.monotonic() - start
install()
cached_spans = guidance.feasible_oligo_spans(template, conditions)
cached_temperatures = guidance.compatible_span_temperatures(template, conditions, cached_spans)
start = time.monotonic()
disk_spans = guidance.feasible_oligo_spans(template, conditions)
disk_temperatures = guidance.compatible_span_temperatures(template, conditions, disk_spans)
cached_seconds = time.monotonic() - start
for arrays in [cached_spans, disk_spans]:
    for left, right in zip([original_spans[0], *original_spans[1]], [arrays[0], *arrays[1]]):
        assert np.array_equal(left, right)
assert np.array_equal(original_temperatures, cached_temperatures, equal_nan=True)
assert np.array_equal(original_temperatures, disk_temperatures, equal_nan=True)
nonempty_counts = {}
for mode in ["direct", "diffusion"]:
    checkpoint = torch.load(root / f"results/checkpoints/{mode}_weak_partial_pilot_seed20261013.pt", weights_only=True)
    model = ConditionalDenoiser(**checkpoint["model_spec"])
    model.load_state_dict(checkpoint["state_dict"])
    outputs = []
    for spans, temperatures in [(original_spans, original_temperatures), (disk_spans, disk_temperatures)]:
        outputs.append(generate(model, template, conditions, budget=5, steps=4, seed=20261013, mode=mode,
            span_allowed=spans[1], span_temperatures=temperatures))
    for fresh, cached in zip(*outputs):
        assert fresh["pair"] == cached["pair"]
        assert np.array_equal(fresh["labels"], cached["labels"])
    nonempty_counts[mode] = sum(row["pair"] is not None for row in outputs[0])
    assert nonempty_counts[mode] > 0
    rule_fresh = guidance.rule_span_sampler(template, conditions, budget=20, seed=20261013,
        spans=original_spans, span_temperatures=original_temperatures)
    rule_cached = guidance.rule_span_sampler(template, conditions, budget=20, seed=20261013,
        spans=disk_spans, span_temperatures=disk_temperatures)
    assert rule_fresh == rule_cached
    assert any(pair is not None for pair in rule_fresh)
write_json(root / "results/checks/span_cache/check.json", {"checked_at": utc_now(), "exact_output_equality": True,
    "cache_sha256": sha256(root / "scripts/span_cache.py"), "array_and_nan_equality": True,
    "model_outputs": "Direct and diffusion candidate positions and every label; rule samples also identical",
    "nonempty_model_candidates": nonempty_counts,
    "fresh_span_seconds": fresh_seconds, "cached_span_seconds": cached_seconds,
    "implementation": str(Path(guidance.__file__).resolve())})
print("Fresh and disk-cached arrays and real model outputs are identical", flush=True)
