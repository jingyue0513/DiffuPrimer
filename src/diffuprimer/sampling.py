from __future__ import annotations

import numpy as np
import torch

from .decode import sample_grammar
from .model import MASK, encode_conditions, encode_dna, reveal_probability
from .sequence import labels_to_pair


@torch.no_grad()
def generate(model, template, conditions, budget=20, steps=16, seed=0,
             mode="diffusion", fixed_forward=None, fixed_reverse=None,
             temperature=1.0, constrained=True, span_allowed=None):
    if mode not in {"diffusion", "direct"} or steps < 1 or budget < 1:
        raise ValueError("Invalid generation mode, step count, or budget")
    model.eval()
    device = next(model.parameters()).device
    rng = np.random.default_rng(seed)
    dna = torch.tensor(encode_dna(template), device=device)[None, :].expand(budget, -1)
    c = torch.tensor(encode_conditions(conditions), device=device)[None, :].expand(budget, -1)
    noisy = np.full((budget, len(template)), MASK, dtype=np.int64)
    if fixed_forward is not None:
        noisy[:, fixed_forward[0]:fixed_forward[1]] = 1
    if fixed_reverse is not None:
        noisy[:, fixed_reverse[0]:fixed_reverse[1]] = 2
    previous = [None] * budget
    failures = [False] * budget
    schedule = [1] if mode == "direct" else list(range(steps, 0, -1))
    for t in schedule:
        time = torch.full((budget,), t / (steps if mode == "diffusion" else 1), device=device)
        logits = model(dna, torch.tensor(noisy, device=device), time, c).cpu().numpy()
        for index in range(budget):
            if failures[index]:
                continue
            if constrained:
                pair = sample_grammar(logits[index], conditions, rng,
                    temperature=temperature, known=noisy[index], fixed_forward=fixed_forward,
                    fixed_reverse=fixed_reverse, previous=previous[index], span_allowed=span_allowed)
                if pair is None:
                    failures[index] = True
                    continue
                clean = pair.labels(len(template))
                previous[index] = pair
            else:
                probabilities = np.exp((logits[index] - logits[index].max(1, keepdims=True)) / temperature)
                probabilities /= probabilities.sum(1, keepdims=True)
                uniforms = rng.random(len(template))
                clean = (uniforms[:, None] > probabilities.cumsum(1)).sum(1).clip(max=2)
            reveal = (noisy[index] == MASK) & (rng.random(len(template)) < reveal_probability(t))
            noisy[index, reveal] = clean[reveal]
    output = []
    for index in range(budget):
        pair = None if failures[index] else labels_to_pair(noisy[index])
        output.append({"proposal_index": index, "pair": pair,
                       "structural_valid": pair is not None,
                       "labels": noisy[index], "mode": mode,
                       "constrained_decoder": constrained})
    return output
