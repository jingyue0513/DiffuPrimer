from __future__ import annotations

import math
from dataclasses import asdict

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from .constraints import Conditions


MASK = 3
DNA_CODES = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
CONDITION_FIELDS = ["primer_min_length", "primer_max_length", "amplicon_min", "amplicon_max",
                    "tm_min", "tm_max", "max_tm_difference", "gc_min", "gc_max"]
CONDITION_SCALES = [25, 25, 350, 350, 65, 65, 5, 100, 100]


def encode_dna(sequence):
    return np.array([DNA_CODES.get(base, 4) for base in sequence.upper()], dtype=np.int64)


def encode_conditions(conditions: Conditions):
    return np.array([getattr(conditions, key) / scale for key, scale in zip(CONDITION_FIELDS, CONDITION_SCALES)], dtype=np.float32)


class ResidualBlock(nn.Module):
    def __init__(self, width, dilation):
        super().__init__()
        self.norm = nn.GroupNorm(8, width)
        self.convolution = nn.Conv1d(width, width, 5, padding=2 * dilation, dilation=dilation)
        self.out = nn.Conv1d(width, width, 1)

    def forward(self, x):
        return x + self.out(F.silu(self.convolution(F.silu(self.norm(x)))))


class ConditionalDenoiser(nn.Module):
    """Fully convolutional x0 predictor. No pretrained weights or sequence generation."""

    def __init__(self, width=96, layers=4, condition_dim=len(CONDITION_FIELDS)):
        super().__init__()
        if width % 8:
            raise ValueError("Width must be a multiple of eight")
        self.spec = {"width": width, "layers": layers, "condition_dim": condition_dim}
        self.dna = nn.Embedding(5, width)
        self.labels = nn.Embedding(4, width)
        self.condition = nn.Sequential(nn.Linear(condition_dim + 4, width), nn.SiLU(), nn.Linear(width, width))
        self.position = nn.Conv1d(2, width, 1)
        self.blocks = nn.ModuleList([ResidualBlock(width, 4 ** i) for i in range(layers)])
        self.head = nn.Sequential(nn.GroupNorm(8, width), nn.SiLU(), nn.Conv1d(width, 3, 1))

    def forward(self, dna, noisy_labels, time_fraction, conditions, valid=None):
        batch, length = dna.shape
        time_fraction = time_fraction.reshape(batch, 1)
        time = torch.cat([time_fraction, time_fraction ** 2,
                          torch.sin(math.pi * time_fraction), torch.cos(math.pi * time_fraction)], dim=1)
        embedding = self.condition(torch.cat([conditions, time], dim=1))
        if valid is None:
            valid = torch.ones_like(dna, dtype=torch.bool)
        actual_lengths = valid.sum(1).clamp(min=1).float()
        relative = torch.arange(length, device=dna.device).float()[None, :] / actual_lengths[:, None]
        position = torch.stack([relative.expand(batch, -1), torch.sin(2 * math.pi * relative)], dim=1)
        x = (self.dna(dna) + self.labels(noisy_labels)).transpose(1, 2)
        x = x + embedding[:, :, None] + self.position(position)
        x = x * valid[:, None, :]
        for block in self.blocks:
            x = block(x) * valid[:, None, :]
        return self.head(x).transpose(1, 2)


def corrupt(labels, t, steps, valid=None, generator=None, fixed=None):
    """Linear absorbing survival alpha_bar(t)=1-t/T."""
    if steps < 1 or torch.any((t < 1) | (t > steps)):
        raise ValueError("t must be in 1..T")
    mask = torch.rand(labels.shape, generator=generator, device=labels.device) < t[:, None] / steps
    if valid is not None:
        mask &= valid
    if fixed is not None:
        mask &= ~fixed
    return torch.where(mask, MASK, labels), mask


def diffusion_loss(logits, clean, masked, t, steps, valid):
    # T/t * average over ALL valid sites of masked-site CE is an unbiased
    # uniformly sampled discrete-time absorbing variational loss, up to constants.
    ce = F.cross_entropy(logits.transpose(1, 2), clean, reduction="none")
    per_example = (ce * masked * valid).sum(1) / valid.sum(1).clamp(min=1)
    return (per_example * steps / t).mean()


def reveal_probability(t):
    if t < 1:
        raise ValueError("t must be positive")
    return 1.0 / t
