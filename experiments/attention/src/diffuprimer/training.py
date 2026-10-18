from __future__ import annotations

import random
from dataclasses import replace

import numpy as np
import torch
from torch.utils.data import Dataset

from .constraints import Conditions
from .model import MASK, encode_conditions, encode_dna
from .sequence import PrimerPair


class PrimerDataset(Dataset):
    def __init__(self, rows, window=512, seed=0, augment=True):
        self.window = window
        self.seed = seed
        self.augment = augment
        self.epoch = 0
        self.groups = {}
        for row in rows:
            pair = PrimerPair(*(row[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
            if pair.amplicon_length <= min(window, 400) and all(18 <= size <= 25 for size in
                [pair.forward_end - pair.forward_start, pair.reverse_end - pair.reverse_start]):
                self.groups.setdefault(row["gene_id"], []).append(row)
        self.genes = sorted(self.groups)

    def __len__(self):
        # Four independent conditional windows per gene per epoch; genes equally weighted.
        return 4 * len(self.genes)

    def __getitem__(self, index):
        rng = random.Random(self.seed + self.epoch * 1000003 + index)
        gene = self.genes[index % len(self.genes)]
        options = self.groups[gene]
        row = options[rng.randrange(len(options))] if self.augment else options[index // len(self.genes) % len(options)]
        template = row["template_sequence"]
        pair = PrimerPair(*(row[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
        low = max(0, pair.reverse_end - self.window)
        high = min(pair.forward_start, max(0, len(template) - self.window))
        offset = rng.randint(low, high) if self.augment else (low + high) // 2
        sequence = template[offset:offset + self.window]
        shifted = PrimerPair(*(x - offset for x in pair.__dict__.values()))
        labels = np.zeros(self.window, dtype=np.int64)
        labels[:len(sequence)] = shifted.labels(len(sequence))
        dna = np.full(self.window, 4, dtype=np.int64)
        dna[:len(sequence)] = encode_dna(sequence)
        valid = np.arange(self.window) < len(sequence)
        length = pair.amplicon_length
        # Conditions contain the observed pair during imitation training. Inference
        # uses user-specified intervals without access to source-pair coordinates.
        ranges = [(60, 120), (80, 200), (150, 300), (200, 400)]
        ranges = [(a, b) for a, b in ranges if a <= length <= b]
        amin, amax = rng.choice(ranges) if self.augment and ranges else (max(36, length - 25), min(400, length + 25))
        conditions = replace(Conditions(), amplicon_min=amin, amplicon_max=amax)
        fixed = np.zeros(self.window, dtype=bool)
        if self.augment and rng.random() < 0.25:
            side = rng.choice([1, 2])
            fixed = labels == side
        return {"dna": torch.from_numpy(dna), "clean": torch.from_numpy(labels),
                "valid": torch.from_numpy(valid), "fixed": torch.from_numpy(fixed),
                "conditions": torch.from_numpy(encode_conditions(conditions))}
