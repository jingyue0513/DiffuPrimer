from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def occurrences(sequence: str, pattern: str) -> list[int]:
    return [match.start() for match in re.finditer(f"(?={re.escape(pattern)})", sequence)]


@dataclass(frozen=True, order=True)
class PrimerPair:
    """Zero-based, half-open coordinates on the forward template strand."""

    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int

    def validate(self, sequence: str):
        if not (0 <= self.forward_start < self.forward_end <= self.reverse_start < self.reverse_end <= len(sequence)):
            raise ValueError("Expected nonoverlapping F then R intervals within the template")

    @property
    def amplicon_length(self):
        return self.reverse_end - self.forward_start

    def sequences(self, template):
        self.validate(template)
        return (template[self.forward_start:self.forward_end],
                reverse_complement(template[self.reverse_start:self.reverse_end]))

    def labels(self, length):
        self.validate("N" * length)
        labels = np.zeros(length, dtype=np.int64)
        labels[self.forward_start:self.forward_end] = 1
        labels[self.reverse_start:self.reverse_end] = 2
        return labels

    def as_dict(self, template):
        forward, reverse = self.sequences(template)
        return {**self.__dict__, "forward_primer": forward, "reverse_primer": reverse,
                "amplicon_length": self.amplicon_length}


def map_primers(template, forward, reverse, expected_length=None):
    template, forward, reverse = template.upper(), forward.upper(), reverse.upper()
    pairs = []
    for left in occurrences(template, forward):
        for right in occurrences(template, reverse_complement(reverse)):
            pair = PrimerPair(left, left + len(forward), right, right + len(reverse))
            if pair.forward_end > pair.reverse_start:
                continue
            if expected_length is None or pair.amplicon_length == int(expected_length):
                pairs.append(pair)
    return pairs


def labels_to_pair(labels):
    labels = np.asarray(labels)
    if np.any((labels < 0) | (labels > 2)):
        return None
    intervals = []
    for value in [1, 2]:
        positions = np.flatnonzero(labels == value)
        if len(positions) == 0 or np.any(np.diff(positions) != 1):
            return None
        intervals.append((int(positions[0]), int(positions[-1]) + 1))
    pair = PrimerPair(*intervals[0], *intervals[1])
    try:
        pair.validate("N" * len(labels))
    except ValueError:
        return None
    return pair
