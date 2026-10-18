import numpy as np
import pytest

from diffuprimer.sequence import PrimerPair, labels_to_pair, map_primers, reverse_complement


def test_reverse_coordinates_and_roundtrip():
    template = "GGAACCTTGCAGTACCGTTACGATGCCTTAGCGATCTAGCTGGACCAT"
    pair = PrimerPair(2, 12, 30, 41)
    forward, reverse = pair.sequences(template)
    assert reverse == reverse_complement(template[30:41])
    assert map_primers(template, forward, reverse, 39) == [pair]
    assert labels_to_pair(pair.labels(len(template))) == pair


def test_ambiguous_locations_remain_explicit():
    assert len(map_primers("AACCGGAACCGGTTTTCCGG", "AACC", "CCGG")) > 1


def test_discontiguous_and_overlapping_masks_rejected():
    assert labels_to_pair([0, 1, 0, 1, 2, 2]) is None
    assert labels_to_pair([2, 2, 0, 1, 1]) is None
    assert labels_to_pair([0, 0, 0]) is None
    with pytest.raises(ValueError):
        PrimerPair(0, 5, 4, 7).validate("A" * 10)


def test_unknown_label_is_not_background():
    assert labels_to_pair(np.array([1, 1, 3, 2, 2])) is None
