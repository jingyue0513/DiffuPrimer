import numpy as np

from diffuprimer.constraints import Conditions
from diffuprimer.decode import grammar_candidates, sample_grammar
from diffuprimer.sequence import PrimerPair


def test_decoder_enforces_orientation_range_and_known_labels():
    c = Conditions(primer_min_length=18, primer_max_length=20, amplicon_min=60, amplicon_max=100)
    logits = np.zeros((128, 3))
    previous = PrimerPair(5, 23, 70, 90)
    known = np.full(128, 3)
    known[10:15] = 1
    known[77:80] = 2
    known[30:60] = 0
    pairs, scores = grammar_candidates(logits, c, known=known, previous=previous)
    assert pairs
    for pair in pairs:
        clean = pair.labels(128)
        assert np.all((known == 3) | (known == clean))
        assert 60 <= pair.amplicon_length <= 100
        assert pair.forward_end <= pair.reverse_start


def test_fixed_primer_is_preserved_even_after_beam_changes():
    c = Conditions(amplicon_min=80, amplicon_max=120)
    rng = np.random.default_rng(2)
    logits = rng.normal(size=(150, 3))
    pair = sample_grammar(logits, c, rng, fixed_forward=(12, 33))
    assert (pair.forward_start, pair.forward_end) == (12, 33)


def test_impossible_condition_has_no_candidate():
    c = Conditions(amplicon_min=80, amplicon_max=100)
    pairs, _ = grammar_candidates(np.zeros((50, 3)), c)
    assert not pairs
