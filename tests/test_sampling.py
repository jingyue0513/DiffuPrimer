import numpy as np
import torch

from diffuprimer.constraints import Conditions
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.sampling import generate


def test_reverse_sampling_preserves_an_exact_fixed_span_and_seed():
    torch.set_num_threads(2)
    torch.manual_seed(4)
    model = ConditionalDenoiser(width=16, layers=2)
    sequence = "ACGT" * 40
    conditions = Conditions(amplicon_min=80, amplicon_max=120)
    first = generate(model, sequence, conditions, budget=3, steps=4, seed=19, fixed_forward=(5, 25))
    second = generate(model, sequence, conditions, budget=3, steps=4, seed=19, fixed_forward=(5, 25))
    for left, right in zip(first, second):
        assert left["pair"] == right["pair"]
        assert np.array_equal(left["labels"], right["labels"])
        assert (left["pair"].forward_start, left["pair"].forward_end) == (5, 25)
        assert 80 <= left["pair"].amplicon_length <= 120
