from pathlib import Path

from diffuprimer.acquisition import parse_detail
from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.sequence import PrimerPair


def test_real_template_primer3_roundtrip_and_shared_constraints():
    path = Path(__file__).resolve().parents[1] / "data/raw/primerbank/6671509a1_detail.html"
    if not path.exists():
        import pytest
        pytest.skip("Local source snapshot not present")
    template = parse_detail(path.read_text())
    c = Conditions()
    pairs, explain = primer3_candidates(template, c, budget=5)
    assert len(pairs) > 0
    for row in pairs:
        pair = PrimerPair(*(row[key] for key in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
        assert pair.as_dict(template)["reverse_primer"] == row["reverse_primer"]
        report = assess_pair(template, pair, c)
        assert report["computational_feasible"], report["failures"]
        assert report["experimental_status"] == "not_tested"
        assert report["specificity_status"] == "not_assessed"
