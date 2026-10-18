from pathlib import Path

from diffuprimer.acquisition import parse_detail, parse_search, parse_validation
from diffuprimer.sequence import map_primers


ROOT = Path(__file__).resolve().parents[1]


def test_real_primerbank_parser_and_coding_coordinate_mapping():
    detail = ROOT / "data/raw/primerbank/6671509a1_detail.html"
    if not detail.exists():
        import pytest
        pytest.skip("Local acquisition calibration snapshot not present")
    records = parse_search(detail.read_text())
    assert len(records) == 1
    row = records[0]
    assert row["primerbank_id"] == "6671509a1"
    assert row["gene_id"] == "11461"
    template = parse_detail(detail.read_text())
    pair, = map_primers(template, row["forward_primer"], row["reverse_primer"], row["amplicon_length"])
    assert (pair.forward_start, pair.forward_end, pair.reverse_start, pair.reverse_end) == (83, 103, 215, 237)


def test_evidence_artifacts_do_not_assign_success():
    html = '<h1>QPCR Validation Results</h1><img src="/display.do?ampPlotId=12"><input name="qpcr" value="1">'
    evidence = parse_validation(html)
    assert evidence["experimental_record_available"]
    assert evidence["amplification_plot_available"]
    assert evidence["experimental_result"] == "unknown"
    assert evidence["evidence_type"] == "unverified"
