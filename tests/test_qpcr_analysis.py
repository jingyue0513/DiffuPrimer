"""Numerical/data-import checks; fixtures are not biological observations."""
import csv
import importlib.util
import math
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location("qpcr_analysis", Path(__file__).resolve().parents[1] / "scripts/analyze_qpcr_records.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def curve_rows(slope=-math.log2(10)):
    return [{"assay_status": "conducted", "raw_data_file": "numerical_fixture", "dilution_level": str(level),
             "input_quantity": 10 ** level, "input_unit": "fixture_units", "cq_value": 35 + slope * level + offset}
            for level in range(5) for offset in [-.1, 0, .1]]


def test_known_doubling_efficiency_and_repeatability():
    samples = [{"assay_status": "conducted", "raw_data_file": "numerical_fixture", "reaction_type": "sample", "cq_value": cq}
               for cq in [20, 20.1, 19.9]]
    result = analysis.analyze_assay(samples, curve_rows())
    assert result["reaction_groups"]["sample"]["mean_cq"] == pytest.approx(20)
    assert result["reaction_groups"]["sample"]["sd_cq"] == pytest.approx(.1)
    assert result["standard_curve"]["slope"] == pytest.approx(-math.log2(10))
    assert result["standard_curve"]["efficiency_percent"] == pytest.approx(100)
    assert result["standard_curve"]["r_squared"] == pytest.approx(1)


@pytest.mark.parametrize("slope", [0, 1])
def test_nonnegative_standard_curve_has_no_efficiency(slope):
    curve = analysis.analyze_assay([], curve_rows(slope))["standard_curve"]
    assert curve["fit_status"] == "nonnegative_slope"
    assert curve["efficiency_percent"] is None
    assert curve["r_squared"] is None or math.isfinite(curve["r_squared"])


def test_empty_records_and_unmeasured_cq():
    result = analysis.analyze_assay([{"assay_status": "not_conducted", "cq_value": None}], [])
    assert result["conducted_reactions"] == 0
    assert result["standard_curve"] is None
    assert analysis.number("Undetermined") is None
    assert analysis.number("0") == 0


def test_mixed_quantity_units_rejected():
    rows = curve_rows()
    rows[0]["input_unit"] = "another_unit"
    with pytest.raises(ValueError, match="different units"):
        analysis.analyze_assay([], rows)


def test_biological_samples_are_not_pooled_as_technical_replicates():
    samples = [{"assay_status": "conducted", "raw_data_file": "numerical_fixture", "reaction_type": "sample",
                "sample_id": source, "cq_value": cq} for source, cq in [("A", 20), ("A", 20.2), ("B", 30), ("B", 30.2)]]
    groups = analysis.analyze_assay(samples, [])["reaction_groups"]
    assert set(groups) == {"sample:A", "sample:B"}
    assert groups["sample:A"]["mean_cq"] == pytest.approx(20.1)
    assert groups["sample:B"]["mean_cq"] == pytest.approx(30.1)
    assert groups["sample:A"]["sd_cq"] == pytest.approx(math.sqrt(.02))


def test_duplicate_wells_and_unconducted_measurement_rejected(tmp_path):
    path = tmp_path / "instrument_records.csv"
    fields = ["plate", "well", "assay_status", "cq", "raw_data_file"]
    row = dict(zip(fields, [1, "A01", "conducted", 20, "numerical_fixture"]))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows([row, row])
    with pytest.raises(ValueError, match="Repeated plate/well"):
        analysis.read_records(path)
    row["assay_status"] = "not_conducted"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerow(row)
    with pytest.raises(ValueError, match="requires assay_status=conducted"):
        analysis.read_records(path)
