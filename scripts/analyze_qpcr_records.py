#!/usr/bin/env python3
"""Analyze measured blind assay records and dilution curves."""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import linregress

from diffuprimer.io import sha256, utc_now, write_json, write_jsonl


def number(value):
    value = str(value or "").strip()
    if value.lower() in {"", "na", "n/a", "undetermined", "not detected"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite measurement")
    return result


def read_records(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    wells = set()
    for row in rows:
        identity = row["plate"], row["well"]
        if identity in wells:
            raise ValueError(f"Repeated plate/well in {path}: {identity}")
        wells.add(identity)
        row["cq_value"] = number(row.get("cq"))
        if row["cq_value"] is not None and row.get("assay_status") != "conducted":
            raise ValueError("A recorded Cq requires assay_status=conducted")
        if row.get("assay_status") == "conducted" and not row.get("raw_data_file", "").strip():
            raise ValueError("Conducted reactions require their instrument export filename")
        if row["cq_value"] is not None and row["cq_value"] <= 0:
            raise ValueError("Recorded Cq must be positive")
    return rows


def analyze_assay(sample_rows, curve_rows):
    observed = [r for r in sample_rows if r.get("assay_status") == "conducted"]
    result = {"conducted_reactions": len(observed), "reaction_groups": {}, "standard_curve": None}
    by_type = defaultdict(list)
    for row in observed:
        by_type[(row["reaction_type"], row.get("sample_id", ""))].append(row)
    sample_counts = {reaction_type: sum(kind == reaction_type for kind, _ in by_type) for reaction_type, _ in by_type}
    for (reaction_type, sample_id), own in sorted(by_type.items()):
        cq = [r["cq_value"] for r in own if r["cq_value"] is not None]
        group_key = reaction_type if sample_counts[reaction_type] == 1 else f"{reaction_type}:{sample_id}"
        result["reaction_groups"][group_key] = {
            "reaction_type": reaction_type, "sample_id": sample_id,
            "conducted": len(own), "numeric_cq": len(cq),
            "mean_cq": float(np.mean(cq)) if cq else None,
            "sd_cq": float(np.std(cq, ddof=1)) if len(cq) >= 2 else None,
            "cq_values": cq,
            "recorded_control_outcomes": [r.get("control_outcome", "") for r in own],
            "recorded_melting_single_peak": [r.get("melting_single_peak", "") for r in own],
            "recorded_gel_single_expected_band": [r.get("gel_single_expected_band", "") for r in own],
            "instrument_exports": sorted({r["raw_data_file"] for r in own})}
    measured = [r for r in curve_rows if r.get("assay_status") == "conducted"]
    if not measured:
        return result
    quantities, units, levels = {}, set(), defaultdict(list)
    for row in measured:
        quantity = number(row.get("input_quantity"))
        unit = row.get("input_unit", "").strip()
        if quantity is None or quantity <= 0 or not unit:
            raise ValueError("Standard curves require positive actual quantities and a unit")
        level = row["dilution_level"]
        if level in quantities and quantities[level] != quantity:
            raise ValueError("Technical replicates at one dilution level have different input quantities")
        quantities[level] = quantity
        units.add(unit)
        levels[level].append(row)
    if len(units) != 1:
        raise ValueError("Standard-curve input quantities use different units")
    points = []
    for level, own in sorted(levels.items(), key=lambda item: quantities[item[0]]):
        cq = [r["cq_value"] for r in own if r["cq_value"] is not None]
        points.append({"dilution_level": level, "input_quantity": quantities[level], "conducted": len(own),
                       "numeric_cq": len(cq), "mean_cq": float(np.mean(cq)) if cq else None,
                       "sd_cq": float(np.std(cq, ddof=1)) if len(cq) >= 2 else None})
    valid = [p for p in points if p["mean_cq"] is not None]
    curve = {"input_unit": next(iter(units)), "levels": points, "conducted_reactions": len(measured),
             "regression": "OLS of mean Cq at each measured level on log10 actual input quantity",
             "fit_status": "insufficient_distinct_quantities", "slope": None, "intercept": None,
             "r_squared": None, "efficiency_percent": None,
             "instrument_exports": sorted({r["raw_data_file"] for r in measured})}
    if len({p["input_quantity"] for p in valid}) >= 3:
        fit = linregress([math.log10(p["input_quantity"]) for p in valid], [p["mean_cq"] for p in valid])
        curve.update(slope=float(fit.slope), intercept=float(fit.intercept),
                     r_squared=float(fit.rvalue ** 2) if math.isfinite(fit.rvalue) else None,
                     slope_standard_error=float(fit.stderr) if math.isfinite(fit.stderr) else None,
                     fit_status="fitted" if fit.slope < 0 else "nonnegative_slope")
        # Cq = intercept + slope * log10(quantity); one cycle multiplies
        # quantity by 10**(-1/slope). Efficiency is that factor minus one.
        if fit.slope < 0:
            exponent = -1 / fit.slope
            curve["efficiency_percent"] = float(math.expm1(exponent * math.log(10)) * 100) if exponent < 300 else None
    result["standard_curve"] = curve
    return result


def main():
    parser = argparse.ArgumentParser(description="Analyze instrument measurements; measured rows use assay_status=conducted and raw_data_file.")
    parser.add_argument("--reactions", required=True)
    parser.add_argument("--curves")
    parser.add_argument("--output", default="results/qpcr_observations")
    args = parser.parse_args()
    reactions_path = Path(args.reactions)
    curves_path = Path(args.curves) if args.curves else None
    reactions = read_records(reactions_path)
    curves = read_records(curves_path) if curves_path else []
    grouped, curve_grouped = defaultdict(list), defaultdict(list)
    for row in reactions:
        grouped[row["blind_id"]].append(row)
    for row in curves:
        curve_grouped[row["blind_id"]].append(row)
    if set(curve_grouped) - set(grouped):
        raise ValueError("Standard curve contains an assay absent from reaction records")
    results = [{"blind_id": blind_id, **analyze_assay(grouped[blind_id], curve_grouped[blind_id])}
               for blind_id in sorted(grouped)]
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    write_jsonl(output / "per_assay.jsonl", results)
    conducted = sum(r["conducted_reactions"] for r in results)
    curves_conducted = sum(r["standard_curve"]["conducted_reactions"] for r in results if r["standard_curve"])
    write_json(output / "status.json", {"at": utc_now(), "stage": "measured_records_analyzed" if conducted or curves_conducted else "awaiting_measured_records",
                                       "planned_assays_in_records": len(results), "conducted_reactions": conducted,
                                       "conducted_curve_reactions": curves_conducted,
                                       "reaction_records_sha256": sha256(reactions_path),
                                       "curve_records_sha256": sha256(curves_path) if curves_path else None,
                                       "per_assay_sha256": sha256(output / "per_assay.jsonl")})
    print(output / "status.json")


if __name__ == "__main__":
    main()
