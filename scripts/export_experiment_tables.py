#!/usr/bin/env python3
"""Export the completed statistical groups and all contrasts as plain CSV."""
import csv
import json
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
output = root / "results/tables"
output.mkdir(parents=True, exist_ok=True)
sources = {
    "selected_primary": root / "results/temperature_tuning/selected_test_summary.json",
    "selected_conditions": root / "results/temperature_tuning/selected_conditions_summary.json",
    "guidance": root / "results/guide_ablation_summary.json",
    "training_duration": root / "experiments/epoch100/results/duration_reference_comparison.json",
    "decoder_ablations": root / "results/complete_secondary_summary.json",
    "exon_association": root / "results/exon_specificity_association.json",
}
artifacts = {}
for label, path in sources.items():
    summary = json.loads(path.read_text())
    for kind, key in [("groups", "groups" if label == "exon_association" else "summaries"),
                      ("contrasts", "paired_comparisons")]:
        rows = summary[key]
        if not rows:
            raise ValueError("Missing completed table rows: " + label + "/" + kind)
        fields = sorted({key for row in rows for key in row})
        destination = output / f"{label}_{kind}.csv"
        with destination.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value
                              for key, value in row.items()} for row in rows)
        artifacts[destination.name] = {"rows": len(rows), "sha256": sha256(destination)}
write_json(output / "provenance.json", {"at": utc_now(), "source_summaries": {str(p): sha256(p) for p in sources.values()},
    "artifacts": artifacts, "policy": "All completed groups and reported contrasts retained; gene and seed denominators included in original statistical fields"})
print(json.dumps(artifacts, indent=2))
