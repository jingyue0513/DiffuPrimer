#!/usr/bin/env python3
"""Package the final paper, figures, statistical records, and source outputs."""
import json
import zipfile
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json

ROOT = Path(__file__).resolve().parents[1]
summary = json.loads((ROOT / "results/final_fixed_seed/study_summary.json").read_text())
assert json.loads((ROOT / "results/final_fixed_seed/completion_check.json").read_text())["passed"]
files = set()


def include(directory, suffixes):
    files.update(p for p in (ROOT / directory).rglob("*") if p.is_file() and p.suffix in suffixes and "__pycache__" not in p.parts)


for name in ["manuscript.tex", "manuscript.pdf", "supplementary.tex", "supplementary.pdf", "final_values.tex", "reported_values.json", "build_manifest.json"]:
    files.add(ROOT / "paper" / name)
for directory, suffixes in [
    ("results/final_fixed_seed", {".json", ".jsonl", ".csv", ".pdf", ".png"}),
    ("results/tables", {".csv"}), ("results/figures/training", {".pdf"}),
    ("results/figures/selected_study", {".pdf"}), ("scripts", {".py", ".cjs"}),
    ("src", {".py"}), ("configs", {".json", ".yaml"}), ("data/splits", {".json", ".jsonl"}),
    ("experiments/final_followup", {".py", ".json", ".jsonl"}),
    ("experiments/joint_tm/src", {".py"}), ("experiments/uniform_pairs/src", {".py"}),
    ("experiments/joint_tm/scripts", {".py"}), ("experiments/uniform_pairs/scripts", {".py"}),
    ("experiments/joint_tm/configs", {".json"}), ("experiments/uniform_pairs/configs", {".json"}),
]:
    include(directory, suffixes)
for source in summary["inputs"]:
    files.update(p for p in Path(source["directory"]).iterdir() if p.suffix in {".json", ".jsonl"})
for name in ["targets_pilot.jsonl", "dataset_partial_pilot.jsonl", "dataset_partial_pilot.manifest.json", "reference_manifest.json", "primerbank_pilot.jsonl"]:
    files.add(ROOT / "data/processed" / name)
for method in ["diffusion", "direct"]:
    files.add(ROOT / f"results/checkpoints/{method}_weak_partial_pilot_seed20261013.pt")
for name in ["README.md", "pyproject.toml", "requirements.lock", "docs/reproduction.md", "docs/data_dictionary.md",
             "results/current_research_status.json", "results/current_research_scope.json", "results/remaining_work.json",
             "results/frozen_protocol.json", "results/computational_completion_check.json",
             "results/temperature_tuning/selected_conditions_summary.json", "results/selected_exon_context_audit.per_candidate.jsonl"]:
    files.add(ROOT / name)
files = sorted(p for p in files if p.is_file())
readme = """DiffuPrimer — formal paper and experimental records, 13 October 2026

paper/manuscript.pdf: main manuscript (11 pages)
paper/supplementary.pdf: supplementary material (5 pages)
results/final_fixed_seed/: full fixed-seed statistics, gene records, final figures and checks
experiments/.../results/evaluations/: actual candidate and assessment files used in final analysis
results/checkpoints/: original master-seed model checkpoints
PACKAGE_SHA256.json: archived-file checksums

The paper uses seed 20261013 and eight experimental questions. New follow-up runs: 18.
Full workspace verification: 264 completed generation/reference-audit runs.
Final analysis: 47 run sources, 235 statistical groups, 2599 paired endpoint comparisons.

To compile the supplied TeX directly, work in paper/ and run pdflatex twice on each
of manuscript.tex and supplementary.tex. Figures and value macros are included.
The analysis provenance retains original absolute workspace paths. The archive
contains the final-analysis outputs and follow-up audit records; full reference
databases and other historical raw outputs remain in the research workspace.
Reference identifiers, official download URLs and checksums accompany the data.
See docs/reproduction.md for the original workspace commands.
"""
archive = ROOT / "paper/DiffuPrimer_paper_and_results.zip"
hashes = {str(p.relative_to(ROOT)): sha256(p) for p in files}
with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in files:
        z.write(p, str(p.relative_to(ROOT)))
    z.writestr("PACKAGE_README.txt", readme)
    z.writestr("PACKAGE_SHA256.json", json.dumps(hashes, indent=2) + "\n")
web_archive = ROOT / "reports/DiffuPrimer_中文研究进展网页.zip"
with zipfile.ZipFile(web_archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in sorted((ROOT / "reports/chinese_research").rglob("*")):
        if p.is_file():
            z.write(p, str(p.relative_to(ROOT / "reports")))
write_json(ROOT / "paper/package_manifest.json", {"created_at": utc_now(), "status": "complete", "paper_archive": str(archive),
    "paper_archive_sha256": sha256(archive), "paper_archive_bytes": archive.stat().st_size, "files": len(files),
    "web_archive": str(web_archive), "web_archive_sha256": sha256(web_archive), "web_archive_bytes": web_archive.stat().st_size})
print(archive)
print(web_archive)
