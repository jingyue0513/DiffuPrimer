#!/usr/bin/env python3
"""Publication vector overview of the implemented design experiment."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from diffuprimer.io import utc_now, write_json


root = Path(__file__).resolve().parents[1]
destination = root / "paper/figures"
fig, ax = plt.subplots(figsize=(11.0, 5.0))
ax.set_xlim(0, 11)
ax.set_ylim(0, 5)
ax.axis("off")
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})


def box(x, y, w, h, title, lines, color):
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.08",
        edgecolor=color, facecolor=color + "12", linewidth=1.3)
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h - 0.22, title, ha="center", va="top", fontsize=11, fontweight="bold", color=color)
    ax.text(x + w / 2, y + h / 2 - 0.12, "\n".join(lines), ha="center", va="center", fontsize=9, linespacing=1.45)


def arrow(start, end, label=None, **kwargs):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13,
        color="#444444", linewidth=1.2, **kwargs))
    if label:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.1, label,
            ha="center", va="bottom", fontsize=8)


box(0.15, 3.0, 2.5, 1.65, "Template and conditions", ["Fixed mouse transcript", "Product length, primer length, Tm, GC", "Optional fixed forward/reverse primer"], "#246c91")
box(3.25, 3.0, 3.0, 1.65, "Conditional denoiser", ["DNA + noisy labels + conditions + time", "96-channel convolutional network", "Predict clean 0 / F / R logits"], "#246c91")
box(7.0, 3.0, 3.75, 1.65, "Constrained pair proposal", ["Two ordered binding intervals", "Single-oligo thermodynamic guidance", "Retain revealed labels and fixed primer"], "#468448")
arrow((2.73, 3.85), (3.15, 3.85))
arrow((6.33, 3.85), (6.9, 3.85))
arrow((8.85, 2.9), (4.75, 2.9), "Reveal masked labels with probability 1/t", connectionstyle="arc3,rad=-0.15")
arrow((4.75, 2.89), (4.75, 3.0))
ax.text(6.8, 2.3, "16 reverse steps; 20 candidate proposals", fontsize=9, ha="center", color="#333333")
box(0.15, 0.25, 2.5, 1.55, "Gene-separated learning", ["350 train / 75 validation genes", "Source positions + computed pairs", "3 training and sampling seeds"], "#246c91")
box(3.25, 0.25, 3.0, 1.55, "Held-out experiments", ["75 test genes; paired target analysis", "Primer3 / direct CNN / rule sampler", "Length conditions and completion"], "#a3642f")
box(7.0, 0.25, 3.75, 1.55, "Shared candidate assessment", ["Pair geometry and thermodynamics", "RNA and genomic BLAST screens", "Unique pairs and amplicon overlap"], "#468448")
arrow((8.85, 2.88), (8.85, 1.9))
arrow((6.9, 1.05), (6.34, 1.05))
arrow((1.4, 1.9), (3.2, 3.0))
fig.tight_layout(pad=0.2)
fig.savefig(destination / "method_overview.pdf")
fig.savefig(destination / "method_overview.png", dpi=220)
plt.close(fig)
write_json(destination / "method_overview.provenance.json", {"generated_at": utc_now(),
    "figure_kind": "schematic of implemented model and frozen evaluation protocol",
    "script": "scripts/plot_overview.py"})
