#!/usr/bin/env python3
"""Publication figure for the mathematical corruption process, not a model result."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np

from diffuprimer.io import read_jsonl, write_json, sha256
from diffuprimer.sequence import PrimerPair


root = Path(__file__).resolve().parents[1]
rows = read_jsonl(root / "data/processed/dataset_partial_pilot.jsonl")
row = next(row for row in rows if row["source_database"] == "PrimerBank" and row["amplicon_length"] <= 300)
pair = PrimerPair(*(row[key] for key in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
offset = max(0, pair.forward_start - 40)
end = min(len(row["template_sequence"]), pair.reverse_end + 40)
clean = pair.labels(len(row["template_sequence"]))[offset:end]
rng = np.random.default_rng(20261013)
uniform = rng.random(len(clean))
times = [0, 4, 8, 12, 16]
masked = np.stack([np.where(uniform < t / 16, 3, clean) for t in times])
colors = ["#f4f0e8", "#df8039", "#147f9d", "#a3aab1"]
fig, axes = plt.subplots(1, 2, figsize=(10, 3.2), gridspec_kw={"width_ratios": [2.2, 1]}, constrained_layout=True)
axes[0].imshow(masked, cmap=ListedColormap(colors), vmin=0, vmax=3, aspect="auto", interpolation="nearest",
               extent=[offset, end, 4.5, -0.5])
axes[0].set_yticks(range(len(times)), [f"t = {t}" for t in times])
axes[0].set_xlabel("Transcript position (zero-based)")
axes[0].set_title(f"Absorbing corruption illustration; PrimerBank {row['primerbank_id']}", fontsize=9)
axes[0].legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(colors, ["background", "F", "R", "MASK"])],
               loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=4, frameon=False, fontsize=8)
t = np.arange(1, 17)
axes[1].plot(np.arange(17), 1 - np.arange(17) / 16, color="#147f9d", label=r"Forward survival $\bar\alpha_t$")
axes[1].plot(t, 1 / t, color="#df8039", label="Backward reveal probability 1/t")
axes[1].set_xlabel("Step t")
axes[1].set_ylabel("Probability")
axes[1].set_ylim(-0.03, 1.03)
axes[1].legend(frameon=False, fontsize=7)
axes[1].spines[["top", "right"]].set_visible(False)
directory = root / "paper/figures"
fig.savefig(directory / "forward_process.pdf", bbox_inches="tight")
fig.savefig(directory / "forward_process.png", dpi=220, bbox_inches="tight")
plt.close(fig)
write_json(directory / "forward_process.provenance.json", {"type": "mathematical illustration, not performance evidence",
    "primerbank_id": row["primerbank_id"], "reference_accession": row["transcript_accession_version"],
    "offset": offset, "end": end, "seed": 20261013, "steps": 16,
    "dataset_sha256": sha256(root / "data/processed/dataset_partial_pilot.jsonl")})
