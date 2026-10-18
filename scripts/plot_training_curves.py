#!/usr/bin/env python3
"""Inspect actual training histories and export separate objective curves."""
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
output = root / "results/figures/training"
output.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
rows, runs, inputs = [], [], []
colors = ["#3d7aa8", "#37926c", "#855ab0"]
for axis, method in zip(axes, ["direct", "diffusion"]):
    for index, seed in enumerate([20261013, 20261014, 20261015]):
        source = root / f"results/{method}_weak_partial_pilot_seed{seed}_training.json"
        record = json.loads(source.read_text())
        history = record["history"]
        if len(history) != 30 or [r["epoch"] for r in history] != list(range(1, 31)):
            raise ValueError("Incomplete primary training history")
        continuation = root / f"experiments/epoch100/results/{method}_weak_partial_pilot_seed{seed}_training.json"
        if continuation.exists():
            extended = json.loads(continuation.read_text())["history"]
            if len(extended) != 100 or any(any(a[key] != b[key] for key in ["epoch", "training_loss", "validation_loss"])
                                           for a, b in zip(history, extended[:30])):
                raise ValueError("Training duration comparison does not share its original optimization path")
            history = extended
            inputs.append({"path": str(continuation), "sha256": sha256(continuation)})
        inputs.append({"path": str(source), "sha256": sha256(source)})
        for row in history:
            if not all(math.isfinite(row[key]) for key in ["training_loss", "validation_loss"]):
                raise ValueError("Nonfinite training objective")
            rows.append({"method": method, "seed": seed, **row})
        epochs = [r["epoch"] for r in history]
        axis.plot(epochs, [r["validation_loss"] for r in history], color=colors[index], label=f"Seed {seed}")
        axis.plot(epochs, [r["training_loss"] for r in history], color=colors[index], alpha=.55, linestyle="--")
        best = min(history, key=lambda r: r["validation_loss"])
        axis.scatter([best["epoch"]], [best["validation_loss"]], color=colors[index], s=24, zorder=4)
        runs.append({"method": method, "seed": seed, "completed_epochs": len(history),
                     "best_validation_epoch": best["epoch"], "best_validation_loss": best["validation_loss"],
                     "last_validation_loss": history[-1]["validation_loss"]})
    axis.axvline(30, color="#999999", linewidth=.8, linestyle=":")
    axis.set(title=method.capitalize(), xlabel="Epoch", ylabel="Cross-entropy" if method == "direct" else "Weighted denoising loss")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(fontsize=7)
fig.legend(handles=[Line2D([], [], color="#555555", label="Validation"),
                    Line2D([], [], color="#555555", linestyle="--", label="Training"),
                    Line2D([], [], color="#555555", marker="o", linestyle="none", label="Best validation epoch")],
           loc="upper center", bbox_to_anchor=(.5, 1.12), ncol=3, frameon=False, fontsize=9)
fig.savefig(output / "training_curves.pdf", bbox_inches="tight")
fig.savefig(output / "training_curves.png", dpi=200, bbox_inches="tight")
plt.close(fig)
with (output / "training_history.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["method", "seed", "epoch", "training_loss", "validation_loss", "elapsed_seconds"])
    writer.writeheader()
    writer.writerows(rows)
write_json(output / "training_checks.json", {"at": utc_now(), "runs": runs, "inputs": inputs,
                                            "losses_finite": True, "hundred_epoch_prefix_equality": True,
                                            "objective_comparison": "Within-model trajectories; direct and diffusion objectives have different scales"})
print(output)
