#!/usr/bin/env python3
"""Inspect completed training histories and duration sensitivity."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
records = []
for mode, ax in zip(["direct", "diffusion"], axes):
    histories = []
    for seed in [20261013, 20261014, 20261015]:
        path = root / f"results/{mode}_weak_partial_pilot_seed{seed}_training.json"
        history = json.loads(path.read_text())["history"]
        assert len(history) == 30 and all(np.isfinite(r["validation_loss"]) and np.isfinite(r["training_loss"]) for r in history)
        checkpoint = root / f"results/checkpoints/{mode}_weak_partial_pilot_seed{seed}.pt"
        selected = torch.load(checkpoint, weights_only=True)["epoch"]
        assert selected == min(history, key=lambda row: row["validation_loss"])["epoch"]
        records.append({"mode": mode, "seed": seed, "epochs": 30, "selected_epoch": selected,
            "history_sha256": sha256(path), "checkpoint_sha256": sha256(checkpoint),
            "selected_validation_loss": history[selected-1]["validation_loss"]})
        histories.append(history)
        ax.plot([r["epoch"] for r in history], [r["validation_loss"] for r in history], color="#9fb3c8", alpha=0.65, lw=1)
    path = root / f"experiments/epoch100/results/{mode}_weak_partial_pilot_seed20261013_training.json"
    history = json.loads(path.read_text())["history"]
    assert len(history) == 100
    assert all(a[key] == b[key] for a,b in zip(histories[0], history[:30]) for key in ["training_loss", "validation_loss"])
    checkpoint = root / f"experiments/epoch100/results/checkpoints/{mode}_weak_partial_pilot_seed20261013.pt"
    selected = torch.load(checkpoint, weights_only=True)["epoch"]
    assert selected == min(history, key=lambda row: row["validation_loss"])["epoch"]
    epochs = [r["epoch"] for r in history]
    ax.plot(epochs, [r["training_loss"] for r in history], color="#6a9d57", lw=1.4, label="Training, 100 epochs")
    ax.plot(epochs, [r["validation_loss"] for r in history], color="#167f9d", lw=1.6, label="Validation, 100 epochs")
    ax.scatter([selected], [history[selected-1]["validation_loss"]], color="#e49b2d", zorder=4, label=f"Selected epoch {selected}")
    ax.axvline(30, ls="--", color="#666666", lw=1)
    ax.set(title=mode.capitalize(), xlabel="Epoch", ylabel="Position-imitation loss")
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(fontsize=7)
    records.append({"mode": mode, "seed": 20261013, "epochs": 100, "selected_epoch": selected,
        "history_sha256": sha256(path), "checkpoint_sha256": sha256(checkpoint),
        "first_30_loss_epochs_exactly_match_primary": True, "selected_validation_loss": history[selected-1]["validation_loss"]})
directory = root / "results/figures"
directory.mkdir(parents=True, exist_ok=True)
fig.savefig(directory / "training_adequacy.pdf")
fig.savefig(directory / "training_adequacy.png", dpi=220)
plt.close(fig)
write_json(root / "results/training_adequacy.json", {"analyzed_at": utc_now(), "runs": records,
    "policy": "Six primary training runs verified complete. Same-seed duration follow-ups preserve the first 30 losses exactly. Best checkpoints selected by validation loss; losses are not compared across objectives."})
print(json.dumps(records, indent=2))
