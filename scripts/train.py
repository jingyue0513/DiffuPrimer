#!/usr/bin/env python3
import argparse
import json
import math
import random
import time
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F
from torch.utils.data import DataLoader

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.model import ConditionalDenoiser, MASK, corrupt, diffusion_loss
from diffuprimer.training import PrimerDataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["diffusion", "direct"], required=True)
    parser.add_argument("--dataset", default="data/processed/dataset_pilot.jsonl")
    parser.add_argument("--seed", type=int, default=20261013)
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--allow-weak", action="store_true", help="Imitate unverified database positions; never label them experimental successes")
    parser.add_argument("--calibration", action="store_true")
    parser.add_argument("--partial-source", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "configs/pilot.json").read_text())
    settings = config["model"]
    torch.set_num_threads(settings["cpu_threads"])
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    torch.use_deterministic_algorithms(True)
    path = root / args.dataset
    manifest = json.loads(path.with_suffix(".manifest.json").read_text())
    if manifest["dataset_sha256"] != sha256(path):
        raise ValueError("Dataset integrity failure")
    if manifest["incomplete_genes"] and not args.calibration and not (args.partial_source and manifest.get("partial_source_frozen")):
        raise ValueError("Incomplete acquisition may only be used for explicitly labeled calibration runs")
    rows = read_jsonl(path)
    allowed = {"experimental_success"}
    if args.allow_weak:
        allowed |= {"unverified", "computational_candidate"}
    rows = [row for row in rows if row["evidence_type"] in allowed]
    train = PrimerDataset([row for row in rows if row["split"] == "train"], settings["window_length"], args.seed, augment=True)
    validation = PrimerDataset([row for row in rows if row["split"] == "validation"], settings["window_length"], args.seed + 1, augment=False)
    if not len(train) or (not len(validation) and not args.calibration):
        raise ValueError("Insufficient eligible training/validation genes for selected evidence regime")
    loader = DataLoader(train, batch_size=settings["batch_size"], shuffle=True, generator=torch.Generator().manual_seed(args.seed))
    val_loader = DataLoader(validation, batch_size=settings["batch_size"])
    model = ConditionalDenoiser(settings["hidden_dim"], settings["layers"])
    optimizer = torch.optim.AdamW(model.parameters(), lr=settings["learning_rate"], weight_decay=0.01)
    steps = settings["diffusion_steps"]
    epochs = args.epochs or settings["epochs"]
    label = "calibration" if args.calibration else ("weak_partial_pilot" if args.partial_source and args.allow_weak else ("weak_pilot" if args.allow_weak else "validated"))
    name = f"{args.mode}_{label}_seed{args.seed}"
    directory = root / "results/checkpoints"
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / f"{name}.pt"
    history, best, started = [], math.inf, time.monotonic()

    def objective(batch, generator=None):
        clean, valid = batch["clean"], batch["valid"]
        if args.mode == "diffusion":
            t = torch.randint(1, steps + 1, (len(clean),), generator=generator)
            noisy, masked = corrupt(clean, t, steps, valid, generator, batch["fixed"])
            logits = model(batch["dna"], noisy, t.float() / steps, batch["conditions"], valid)
            return diffusion_loss(logits, clean, masked, t, steps, valid)
        noisy = torch.full_like(clean, MASK)
        noisy[batch["fixed"]] = clean[batch["fixed"]]
        logits = model(batch["dna"], noisy, torch.ones(len(clean)), batch["conditions"], valid)
        ce = F.cross_entropy(logits.transpose(1, 2), clean, reduction="none")
        supervised = valid & ~batch["fixed"]
        return (ce * supervised).sum() / supervised.sum().clamp(min=1)

    for epoch in range(epochs):
        train.epoch = epoch
        model.train()
        losses = []
        for batch in loader:
            optimizer.zero_grad(set_to_none=True)
            loss = objective(batch)
            if not torch.isfinite(loss):
                raise ValueError("Nonfinite loss")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            losses.append(loss.item())
        model.eval()
        val_losses = []
        generator = torch.Generator().manual_seed(args.seed + 1000)
        with torch.no_grad():
            for batch in val_loader:
                val_losses.append(objective(batch, generator).item())
        train_loss = float(np.mean(losses))
        val_loss = float(np.mean(val_losses)) if val_losses else None
        item = {"epoch": epoch + 1, "training_loss": train_loss, "validation_loss": val_loss,
                "elapsed_seconds": time.monotonic() - started}
        history.append(item)
        print(json.dumps({"run": name, **item}), flush=True)
        comparison = val_loss if val_loss is not None else train_loss
        if comparison < best:
            best = comparison
            torch.save({"state_dict": model.state_dict(), "model_spec": model.spec,
                        "mode": args.mode, "epoch": epoch + 1, "seed": args.seed, "steps": steps,
                        "dataset_sha256": sha256(path), "split_sha256": manifest["split_sha256"],
                        "stage": label, "allowed_evidence": sorted(allowed)}, destination)
        write_json(root / f"results/{name}_training.json", {
            "run": name, "started_at": utc_now(), "history": history, "train_genes": len(train.genes),
            "validation_genes": len(validation.genes), "parameters": sum(p.numel() for p in model.parameters()),
            "config": config, "arguments": vars(args), "checkpoint": str(destination),
            "interpretation": "Position-imitation objective; no measured qPCR efficacy is inferred from this loss."})
    print(f"Checkpoint: {destination}", flush=True)


if __name__ == "__main__":
    main()
