#!/usr/bin/env python3
"""Run native pretrained PrimerCast quick-design on frozen transcript targets.

The upstream candidate search, DNA-duplex energy calculation, Bowtie2 mapping,
feature construction, neural scoring and native postfilters execute unchanged.
Only supported task bounds are supplied through its public parameter interface.
"""
import argparse
import contextlib
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).parent / "vendor/PrimerCast"
sys.path[:0] = [str(ROOT / "src"), str(VENDOR / "src")]
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import numpy as np
import pandas as pd
import torch
from primercast.commands import quick_design
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import summarize_candidates
from diffuprimer.sequence import PrimerPair, reverse_complement


def occurrences(template, sequence):
    start = template.find(sequence)
    while start >= 0:
        yield start
        start = template.find(sequence, start + 1)


def map_pair(template, forward, reverse, conditions):
    # The common endpoint counts unique synthesis-sequence pairs; when a sequence
    # has repeated perfect sites, choose a compatible coordinate witness.
    matches = []
    for f in occurrences(template, forward):
        for r in occurrences(template, reverse_complement(reverse)):
            if f + len(forward) <= r and conditions.amplicon_min <= r + len(reverse) - f <= conditions.amplicon_max:
                matches.append(PrimerPair(f, f + len(forward), r, r + len(reverse)))
    return matches[0] if matches else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--split", choices=["validation", "test", "external_test"], default="test")
    p.add_argument("--targets-path",default="data/processed/targets_pilot.jsonl")
    p.add_argument("--output-root",default="experiments/attention/results/evaluations")
    p.add_argument("--amplicon-min", type=int, default=80)
    p.add_argument("--amplicon-max", type=int, default=200)
    p.add_argument("--budget", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261013)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args()
    config = json.loads((ROOT / "configs/pilot.json").read_text())
    data = dict(config["conditions"], amplicon_min=a.amplicon_min, amplicon_max=a.amplicon_max)
    conditions = Conditions(**data)
    target_path=ROOT / a.targets_path
    targets = [t for t in read_jsonl(target_path) if t["split"] == a.split]
    targets.sort(key=lambda t: int(hashlib.sha256(f"{t['gene_id']}:{a.seed}".encode()).hexdigest()[:8], 16))
    name = f"primercast_design_{a.amplicon_min}-{a.amplicon_max}_{a.split}_seed{a.seed}"
    out = ROOT / a.output_root / name
    out.mkdir(parents=True, exist_ok=True)
    params = {"MAX_PRIMER_CANDIDATES":10000, "TM_MIN":conditions.tm_min, "TM_MAX":conditions.tm_max,
        "GC_MAX":conditions.gc_max, "DG_MIN":-6, "PRIMER_LEN_MIN":conditions.primer_min_length,
        "PRIMER_LEN_MAX":conditions.primer_max_length, "AMPLEN_MIN":conditions.amplicon_min,
        "AMPLEN_MAX":conditions.amplicon_max, "QUICK_COV_MIN":1.0, "QUICK_ACT_MIN":1.0,
        "QUICK_MIN_PAIRS":a.budget, "OFFLEN_MIN":36, "OFFLEN_MAX":1000, "NUM_TOP_SENSITIVITY":200}
    param_path = out / "native_params.txt"
    param_path.write_text("".join(f"{k} = {v}\n" for k, v in params.items()))
    snapshot = json.loads((VENDOR / "source_snapshot.json").read_text())
    manifest = {"method":"primercast", "stage":"attention_external_comparison" if a.split=="external_test" else "attention_comparison", "task":"design", "arguments":vars(a),
        "conditions":data,"split":a.split, "upstream_repository":snapshot["repository"], "upstream_commit":snapshot["commit"],
        "upstream_source_sha256":sha256(VENDOR / "src/primercast/commands/quick_design.py"),
        "upstream_classifier_sha256":sha256(VENDOR / "src/primercast/data/combined_classifier.pth"),
        "upstream_regressor_sha256":sha256(VENDOR / "src/primercast/data/combined_regressor.pth"),
        "adapter_sha256":sha256(Path(__file__)), "native_parameters":params,
        "target_sha256":sha256(target_path),
        "expected_targets":[t["gene_id"] for t in targets],
        "candidate_accounting":"Top K native on-target quick-design candidates ranked by native score; shortages and shared-assessment failures count against K",
        "supported_controls":"Length/Tm/GC-max via native interface; additional common bounds evaluated after selection",
        "native_reference_search":"Bowtie2 against the intended transcript for native ML scoring; common full RNA/genome screening downstream",
        "device":"cpu", "physical_measurement":"No new assay measurements; native pretrained model scores",
        "fasta_headers":"RefSeq accession versions, preserving string IDs during native pandas parsing",
        "rna_duplex_version":subprocess.check_output(["RNAduplex","--version"],text=True).strip(),
        "bowtie2_version":subprocess.check_output(["bowtie2","--version"],text=True).splitlines()[0]}
    write_json(out / "manifest.json", manifest)
    ledger = read_jsonl(out / "per_gene.jsonl") if (out / "per_gene.jsonl").exists() else []
    all_candidates = read_jsonl(out / "candidates.jsonl") if (out / "candidates.jsonl").exists() else []
    completed = {r["gene_id"] for r in ledger}
    torch.set_num_threads(a.threads)
    for target in targets:
        gene = target["gene_id"]
        if gene in completed:
            continue
        started = time.monotonic()
        seed = int(hashlib.sha256(f"{gene}:{a.seed}".encode()).hexdigest()[:8],16)
        random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
        own = out / "native" / gene
        own.mkdir(parents=True, exist_ok=True)
        fasta = own / "target.fasta"
        fasta.write_text(f">{target['transcript_accession_version']}\n{target['template_sequence']}\n")
        csv_path = own / "ranked.csv"
        namespace = argparse.Namespace(target=str(fasta),selected=str(fasta),out=str(csv_path),
            param_file=str(param_path),name=gene,msa=None,target_window=None,tmpdir=str(own),
            init_fa=None,init_feat=None,probe_mode=False,probe_fa=None,probe_feat=None,probe_csv=None,
            probe_pair_csv=None,threads=a.threads)
        with (own / "native.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            quick_design.run(namespace)
        native = pd.read_csv(csv_path) if csv_path.exists() and csv_path.stat().st_size else pd.DataFrame()
        if not native.empty:
            native = native.sort_values("score",ascending=False,kind="stable").head(a.budget)
        candidates = []
        for i, native_row in enumerate(native.to_dict("records")):
            forward, reverse = native_row["fwd_seq"], native_row["rev_seq"]
            pair = map_pair(target["template_sequence"],forward,reverse,conditions)
            row = {"gene_id":gene,"transcript_accession_version":target["transcript_accession_version"],
                "method":"primercast","split":a.split,"proposal_index":i,"candidate_id":f"{name}:{gene}:{i}",
                "structural_valid":pair is not None,"native_score":float(native_row["score"]),
                "native_coverage":float(native_row["coverage"]),"native_predicted_activity":float(native_row["activity"]),
                "requested_amplicon_min":conditions.amplicon_min,"requested_amplicon_max":conditions.amplicon_max}
            if pair is not None:
                row.update(pair.as_dict(target["template_sequence"]))
                assert row["forward_primer"]==forward and row["reverse_primer"]==reverse
                row["assessment"] = assess_pair(target["template_sequence"],pair,conditions)
            else:
                row.update(forward_primer=forward,reverse_primer=reverse,
                    assessment={"computational_feasible":False,"failures":["no_compatible_template_pair"]})
            candidates.append(row)
        row = {"gene_id":gene,"method":"primercast","split":a.split,"task":"design","conditions":data,
            **summarize_candidates(candidates,a.budget),"elapsed_seconds":time.monotonic()-started,
            "native_returned_pairs_before_budget":sum(1 for _ in csv_path.open())-1 if csv_path.exists() and csv_path.stat().st_size else 0}
        ledger.append(row); all_candidates.extend(candidates)
        write_jsonl(out / "per_gene.jsonl",ledger);write_jsonl(out / "candidates.jsonl",all_candidates)
        print(json.dumps({"run":name,"completed_genes":len(ledger),"gene":gene,
            "unique_feasible_pairs":row["unique_feasible_pairs"],"seconds":row["elapsed_seconds"]}),flush=True)
    write_json(out / "completion.json",{"completed_at":utc_now(),"targets":len(ledger),
        "all_expected_targets_evaluated":len(ledger)==len(targets),
        "per_gene_sha256":sha256(out / "per_gene.jsonl"),"candidates_sha256":sha256(out / "candidates.jsonl")})


if __name__ == "__main__":
    main()
