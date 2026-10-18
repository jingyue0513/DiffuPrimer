#!/usr/bin/env python3
"""Evaluate native external primer tools, without loading any learned generator."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT / "src"))
from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import summarize_candidates
from diffuprimer.sequence import PrimerPair


def run_primer3(args):
    target_path=ROOT / args.targets_path
    targets=[r for r in read_jsonl(target_path) if r["split"]=="external_test"]
    targets.sort(key=lambda t:int(hashlib.sha256(f"{t['gene_id']}:{args.seed}".encode()).hexdigest()[:8],16))
    cdata=dict(json.loads((ROOT / "configs/pilot.json").read_text())["conditions"],amplicon_min=80,amplicon_max=200)
    conditions=Conditions(**cdata)
    name=f"primer3_design_80-200_external_test_seed{args.seed}"
    out=ROOT / args.output_root / name
    out.mkdir(parents=True,exist_ok=True)
    write_json(out / "manifest.json",{"method":"primer3","stage":"attention_external_comparison","split":"external_test",
        "task":"design","arguments":vars(args),"conditions":cdata,"target_sha256":sha256(target_path),
        "expected_targets":[r["gene_id"] for r in targets],"baseline_sha256":sha256(ROOT / "src/diffuprimer/baselines.py"),
        "driver_sha256":sha256(Path(__file__)),"candidate_accounting":"K=20 native returned pairs; shortages and duplicates remain in the budget",
        "external_cohort_manifest_sha256":sha256(ROOT / "experiments/attention_external/manifest.json")})
    ledger=[];all_candidates=[]
    for target in targets:
        started=time.monotonic();template=target["template_sequence"];gene=target["gene_id"]
        candidates,explain=primer3_candidates(template,conditions,args.budget)
        for i,row in enumerate(candidates):
            row.update(gene_id=gene,method="primer3",split="external_test",proposal_index=i,
                candidate_id=f"{name}:{gene}:{i}",transcript_accession_version=target["transcript_accession_version"],
                structural_valid=True,requested_amplicon_min=80,requested_amplicon_max=200)
            pair=PrimerPair(*(row[k] for k in ["forward_start","forward_end","reverse_start","reverse_end"]))
            row["assessment"]=assess_pair(template,pair,conditions)
        summary={"gene_id":gene,"method":"primer3","split":"external_test","task":"design","conditions":cdata,
            **summarize_candidates(candidates,args.budget),"elapsed_seconds":time.monotonic()-started,"primer3_explain":explain}
        ledger.append(summary);all_candidates.extend(candidates)
        write_jsonl(out / "per_gene.jsonl",ledger);write_jsonl(out / "candidates.jsonl",all_candidates)
        print(json.dumps({"method":"primer3","completed_genes":len(ledger),"gene":gene,"unique_feasible_pairs":summary["unique_feasible_pairs"]}),flush=True)
    write_json(out / "completion.json",{"completed_at":utc_now(),"targets":len(ledger),"all_expected_targets_evaluated":len(ledger)==len(targets),
        "split":"external_test","per_gene_sha256":sha256(out / "per_gene.jsonl"),"candidates_sha256":sha256(out / "candidates.jsonl")})


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--method",choices=["primer3","qprimer","primercast"],required=True)
    p.add_argument("--targets-path",default="data/processed/attention_external_targets.jsonl")
    p.add_argument("--output-root",default="experiments/attention_external/results/evaluations")
    p.add_argument("--seed",type=int,default=20261013)
    p.add_argument("--budget",type=int,default=20)
    args=p.parse_args()
    if args.method=="primer3":run_primer3(args);return
    python=ROOT / (".venv_primercast/bin/python" if args.method=="primercast" else ".venv/bin/python")
    script=ROOT / f"experiments/attention/baselines/run_{args.method}.py"
    command=[str(python),str(script),"--targets-path",args.targets_path,"--output-root",args.output_root,
        "--split","external_test","--budget",str(args.budget),"--seed",str(args.seed)]
    env=dict(os.environ)
    if args.method=="primercast":
        env["PATH"]=str(ROOT / ".venv_primercast/bin")+os.pathsep+env["PATH"]
        env["CUDA_VISIBLE_DEVICES"]=""
    subprocess.run(command,check=True,env=env)


if __name__=="__main__":main()
