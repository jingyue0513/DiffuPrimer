#!/usr/bin/env python3
"""Read-only verification of reported cohorts, sequences, anchors and source hashes."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from diffuprimer.io import read_jsonl,sha256,utc_now
from diffuprimer.sequence import reverse_complement

spec=importlib.util.spec_from_file_location("attention_analysis",Path(__file__).with_name("analyze.py"))
analysis=importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",default="experiments/attention/results/study_summary.json")
    parser.add_argument("--allow-partial",action="store_true")
    parser.add_argument("--protocol",help="Optionally require the declared final experiment protocol")
    args=parser.parse_args()
    summary_path=ROOT/args.summary
    summary_bytes=summary_path.read_bytes()
    summary=json.loads(summary_bytes)
    summary_sha=hashlib.sha256(summary_bytes).hexdigest()
    if summary["status"]!="complete" and not args.allow_partial:
        parser.error("Summary is incomplete; --allow-partial verifies completed sources only")
    pending_protocol=[]
    if args.protocol:
        protocol_path=ROOT/args.protocol
        protocol=json.loads(protocol_path.read_text())
        planned=analysis.planned_names(protocol,summary.get("evaluation_split","test"))
        present={s.get("name",Path(s["directory"]).name) for s in summary["inputs"]+summary.get("validation_inputs",[])}
        pending_protocol=sorted(planned-present)
        if not args.allow_partial and (pending_protocol or summary.get("protocol_sha256")!=sha256(protocol_path)):
            raise ValueError("Summary does not satisfy the declared final protocol")
    targets_path=Path(summary.get("targets_path",ROOT/"data/processed/targets_pilot.jsonl"))
    if summary.get("targets_sha256") and sha256(targets_path)!=summary["targets_sha256"]:
        raise ValueError("Frozen targets changed after analysis")
    targets=read_jsonl(targets_path)
    split=summary.get("evaluation_split","test")
    expected={t["gene_id"]:t for t in targets if t["split"]==split}
    if len(expected)!=summary["total_test_genes"]:
        raise ValueError("Summary gene denominator disagrees with its frozen cohort")
    if split=="external_test":
        original={r["gene_id"] for r in read_jsonl(ROOT/"data/processed/targets_pilot.jsonl")}
        if summary.get("old_sources_included") or set(expected)&original:
            raise ValueError("Independent extra genes have been pooled with original study genes")
    records=read_jsonl(summary_path.with_suffix(".per_gene.jsonl"))
    grouped=defaultdict(list)
    for record in records:
        if record["gene_id"] not in expected:
            raise ValueError("Reported record belongs to a foreign cohort")
        grouped[record["context"],record["method"],record["variant"],record["budget"]].append(record)
    group_values=0
    for g in summary["summaries"]:
        own=grouped[g["context"],g["method"],g["variant"],g["budget"]]
        if len(own)!=len(expected) or {r["gene_id"] for r in own}!=set(expected):
            raise ValueError("A reported group does not retain the complete target cohort")
        for endpoint in analysis.ENDPOINTS:
            values=[r[endpoint] for r in own if r[endpoint] is not None]
            mean=sum(values)/len(values) if values else None
            if mean is None:
                assert g[endpoint] is None
            else:
                assert math.isclose(mean,g[endpoint],abs_tol=1e-12,rel_tol=1e-12)
            assert len(values)==g[endpoint+"_n"]
            group_values+=1
    anchors={r["gene_id"]:r for r in read_jsonl(ROOT/"experiments/final_followup/redesign_tasks.jsonl")}
    sequence_pairs,anchored_pairs,returned_candidates,feasible_pairs,actual_anchor_violations=0,0,0,0,0
    validated_sources=0
    for source in summary["inputs"]:
        directory=Path(source["directory"])
        manifest,candidates,ledger,screens,hashes=analysis.load_source(directory,expected,summary["seed"],pinned=source,split=split)
        task=manifest.get("task")
        side=manifest.get("fixed_side") or manifest["arguments"].get("fixed_side") or manifest["arguments"].get("repair")
        anchored=task in ("completion","redesign") or side is not None
        per_target=defaultdict(int)
        for c in candidates:
            per_target[c["gene_id"]]+=1
            returned_candidates+=1
            if not c["structural_valid"]:
                assert not c["assessment"]["computational_feasible"]
                continue
            template=expected[c["gene_id"]]["template_sequence"]
            sF,eF,sR,eR=(c[k] for k in ("forward_start","forward_end","reverse_start","reverse_end"))
            if not 0<=sF<eF<=sR<eR<=len(template):
                raise ValueError("Structurally valid candidate has incorrect transcript coordinates")
            if c["forward_primer"]!=template[sF:eF] or c["reverse_primer"]!=reverse_complement(template[sR:eR]):
                raise ValueError("Synthesized F/R sequences disagree with fixed-template binding intervals")
            if c["amplicon_length"]!=eR-sF:
                raise ValueError("Product length disagrees with outer binding coordinates")
            sequence_pairs+=1
            if c["assessment"]["computational_feasible"]:
                conditions=ledger[c["gene_id"]].get("conditions",manifest["conditions"])
                lo=c.get("requested_amplicon_min",conditions["amplicon_min"])
                hi=c.get("requested_amplicon_max",conditions["amplicon_max"])
                assert lo<=eR-sF<=hi
                assert conditions["primer_min_length"]<=eF-sF<=conditions["primer_max_length"]
                assert conditions["primer_min_length"]<=eR-sR<=conditions["primer_max_length"]
                feasible_pairs+=1
            if anchored:
                supplied=anchors[c["gene_id"]]["anchor"]
                actual=(sF,eF)==(supplied["forward_start"],supplied["forward_end"]) if side=="forward" else (sR,eR)==(supplied["reverse_start"],supplied["reverse_end"])
                sequence_key="forward_primer" if side=="forward" else "reverse_primer"
                actual=actual and c[sequence_key]==supplied[sequence_key]
                if c.get("anchor_preserved",actual)!=actual:
                    raise ValueError("Declared anchor preservation disagrees with exact coordinates/sequence")
                if not actual:
                    actual_anchor_violations+=1
                    assert not c["assessment"]["computational_feasible"]
                anchored_pairs+=1
        if any(n>20 for n in per_target.values()):
            raise ValueError("A source returned more than the common proposal budget")
        validated_sources+=1
    # Check that all reported paired units are genuine common genes within this cohort.
    paired_values=0
    for contrast in summary["paired_comparisons"]:
        key=(contrast["context"],contrast["left_method"],contrast["left_variant"],contrast["budget"])
        other=(contrast["context"],contrast["right_method"],contrast["right_variant"],contrast["budget"])
        left={r["gene_id"]:r for r in grouped[key]};right={r["gene_id"]:r for r in grouped[other]}
        endpoint=contrast["endpoint"]
        genes=sorted(g for g in left.keys()&right.keys() if left[g][endpoint] is not None and right[g][endpoint] is not None)
        assert genes==contrast["joint_gene_ids"] and len(genes)==contrast["n"]
        mean=sum(left[g][endpoint]-right[g][endpoint] for g in genes)/len(genes)
        assert math.isclose(mean,contrast["mean_difference"],abs_tol=1e-12,rel_tol=1e-12)
        paired_values+=1
    print(json.dumps({"at":utc_now(),"passed":True,"status":summary["status"],"summary_sha256":summary_sha,
        "pending_declared_protocol_sources":pending_protocol,
        "evaluation_split":split,"total_genes":len(expected),"sources_checked":validated_sources,
        "group_endpoints_checked":group_values,"paired_estimates_checked":paired_values,"returned_candidates":returned_candidates,
        "valid_F_R_sequence_pairs_checked":sequence_pairs,"feasible_geometry_pairs_checked":feasible_pairs,
        "exact_anchor_pairs_checked":anchored_pairs,"anchor_violations_among_returned_candidates":actual_anchor_violations,
        "mode":"read_only; no evaluation, candidate, checkpoint or summary files modified"},indent=2))


if __name__=="__main__":main()
