#!/usr/bin/env python3
"""Map independent public literature oligos to frozen transcript references.

This evaluates sequence provenance and coverage, without treating reconstructed
combinations of literature strings as verified assays or as PICKR ranked pairs.
"""
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT / "src"))
import ahocorasick
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.sequence import reverse_complement
from diffuprimer.constraints import Conditions, oligo_properties


def main():
    source = ROOT / "data/raw/pickr/raw_primers_3.csv.gz"
    out = ROOT / "experiments/attention/results/pickr_literature"
    out.mkdir(parents=True,exist_ok=True)
    counts = Counter(); raw_rows = eligible_rows = 0
    with gzip.open(source,"rt") as handle:
        for row in csv.DictReader(handle):
            raw_rows += 1
            sequence = row["Sequence"].strip().upper()
            if 18 <= len(sequence) <= 25 and set(sequence) <= set("ACGT"):
                counts[sequence] += 1; eligible_rows += 1
    retained = {seq:n for seq,n in counts.items() if n >= 3}
    print(json.dumps({"raw_rows":raw_rows,"eligible_sequences":len(counts),"repeated_sequences":len(retained)}),flush=True)
    tree = ahocorasick.Automaton()
    patterns = defaultdict(list)
    for sequence in retained:
        patterns[sequence].append((sequence,"forward"))
        patterns[reverse_complement(sequence)].append((sequence,"reverse"))
    for pattern,records in patterns.items():
        tree.add_word(pattern,(len(pattern),records))
    tree.make_automaton()
    target_path = ROOT / "data/processed/targets_pilot.jsonl"
    targets = read_jsonl(target_path)
    mapped = []; matched_sequences = set()
    for target in targets:
        for end,(length,records) in tree.iter(target["template_sequence"]):
            for sequence,role in records:
                mapped.append({"gene_id":target["gene_id"],"split":target["split"],
                    "transcript_accession_version":target["transcript_accession_version"],
                    "sequence":sequence,"role":role,"start":end-length+1,"end":end+1,
                    "source_occurrences":retained[sequence]})
                matched_sequences.add(sequence)
    citations = defaultdict(set)
    with gzip.open(source,"rt") as handle:
        for row in csv.DictReader(handle):
            sequence = row["Sequence"].strip().upper()
            if sequence in matched_sequences and row["PMCID"]:
                citations[sequence].add(row["PMCID"])
    for row in mapped:
        row["source_pmcids"] = sorted(citations[row["sequence"]])
        row["distinct_source_publications"] = len(row["source_pmcids"])
        prop = oligo_properties(row["sequence"],Conditions())
        row["individual_oligo_feasible"] = (
            35 <= prop["gc_percent"] <= 65 and 57 <= prop["tm"] <= 63
            and prop["max_homopolymer"] <= 5 and prop["hairpin_tm"] <= 47 and prop["homodimer_tm"] <= 47)
    write_jsonl(out / "mapped_oligos.jsonl",mapped)
    summaries = []
    for target in targets:
        own = [row for row in mapped if row["gene_id"]==target["gene_id"]]
        by_role = {role:{row["sequence"] for row in own if row["role"]==role} for role in ["forward","reverse"]}
        summaries.append({"gene_id":target["gene_id"],"split":target["split"],
            "mapped_unique_oligos":len({row["sequence"] for row in own}),
            "forward_unique_oligos":len(by_role["forward"]),"reverse_unique_oligos":len(by_role["reverse"]),
            "has_both_orientations":bool(by_role["forward"] and by_role["reverse"]),
            "individual_feasible_unique_oligos":len({row["sequence"] for row in own if row["individual_oligo_feasible"]})})
    write_jsonl(out / "per_gene.jsonl",summaries)
    summary = {"completed_at":utc_now(),"source_sha256":sha256(source),"source_snapshot":json.loads((source.parent / "source_snapshot.json").read_text()),
        "targets_sha256":sha256(target_path),"mapper_sha256":sha256(Path(__file__)),
        "raw_rows":raw_rows,"eligible_rows_18_to_25_nt_ACGT":eligible_rows,"eligible_unique_sequences":len(counts),
        "sequences_with_at_least_three_source_occurrences":len(retained),
        "matched_unique_literature_sequences":len(matched_sequences),"mapped_binding_intervals":len(mapped),
        "mapped_oligos_sha256":sha256(out / "mapped_oligos.jsonl"),
        "per_split":{split:{"genes":sum(r["split"]==split for r in summaries),
            "genes_with_literature_oligos":sum(r["split"]==split and r["mapped_unique_oligos"]>0 for r in summaries),
            "genes_with_individual_feasible_literature_oligos":sum(r["split"]==split and r["individual_feasible_unique_oligos"]>0 for r in summaries),
            "genes_with_both_orientations":sum(r["split"]==split and r["has_both_orientations"] for r in summaries)} for split in ["train","validation","test"]},
        "usage":"Independent literature-source mapping only; no data supplied to model training; no ranked-pair reconstruction or experimental-outcome labels"}
    write_json(out / "summary.json",summary)
    print(json.dumps(summary),flush=True)


if __name__=="__main__":
    main()
