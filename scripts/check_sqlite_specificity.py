#!/usr/bin/env python3
"""Verify database decisions against complete synthetic alignment enumeration."""
import json
import random
from pathlib import Path

from bowtie_audit import parse_sam
from bowtie_sqlite_cache import AlignmentCache
from diffuprimer.io import sha256, utc_now, write_json
from diffuprimer.specificity import potential_amplicons, classify_transcript_products
from diffuprimer.sequence import reverse_complement


root = Path(__file__).resolve().parents[1]
source = root / "results/checks/bowtie_exhaustive"
queries, name, parts = {}, None, []
for line in (source / "queries.fasta").read_text().splitlines():
    if line.startswith(">"):
        if name is not None:
            queries[name] = "".join(parts)
        name, parts = line[1:], []
    else:
        parts.append(line)
queries[name] = "".join(parts)
directory = root / "results/checks/sqlite_specificity"
directory.mkdir(parents=True, exist_ok=True)
path = directory / (sha256(root / "scripts/bowtie_sqlite_cache.py")[:16] + ".sqlite")
cache = AlignmentCache(path)
cache.ingest(source / "hits.sam", queries, 0)
cache.ingest(source / "hits.sam", queries, 1)
registered = cache.register(queries.values())
hits = parse_sam(source / "hits.sam", queries)
comparisons = products = 0
mapping = {"random": "1"}
for f, forward in queries.items():
    for r, reverse in queries.items():
        expected = potential_amplicons(hits.get(f, []),hits.get(r, []))
        expected_roles = {(p["subject"],p["left_primer"],p["right_primer"]) for p in expected}
        for reference in [0,1]:
            observed = cache.products(registered[forward][0],registered[reverse][0],reference)
            assert {(p["subject"],p["left_primer"],p["right_primer"]) for p in observed} == expected_roles
            assert all(p in expected for p in observed)
            assert classify_transcript_products(expected,"1",mapping)["gene_specificity_screen_pass"] == classify_transcript_products(observed,"1",mapping)["gene_specificity_screen_pass"]
            products += len(observed)
            comparisons += 1
assert products > 0
random_source = directory / "random_and_boundary.sam"
rng = random.Random(20261013)
random_queries, random_hits, sam = {}, {}, []
for index in range(24):
    sequence = "".join(rng.choice("ACGT") for _ in range(18+index%8))
    qid = f"random{index}"
    random_queries[qid] = sequence
    random_hits[qid] = []
    for item in range(75):
        strand = rng.choice(["+","-"])
        start = rng.choice([0,18,36-len(sequence),1000-len(sequence),1001-len(sequence),rng.randrange(3500)])
        subject = rng.choice(["on","off","unknown"])
        hit = {"sseqid":subject,"strand":strand,"start":start,"end":start+len(sequence)}
        random_hits[qid].append(hit)
        flag = 16 if strand == "-" else 0
        oriented = reverse_complement(sequence) if flag else sequence
        line = f"{qid}\t{flag}\t{subject}\t{start+1}\t255\t{len(sequence)}M\t*\t0\t0\t{oriented}\t*\tMD:Z:{len(sequence)}\n"
        sam.append(line)
        if item % 5 == 0:
            sam.append(line)
random_source.write_text("".join(sam))
random_cache = AlignmentCache(directory / (sha256(random_source)[:16] + "_random.sqlite"))
for reference in [0,1]:
    random_cache.ingest(random_source,random_queries,reference)
registered = random_cache.register(random_queries.values())
random_comparisons = 0
for index in range(120):
    f,r = rng.sample(list(random_queries),2)
    if index % 5 == 0:
        r=f
    expected = potential_amplicons(random_hits[f],random_hits[r])
    expected_roles = {(p['subject'],p['left_primer'],p['right_primer']) for p in expected}
    for reference in [0,1]:
        observed = random_cache.products(registered[random_queries[f]][0],registered[random_queries[r]][0],reference)
        assert {(p['subject'],p['left_primer'],p['right_primer']) for p in observed} == expected_roles
        assert all(p in expected for p in observed)
        assert classify_transcript_products(expected,'1',{'on':'1','off':'2'})['gene_specificity_screen_pass'] == classify_transcript_products(observed,'1',{'on':'1','off':'2'})['gene_specificity_screen_pass']
        random_comparisons += 1
write_json(directory / "check.json", {"checked_at": utc_now(), "query_pairs_and_references": comparisons,
    "additional_random_and_boundary_comparisons":random_comparisons,
    "cases":["duplicate sites","same oligo","multiple and unknown subjects","variable oligo lengths","36/1000 inclusive boundaries","all four primer role combinations"],
    "verified_coordinate_witnesses": products, "exact_subject_role_coverage": True,
    "exact_gene_specificity_decisions": True, "exact_genomic_existence_decisions": True,
    "cache_sha256": sha256(root / "scripts/bowtie_sqlite_cache.py"), "parser_sha256": sha256(root / "scripts/bowtie_audit.py"),
    "evidence_format": "One actual product coordinate per subject and primer-role combination; all binding sites retained in the database"})
print(f"SQLite decisions match complete enumeration: {comparisons} pair/reference comparisons",flush=True)
