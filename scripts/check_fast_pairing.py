#!/usr/bin/env python3
"""Verify ordered output equality against the frozen pairing kernel."""
import random
import time
from pathlib import Path

from amplicon_pairing_fast import potential_amplicons as fast
from amplicon_witnesses import potential_amplicon_witnesses as witnesses
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json
from diffuprimer.specificity import parse_blast_tsv, potential_amplicons as original


root = Path(__file__).resolve().parents[1]
rng = random.Random(20261013)


def compare_witnesses(pools, minimum=36, maximum=1000):
    full = original(*pools, minimum=minimum, maximum=maximum)
    compact = witnesses(*pools, minimum=minimum, maximum=maximum)
    roles = lambda rows: {(r["subject"], r["left_primer"], r["right_primer"]) for r in rows}
    assert roles(full) == roles(compact)
    assert all(row in full for row in compact)


for i in range(200):
    pools = []
    for side in range(2):
        hits = []
        for j in range(rng.randrange(0, 150)):
            start = rng.randrange(0, 3000)
            hits.append({"sseqid": rng.choice(["A", "B", "C"]), "strand": rng.choice(["+", "-"]),
                         "start": start, "end": start + rng.randrange(18, 26)})
        if hits and i % 3 == 0:
            hits += hits[:5]
        rng.shuffle(hits)
        pools.append(hits)
    if i % 5 == 0:
        pools[1] = pools[0]
    minimum, maximum = [(36, 1000), (80, 200), (36,36), (999,1000)][i % 4]
    assert original(*pools, minimum=minimum, maximum=maximum) == fast(*pools, minimum=minimum, maximum=maximum)
    compare_witnesses(pools, minimum, maximum)
hits = [{"sseqid": "large", "strand": "+" if i % 2 else "-", "start": i*1000, "end": i*1000+20} for i in range(3000)]
start = time.monotonic()
expected = original(hits, hits)
old_seconds = time.monotonic() - start
start = time.monotonic()
observed = fast(hits, hits)
new_seconds = time.monotonic() - start
assert observed == expected
calibration = root / "results/evaluations/primer3_calibration_validation_seed20261013_amp80-200"
queries, qid, sequence = {}, None, []
for line in (calibration / "specificity_queries.fasta").read_text().splitlines():
    if line.startswith(">"):
        if qid is not None:
            queries["".join(sequence)] = qid
        qid, sequence = line[1:], []
    else:
        sequence.append(line)
queries["".join(sequence)] = qid
real_cases = 0
for name in ["transcriptome", "genome"]:
    hits = parse_blast_tsv(calibration / f"blast_{name}.tsv")
    for row in read_jsonl(calibration / "candidates.jsonl"):
        if not row["assessment"]["computational_feasible"]:
            continue
        f, r = queries[row["forward_primer"]], queries[row["reverse_primer"]]
        assert original(hits.get(f, []), hits.get(r, [])) == fast(hits.get(f, []), hits.get(r, []))
        compare_witnesses((hits.get(f, []), hits.get(r, [])))
        real_cases += 1
write_json(root / "results/checks/fast_pairing/check.json", {"checked_at": utc_now(),
    "exact_ordered_output_equality": True, "random_cases": 200, "large_case_hits_per_oligo": len(hits),
    "original_seconds": old_seconds, "indexed_seconds": new_seconds,
    "real_reference_candidate_reference_comparisons": real_cases,
    "original_sha256": sha256(root / "src/diffuprimer/specificity.py"),
    "indexed_sha256": sha256(root / "scripts/amplicon_pairing_fast.py"),
    "witness_sha256": sha256(root / "scripts/amplicon_witnesses.py"),
    "witnesses_preserve_all_subject_and_role_existence_decisions": True,
    "cases": ["F/R,R/F,F/F,R/R", "same oligo on both sides", "duplicate and unordered hits", "multiple subjects", "inclusive length boundaries", "empty hits"]})
print("Distance-indexed pairing preserves exact ordered products", flush=True)
