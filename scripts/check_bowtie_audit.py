#!/usr/bin/env python3
"""Compare the actual aligner/parser against exhaustive synthetic scanning."""
import random
import subprocess
from pathlib import Path

from bowtie_audit import alignment_command, parse_sam
from diffuprimer.io import sha256, utc_now, write_json
from diffuprimer.sequence import reverse_complement


root = Path(__file__).resolve().parents[1]
directory = root / "results/checks/bowtie_exhaustive"
directory.mkdir(parents=True, exist_ok=True)
rng = random.Random(20261013)
queries = {"P18": "ACGTGCTAGCTACGATGC", "P20": "CGTACAGCTGCATGACTGAC", "P25": "GCTACGTAGCATGCATGCTAGCAGT"}
references = {"random": "".join(rng.choices("ACGT", k=2000))}
references["softmasked_plus"] = "gggg" + queries["P18"].lower() + "cccc"
references["softmasked_minus"] = "gggg" + reverse_complement(queries["P18"]).lower() + "cccc"
for qid, sequence in queries.items():
    patterns = [[], [0], [6, 12], [len(sequence)-1], [0, 6, 12]]
    for number, positions in enumerate(patterns):
        changed = list(sequence)
        for pos in positions:
            changed[pos] = next(base for base in "ACGT" if base != changed[pos])
        mutated = "".join(changed)
        for strand, site in [("plus", mutated), ("minus", reverse_complement(mutated))]:
            references[f"{qid}_{number}_{strand}"] = "G" * 11 + site + "C" * 13
for i in range(20):
    length = rng.choice([18, 19, 20, 21, 25])
    start = rng.randrange(0, 2000-length)
    sequence = references["random"][start:start+length]
    if i % 2:
        sequence = reverse_complement(sequence)
    if i % 3:
        positions = rng.sample(range(length-5), i % 3)
        sequence = list(sequence)
        for pos in positions:
            sequence[pos] = next(base for base in "ACGT" if base != sequence[pos])
        sequence = "".join(sequence)
    queries[f"random{i}"] = sequence
reference = directory / "reference.fasta"
query = directory / "queries.fasta"
reference.write_text("".join(f">{name}\n{seq}\n" for name, seq in references.items()))
query.write_text("".join(f">{name}\n{seq}\n" for name, seq in queries.items()))
binary = root / "tools/bowtie-1.3.1-linux-x86_64/bowtie"
with (directory / "build.log").open("w") as handle:
    subprocess.run([str(binary.parent / "bowtie-build"), str(reference), str(directory / "index")],
        check=True, stdout=handle, stderr=subprocess.STDOUT)
output = directory / "hits.sam"
command = alignment_command(binary, directory / "index", query, output, threads=1)
with (directory / "align.log").open("w") as handle:
    subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
observed = parse_sam(output, queries)
expected = set()
for qid, primer in queries.items():
    for name, seq in references.items():
        for start in range(len(seq)-len(primer)+1):
            site = seq[start:start+len(primer)].upper()
            for strand, aligned in [("+", site), ("-", reverse_complement(site))]:
                positions = [i for i, (a,b) in enumerate(zip(primer, aligned)) if a != b]
                if len(positions) <= 2 and all(i < len(primer)-5 for i in positions):
                    expected.add((qid, name, strand, start, start+len(primer), tuple(positions)))
actual = {(qid, hit["sseqid"], hit["strand"], hit["start"], hit["end"], tuple(hit["mismatch_positions_original"]))
    for qid, hits in observed.items() for hit in hits}
if actual != expected:
    write_json(directory / "failure.json", {"missing": sorted(expected-actual), "extra": sorted(actual-expected)})
    raise ValueError("Actual aligner differs from exhaustive synthetic reference scanning")
assert ("P18", "P18_2_plus", "+", 11, 29, (6,12)) in actual
assert ("P18", "P18_2_minus", "-", 11, 29, (6,12)) in actual
write_json(directory / "check.json", {"passed_at": utc_now(), "queries": len(queries), "references": len(references),
    "expected_hits": len(expected), "observed_hits": len(actual), "exact_set_equality": True,
    "cases": ["both strands", "5prime terminal mismatch", "two dispersed mismatches without a seven-base exact seed",
              "3prime terminal mismatch rejected", "three mismatches rejected", "18/20/25-base oligos", "all mismatch strata", "softmasked lowercase reference bases"],
    "parser_sha256": sha256(root / "scripts/bowtie_audit.py"), "command": command})
print(f"Exact agreement with exhaustive scanning: {len(actual)} hits", flush=True)
