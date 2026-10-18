from __future__ import annotations

from collections import defaultdict
import re


def normalize_subject_id(identifier):
    match = re.search(r"(?:^|\|)([NX][MRCGW]_\d+\.\d+)(?:\||$)", identifier)
    return match.group(1) if match else identifier


def qualifying_hit(row, max_mismatches=2, perfect_3prime=5):
    """A declared limited screen, not a surrogate for experimental specificity."""
    if row["qstart"] != 1 or row["qend"] != row["qlen"] or row["gaps"] != 0:
        return False
    if row["mismatch"] > max_mismatches:
        return False
    return row["qseq"][-perfect_3prime:] == row["sseq"][-perfect_3prime:]


def parse_blast_tsv(path):
    fields = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length",
              "mismatch", "gaps", "qlen", "qseq", "sseq", "bitscore", "evalue"]
    numeric = {"qstart", "qend", "sstart", "send", "length", "mismatch", "gaps", "qlen"}
    result = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            values = line.rstrip().split("\t")
            if len(values) != len(fields):
                raise ValueError("Unexpected BLAST format")
            row = dict(zip(fields, values))
            row["raw_sseqid"] = row["sseqid"]
            row["sseqid"] = normalize_subject_id(row["sseqid"])
            for key in numeric:
                row[key] = int(row[key])
            row["bitscore"] = float(row["bitscore"])
            row["evalue"] = float(row["evalue"])
            if qualifying_hit(row):
                row["strand"] = "+" if row["sstart"] < row["send"] else "-"
                row["start"] = min(row["sstart"], row["send"]) - 1
                row["end"] = max(row["sstart"], row["send"])
                result[row["qseqid"]].append(row)
    return result


def potential_amplicons(forward_hits, reverse_hits, minimum=36, maximum=1000):
    """Pair inward-facing binding sites, including F/F and R/R products."""
    by_subject = defaultdict(list)
    for label, hits in [("F", forward_hits), ("R", reverse_hits)]:
        for hit in hits:
            by_subject[hit["sseqid"]].append({**hit, "primer_role": label})
    products = {}
    for subject, hits in by_subject.items():
        plus = [hit for hit in hits if hit["strand"] == "+"]
        minus = [hit for hit in hits if hit["strand"] == "-"]
        for left in plus:
            for right in minus:
                product_length = right["end"] - left["start"]
                if left["end"] <= right["start"] and minimum <= product_length <= maximum:
                    key = (subject, left["start"], right["end"], left["primer_role"], right["primer_role"])
                    products[key] = {"subject": subject, "start": left["start"], "end": right["end"],
                                     "amplicon_length": product_length,
                                     "left_primer": left["primer_role"], "right_primer": right["primer_role"]}
    return list(products.values())


def classify_transcript_products(products, target_gene, transcript_gene_map):
    off_target, unknown, on_target = [], [], []
    for product in products:
        gene = transcript_gene_map.get(product["subject"])
        item = {**product, "gene_id": gene}
        if gene is None:
            unknown.append(item)
        elif gene != target_gene:
            off_target.append(item)
        else:
            on_target.append(item)
    return {"on_target_products": on_target, "off_target_products": off_target,
            "unresolved_products": unknown,
            "gene_specificity_screen_pass": bool(on_target) and not off_target and not unknown,
            "scope": "full-length ungapped hits, <=2 mismatches, perfect terminal 5 bases; potential products 36..1000 bp"}
