"""Distance-indexed pairing, preserving the original product values and order."""
from bisect import bisect_left, bisect_right
from collections import defaultdict


def potential_amplicons(forward_hits, reverse_hits, minimum=36, maximum=1000):
    by_subject = defaultdict(list)
    for label, hits in [("F", forward_hits), ("R", reverse_hits)]:
        for hit in hits:
            by_subject[hit["sseqid"]].append({**hit, "primer_role": label})
    products = {}
    for subject, hits in by_subject.items():
        plus = [hit for hit in hits if hit["strand"] == "+"]
        minus = [hit for hit in hits if hit["strand"] == "-"]
        ranked = sorted(enumerate(minus), key=lambda item: item[1]["end"])
        ends = [hit["end"] for _, hit in ranked]
        for left in plus:
            lower = bisect_left(ends, left["start"] + minimum)
            upper = bisect_right(ends, left["start"] + maximum)
            for _, right in sorted(ranked[lower:upper], key=lambda item: item[0]):
                if left["end"] <= right["start"]:
                    key = (subject, left["start"], right["end"], left["primer_role"], right["primer_role"])
                    products[key] = {"subject": subject, "start": left["start"], "end": right["end"],
                        "amplicon_length": right["end"] - left["start"],
                        "left_primer": left["primer_role"], "right_primer": right["primer_role"]}
    return list(products.values())
