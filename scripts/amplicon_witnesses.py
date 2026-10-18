"""Exact product-existence decisions with one coordinate witness per subject/role."""
from bisect import bisect_left, bisect_right
from collections import defaultdict


def potential_amplicon_witnesses(forward_hits, reverse_hits, minimum=36, maximum=1000):
    grouped = defaultdict(lambda: defaultdict(list))
    for role, hits in [("F", forward_hits), ("R", reverse_hits)]:
        for hit in hits:
            grouped[hit["sseqid"]][(role, hit["strand"])].append(hit)
    products = []
    for subject, groups in grouped.items():
        for right_role in ["F", "R"]:
            right_hits = sorted(groups[(right_role, "-")], key=lambda hit: hit["end"])
            ends = [hit["end"] for hit in right_hits]
            if not ends:
                continue
            for left_role in ["F", "R"]:
                found = False
                for left in groups[(left_role, "+")]:
                    lower = bisect_left(ends, left["start"] + minimum)
                    upper = bisect_right(ends, left["start"] + maximum)
                    for index in range(lower, upper):
                        right = right_hits[index]
                        if left["end"] <= right["start"]:
                            products.append({"subject": subject, "start": left["start"], "end": right["end"],
                                "amplicon_length": right["end"] - left["start"],
                                "left_primer": left_role, "right_primer": right_role})
                            found = True
                            break
                    if found:
                        break
    return products
