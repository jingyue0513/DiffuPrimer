"""Cache binding-site indexes and retain every subject/primer-role decision."""
from bisect import bisect_left, bisect_right
from collections import defaultdict


sources, parts = {}, {}


def prepared(hits):
    key = id(hits)
    if key in sources:
        original, groups = sources[key]
        if original is not hits:
            raise ValueError("Binding-site cache identity changed")
        return groups
    grouped = defaultdict(lambda: {"+": [], "-": []})
    for hit in hits:
        grouped[hit["sseqid"]][hit["strand"]].append(hit)
    result = {}
    for subject, group in grouped.items():
        plus = sorted(group["+"], key=lambda h: h["start"])
        minus = sorted(group["-"], key=lambda h: h["end"])
        result[subject] = (plus, [h["start"] for h in plus], minus, [h["end"] for h in minus])
    sources[key] = hits, result
    return result


def directional(left_hits, right_hits, minimum, maximum):
    key = id(left_hits), id(right_hits), minimum, maximum
    if key in parts:
        return parts[key]
    left, right = prepared(left_hits), prepared(right_hits)
    result = []
    for subject in sorted(left.keys() & right.keys()):
        plus, starts, _, _ = left[subject]
        _, _, minus, ends = right[subject]
        if not plus or not minus:
            continue
        witness = None
        if len(plus) <= len(minus):
            for first in plus:
                lower = bisect_left(ends, first["start"] + minimum)
                upper = bisect_right(ends, first["start"] + maximum)
                for index in range(lower, upper):
                    last = minus[index]
                    if first["end"] <= last["start"]:
                        witness = first, last
                        break
                if witness:
                    break
        else:
            for last in minus:
                lower = bisect_left(starts, last["end"] - maximum)
                upper = bisect_right(starts, last["end"] - minimum)
                for index in range(lower, upper):
                    first = plus[index]
                    if first["end"] <= last["start"]:
                        witness = first, last
                        break
                if witness:
                    break
        if witness:
            first, last = witness
            result.append({"subject": subject, "start": first["start"], "end": last["end"],
                           "amplicon_length": last["end"] - first["start"]})
    parts[key] = result
    return result


def potential_amplicons(forward_hits, reverse_hits, minimum=36, maximum=1000):
    result = []
    for left_role, left in [("F", forward_hits), ("R", reverse_hits)]:
        for right_role, right in [("F", forward_hits), ("R", reverse_hits)]:
            result.extend({**row, "left_primer": left_role, "right_primer": right_role}
                          for row in directional(left, right, minimum, maximum))
    return result
