"""Full-length ungapped short-oligo alignment with explicit terminal matching."""
from collections import defaultdict
import re

from diffuprimer.sequence import reverse_complement
from diffuprimer.specificity import normalize_subject_id


def parse_sam_line(line, queries, max_mismatches=2, perfect_3prime=5):
        if line.startswith("@"):
            return None
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            raise ValueError("Malformed SAM alignment")
        qid, flag, subject, position, _, cigar = fields[:6]
        flag = int(flag)
        if flag & 4:
            return None
        original = queries[qid]
        reverse = bool(flag & 16)
        oriented = reverse_complement(original) if reverse else original
        if fields[9] != oriented or cigar != str(len(original)) + "M":
            raise ValueError("Audit requires full-length ungapped, untrimmed alignments")
        tags = {item.split(":", 2)[0]: item.split(":", 2)[2] for item in fields[11:]}
        md = tags.get("MD")
        if md is None or "^" in md:
            raise ValueError("Missing or gapped mismatch descriptor")
        offset, mismatches = 0, []
        for token in re.findall(r"\d+|[A-Z]", md):
            if token.isdigit():
                offset += int(token)
            else:
                mismatches.append(len(original) - 1 - offset if reverse else offset)
                offset += 1
        if offset != len(original):
            raise ValueError("Mismatch descriptor does not span the full oligo")
        if len(mismatches) > max_mismatches or any(p >= len(original) - perfect_3prime for p in mismatches):
            return None
        start = int(position) - 1
        return qid, {"qseqid": qid, "sseqid": normalize_subject_id(subject),
            "strand": "-" if reverse else "+", "start": start, "end": start + len(original),
            "mismatch": len(mismatches), "mismatch_positions_original": sorted(mismatches)}


def parse_sam(path, queries, max_mismatches=2, perfect_3prime=5):
    hits = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            parsed = parse_sam_line(line, queries, max_mismatches, perfect_3prime)
            if parsed is not None:
                qid, hit = parsed
                hits[qid].append(hit)
    return hits


def alignment_command(binary, index, queries, output, threads=4):
    # All mismatch strata and both orientations, with no hit-count ceiling.
    return [str(binary), "-v", "2", "-a", "--best", "--tryhard", "-f", "-S",
            "-p", str(threads), str(index), str(queries), str(output)]
