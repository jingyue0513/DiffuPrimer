from __future__ import annotations

import gzip
import hashlib
import re
import time
import threading
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import requests
from bs4 import BeautifulSoup
from Bio import SeqIO

from .io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from .sequence import map_primers


BASE = "https://pga.mgh.harvard.edu"
REFSEQ_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/"


class SourceAccessBlocked(RuntimeError):
    pass


class RateLimiter:
    def __init__(self, delay=1.0):
        self.delay = delay
        self.last_request = 0.0
        self.lock = threading.Lock()
        self.stopped = threading.Event()

    def wait(self):
        with self.lock:
            if self.stopped.is_set():
                raise SourceAccessBlocked("Source denied or limited access; no further requests are permitted in this run")
            time.sleep(max(0, self.delay - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()


class CachedHTTP:
    def __init__(self, cache, delay=1.0, limiter=None):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.delay = delay
        self.last_request = 0.0
        self.limiter = limiter or RateLimiter(delay)
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "DiffuPrimer academic research pilot/0.1 (cached, rate-limited requests)"

    def fetch(self, url, data=None):
        key = hashlib.sha256((url + repr(sorted((data or {}).items()))).encode()).hexdigest()
        path = self.cache / f"{key}.html"
        if path.exists() and path.with_suffix(".json").exists():
            return path.read_text(), str(path)
        for attempt in range(3):
            self.limiter.wait()
            try:
                response = self.session.post(url, data=data, timeout=40) if data is not None else self.session.get(url, timeout=40)
                if response.status_code in {401, 403, 429}:
                    self.limiter.stopped.set()
                    raise SourceAccessBlocked(f"HTTP {response.status_code} from {url}; stopped source requests")
                response.raise_for_status()
                if len(response.content) < 100:
                    raise ValueError("Empty or truncated server response")
                path.write_text(response.text)
                write_json(path.with_suffix(".json"), {"url": url, "post_data": data,
                    "retrieved_at": utc_now(), "sha256": sha256(path), "status_code": response.status_code})
                return response.text, str(path)
            except (requests.RequestException, ValueError):
                if attempt == 2:
                    raise
                time.sleep(2 ** (attempt + 1))


def table_fields(soup):
    fields = {}
    for row in soup.select("tr"):
        cells = row.find_all("td", recursive=False)
        if len(cells) >= 2:
            fields[cells[0].get_text(" ", strip=True)] = [cell.get_text(" ", strip=True) for cell in cells[1:]]
    return fields


def parse_search(html):
    soup = BeautifulSoup(html, "lxml")
    fields = table_fields(soup)
    gene = fields.get("NCBI GeneID", [None])[0]
    accession = fields.get("GenBank Accession", [None])[0]
    rows = []
    for table in soup.find_all("table"):
        values = table_fields(table)
        if not {"PrimerBank ID", "Forward Primer", "Reverse Primer", "Amplicon Size"} <= values.keys():
            continue
        # Only the innermost pair table; outer layout tables contain several pairs.
        if table.find("table") is not None:
            continue
        pid = values["PrimerBank ID"][0]
        validation = table.find("a", href=re.compile("displayResult"))
        rows.append({"gene_id": gene, "source_transcript_accession": accession,
                     "primerbank_id": pid, "forward_primer": values["Forward Primer"][0].upper(),
                     "reverse_primer": values["Reverse Primer"][0].upper(),
                     "amplicon_length": int(values["Amplicon Size"][0]),
                     "source_forward_location": values["Forward Primer"][-1],
                     "source_reverse_location": values["Reverse Primer"][-1],
                     "validation_url": BASE + "/rtpcr/displayResult.do?v=2&primerPairId=" + pid if validation else None})
    return rows


def parse_detail(html):
    soup = BeautifulSoup(html, "lxml")
    pre = soup.find("pre")
    if pre is None:
        raise ValueError("Source coding template missing")
    sequence = re.sub(r"[^ACGTNacgtn]", "", pre.get_text()).upper()
    if len(sequence) < 40:
        raise ValueError("Source coding template too short")
    return sequence


def parse_validation(html):
    """Extract observable evidence. A linked result page does NOT imply success."""
    soup = BeautifulSoup(html, "lxml")
    text = soup.get_text(" ", strip=True)
    source_version = re.search(r"ACCESSION\|([NX][MR]_\d+\.\d+)\|", text)
    image_sources = [x.get("src", "") for x in soup.find_all("img")]
    match = re.search(r"Percent Identity:\s*([\d.]+).*?Match Length:\s*(\d+)", text)
    return {
        "experimental_record_available": "QPCR Validation Results" in text,
        "amplification_plot_available": any("ampPlotId=" in x for x in image_sources),
        "dissociation_curve_available": any("dissCurveId=" in x for x in image_sources),
        "gel_image_available": any("egelId=" in x or "gelId=" in x for x in image_sources),
        "source_validation_accession_version": source_version.group(1) if source_version else None,
        "reported_blast_identity": float(match.group(1)) if match else None,
        "reported_blast_match_length": int(match.group(2)) if match else None,
        "experimental_result": "unknown", "evidence_type": "unverified",
        "label_reason": "No explicit assay outcome extracted; plots and links are evidence artifacts, not success labels.",
    }


def reference_gene_catalog(gff_path):
    genes, accessions = {}, {}
    with gzip.open(gff_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "mRNA"}:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            ids = re.search(r"(?:^|,)GeneID:(\d+)(?:,|$)", attrs.get("Dbxref", ""))
            if ids is None:
                continue
            gene_id = ids.group(1)
            if fields[2] == "gene" and attrs.get("gene_biotype") == "protein_coding":
                genes[gene_id] = {"gene_id": gene_id, "gene_symbol": attrs.get("Name"), "reference_assembly": "GCF_000001635.27"}
            elif fields[2] == "mRNA" and attrs.get("transcript_id", "").startswith("NM_"):
                accessions.setdefault(gene_id, []).append(attrs["transcript_id"])
    return [{**genes[gid], "reference_transcripts": sorted(set(accessions[gid]))}
            for gid in sorted(genes, key=int) if gid in accessions]


def load_reference_sequences(fasta_path, accessions):
    selected = {}
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in accessions:
                selected[record.id] = str(record.seq).upper()
    return selected


def remap_record(record, source_template, references):
    source_matches = map_primers(source_template, record["forward_primer"], record["reverse_primer"], record["amplicon_length"])
    if len(source_matches) != 1:
        return {**record, "mapping_status": "ambiguous_source" if source_matches else "unmapped_source"}
    preferred_base = record.get("source_transcript_accession", "")
    candidates = []
    for accession, sequence in references.items():
        pairs = map_primers(sequence, record["forward_primer"], record["reverse_primer"], record["amplicon_length"])
        if len(pairs) == 1:
            candidates.append((not accession.startswith(preferred_base + "."), accession, sequence, pairs[0]))
    if not candidates:
        return {**record, "mapping_status": "unmapped_current_reference"}
    _, accession, sequence, pair = sorted(candidates)[0]
    return {**record, "mapping_status": "unique_exact_pair_remap",
            "template_sequence": sequence, "transcript_accession_version": accession,
            **pair.as_dict(sequence), "source_template_sequence": source_template,
            "source_template_sha256": hashlib.sha256(source_template.encode()).hexdigest(),
            "source_template_exact_subsequence": source_template in sequence,
            "source_pair_coordinates": source_matches[0].__dict__,
            "reference_assembly": "GCF_000001635.27", "database_version": "live snapshot with retrieval timestamp and hash",
            "annotation_version": "GCF_000001635.27-RS_2024_02"}


def acquire_gene(gene, http, reference_sequences):
    html, search_path = http.fetch(BASE + "/cgi-bin/primerbank/new_search2.cgi",
        {"selectBox": "NCBI Gene Symbol", "species": "Mouse", "searchBox": gene["gene_symbol"], "Submit": "Submit"})
    output = []
    for record in parse_search(html):
        record.update({"species": "Mus musculus", "requested_gene_id": gene["gene_id"],
                       "source_database": "PrimerBank", "search_snapshot": search_path})
        if record["gene_id"] != gene["gene_id"]:
            record.update({"mapping_status": "gene_resolution_mismatch", "evidence_type": "unverified", "experimental_result": "unknown"})
            output.append(record)
            continue
        url = BASE + "/cgi-bin/primerbank/new_displayDetail2.cgi?primerID=" + record["primerbank_id"]
        detail, detail_path = http.fetch(url)
        record.update({"source_url": url, "detail_snapshot": detail_path,
                       "evidence_type": "unverified", "experimental_result": "unknown",
                       "experimental_record_available": False, "failure_reason": None})
        if record["validation_url"]:
            validation, path = http.fetch(record["validation_url"])
            record.update(parse_validation(validation))
            record["validation_snapshot"] = path
        subset = {acc: reference_sequences[acc] for acc in gene["reference_transcripts"] if acc in reference_sequences}
        record = remap_record(record, parse_detail(detail), subset)
        output.append(record)
    return output
