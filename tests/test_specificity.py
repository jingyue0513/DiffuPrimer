from diffuprimer.specificity import classify_transcript_products, potential_amplicons, qualifying_hit, normalize_subject_id


def test_only_inward_pairs_generate_products_and_isoforms_are_gene_grouped():
    f = {"sseqid": "NM_1.1", "strand": "+", "start": 10, "end": 30}
    r = {"sseqid": "NM_1.1", "strand": "-", "start": 80, "end": 100}
    assert potential_amplicons([f], [dict(r, strand="+")]) == []
    products = potential_amplicons([f], [r])
    assert len(products) == 1
    assert products[0]["amplicon_length"] == 90
    assert classify_transcript_products(products, "1", {"NM_1.1": "1"})["gene_specificity_screen_pass"]
    assert not classify_transcript_products(products, "2", {"NM_1.1": "1"})["gene_specificity_screen_pass"]
    assert not classify_transcript_products(products, "1", {})["gene_specificity_screen_pass"]


def test_terminal_mismatch_and_partial_hits_are_not_silently_exact():
    hit = {"qstart": 1, "qend": 20, "qlen": 20, "gaps": 0, "mismatch": 1,
           "qseq": "A" * 20, "sseq": "C" + "A" * 19}
    assert qualifying_hit(hit)
    assert not qualifying_hit(dict(hit, sseq="A" * 19 + "C"))
    assert not qualifying_hit(dict(hit, qstart=2))


def test_same_primer_pair_products_are_screened():
    plus = {"sseqid": "offtarget", "strand": "+", "start": 1, "end": 21}
    minus = {"sseqid": "offtarget", "strand": "-", "start": 70, "end": 90}
    products = potential_amplicons([plus, minus], [])
    assert products[0]["left_primer"] == products[0]["right_primer"] == "F"


def test_blast_reference_accessions_match_annotation_keys():
    assert normalize_subject_id("ref|NM_133883.2|") == "NM_133883.2"
    assert normalize_subject_id("ref|NC_000067.7|") == "NC_000067.7"
    assert normalize_subject_id("NM_133883.2") == "NM_133883.2"
    assert normalize_subject_id("custom_contig") == "custom_contig"
