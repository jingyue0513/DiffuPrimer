from diffuprimer.splits import assign_cluster_splits


def test_similarity_chain_cannot_leak_and_is_reproducible():
    genes = [str(i) for i in range(40)]
    edges = [("1", "2"), ("2", "3"), ("5", "6")]
    first = assign_cluster_splits(genes, edges, seed=18)
    assert first == assign_cluster_splits(genes, reversed(edges), seed=18)
    assignments = {x["gene_id"]: x for x in first}
    assert len(first) == 40
    assert assignments["1"]["split"] == assignments["3"]["split"]
    assert assignments["1"]["cluster_id"] == assignments["3"]["cluster_id"]
    assert set(x["split"] for x in first) == {"train", "validation", "test"}
