from __future__ import annotations

import random


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def assign_cluster_splits(genes, edges, fractions=(0.7, 0.15, 0.15), seed=20261013):
    if len(fractions) != 3 or abs(sum(fractions) - 1) > 1e-8 or min(fractions) <= 0:
        raise ValueError("Require three positive fractions summing to one")
    uf = UnionFind(genes)
    for a, b in edges:
        uf.union(a, b)
    groups = {}
    for gene in sorted(genes):
        groups.setdefault(uf.find(gene), []).append(gene)
    clusters = list(groups.values())
    random.Random(seed).shuffle(clusters)
    clusters.sort(key=len, reverse=True)
    targets = [fraction * len(genes) for fraction in fractions]
    counts = [0, 0, 0]
    output = []
    for cluster in clusters:
        # Place a whole component into the partition with most remaining capacity.
        choice = max(range(3), key=lambda i: targets[i] - counts[i])
        counts[choice] += len(cluster)
        name = ("train", "validation", "test")[choice]
        for gene in cluster:
            output.append({"gene_id": gene, "cluster_id": min(cluster), "split": name})
    assignment = {x["gene_id"]: x["split"] for x in output}
    assert all(assignment[a] == assignment[b] for a, b in edges)
    return sorted(output, key=lambda x: int(x["gene_id"]))
