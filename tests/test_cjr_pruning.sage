"""Compare public graph families with a deliberately unpruned small search."""

from itertools import product
from log_glsm.infinity.graphs import (
    BipartiteLocalizationGraph,
    PlaneCubicGraphEnumerator,
    _connected_bipartite,
    _positive_compositions,
    _weak_compositions,
)


def reference_graphs(g, n, d):
    # The former, larger bound is independent of the new contact budget.
    bound = 4*g + 2*n + d - 1
    found = {}
    for graph in PlaneCubicGraphEnumerator(g, n, d)._one_level_graphs():
        found[graph.canonical_key()] = graph
    for ni in range(1, g + 1):
        for nz in range(1, bound - ni + 1):
            nv = nz + ni
            for ne in range(nv - 1, nv + g - ni):
                betti = ne - nv + 1
                # Enumerate all endpoint multiplicities, without orbit or
                # leaf pruning; validation is the reference specification.
                for counts in _weak_compositions(ne, nz * ni):
                    pairs = tuple((z, i) for z in range(nz) for i in range(ni)
                                  for _ in range(counts[z*ni+i]))
                    if not _connected_bipartite(nz, ni, pairs):
                        continue
                    for extra in _weak_compositions(g-betti-ni, nv):
                        zg = extra[:nz]
                        ig = tuple(1+x for x in extra[nz:])
                        for degrees in _weak_compositions(d, nv):
                            choices = []
                            for i in range(ni):
                                val = sum(j == i for _, j in pairs)
                                choices.append(tuple(_positive_compositions(
                                    2*ig[i]-2+val-3*degrees[nz+i], val)))
                            for local in product(*choices):
                                offsets = [0]*ni
                                edges = []
                                for z, i in pairs:
                                    edges.append((z, i, local[i][offsets[i]]))
                                    offsets[i] += 1
                                for marks in product(range(nz), repeat=n):
                                    graph = BipartiteLocalizationGraph(
                                        zg, degrees[:nz], ig, degrees[nz:],
                                        marks, tuple(edges))
                                    if graph.is_valid():
                                        found[graph.canonical_key()] = graph
    return found


for data in ((1, 1, 0), (1, 0, 1), (2, 0, 0)):
    reference = reference_graphs(*data)
    actual = PlaneCubicGraphEnumerator(*data).graphs()
    assert {g.canonical_key(): g.automorphism_order() for g in actual} == {
        key: g.automorphism_order() for key, g in reference.items()}
    assert all(g.total_genus() == data[0] and g.total_degree() == data[2]
               for g in actual)

print("all independent CJR pruning comparisons passed")
