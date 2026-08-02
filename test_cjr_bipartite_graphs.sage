"""Regression tests for the finite plane-cubic CJR graph enumerator."""

load("cjr_bipartite_graphs.sage")


def run_tests():
    assert list(_weak_compositions(3, 2)) == [
        (0, 3), (1, 2), (2, 1), (3, 0)
    ]
    assert list(_positive_compositions(4, 2)) == [
        (1, 3), (2, 2), (3, 1)
    ]

    # Figure 4 of CJR I: the six degree-zero genus-two graph types without
    # ordinary legs.  Their automorphism orders are 1,2,2,2,1,1.
    genus_two = PlaneCubicGraphEnumerator(2, 0, 0)
    graphs = genus_two.graphs()
    assert genus_two.mixed_vertex_bound() == 7
    assert len(graphs) == 6
    assert sorted(graph.automorphism_order() for graph in graphs) == [1, 1, 1, 2, 2, 2]
    assert sum(graph.infinity_count == 0 for graph in graphs) == 1
    assert sum(graph.infinity_count == 2 for graph in graphs) == 1
    assert all(graph.is_valid() for graph in graphs)
    assert all(graph.total_genus() == 2 for graph in graphs)
    assert all(graph.total_degree() == 0 for graph in graphs)
    assert len({graph.canonical_key() for graph in graphs}) == 6

    contact_profiles = sorted(
        tuple(sorted(c for _, _, c in graph.edges)) for graph in graphs
    )
    assert contact_profiles == [(), (1,), (1, 1), (1, 1), (2, 2), (3,)]
    two_infinity = [graph for graph in graphs
                    if graph.infinity_count == 2][0]
    two_tails = [graph for graph in graphs
                 if graph.zero_count == 2][0]
    assert (two_infinity.vertex_automorphism_order(),
            two_infinity.edge_automorphism_order()) == (2, 1)
    assert (two_tails.vertex_automorphism_order(),
            two_tails.edge_automorphism_order()) == (2, 1)

    # The marked genus-one sector consists of the pure zero graph and the
    # contact-one marked-unstable/infinity graph.
    genus_one_marked = enumerate_plane_cubic_graphs(1, 1, 0)
    assert len(genus_one_marked) == 2
    mixed = [graph for graph in genus_one_marked if graph.infinity_count][0]
    assert mixed.zero_vertex_type(0) == "marked"
    assert mixed.infinity_balance(0) == (1, 1)
    assert mixed.infinity_contact_profile(0) == (-1,)
    assert mixed.infinity_vertex_data(0) == (1, 0, (-1,))
    assert mixed.localization_weight() == 1

    # CJR III, Example 9.1 and Figure 1: with two compact markings there are
    # exactly three types (pure zero, one stable zero bridge, two marked tails).
    genus_one_two_markings = enumerate_plane_cubic_graphs(1, 2, 0)
    assert len(genus_one_two_markings) == 3
    assert sorted(graph.zero_count for graph in genus_one_two_markings) == [1, 1, 2]

    # Positive ambient degree may lie at zero level.  The degree-one,
    # genus-one unmarked sector has its one-level and mixed representatives.
    degree_one = enumerate_plane_cubic_graphs(1, 0, 1)
    assert len(degree_one) == 2
    assert all(graph.total_degree() == 1 for graph in degree_one)

    # Canonical labeling forgets the chosen names of vertices but remembers
    # labeled markings, contacts, decorations, and parallel-edge symmetry.
    original = BipartiteLocalizationGraph(
        (0, 0), (0, 0), (2,), (0,), (0,),
        ((0, 0, 2), (1, 0, 2)),
    )
    swapped = BipartiteLocalizationGraph(
        (0, 0), (0, 0), (2,), (0,), (1,),
        ((0, 0, 2), (1, 0, 2)),
    )
    assert original.is_valid() and swapped.is_valid()
    assert original.canonical_key() == swapped.canonical_key()
    assert original == swapped
    assert hash(original) == hash(swapped)
    # The marking distinguishes the leaves, so the automorphism is trivial.
    assert original.automorphism_order() == 1

    double_edge = BipartiteLocalizationGraph(
        (0,), (0,), (1,), (0,), (),
        ((0, 0, 1), (0, 0, 1)),
    )
    assert double_edge.is_valid()
    assert double_edge.zero_vertex_type(0) == "nodal"
    assert double_edge.edge_automorphism_order() == 2
    assert double_edge.vertex_automorphism_order() == 1
    assert double_edge.automorphism_order() == 2

    # CJR nonspecial zero-level tails require contact strictly greater than 1.
    forbidden_contact_one_tail = BipartiteLocalizationGraph(
        (0,), (0,), (1,), (0,), (), ((0, 0, 1),)
    )
    assert not forbidden_contact_one_tail.is_valid()
    assert any("contact > 1" in error
               for error in forbidden_contact_one_tail.validation_errors())

    unbalanced = BipartiteLocalizationGraph(
        (0,), (1,), (1,), (0,), (), ((0, 0, 2),)
    )
    assert not unbalanced.is_valid()
    assert any("violates balance" in error
               for error in unbalanced.validation_errors())

    malformed = BipartiteLocalizationGraph(
        (0,), (), (1,), (0,), (2,), ((3, 4, 1),)
    )
    assert not malformed.is_valid()
    assert "zero genus/degree arrays have different lengths" in \
        malformed.validation_errors()

    # A somewhat larger search exercises labeled markings, cycles, multiple
    # infinity vertices, and canonical deduplication.
    genus_two_marked = enumerate_plane_cubic_graphs(2, 1, 0)
    assert len(genus_two_marked) == 11
    assert len({graph.canonical_key() for graph in genus_two_marked}) == 11
    assert all(graph.is_valid() for graph in genus_two_marked)
    dependencies = sorted(set(
        data for graph in genus_two_marked
        for data in graph.all_infinity_vertex_data()
    ))
    assert dependencies == [
        (1, 0, (-1,)),
        (1, 0, (-1, -1)),
        (1, 0, (-1, -1, -1)),
        (2, 0, (-3,)),
        (2, 0, (-3, -1)),
        (2, 0, (-2, -2)),
        (2, 0, (-2, -2, -1)),
    ]


run_tests()
print("all CJR bipartite-graph enumerator tests passed")
