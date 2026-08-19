"""Regression tests for the finite plane-cubic CJR graph enumerator."""

load("cjr_bipartite_graphs.sage")


def run_tests():
    # A constant genus-one connected stable-map probe needs a special point.
    # The global guard runs before CJR's vertex classification, which is
    # stated only inside globally X-stable data.
    assert PlaneCubicGraphEnumerator(1, 0, 0).graphs() == tuple()

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


def _brute_force_endpoint_multigraphs(zero_count, infinity_count, edge_count):
    """The pre-canonicalization search: every labeled endpoint multigraph."""
    cell_count = zero_count * infinity_count
    for multiplicities in _weak_compositions(edge_count, cell_count):
        pairs = []
        for flat_index, multiplicity in enumerate(multiplicities):
            zero = flat_index // infinity_count
            infinity = flat_index % infinity_count
            pairs.extend([(zero, infinity)] * multiplicity)
        if _connected_bipartite(zero_count, infinity_count, pairs):
            yield tuple(pairs)


def test_canonical_incidence_matrices_lose_no_orbit():
    r"""Row canonicalization must keep one representative of every orbit.

    The generator emits incidence matrices only up to relabeling of the zero
    vertices.  Sorting each brute-force matrix's rows must land exactly on
    the emitted set, which is what makes the reduction lossless once every
    decoration is enumerated downstream.
    """
    for zero_count in range(1, 5):
        for infinity_count in range(1, 3):
            for edge_count in range(1, 6):
                enumerator = PlaneCubicGraphEnumerator(2, 0, 0)
                emitted = set(enumerator._canonical_incidence_matrices(
                    zero_count, infinity_count, edge_count
                ))
                # Every emitted matrix is in canonical order.
                for matrix in emitted:
                    keys = [(sum(row), row) for row in matrix]
                    assert keys == sorted(keys, reverse=True)
                    assert all(sum(row) >= 1 for row in matrix)

                brute = set()
                for pairs in _brute_force_endpoint_multigraphs(
                        zero_count, infinity_count, edge_count):
                    matrix = [[0] * infinity_count
                              for _ in range(zero_count)]
                    for zero, infinity in pairs:
                        matrix[zero][infinity] += 1
                    rows = [tuple(Integer(x) for x in row) for row in matrix]
                    rows.sort(key=lambda row: (sum(row), row), reverse=True)
                    brute.add(tuple(rows))
                # A disconnected brute-force matrix is never emitted, and the
                # emitted set applies its own connectivity filter later, so
                # compare on the connected orbits only.
                assert brute.issubset(emitted), (
                    zero_count, infinity_count, edge_count,
                    sorted(brute - emitted),
                )


def test_fast_zero_vertex_types_match_graph_classification():
    """The inner-loop classifier must exactly match the public graph API."""
    examples = (
        ((0,), (0,), ((0, 0),), ()),
        ((0,), (0,), ((0, 0),), (0,)),
        ((0,), (0,), ((0, 0), (0, 0)), ()),
        ((0,), (0,), ((0, 0),), (0, 0)),
        ((1,), (0,), ((0, 0),), ()),
        ((0,), (2,), ((0, 0),), ()),
        ((0, 0, 1), (0, 0, 0),
         ((0, 0), (1, 0), (1, 0), (2, 0)), (0,)),
    )
    for zero_genera, zero_degrees, pairs, marking_vertices in examples:
        fast = PlaneCubicGraphEnumerator._zero_vertex_types(
            zero_genera, zero_degrees, pairs, marking_vertices
        )
        graph = BipartiteLocalizationGraph(
            zero_genera, zero_degrees, (3,), (0,), marking_vertices,
            tuple((zero, infinity, 1) for zero, infinity in pairs),
        )
        public = tuple(
            graph.zero_vertex_type(zero)
            for zero in range(graph.zero_count)
        )
        assert fast == public, (fast, public)


def test_canonical_marking_decorations():
    """Interchangeable zero vertices use restricted-growth assignments."""
    equivalent = tuple(
        PlaneCubicGraphEnumerator._canonical_marking_decorations(
            2, (0, 0), (0, 0), ((0, 0), (1, 0))
        )
    )
    assert equivalent == ((0, 0), (0, 1))

    distinguished = tuple(
        PlaneCubicGraphEnumerator._canonical_marking_decorations(
            2, (0, 1), (0, 0), ((0, 0), (1, 0))
        )
    )
    assert set(distinguished) == {(0, 0), (0, 1), (1, 0), (1, 1)}


def test_genus_three_enumeration():
    """Genus three is reachable and exposes the positive-degree sector."""
    unmarked = enumerate_plane_cubic_graphs(3, 0, 0)
    assert len(unmarked) == 32
    assert all(graph.is_valid() for graph in unmarked)
    assert len({graph.canonical_key() for graph in unmarked}) == 32

    marked = enumerate_plane_cubic_graphs(3, 1, 0)
    assert len(marked) == 85
    assert all(graph.is_valid() for graph in marked)
    assert len({graph.canonical_key() for graph in marked}) == 85

    dependencies = sorted(set(
        data for graph in marked
        for data in graph.all_infinity_vertex_data()
    ))
    # Genus three is the first genus admitting a positive infinity degree,
    # but a degree-zero probe cannot carry one: 3*D <= 2g-2 needs the degree
    # to come from somewhere.
    assert all(degree == 0 for _, degree, _ in dependencies)
    genus_three = [item for item in dependencies if item[0] == 3]
    assert (3, 0, (-5,)) in genus_three
    assert (3, 0, (-4, -2)) in genus_three
    assert (3, 0, (-3, -3)) in genus_three
    # Every occurring profile satisfies the balance equation.
    for genus, degree, contacts in dependencies:
        assert 3 * degree - (2 * genus - 2) == sum(c + 1 for c in contacts)


run_tests()
test_canonical_incidence_matrices_lose_no_orbit()
test_fast_zero_vertex_types_match_graph_classification()
test_canonical_marking_decorations()
test_genus_three_enumeration()
print("all CJR bipartite-graph enumerator tests passed")
