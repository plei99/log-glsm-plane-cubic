"""Regression tests for full base-torus localization of zero vertices."""

load("o3_fixed_locus_graphs.sage")

import os
import tempfile


def run_tests():
    rings = PlaneCubicCoefficientRing()
    lambdas = rings.base_weights()

    divisor_lift = TwistedInsertion.H_vanishing_at(1)
    assert divisor_lift.restriction(1, lambdas) == 0
    assert divisor_lift.restriction(0, lambdas) == lambdas[0] - lambdas[1]
    point_lift = TwistedInsertion.point_supported_at(2)
    assert point_lift.allowed_fixed_points() == (2,)
    assert point_lift.restriction(0, lambdas) == 0
    assert point_lift.restriction(2, lambdas) == (
        (lambdas[2] - lambdas[0]) * (lambdas[2] - lambdas[1])
    )

    degree_one = P2FixedLocusGraphEnumerator(0, 0, 1)
    assert degree_one.graph_count() == 3
    assert all(graph.automorphism_order() == 1
               for graph in degree_one.graphs())
    assert all(graph.deck_order() == 1 for graph in degree_one.graphs())

    parallel = [
        graph for graph in P2FixedLocusGraphEnumerator(1, 0, 2).graphs()
        if graph.vertex_count == 2 and graph.edge_count == 2
    ]
    assert len(parallel) == 3
    assert all(graph.automorphism_order() == 2 for graph in parallel)

    ordinary = P2FixedLocusEvaluator(
        rings, include_twist=False, lift_strategy="sparse"
    )
    numeric = P2FixedLocusEvaluator(
        rings, include_twist=False, lift_strategy="sparse",
        base_weight_specialization=(0, 1, 3),
    )
    constant_request = TwistedZeroVertexRequest(
        0, 0, (TwistedInsertion(0),) * 3
    )
    assert numeric.evaluate(constant_request) == rings.specialize_base_weights(
        ordinary.evaluate(constant_request), (0, 1, 3)
    )
    nonequivariant = P2FixedLocusEvaluator(
        rings, include_twist=True, lift_strategy="sparse",
        base_weight_specialization="nonequivariant",
    )
    epsilon = rings.full_lambda0
    assert nonequivariant.lambdas == (0, epsilon, epsilon ** 2)

    # The former linear ray epsilon*(0,1,3) makes the two flag weights at
    # the nodal colour-1 vertex cancel: epsilon + (-2*epsilon)/2 = 0.
    # This concrete degree-three locus guards the nonresonant curved path.
    resonant_request = TwistedZeroVertexRequest(
        0, 3, (TwistedInsertion(0),)
    )
    resonant_graph = next(
        graph for graph in nonequivariant.fixed_graphs(resonant_request)
        if graph.vertex_colours == (0, 2, 1)
        and graph.edges == ((0, 2, 1), (1, 2, 2))
    )
    nodal_vertex = 2
    flag_sum = sum(
        resonant_graph.flag_weight(
            edge_index, nodal_vertex, nonequivariant.lambdas
        )
        for edge_index in resonant_graph.incident_edges(nodal_vertex)
    )
    assert flag_sum == epsilon * (3 - epsilon) / 2
    assert flag_sum != 0
    standard = P2FixedLocusEvaluator(
        rings, include_twist=False, lift_strategy="standard"
    )
    line = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(2), TwistedInsertion(2),
    ))
    same_fixed_loci = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(2, scale=2, label="scaled"),
        TwistedInsertion(2, label="renamed"),
    ))
    assert ordinary.fixed_graphs(same_fixed_loci) is ordinary.fixed_graphs(line)
    assert ordinary.graph_count(line) == 1
    assert standard.graph_count(line) == 12
    assert ordinary.evaluate(line) == 1
    assert standard.evaluate(line) == 1
    cache_entries = ordinary.cache_info()["entries"]
    assert ordinary.evaluate(same_fixed_loci) == 2
    assert ordinary.cache_info()["entries"] == cache_entries
    assert ordinary.cache_info()["hits"] >= 1

    with tempfile.TemporaryDirectory() as directory:
        cache_path = os.path.join(directory, "zero-vertices.json")
        persistent = P2FixedLocusEvaluator(
            rings, include_twist=False, lift_strategy="sparse",
            cache_path=cache_path, autosave=True,
        )
        assert persistent.evaluate(line) == 1
        assert os.path.exists(cache_path)
        assert persistent.cache_info()["misses"] == 1
        restored = P2FixedLocusEvaluator(
            rings, include_twist=False, lift_strategy="sparse",
            cache_path=cache_path,
        )
        assert restored.is_cached(line)
        assert restored.evaluate(line) == 1
        assert restored.cache_info()["hits"] == 1

        sqlite_path = os.path.join(directory, "zero-vertices.sqlite")
        sqlite_cache = P2FixedLocusEvaluator(
            rings, include_twist=False, lift_strategy="sparse",
            cache_path=sqlite_path, autosave=True,
        )
        assert sqlite_cache.evaluate(line) == 1
        assert sqlite_cache.cache_info()["cache_backend"] == "sqlite"
        assert sqlite_cache.cache_info()["persistent_entries"] == 1
        sqlite_cache._sqlite_store.close()

        sqlite_restored = P2FixedLocusEvaluator(
            rings, include_twist=False, lift_strategy="sparse",
            cache_path=sqlite_path,
        )
        relabelled_line = TwistedZeroVertexRequest(0, 1, (
            TwistedInsertion(2, label="first"),
            TwistedInsertion(2, label="second"),
        ))
        assert sqlite_restored.evaluate(relabelled_line) == 1
        assert sqlite_restored.cache_info()["hits"] == 1
        sqlite_restored._sqlite_store.close()

    line_divisor = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(2), TwistedInsertion(2), TwistedInsertion(1),
    ))
    planned_line_divisor = ordinary.planned_insertions(line_divisor)
    assert tuple(item.lift_kind for item in planned_line_divisor) == (
        "vanish", "support", "support"
    )
    assert ordinary.graph_count(line_divisor) == 1
    assert ordinary.evaluate(line_divisor) == 1

    # String and dilaton equations test the sign of psi on an unstable
    # marked endpoint as well as the stable flag expansions.
    string_request = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(0), TwistedInsertion(2, 1), TwistedInsertion(2),
    ))
    assert ordinary.evaluate(string_request) == 1
    dilaton_request = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(0, 1), TwistedInsertion(2), TwistedInsertion(2),
    ))
    assert ordinary.evaluate(dilaton_request) == 0

    conic = TwistedZeroVertexRequest(
        0, 2, tuple(TwistedInsertion(2) for _ in range(5))
    )
    assert P2FixedLocusGraphEnumerator(0, 5, 2).graph_count() == 1557
    assert ordinary.graph_count(conic) == 3
    assert standard.graph_count(conic) == 1557
    assert ordinary.evaluate(conic) == 1
    assert standard.evaluate(conic) == 1

    # Mapping-to-a-point: int_[Mbar_1,1 x P2] H e(E^* tensor TP2)=-1/8.
    constant_genus_one = TwistedZeroVertexRequest(
        1, 0, (TwistedInsertion(1),)
    )
    localized_constant_genus_one = sum(
        ordinary.graph_contribution(
            graph, ordinary.planned_insertions(constant_genus_one)
        )
        for graph in ordinary.fixed_graphs(constant_genus_one)
    )
    assert ordinary.evaluate(constant_genus_one) == -QQ(1) / 8
    assert ordinary.evaluate(constant_genus_one) \
        == localized_constant_genus_one

    twisted = P2FixedLocusEvaluator(rings, include_twist=True)
    constant_rational = TwistedZeroVertexRequest(0, 0, (
        TwistedInsertion(1), TwistedInsertion(1), TwistedInsertion(0),
    ))
    value = twisted.evaluate(constant_rational)
    assert value(0, 0, 0, rings.t) == rings.t

    cubic_rational = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(1), TwistedInsertion(0), TwistedInsertion(0),
    ))
    assert twisted.evaluate(cubic_rational) == 0

    full = FullTwistedZeroVertexBackend(rings)
    assert full.is_cached(constant_rational)
    assert full.evaluate(constant_rational) == twisted.evaluate(
        constant_rational
    )
    assert full.evaluate(constant_genus_one) == twisted.evaluate(
        constant_genus_one
    )
    genus_two = TwistedZeroVertexRequest(2, 0, ())
    localized_genus_two = sum(
        twisted.graph_contribution(
            graph, twisted.planned_insertions(genus_two)
        )
        for graph in twisted.fixed_graphs(genus_two)
    )
    assert full.evaluate(genus_two) != 0
    assert full.evaluate(genus_two) == localized_genus_two
    assert "constant-map" in full.provenance(genus_two)

    # The all-genus constant-map path uses the closed top-Hodge-triple
    # formula in the first previously pathological genus-three case.
    genus_three = TwistedZeroVertexRequest(3, 0, ())
    localized_genus_three = sum(
        twisted.graph_contribution(
            graph, twisted.planned_insertions(genus_three)
        )
        for graph in twisted.fixed_graphs(genus_three)
    )
    assert full.evaluate(genus_three) == localized_genus_three
    assert full.cache_info()["hodge_cache"]["lambda_product_entries"] == 0


run_tests()
print("all full O(3)-twisted fixed-locus tests passed")
