"""Regression tests for full base-torus localization of zero vertices."""

load("o3_fixed_locus_graphs.sage")


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
    standard = P2FixedLocusEvaluator(
        rings, include_twist=False, lift_strategy="standard"
    )
    line = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(2), TwistedInsertion(2),
    ))
    assert ordinary.graph_count(line) == 1
    assert standard.graph_count(line) == 12
    assert ordinary.evaluate(line) == 1
    assert standard.evaluate(line) == 1

    line_divisor = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(2), TwistedInsertion(2), TwistedInsertion(1),
    ))
    planned_line_divisor = ordinary.planned_insertions(line_divisor)
    assert tuple(item.lift_kind for item in planned_line_divisor) == (
        "support", "support", "vanish"
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
    assert ordinary.evaluate(constant_genus_one) == -QQ(1) / 8

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
    assert full.evaluate(constant_genus_one) == twisted.evaluate(
        constant_genus_one
    )
    genus_two = TwistedZeroVertexRequest(2, 0, ())
    assert full.evaluate(genus_two) != 0
    assert "admcycles" in full.provenance(genus_two)


run_tests()
print("all full O(3)-twisted fixed-locus tests passed")
