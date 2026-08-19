"""Tests for exact per-graph CJR contribution compilation."""

load("cjr_graph_contributions.sage")


def run_tests():
    rings = PlaneCubicCoefficientRing(12)
    compiler = PlaneCubicGraphContributionCompiler(rings)
    assert compiler.base_weight_specialization == "nonequivariant"
    try:
        compiler.compile_probe(ProbeSpec(1, 0, ()))
        raise AssertionError("an unmarked constant genus-one probe is unstable")
    except ValueError:
        pass

    # Auxiliary base-torus weights are specialized only after the complete
    # fixed-locus sum.  CJR coefficient extraction is at t=infinity, hence in
    # QQ((t^(-1))).
    l0, l1, l2 = rings.base_weights()
    t = rings.t
    rational = (l0 + l1 * t + t**2) / (l2 + t)
    specialized = rings.specialize_base_weights(rational, (0, 1, 3))
    nonequivariant = rings.nonequivariant_base_limit(rational)
    for power in range(4):
        assert compiler.laurent_coefficient(rational, power) \
            == rings.laurent_coefficient_at_infinity(nonequivariant, power)

    numeric_compiler = PlaneCubicGraphContributionCompiler(
        rings, base_weight_specialization=(0, 1, 3)
    )
    for power in range(4):
        assert numeric_compiler.laurent_coefficient(rational, power) \
            == rings.laurent_coefficient_at_infinity(specialized, power)

    symbolic_compiler = PlaneCubicGraphContributionCompiler(
        rings, base_weight_specialization=None
    )
    assert symbolic_compiler.laurent_coefficient(rational, 0) \
        == rings.laurent_coefficient_at_infinity(rational, 0)

    genus_one_probe = ProbeSpec.stationary(1, 0, (0,))
    cached_genus_one_graphs = compiler._graphs_for_probe(genus_one_probe)
    alternate_genus_one_probe = ProbeSpec(
        1, 0, (EllipticInsertion("unit", 1),)
    )
    assert compiler._graphs_for_probe(alternate_genus_one_probe) \
        is cached_genus_one_graphs

    stabilized_descendant = ProbeSpec.stationary(2, 0, (2,))
    stabilized_compilation = compiler.compile_probe(stabilized_descendant)
    assert stabilized_compilation.is_complete
    stabilized_vertex = EffectiveVertex(
        2, 0, (-2, -2, -1), psi_min=0,
        insertions=(0, 0, 1), contact_psi=(0, 0, 2),
    )
    assert stabilized_compilation.polynomial.terms == {
        (stabilized_vertex,): rings.full(-QQ(1) / 6)
    }
    assert stabilized_vertex.is_dimension_zero()

    # The same discrete descendant with CJR's log-domain convention retains
    # equation (8.21) on marked zero tails and therefore has no contact psi.
    log_descendant = ProbeSpec(
        2, 0, (EllipticInsertion("point", 2),),
        psi_convention="log",
    )
    log_compilation = compiler.compile_probe(log_descendant)
    assert log_compilation.is_complete
    assert log_compilation.polynomial.terms != \
        stabilized_compilation.polynomial.terms
    assert all(
        not vertex.contact_psi
        for factors in log_compilation.polynomial.terms
        for vertex in factors
    )
    genus_one_graphs = PlaneCubicGraphEnumerator(1, 1, 0).graphs()
    mixed = [graph for graph in genus_one_graphs if graph.infinity_count][0]
    contribution = compiler.compile_graph(genus_one_probe, mixed)
    assert contribution.supported
    basic = EffectiveVertex(1, 0, (-1,), insertions=(1,))
    assert contribution.polynomial.terms == {(basic,): rings.full(-QQ(1)/3)}

    # Equal infinity insertion states are combined before EffectiveVertex
    # construction.  This sparse aggregation prevents the high-degree flag
    # Cartesian product from repeating identical full-field coercions.
    synthetic_options = (((
        rings.full(2), {0: (ZZ.one(), ZZ.zero())}
    ), (
        rings.full(3), {0: (ZZ.one(), ZZ.zero())}
    )),)
    combined = compiler._combined_zero_states(mixed, synthetic_options)
    assert combined == ((
        rings.full(5) * rings.full(mixed.localization_weight()),
        {0: (ZZ.one(), ZZ.zero())},
    ),)

    pure_zero = [graph for graph in genus_one_graphs if not graph.infinity_count][0]
    pure_zero_contribution = compiler.compile_graph(genus_one_probe, pure_zero)
    assert pure_zero_contribution.supported
    assert tuple() in pure_zero_contribution.polynomial.terms
    assert pure_zero_contribution.polynomial.terms[tuple()] != 0
    assert compiler.laurent_coefficient(
        pure_zero_contribution.polynomial.terms[tuple()], 0
    ) == -QQ(1) / 24

    # The old partial backend remains injectable and keeps its exact failure
    # diagnostics; it is no longer the compiler default.
    partial_compiler = PlaneCubicGraphContributionCompiler(
        rings, TwistedZeroVertexBackend(rings)
    )
    unsupported = partial_compiler.compile_graph(genus_one_probe, pure_zero)
    assert not unsupported.supported
    assert "no exact individual O(3)-twisted" in unsupported.unsupported[0]

    genus_two_probe = ProbeSpec(2, 0, ())
    genus_two_graphs = PlaneCubicGraphEnumerator(2, 0, 0).graphs()
    two_tail_star = [
        graph for graph in genus_two_graphs
        if graph.zero_count == 2 and graph.infinity_count == 1
    ][0]
    star_contribution = compiler.compile_graph(genus_two_probe, two_tail_star)
    assert star_contribution.supported
    basic_g2 = EffectiveVertex(
        2, 0, (-2, -2), psi_min=2, insertions=(0, 0)
    )
    assert basic_g2.is_dimension_zero()
    assert tuple(star_contribution.polynomial.terms) == ((basic_g2,),)
    assert star_contribution.polynomial.terms[(basic_g2,)] == -QQ(1)/(2*t**2)

    compilation = compiler.compile_probe(genus_two_probe)
    assert len(compilation.contributions) == 6
    assert not compilation.unsupported
    assert all(item.supported for item in compilation.contributions)
    assert all(
        vertex.is_dimension_zero()
        for contribution in compilation.contributions
        for factors in contribution.polynomial.terms
        for vertex in factors
    )
    report = compilation.report()
    assert report["graph_count"] == 6
    assert report["supported_graphs"] == 6
    assert report["unsupported_graphs"] == 0

    partial_compilation = partial_compiler.compile_probe(genus_two_probe)
    assert len(partial_compilation.unsupported) == 2


run_tests()
print("all CJR graph-contribution compiler tests passed")
