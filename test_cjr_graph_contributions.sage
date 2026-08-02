"""Tests for exact per-graph CJR contribution compilation."""

load("cjr_graph_contributions.sage")


def run_tests():
    rings = PlaneCubicCoefficientRing(12)
    compiler = PlaneCubicGraphContributionCompiler(rings)

    genus_one_probe = ProbeSpec.stationary(1, 0, (0,))
    genus_one_graphs = PlaneCubicGraphEnumerator(1, 1, 0).graphs()
    mixed = [graph for graph in genus_one_graphs if graph.infinity_count][0]
    contribution = compiler.compile_graph(genus_one_probe, mixed)
    assert contribution.supported
    basic = EffectiveVertex(1, 0, (-1,), insertions=(1,))
    assert contribution.polynomial.terms == {(basic,): rings.full(-QQ(1)/3)}

    pure_zero = [graph for graph in genus_one_graphs if not graph.infinity_count][0]
    pure_zero_contribution = compiler.compile_graph(genus_one_probe, pure_zero)
    assert pure_zero_contribution.supported
    assert tuple() in pure_zero_contribution.polynomial.terms
    assert pure_zero_contribution.polynomial.terms[tuple()] != 0

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
    basic_g2 = EffectiveVertex(2, 0, (-2, -2), insertions=(0, 0))
    assert star_contribution.polynomial.terms[(basic_g2,)] == -QQ(1)/2

    compilation = compiler.compile_probe(genus_two_probe)
    assert len(compilation.contributions) == 6
    assert not compilation.unsupported
    assert all(item.supported for item in compilation.contributions)
    report = compilation.report()
    assert report["graph_count"] == 6
    assert report["supported_graphs"] == 6
    assert report["unsupported_graphs"] == 0

    partial_compilation = partial_compiler.compile_probe(genus_two_probe)
    assert len(partial_compilation.unsupported) == 2


run_tests()
print("all CJR graph-contribution compiler tests passed")
