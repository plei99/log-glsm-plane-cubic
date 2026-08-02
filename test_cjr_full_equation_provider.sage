"""Tests for rank-aware probes and the full equation-provider shell."""

load("cjr_full_equation_provider.sage")


def run_tests():
    x = EffectiveVertex(2, 0, (-3,), label="x")
    y = EffectiveVertex(2, 0, (-2, -2), label="y")
    probe = ProbeSpec.stationary(2, 0, (2,))
    relation1 = CompiledLocalizationRelation(
        probe, 0, QQ, 5, terms=((2, (x,)), (1, (y,)))
    )
    relation2 = CompiledLocalizationRelation(
        probe, -1, QQ, 0, terms=((1, (x,)), (3, (y,)))
    )
    dependent = CompiledLocalizationRelation(
        probe, -2, QQ, 0, terms=((4, (x,)), (2, (y,)))
    )
    factory = ProbeFactory(QQ)
    selection = factory.select_relations(
        (x, y), (relation1, dependent, relation2)
    )
    assert selection.is_full_rank
    assert selection.rank == 2
    assert len(selection.relations) == 2
    assert selection.kernel_basis() == tuple()

    candidates = factory.stationary_candidates(2, 0, 3)
    assert all(candidate.is_dimension_zero() for candidate in candidates)
    assert {candidate.marking_count for candidate in candidates} == {1, 2, 3}

    provider = PlaneCubicFullEquationProvider(laurent_precision=12)
    genus_one = ProbeSpec.stationary(1, 0, (0,))
    complete = provider.relation(genus_one, 0, require_complete=False)
    assert complete.is_complete
    report = provider.relation_report(genus_one)
    assert report["complete"]
    assert report["compilation"]["graph_count"] == 2
    assert provider.relation(genus_one, 0, require_complete=True).is_complete

    partial_rings = PlaneCubicCoefficientRing(12)
    partial_provider = PlaneCubicFullEquationProvider(
        rings=partial_rings,
        twisted_backend=TwistedZeroVertexBackend(partial_rings),
    )
    incomplete = partial_provider.relation(
        genus_one, 0, require_complete=False
    )
    assert not incomplete.is_complete
    try:
        partial_provider.relation(genus_one, 0, require_complete=True)
        raise AssertionError("an injected partial backend must still block")
    except UnsupportedGeometryError:
        pass

    result = genus_two_resummed_end_to_end(6)
    assert result["b2"] == -QQ(1)/24
    assert result["b4_aggregate"] == QQ(1)/2880
    assert result["coefficients"][:3] == (
        QQ(7)/5760, QQ(1)/24, QQ(9)/8
    )
    assert result["known_series"] == result["reconstructed_series"]
    report = genus_two_end_to_end_report(2)
    assert report["values"][2]["ambient_degree"] == 6
    assert report["values"][2]["value"] == "9/8"
    assert report["reconstruction_check"]

    # Generic-field DP support.
    polynomial = PolynomialRing(QQ, "u")
    field = polynomial.fraction_field()
    u = field.gen()

    def field_provider(vertex):
        if vertex == x:
            return LocalizationEquation(
                x, u, terms=((u, (x,)),), coefficient_field=field
            )
        raise KeyError(vertex)

    field_solver = InfinityVertexDP(field_provider, coefficient_field=field)
    assert field_solver.solve(x) == 1


run_tests()
print("all full CJR equation-provider tests passed")
