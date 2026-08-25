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

    basis = ExactIncrementalRowBasis(QQ, 3)
    assert basis.add((1, 2, 3))
    assert not basis.add((2, 4, 6))
    assert basis.add((0, 1, 1))
    assert basis.rank == 2

    hybrid_provider = PlaneCubicFullEquationProvider(
        zero_vertex_strategy="hybrid"
    )
    assert hybrid_provider.zero_vertex_strategy == "hybrid"
    assert isinstance(
        hybrid_provider.twisted_backend, HybridTwistedZeroVertexBackend
    )
    assert isinstance(
        hybrid_provider.twisted_backend.givental,
        O3CalibratedGiventalBackend,
    )
    routed_positive = TwistedZeroVertexRequest(
        1, 1, (TwistedInsertion(1, 1),)
    )
    assert hybrid_provider.twisted_backend.evaluate(routed_positive) \
        == FullTwistedZeroVertexBackend().evaluate(routed_positive)

    # An unresolved term not containing a block target still contributes to
    # the right-hand side and must keep the row from being selected.
    unresolved = EffectiveVertex(1, 0, (-1,), label="unresolved")
    unsafe = CompiledLocalizationRelation(
        probe, 0, QQ, 1,
        terms=((1, (x,)), (1, (unresolved,))),
    )
    unsafe_selection = factory.select_relations((x,), (unsafe,))
    assert not unsafe_selection.is_full_rank
    assert "unresolved lower factors" in unsafe_selection.rejected[0][1]

    candidates = factory.stationary_candidates(2, 0, 3)
    assert all(candidate.is_dimension_zero() for candidate in candidates)
    assert {candidate.marking_count for candidate in candidates} == {1, 2, 3}

    one_point = ProbeSpec.stationary(2, 0, (2,))
    string_relative, dilaton_relative = factory.unit_relatives(one_point)
    assert string_relative.is_dimension_zero()
    assert dilaton_relative.is_dimension_zero()
    assert tuple((item.kind, item.psi_power)
                 for item in string_relative.insertions) == (
        ("unit", 0), ("point", 3)
    )
    assert tuple((item.kind, item.psi_power)
                 for item in dilaton_relative.insertions) == (
        ("unit", 1), ("point", 2)
    )
    known_backend = EllipticProbeValueBackend()
    assert known_backend.value(string_relative) == QQ(7) / 5760
    assert known_backend.value(dilaton_relative) == QQ(7) / 1920

    mixed = factory.mixed_unit_candidates(
        3, 0, max_point_markings=1, max_unit_insertions=2
    )
    assert len(mixed) == 5
    assert all(candidate.is_dimension_zero() for candidate in mixed)
    assert sum(
        sum(item.kind == "unit" for item in candidate.insertions) == 2
        for candidate in mixed
    ) == 3
    assert all(
        item.psi_power <= 1
        for candidate in mixed for item in candidate.insertions
        if item.kind == "unit"
    )
    assert all(known_backend.value(candidate) in QQ for candidate in mixed)

    chow_probes = factory.chow_candidates(
        3, 0, max_markings=1, max_unit_insertions=1,
        primary_only=True,
    )
    assert len(chow_probes) == 3
    assert all(probe.dimension_defect() > 0 for probe in chow_probes)
    assert any(not probe.insertions for probe in chow_probes)

    provider = PlaneCubicFullEquationProvider(laurent_precision=12)
    genus_one = ProbeSpec.stationary(1, 0, (0,))
    complete = provider.relation(genus_one, 0, require_complete=False)
    assert complete.is_complete
    report = provider.relation_report(genus_one)
    assert report["complete"]
    assert report["compilation"]["graph_count"] == 2
    assert provider.relation(genus_one, 0, require_complete=True).is_complete

    stabilized_relation = provider.relation(
        ProbeSpec.stationary(2, 0, (2,)), 0
    )
    assert stabilized_relation.is_complete
    assert stabilized_relation.known_gw == -QQ(7) / 5760
    assert len(stabilized_relation.terms) == 1
    stabilized_factor = stabilized_relation.terms[0][1][0]
    assert stabilized_factor.contact_psi == (0, 0, 2)

    # CJR I's virtual stabilization-boundary comparison supplies the compact
    # value for the local log-domain descendant without identifying the two
    # cotangent-line classes.
    log_descendant = ProbeSpec(
        2, 0, (EllipticInsertion("point", 2),),
        psi_convention="log",
    )
    assert provider.relation(log_descendant, -1).is_complete
    log_relation = provider.relation(log_descendant, 0)
    assert log_relation.is_complete
    assert log_relation.known_gw == -QQ(7) / 5760
    assert all(
        not factor.contact_psi
        for coefficient, factors in log_relation.terms
        for factor in factors
    )
    safe_candidates = provider.candidate_relations(
        2, 0, max_markings=1, t_powers=(0, -1)
    )
    assert len(safe_candidates) == 2
    assert {relation.t_power for relation in safe_candidates} == {0, -1}
    assert all(
        relation.probe.psi_convention == "stabilized"
        for relation in safe_candidates
    )

    primary_class_probe = ProbeSpec(
        3, 0, (EllipticInsertion("point", 0),),
        label="primary Chow-class test",
        psi_convention="stabilized",
    )
    # The first numerical coefficient of this defect-four Chow probe is
    # t^-4.  Higher powers remain positive-dimensional and are rejected.
    primary_class_relation = provider.relation(primary_class_probe, -4)
    primary_target = EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(1,)
    )
    assert any(
        factors == (primary_target,) and coefficient == -QQ(625) / 72
        for coefficient, factors in primary_class_relation.terms
    )
    assert primary_class_relation.known_gw == 0
    for power in (-3, 1):
        try:
            provider.relation(primary_class_probe, power)
            raise AssertionError(
                "positive-dimensional Chow coefficients need a test class"
            )
        except UnsupportedGeometryError:
            pass
    try:
        provider.relation(genus_one, 1)
        raise AssertionError("positive t powers are not scalar rows")
    except UnsupportedGeometryError:
        pass

    primary_chow_rows = provider.candidate_relations(
        3, 0, max_markings=1,
        t_powers=(-1, -2, -3, -4),
        include_chow_relations=True,
        max_chow_unit_insertions=1,
        chow_primary_only=True,
    )
    assert all(
        relation.probe.dimension_defect() + relation.t_power <= 0
        for relation in primary_chow_rows
    )
    assert any(
        relation.probe.insertions == primary_class_probe.insertions
        and relation.probe.dimension_defect() == 4
        and relation.t_power == -4
        for relation in primary_chow_rows
    )

    effective_rows = provider.candidate_relations(
        3, 0, max_markings=1, t_powers=(-4,),
        include_chow_relations=True,
        max_chow_unit_insertions=0,
        chow_primary_only=True,
        effective_basis_only=True,
    )
    assert len(effective_rows) == 3  # descendant, unmarked, primary point
    assert sum(relation.probe.has_descendants
               for relation in effective_rows) == 1
    assert all(
        relation.probe.psi_convention == "log"
        for relation in effective_rows if relation.probe.has_descendants
    )
    assert all(
        not factor.contact_psi
        for relation in effective_rows
        for _, factors in relation.terms
        for factor in factors
    )

    # Regression for the apparent genus-one contradiction created by
    # scalarizing positive-dimensional coefficients.  The defect-two
    # two-unit probe becomes numerical only at t^-2, where its complete row
    # is V(psi*1) - 3 V(H) - V(psi^2*1,1) = 0.
    genus_one_two_units = ProbeSpec(
        1, 0,
        (EllipticInsertion("unit", 0), EllipticInsertion("unit", 0)),
        psi_convention="stabilized",
    )
    genus_one_chow = provider.relation(genus_one_two_units, -2)
    expected_terms = {
        EffectiveVertex(
            1, 0, (-1,), psi_min=1, insertions=(0,)
        ): QQ.one(),
        EffectiveVertex(
            1, 0, (-1,), psi_min=0, insertions=(1,)
        ): -QQ(3),
        EffectiveVertex(
            1, 0, (-1, -1), psi_min=2, insertions=(0, 0)
        ): -QQ.one(),
    }
    assert {
        factors[0]: coefficient
        for coefficient, factors in genus_one_chow.terms
        if len(factors) == 1
    } == expected_terms
    assert genus_one_chow.known_gw == 0
    genus_one_values = {
        tuple(expected_terms)[0]: QQ(1) / 8,
        tuple(expected_terms)[1]: QQ.zero(),
        tuple(expected_terms)[2]: QQ(1) / 8,
    }
    assert sum(
        coefficient * genus_one_values[factors[0]]
        for coefficient, factors in genus_one_chow.terms
    ) == 0

    # Equivariant Euler powers do not lower the largest possible flag-psi
    # power.  The compiler must truncate by the ambient stable-map dimension,
    # or it misses a term already in the genus-one ambient-degree-one row.
    degree_one = provider.relation(
        ProbeSpec.stationary(1, 1, (0,)), 0
    )
    degree_one_terms = {
        factors[0]: coefficient
        for coefficient, factors in degree_one.terms
        if len(factors) == 1
    }
    genus_one_psi = EffectiveVertex(
        1, 0, (-1,), psi_min=1, insertions=(0,)
    )
    assert degree_one.twisted_zero_level == QQ(5) / 24
    assert degree_one_terms[genus_one_psi] == -QQ(5) / 3
    # CJR, Punctured logarithmic R-maps, (9.17) and (9.21), gives
    # V(psi_min) = 1/8; all other terms in this row vanish by the
    # genus-one string/divisor axioms.
    assert degree_one.twisted_zero_level \
        + degree_one_terms[genus_one_psi] * QQ(1) / 8 == 0

    try:
        genus_two_profile_split_rank_witness(provider)
        raise AssertionError("the obsolete mixed-psi witness must be rejected")
    except UnsupportedGeometryError:
        pass

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


def test_chow_family_covers_the_unmarked_probe():
    r"""The zero-marking probe must stay in the graded Chow family.

    ``stationary_candidates`` deliberately starts at one marking, because a
    zero-marking probe has no dimension-zero stationary value.  Its scalar
    rows live at ``t^(-(2g-2))`` and are supplied by ``chow_candidates``;
    this test pins that coverage so a refactor cannot silently drop the
    unmarked rows again.
    """
    factory = ProbeFactory(QQ)
    chow = CJRInfinityChowBackend()
    for genus in (2, 3):
        candidates = factory.chow_candidates(
            genus, 0, max_markings=1, primary_only=True
        )
        unmarked = [probe for probe in candidates
                    if probe.marking_count == 0]
        assert len(unmarked) == 1, candidates
        probe = unmarked[0]
        defect = probe.dimension_defect()
        assert defect == 2 * genus - 2
        powers = chow.scalar_powers(
            probe, tuple(range(0, -defect - 2, -1))
        )
        assert powers == tuple(range(-defect, -defect - 2, -1)), powers
    # Genus one has no stable unmarked probe: defect zero means the t^0
    # coefficient is the honest stationary value, handled elsewhere.
    genus_one = factory.chow_candidates(1, 0, max_markings=1,
                                        primary_only=True)
    assert all(probe.marking_count >= 1 or probe.dimension_defect() > 0
               for probe in genus_one)


run_tests()
test_chow_family_covers_the_unmarked_probe()
print("all full CJR equation-provider tests passed")
