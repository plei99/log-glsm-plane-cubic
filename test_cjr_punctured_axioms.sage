r"""Tests for the CJR unit-removal axioms in the plane-cubic basis."""

load("cjr_punctured_axioms.sage")
load("log_glsm_infinity_orchestrator.sage")


def coefficient_map(relation):
    answer = {}
    for coefficient, factors in relation.terms:
        assert len(factors) == 1
        answer[factors[0]] = answer.get(factors[0], 0) + coefficient
    return answer


def residual(relation, values):
    return relation.known_gw - relation.twisted_zero_level - sum(
        coefficient * prod(values[factor] for factor in factors)
        for coefficient, factors in relation.terms
    )


def run_tests():
    axioms = PlaneCubicPuncturedAxioms()

    base_psi = EffectiveVertex(
        2, 0, (-3,), psi_min=1, insertions=(0,)
    )
    base_h = EffectiveVertex(
        2, 0, (-3,), psi_min=0, insertions=(1,)
    )
    string_vertex = EffectiveVertex(
        2, 0, (-3, -1), psi_min=2, insertions=(0, 0)
    )
    string = axioms.string_relation(string_vertex, 1)
    assert coefficient_map(string) == {
        string_vertex: 1,
        base_psi: -3,
        base_h: 9,
    }

    divisor_vertex = EffectiveVertex(
        2, 0, (-3, -1), psi_min=1, insertions=(0, 1)
    )
    divisor = axioms.divisor_relation(divisor_vertex, 1)
    assert coefficient_map(divisor) == {
        divisor_vertex: 1,
        base_h: -3,
    }

    dilaton = axioms.dilaton_relation(string_vertex, 1)
    assert coefficient_map(dilaton) == {
        string_vertex: 1,
        divisor_vertex: 3,
        base_psi: -3,
    }

    # The dilaton row is string + 3*divisor after using plane-cubic balance.
    combined = coefficient_map(string)
    for vertex, coefficient in coefficient_map(divisor).items():
        combined[vertex] = combined.get(vertex, 0) + 3 * coefficient
        if not combined[vertex]:
            del combined[vertex]
    assert combined == coefficient_map(dilaton)

    fundamental = EffectiveVertex(
        2, 0, (-3, -1), psi_min=0, insertions=(2, 0)
    )
    assert coefficient_map(axioms.string_relation(fundamental, 1)) == {
        fundamental: 1,
    }
    assert axioms._raise_insertion(base_h, 0, 2) is None  # H^3=0.

    # Relation metadata uses the same ambient P^2 degree as its effective
    # vertices.  It must not apply the elliptic intrinsic-to-ambient factor a
    # second time in the positive infinity-degree sector.
    positive_degree = EffectiveVertex(
        3, 1, (-2, -1), psi_min=3, insertions=(2, 0)
    )
    positive_degree_string = axioms.string_relation(positive_degree, 1)
    assert positive_degree_string.probe.ambient_degree == 1

    for relation in (string, divisor, dilaton):
        assert all(
            factor.is_balanced() and factor.is_dimension_zero()
            for coefficient, factors in relation.terms
            for factor in factors
        )

    genus_one_h = EffectiveVertex(
        1, 0, (-1,), psi_min=0, insertions=(1,)
    )
    genus_one_psi = EffectiveVertex(
        1, 0, (-1,), psi_min=1, insertions=(0,)
    )
    assert axioms.known_value(genus_one_h) == 0
    assert axioms.known_value(genus_one_psi) == QQ(1) / 8

    contact_descendant = EffectiveVertex(
        2, 0, (-3,), psi_min=0, insertions=(0,), contact_psi=(1,)
    )
    assert axioms.relations_for(contact_descendant) == tuple()
    assert axioms.known_value(contact_descendant) is None

    # Repeated unit removal computes a higher-valence genus-one invariant
    # before any compact probe is compiled.
    genus_one_two_units = EffectiveVertex(
        1, 0, (-1, -1), psi_min=2, insertions=(0, 0)
    )
    config = InfinityOrchestrationConfig(
        max_markings=1,
        t_powers=(0,),
        include_unit_relatives=False,
        max_genus=1,
        max_ambient_degree=0,
    )
    orchestrator = InfinityVertexOrchestrator(
        object(),
        (genus_one_two_units,),
        config=config,
        coefficient_field=QQ,
        candidate_provider=lambda *args: tuple(),
    )
    report = orchestrator.run()
    assert report["status"] == "complete"
    assert report["stage"] == -1
    assert orchestrator.solver.values[genus_one_two_units] == QQ(1) / 8
    assert all(
        not residual(orchestrator.relations[key], orchestrator.solver.values)
        for key in orchestrator.active_relation_keys
    )

    # The full CJR localization provider uses the t=infinity expansion and the
    # non-equivariant auxiliary-base limit.  Its genus-one row must therefore
    # agree with the closed punctured value, and every negative genus-two row
    # must vanish after the genus-one CJR reductions.
    full_provider = PlaneCubicFullEquationProvider(laurent_precision=8)
    genus_one_probe = ProbeSpec.stationary(
        1, 0, (0,), label="stationary candidate"
    )
    genus_one_row = full_provider.relation(genus_one_probe, 0)
    assert residual(genus_one_row, {genus_one_h: 0}) == 0
    assert full_provider.relation(genus_one_probe, -1).twisted_zero_level == 0

    genus_two_probe = ProbeSpec(
        2, 0, (EllipticInsertion("point", 2),),
        label="log-domain descendant diagnostic",
        psi_convention="log",
    )
    negative_rows = tuple(
        full_provider.relation(genus_two_probe, power)
        for power in (-1, -2, -3, -4)
    )
    genus_one_factors = tuple(sorted(set(
        factor
        for row in negative_rows
        for coefficient, factors in row.terms
        for factor in factors
    ), key=lambda item: item.order_key()))
    if genus_one_factors:
        lower = InfinityVertexOrchestrator(
            object(),
            genus_one_factors,
            config=InfinityOrchestrationConfig(
                max_markings=1,
                t_powers=(0,),
                include_unit_relatives=False,
                max_genus=1,
                max_ambient_degree=0,
            ),
            coefficient_field=QQ,
            candidate_provider=lambda *args: tuple(),
        )
        assert lower.run()["status"] == "complete"
        lower_values = lower.solver.values
    else:
        lower_values = {}
    assert all(not residual(row, lower_values) for row in negative_rows)

    # Contradictory rows must not disappear after all variables have values.
    bad = CompiledLocalizationRelation(
        ProbeSpec.stationary(1, 0, (0,), label="contradiction"),
        0, QQ, 1, terms=((1, (genus_one_h,)),),
    )
    broken = InfinityVertexOrchestrator(
        object(),
        (genus_one_h,),
        initial_values={genus_one_h: 0},
        config=InfinityOrchestrationConfig(
            max_markings=1,
            t_powers=(0,),
            include_unit_relatives=False,
            include_punctured_axioms=False,
            max_genus=1,
            max_ambient_degree=0,
        ),
        coefficient_field=QQ,
        candidate_provider=lambda *args: (bad,),
    )
    broken.relations[broken._relation_key(bad)] = bad
    broken._activate_incident_relations()
    try:
        broken.solve_next_block()
        raise AssertionError("a nonzero solved residual must be rejected")
    except InconsistentLocalizationRelationError:
        pass


def test_genus_one_closed_form_matches_cjr_ii():
    r"""CJR II (9.17) and (9.21) must reproduce the hard-coded base values.

    This pins three conventions at once: the sign of ``psi_min``, the ``3H^2``
    normalization of the reduced cycle, and the ambient-insertion basis.  It
    is a genuine two-sided check -- the divisor value discriminates classes
    that the ``psi_min`` value alone cannot.
    """
    axioms = PlaneCubicPuncturedAxioms()
    closed = axioms.genus_one_closed_form()

    ring = closed["reduced_cycle"].parent()
    lam, hyperplane = ring.gens()
    assert closed["reduced_cycle"] == 3 * hyperplane**2
    assert closed["psi_min"] == lam - 3 * hyperplane

    psi_vertex = EffectiveVertex(1, 0, (-1,), psi_min=1, insertions=(0,))
    divisor_vertex = EffectiveVertex(1, 0, (-1,), psi_min=0, insertions=(1,))
    assert psi_vertex.is_dimension_zero()
    assert divisor_vertex.is_dimension_zero()

    assert closed["psi_min_value"] == QQ(1) / 8
    assert closed["divisor_value"] == QQ(0)
    # The derivation and the tabulated CJR values must agree.
    assert axioms.known_value(psi_vertex) == closed["psi_min_value"]
    assert axioms.known_value(divisor_vertex) == closed["divisor_value"]


def test_cycle_unit_axiom_collapses_to_the_string_equation():
    r"""CJR II (8.8) yields no rows beyond (8.3) on dimension-zero data.

    Recording this as a test stops the cycle-valued unit axiom from being
    re-proposed as a source of extra independent relations.
    """
    axioms = PlaneCubicPuncturedAxioms()
    # red dim = n + 3D = 2, so psi_min = 2 makes the cycle a point class.
    vertex = EffectiveVertex(
        2, 0, (-3, -1), psi_min=2, insertions=(0, 0)
    )
    assert vertex.is_balanced()
    assert vertex.is_dimension_zero()
    unit_index = axioms._unit_positions(vertex, 0)[0]

    assert axioms.cycle_unit_axiom_test_classes(vertex, unit_index) == (1,)
    collapsed = axioms.cycle_unit_axiom_relation(vertex, unit_index)
    string = axioms.string_relation(vertex, unit_index)
    assert axioms._relation_signature(collapsed) == \
        axioms._relation_signature(string)

    # On a positive-dimensional effective cycle the axiom is not numerical,
    # and the module says so instead of inventing a row.
    higher = EffectiveVertex(
        2, 0, (-3, -1), psi_min=1, insertions=(0, 0)
    )
    assert not higher.is_dimension_zero()
    try:
        axioms.cycle_unit_axiom_test_classes(
            higher, axioms._unit_positions(higher, 0)[0]
        )
        raise AssertionError("a positive-dimensional cycle must be rejected")
    except UnsupportedGeometryError:
        pass


run_tests()
test_genus_one_closed_form_matches_cjr_ii()
test_cycle_unit_axiom_collapses_to_the_string_equation()
print("all CJR punctured-axiom tests passed")
