"""Tests for the CJR-I stabilization-boundary comparison."""

load("cjr_stabilization_boundary.sage")
load("cjr_full_equation_provider.sage")


def run_tests():
    backend = EllipticProbeValueBackend()
    comparison = StabilizationBoundaryComparison(backend)

    stabilized = ProbeSpec.stationary(
        3, 0, (4,), label="genus-three one-point descendant"
    )
    log_probe = comparison.log_probe(stabilized)
    assert log_probe.psi_convention == "log"
    assert comparison.stabilized_probe(log_probe) == stabilized

    # The two classes remain convention-separated.  Only the virtual
    # comparison adapter is allowed to transfer the known number.
    try:
        backend.reduced_log_glsm_value(log_probe)
        raise AssertionError("the raw known backend must reject log psi")
    except UnsupportedGeometryError:
        pass
    assert comparison.reduced_value(log_probe) == -QQ(31) / 967680

    provider = PlaneCubicFullEquationProvider(
        laurent_precision=8,
        zero_vertex_cache="tmp/test-stabilization-boundary-zero.sqlite",
    )
    relation = provider.boundary_compared_relation(stabilized, 0)
    assert relation.probe == log_probe
    target = EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(1,)
    )
    assert relation.known_gw == -QQ(31) / 967680
    assert any(
        factors == (target,) and coefficient == -QQ(625) / 72
        for coefficient, factors in relation.terms
    )
    assert all(
        not factor.contact_psi
        for coefficient, factors in relation.terms
        for factor in factors
    )

    report = provider.relation_report(log_probe)
    comparison_report = report["stabilization_boundary_comparison"]
    assert comparison_report["log_probe"]["psi_convention"] == "log"
    assert comparison_report["stabilized_probe"]["psi_convention"] \
        == "stabilized"
    assert comparison_report["reduced_value"] == "-31/967680"

    # Effective-basis generation now includes the boundary-compared
    # descendant row and never emits an ordinary contact cotangent class.
    rows = provider.candidate_relations(
        3, 0, max_markings=1, t_powers=(0,),
        include_chow_relations=True,
        max_chow_unit_insertions=0,
        chow_primary_only=True,
        effective_basis_only=True,
    )
    assert len(rows) == 1
    assert rows[0].probe.psi_convention == "log"
    assert rows[0].probe.stationary_descendants() == (4,)
    assert all(
        not factor.contact_psi
        for row in rows for coefficient, factors in row.terms
        for factor in factors
    )


run_tests()
print("all stabilization-boundary comparison tests passed")
