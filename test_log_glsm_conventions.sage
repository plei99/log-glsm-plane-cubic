"""Tests for shared plane-cubic log-GLSM conventions."""

load("log_glsm_conventions.sage")


def run_tests():
    assert specialized_zero_vertex_cache_path(
        "zero.sqlite", "nonequivariant"
    ) == "zero.weights-nonequivariant-quadratic.sqlite"
    assert specialized_zero_vertex_cache_path(
        "zero.weights-nonequivariant-quadratic.sqlite", "nonequivariant"
    ) == "zero.weights-nonequivariant-quadratic.sqlite"

    point = EllipticInsertion("point", 2)
    assert point.total_codimension == 3
    assert point.ambient_class() == (QQ(1) / 3, 1)
    assert EllipticInsertion("unit", 0).ambient_class() == (1, 0)
    assert EllipticInsertion.from_record(point.to_record()) == point

    probe = ProbeSpec.stationary(2, 9, (2,), label="tau2")
    assert probe.intrinsic_degree == 3
    assert probe.virtual_dimension == 3
    assert probe.insertion_codimension == 3
    assert probe.is_dimension_zero()
    assert probe.stationary_descendants() == (2,)
    assert probe.psi_convention == "stabilized"
    assert probe.has_descendants
    assert probe.is_stable_map_type()
    assert ProbeSpec.from_record(probe.to_record()) == probe

    assert not ProbeSpec(1, 0, ()).is_stable_map_type()
    assert ProbeSpec(1, 0, (EllipticInsertion("unit"),)).is_stable_map_type()
    assert ProbeSpec(0, 1, ()).is_stable_map_type()

    log_probe = ProbeSpec(
        2, 9, (EllipticInsertion("point", 2),),
        label="tau2",
        psi_convention="log",
    )
    assert log_probe != probe
    assert probe.with_psi_convention("log") == log_probe
    assert ProbeSpec.from_record(log_probe.to_record()) == log_probe
    try:
        ProbeSpec(1, 0, (), psi_convention="ambiguous")
        raise AssertionError("an unknown psi convention must be rejected")
    except ValueError:
        pass

    assert PlaneCubicDimension.plane_virtual_dimension(2, 4, 3) == 16
    assert PlaneCubicDimension.o3_euler_virtual_rank(2, 4) == 11
    assert PlaneCubicDimension.twisted_zero_virtual_dimension(2, 4, 3) == 5
    assert PlaneCubicDimension.elliptic_virtual_dimension(2, 3) == 5
    assert PlaneCubicDimension.infinity_balance_defect(2, 0, (-3,)) == 0
    assert PlaneCubicDimension.infinity_balance_defect(2, 0, (-2, -2)) == 0
    assert PlaneCubicDimension.infinity_reduced_virtual_dimension(
        3, 0, (-5,)
    ) == 1
    assert PlaneCubicDimension.infinity_required_psi_min_power(
        3, 0, (-5,), 0
    ) == 1
    assert PlaneCubicDimension.infinity_required_psi_min_power(
        3, 0, (-5,), 1
    ) == 0
    assert PlaneCubicDimension.infinity_required_psi_min_power(
        3, 0, (-5,), 2
    ) == -1
    assert PlaneCubicDimension.infinity_dimension_defect(
        3, 0, (-5,), psi_min=0, insertion_codimension=1
    ) == 0
    assert PlaneCubicDimension.psi_min_power_bound(2, 2) == 2
    assert PlaneCubicDimension.psi_min_power_bound(
        3, 1, insertion_codimension=1, ambient_degree=1
    ) == 3

    rings = PlaneCubicCoefficientRing(20)
    l0, l1, l2 = rings.base_weights()
    t = rings.t
    value = 1 / (t - l0)
    series = rings.to_laurent(value, 5)
    assert series[0] == -1 / rings.base_field.gen(0)
    for power in range(5):
        assert rings.laurent_coefficient(value, power) == series[power]

    # Check a genuine pole, including a denominator with nonzero t-valuation.
    pole = (l0 + l1 * t + t**3) / (t**2 * (l2 + t + t**2))
    pole_series = rings.to_laurent(pole, 12)
    for power in range(-3, 7):
        assert rings.laurent_coefficient(pole, power) == pole_series[power]
    assert rings.laurent_coefficient(rings.full(0), -100) == 0

    # CJR localization expands at t=infinity, equivalently at u=1/t near 0.
    at_infinity = (t**2 + 2*t + 3) / (t - 1)
    assert rings.laurent_coefficient_at_infinity(at_infinity, 1) == 1
    assert rings.laurent_coefficient_at_infinity(at_infinity, 0) == 3
    assert rings.laurent_coefficient_at_infinity(at_infinity, -1) == 6
    assert rings.laurent_coefficient_at_infinity(at_infinity, -2) == 6
    assert rings.laurent_coefficient_at_infinity(at_infinity, 2) == 0
    assert rings.laurent_coefficient_at_infinity(
        rings.full(0), -100
    ) == 0
    assert rings.laurent_coefficient_at_infinity(1 / (t - l0), -1) == 1
    assert rings.laurent_coefficient_at_infinity(1 / (t - l0), -2) == l0
    specialized = rings.specialize_base_weights(1 / (t - l1), (0, 1, 3))
    assert specialized == 1 / (t - 1)

    assert plane_cubic_reduced_sign(2, 0) == -1
    assert plane_cubic_reduced_sign(2, 1) == 1


run_tests()
print("all log-GLSM convention tests passed")
