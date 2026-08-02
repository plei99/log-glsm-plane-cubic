"""Tests for shared plane-cubic log-GLSM conventions."""

load("log_glsm_conventions.sage")


def run_tests():
    point = EllipticInsertion("point", 2)
    assert point.total_codimension == 3
    assert point.ambient_class() == (QQ(1) / 3, 1)
    assert EllipticInsertion("unit", 0).ambient_class() == (1, 0)

    probe = ProbeSpec.stationary(2, 9, (2,), label="tau2")
    assert probe.intrinsic_degree == 3
    assert probe.virtual_dimension == 3
    assert probe.insertion_codimension == 3
    assert probe.is_dimension_zero()
    assert probe.stationary_descendants() == (2,)

    assert PlaneCubicDimension.plane_virtual_dimension(2, 4, 3) == 16
    assert PlaneCubicDimension.o3_euler_virtual_rank(2, 4) == 11
    assert PlaneCubicDimension.twisted_zero_virtual_dimension(2, 4, 3) == 5
    assert PlaneCubicDimension.elliptic_virtual_dimension(2, 3) == 5
    assert PlaneCubicDimension.infinity_balance_defect(2, 0, (-3,)) == 0
    assert PlaneCubicDimension.infinity_balance_defect(2, 0, (-2, -2)) == 0
    assert PlaneCubicDimension.psi_min_power_bound(2, 2) == 5

    rings = PlaneCubicCoefficientRing(8)
    l0, l1, l2 = rings.base_weights()
    t = rings.t
    value = 1 / (t - l0)
    series = rings.to_laurent(value, 5)
    assert series[0] == -1 / rings.base_field.gen(0)
    specialized = rings.specialize_base_weights(1 / (t - l1), (0, 1, 3))
    assert specialized == 1 / (t - 1)

    assert plane_cubic_reduced_sign(2, 0) == -1
    assert plane_cubic_reduced_sign(2, 1) == 1


run_tests()
print("all log-GLSM convention tests passed")
