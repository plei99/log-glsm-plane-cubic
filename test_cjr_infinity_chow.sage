"""Tests for residual Chow-degree tracking in CJR localization."""

load("cjr_infinity_chow.sage")


def run_tests():
    backend = CJRInfinityChowBackend()

    numerical = ProbeSpec.stationary(2, 0, (2,))
    assert numerical.dimension_defect() == 0
    assert backend.degree(numerical, 0).is_zero_cycle
    assert backend.degree(numerical, -1).is_forced_zero
    assert backend.scalar_powers(numerical, (2, 1, 0, -1, -2)) \
        == (0, -1, -2)
    try:
        backend.require_scalar(numerical, 1)
        raise AssertionError("positive powers need a Chow-class pairing")
    except UnsupportedGeometryError:
        pass

    primary = ProbeSpec(
        3, 0, (EllipticInsertion("point", 0),),
        psi_convention="stabilized",
    )
    assert primary.dimension_defect() == 4
    assert backend.numerical_t_power(primary) == -4
    assert backend.degree(primary, -3).residual_dimension == 1
    assert backend.degree(primary, -4).is_zero_cycle
    assert backend.degree(primary, -5).is_forced_zero
    assert backend.scalar_powers(primary, (-1, -2, -3, -4, -5)) \
        == (-4, -5)


run_tests()
print("all CJR infinity-Chow backend tests passed")
