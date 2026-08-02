"""Tests for exact psi and lambda-g intersection backends."""

load("hodge_integrals.sage")


def run_tests():
    psi = PsiIntersectionBackend()
    assert psi.integral(0, (0, 0, 0)) == 1
    assert psi.integral(0, (1, 0, 0, 0)) == 1
    assert psi.integral(0, (2, 0, 0, 0, 0)) == 1
    assert psi.integral(0, (1, 1, 0, 0, 0)) == 2
    assert psi.integral(1, (1,)) == QQ(1) / 24
    assert psi.integral(1, (2, 0)) == QQ(1) / 24
    assert psi.integral(2, (4,)) == QQ(1) / 1152
    assert psi.integral(2, (3, 2)) == QQ(29) / 5760
    assert psi.integral(2, (1,)) == 0

    hodge = HodgeIntegralBackend()
    assert hodge.lambda_g_constant(1) == QQ(1) / 24
    assert hodge.lambda_g_constant(2) == QQ(7) / 5760
    assert hodge.integral(1, (0,), (1,)) == QQ(1) / 24
    assert hodge.integral(2, (2,), (2,)) == QQ(7) / 5760
    assert hodge.integral(2, (2, 1), (2,)) == QQ(7) / 1920
    assert hodge.integral(2, (1, 1), (2,)) == 0

    try:
        hodge.integral(2, (3,), (1,))
        raise AssertionError("unsupported lambda_1 integral should be explicit")
    except UnsupportedGeometryError:
        pass

    adm = AdmcyclesHodgeIntegralBackend()
    assert adm.integral(2, (3,), (1,)) == QQ(1) / 480
    assert adm.integral(2, (2,), (2,)) == QQ(7) / 5760
    assert adm.integral(2, (2, 1, 1), (1, 1)) == QQ(7) / 240


run_tests()
print("all Hodge-integral backend tests passed")
