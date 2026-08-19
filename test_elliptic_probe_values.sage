"""Tests for known elliptic ProbeSpec values."""

load("elliptic_probe_values.sage")


def run_tests():
    backend = EllipticProbeValueBackend()
    try:
        backend.value(ProbeSpec(1, 0, ()))
        raise AssertionError("an unmarked constant genus-one probe is unstable")
    except UnsupportedGeometryError:
        pass
    expected = {
        0: QQ(7)/5760,
        3: QQ(1)/24,
        6: QQ(9)/8,
    }
    for ambient_degree, value in expected.items():
        probe = ProbeSpec.stationary(2, ambient_degree, (2,))
        assert backend.value(probe) == value
        assert backend.reduced_stabilized_value(probe) \
            == plane_cubic_reduced_sign(2, ambient_degree) * value

    genus_one = ProbeSpec.stationary(1, 0, (0,))
    assert backend.value(genus_one) == -QQ(1)/24
    assert backend.reduced_log_glsm_value(genus_one) == -QQ(1)/24

    # String: <tau_0(1) tau_1(pt)>_(1,0) = <tau_0(pt)>_(1,0).
    string_probe = ProbeSpec(1, 0, (
        EllipticInsertion("unit", 0),
        EllipticInsertion("point", 1),
    ))
    assert string_probe.is_dimension_zero()
    assert backend.value(string_probe) == -QQ(1)/24

    # Dilaton: <tau_1(1) tau_0(pt)>_1 = 1*<tau_0(pt)>_1.
    dilaton_probe = ProbeSpec(1, 0, (
        EllipticInsertion("unit", 1),
        EllipticInsertion("point", 0),
    ))
    assert backend.value(dilaton_probe) == -QQ(1)/24

    # Constant maps pull every point insertion back from the same elliptic
    # curve, hence any degree-zero stationary invariant with at least two
    # point classes vanishes without a Bloch--Okounkov expansion.
    many_points = ProbeSpec.stationary(3, 0, (1, 1, 1, 1))
    assert many_points.is_dimension_zero()
    assert backend.value(many_points) == 0

    assert backend.value(ProbeSpec.stationary(2, 0, (1,))) == 0

    for ambient_degree in (1, 2, 4, 5):
        empty_probe = ProbeSpec.stationary(2, ambient_degree, (2,))
        assert backend.value(empty_probe) == 0
        assert backend.reduced_stabilized_value(empty_probe) == 0
        assert backend.reduced_log_glsm_value(empty_probe) == 0

    # CJR's graph formula inserts the cotangent line before stabilization.
    # Bloch--Okounkov supplies the stable-map cotangent line after
    # stabilization, so the former is not a valid known side for descendants.
    log_descendant = ProbeSpec(
        2, 0, (EllipticInsertion("point", 2),),
        psi_convention="log",
    )
    try:
        backend.value(log_descendant)
        raise AssertionError("log-domain descendants need a separate backend")
    except UnsupportedGeometryError:
        pass
    stabilized = ProbeSpec.stationary(2, 0, (2,))
    assert backend.reduced_log_glsm_value(stabilized) == -QQ(7) / 5760


run_tests()
print("all elliptic ProbeSpec value tests passed")
