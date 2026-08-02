"""Tests for known elliptic ProbeSpec values."""

load("elliptic_probe_values.sage")


def run_tests():
    backend = EllipticProbeValueBackend()
    expected = {
        0: QQ(7)/5760,
        3: QQ(1)/24,
        6: QQ(9)/8,
    }
    for ambient_degree, value in expected.items():
        probe = ProbeSpec.stationary(2, ambient_degree, (2,))
        assert backend.value(probe) == value

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

    assert backend.value(ProbeSpec.stationary(2, 0, (1,))) == 0


run_tests()
print("all elliptic ProbeSpec value tests passed")
