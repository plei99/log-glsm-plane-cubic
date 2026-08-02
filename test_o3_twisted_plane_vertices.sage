"""Tests for exact supported O(3)-twisted plane backends."""

load("o3_twisted_plane_vertices.sage")


def run_tests():
    rings = PlaneCubicCoefficientRing()
    backend = TwistedZeroVertexBackend(rings)
    request_t = TwistedZeroVertexRequest(0, 0, (
        TwistedInsertion(1, 0),
        TwistedInsertion(1, 0),
        TwistedInsertion(0, 0),
    ))
    assert backend.evaluate(request_t) == rings.t

    request_minus_three = TwistedZeroVertexRequest(0, 0, (
        TwistedInsertion(1, 0),
        TwistedInsertion(0, 0),
        TwistedInsertion(0, 0),
    ))
    assert backend.evaluate(request_minus_three) == -3

    request_psi = TwistedZeroVertexRequest(0, 0, (
        TwistedInsertion(1, 1),
        TwistedInsertion(1, 0),
        TwistedInsertion(0, 0),
        TwistedInsertion(0, 0),
    ))
    assert backend.evaluate(request_psi) == rings.t

    unsupported = TwistedZeroVertexRequest(1, 0, (TwistedInsertion(1, 0),))
    try:
        backend.evaluate(unsupported)
        raise AssertionError("unsupported higher-genus vertex must not become zero")
    except UnsupportedGeometryError:
        pass
    backend.register(unsupported, QQ(1)/24, "test value")
    assert backend.evaluate(unsupported) == QQ(1)/24
    assert backend.provenance(unsupported) == "test value"

    i_backend = O3TwistedIBackend()
    assert i_backend.restriction(0, 0) == 1
    degree_one = i_backend.restriction(0, 1)
    l0, l1, l2 = i_backend.lambdas
    assert degree_one == prod(
        i_backend.s + 3*l0 + m*i_backend.z for m in range(1, 4)
    ) / prod(
        l0 - other + i_backend.z for other in (l0, l1, l2)
    )
    assert len(i_backend.section_weights_on_line(0, 1, 2)) == 7

    resummed = ResummedO3TwistedTheory(4)
    assert resummed.primitive_block(1)[1] == 1
    assert resummed.primitive_block(2)[1] == QQ(1)/12
    assert resummed.genus_two_zero_level()[1] == QQ(1)/12


run_tests()
print("all O(3)-twisted plane backend tests passed")
