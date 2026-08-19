"""Tests for exact psi and lambda-g intersection backends."""

load("hodge_integrals.sage")

import os
import tempfile


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
    assert hodge.top_hodge_triple_constant(2) == QQ(1) / 5760
    assert hodge.top_hodge_triple_constant(3) == QQ(1) / 1451520
    assert hodge.integral(1, (0,), (1,)) == QQ(1) / 24
    assert hodge.integral(2, (2,), (2,)) == QQ(7) / 5760
    assert hodge.integral(2, (2, 1), (2,)) == QQ(7) / 1920
    assert hodge.integral(2, (1, 1), (2,)) == 0
    # Mumford's relation removes repeated low-genus Hodge factors before an
    # admcycles class is ever constructed.
    repeated = HodgeIntegralRequest(2, (2, 1, 1), (1, 1))
    factor, normalized = repeated.low_genus_normal_form()
    assert factor == 2
    assert normalized.lambda_indices == (2,)
    assert hodge.integral(repeated) == QQ(7) / 240
    assert HodgeIntegralRequest(1, (0,), (1, 1)).low_genus_normal_form()[0] == 0
    genus_three = HodgeIntegralRequest(3, (), (1, 1, 2, 2))
    factor, normalized = genus_three.low_genus_normal_form()
    assert factor == 4
    assert normalized.lambda_indices == (1, 2, 3)
    assert HodgeIntegralRequest(3, (), (3, 3)).low_genus_normal_form()[0] == 0
    assert hodge.integral(3, (), (1, 2, 3)) == QQ(1) / 1451520
    assert hodge.integral(3, (1,), (1, 2, 3)) == QQ(1) / 362880
    assert hodge.integral(3, (2, 0), (1, 2, 3)) == QQ(1) / 362880

    try:
        hodge.integral(2, (3,), (1,))
        raise AssertionError("unsupported lambda_1 integral should be explicit")
    except UnsupportedGeometryError:
        pass

    adm = AdmcyclesHodgeIntegralBackend()
    # Do not forget into the unstable (g,n)=(1,0) space.
    assert adm.integral(1, (0,), (1,)) == QQ(1) / 24
    assert adm.integral(2, (3,), (1,)) == QQ(1) / 480
    # String/dilaton reduction applies to every Hodge monomial before the
    # general admcycles path.  These identities also guard the genus-three
    # performance fix used by the infinity-vertex computation.
    assert adm.integral(2, (4, 0), (1,)) == QQ(1) / 480
    assert adm.integral(2, (3, 1), (1,)) == 3 * QQ(1) / 480
    genus_three_small = adm.integral(3, (2,), (2, 3))
    assert adm.integral(3, (4, 0, 0), (2, 3)) == genus_three_small
    assert adm.integral(2, (2,), (2,)) == QQ(7) / 5760
    assert adm.integral(2, (2, 1, 1), (1, 1)) == QQ(7) / 240
    # Marking permutations represent the same integral and share the cache.
    permuted = HodgeIntegralRequest(2, (1, 2, 1), (1, 1))
    canonical = HodgeIntegralRequest(2, (2, 1, 1), (1, 1))
    assert permuted == canonical
    cache_size = len(adm._cache)
    assert adm.integral(permuted) == QQ(7) / 240
    assert len(adm._cache) == cache_size

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "hodge.sqlite")
        persistent = AdmcyclesHodgeIntegralBackend(cache_path=path)
        request = HodgeIntegralRequest(2, (3,), (1,))
        assert persistent.integral(request) == QQ(1) / 480
        assert persistent.cache_info()["persistent_writes"] == 1
        diagnostic = HodgeIntegralRequest(3, (4,), (2,))
        diagnostic_key = persistent.store.begin_evaluation(diagnostic)
        diagnostic_value = persistent.store.connection.execute(
            "SELECT value FROM metadata WHERE key = ?", (diagnostic_key,)
        ).fetchone()[0]
        assert json.loads(diagnostic_value) == diagnostic.to_record()
        persistent.store.finish_evaluation(diagnostic_key)
        assert persistent.store.connection.execute(
            "SELECT value FROM metadata WHERE key = ?", (diagnostic_key,)
        ).fetchone() is None
        persistent.store.close()

        restored = AdmcyclesHodgeIntegralBackend(cache_path=path)
        assert restored.integral(request) == QQ(1) / 480
        assert restored.cache_info()["persistent_hits"] == 1
        assert restored.cache_info()["persistent_entries"] == 1
        restored.store.close()


run_tests()
print("all Hodge-integral backend tests passed")
