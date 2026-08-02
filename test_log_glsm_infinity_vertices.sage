"""Exact checks for the genus-two log-GLSM infinity-vertex inversion."""

load("log_glsm_infinity_vertices.sage")


def run_tests():
    polynomial, gw = connected_stationary_qseries([2], 6)
    P = polynomial.parent()
    E2, E4, E6 = P.gens()
    assert polynomial == E2^2 / 2 + E4 / 12
    assert tuple(gw[d] for d in range(7)) == (
        QQ(7) / 5760,
        QQ(1) / 24,
        QQ(9) / 8,
        QQ(31) / 6,
        QQ(343) / 24,
        QQ(117) / 4,
        QQ(111) / 2,
    )

    assert infinity_contact_profiles(0, 3) == tuple()
    assert infinity_contact_profiles(1, 1) == ((-1,),)
    assert infinity_contact_profiles(1, 3) == ((-1, -1, -1),)
    assert infinity_contact_profiles(2, 1) == ((-3,),)
    assert infinity_contact_profiles(2, 2) == ((-3, -1), (-2, -2))

    result = reconstruct_genus_two_infinity_vertices(6)
    assert result["known_genus_one_polynomial"] == E2
    assert result["genus_one_localization_residual"][0] == -QQ(1) / 24
    assert result["genus_one_basic_contact_minus_one"] == -QQ(1) / 24
    assert result["genus_two_primitive_profile_combination"] == QQ(1) / 2880
    assert result["two_genus_one_gluing"] == QQ(1) / 1152
    assert result["assembled_degree_zero_infinity"] == QQ(7) / 5760
    assert all(result["checks"].values())

    residual = result["localization_residual"]
    assert tuple(residual[d] for d in range(7)) == (
        QQ(7) / 5760,
        -QQ(1) / 24,
        -QQ(1) / 8,
        -QQ(1) / 6,
        -QQ(7) / 24,
        -QQ(1) / 4,
        -QQ(1) / 2,
    )

    report = infinity_vertex_report(3)
    assert report["known_genus_one_polynomial"] == "E2"
    assert report["infinity_vertices"]["genus_one_basic_contact_minus_one"] == "-1/24"
    assert report["infinity_vertices"]["genus_two_primitive_profile_combination"] == "1/2880"
    assert all(report["checks"].values())


run_tests()
print("all log-GLSM infinity-vertex tests passed")
