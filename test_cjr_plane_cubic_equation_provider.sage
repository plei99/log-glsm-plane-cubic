"""Tests for the concrete one-point CJR equation provider."""

load("cjr_plane_cubic_equation_provider.sage")


def run_tests():
    basic_g2 = basic_minus_two_vertex(2, 0)
    basic_g3_degree_zero = basic_minus_two_vertex(3, 0)
    basic_g3_degree_one = basic_minus_two_vertex(3, 1)
    assert basic_g2.contacts == (-2, -2)
    assert basic_g3_degree_zero.contacts == (-2, -2, -2, -2)
    assert basic_g3_degree_one.contacts == (-2,)
    assert basic_star_diagonal_coefficient(basic_g2) == QQ(1) / 2
    assert basic_star_diagonal_coefficient(basic_g3_degree_zero) == QQ(1) / 24
    assert basic_star_diagonal_coefficient(basic_g3_degree_one) == 1
    w = SR.var("w")
    assert unstable_leaf_edge_pair_factor(1, w) == w
    assert unstable_leaf_edge_pair_factor(2, w) == 1
    assert unstable_leaf_edge_pair_factor(3, w) == QQ(3) / (2 * w)

    provider = PlaneCubicOnePointEquationProvider(
        max_genus=4, verification_degree=6
    )
    graph_enumerator = provider.graph_enumerator(1, 0)
    assert graph_enumerator.graph_count() == 2
    assert tuple(graph.canonical_key()
                 for graph in provider.localization_graphs(1, 0)) == tuple(
        graph.canonical_key() for graph in graph_enumerator.graphs()
    )
    assert provider.combinatorial_infinity_dependencies(1, 0) == (
        (1, 0, (-1,)),
    )
    solver, values = provider.solve()

    expected = {
        1: -QQ(1) / 24,
        2: QQ(1) / 2880,
        3: -QQ(1) / 181440,
        4: QQ(1) / 9676800,
    }
    for genus, value in expected.items():
        vertex = provider.vertices[genus]
        assert vertex.is_balanced()
        assert values[vertex] == value
        assert values[vertex] == provider.expected_bernoulli_value(genus)

    # Integer partitions index the localization monomials.
    assert len(provider.equation_for_genus(1).terms) == 1
    assert len(provider.equation_for_genus(2).terms) == 2
    assert len(provider.equation_for_genus(3).terms) == 3
    assert len(provider.equation_for_genus(4).terms) == 5

    genus_three = provider.equation_for_genus(3)
    b1 = provider.vertices[1]
    b2 = provider.vertices[2]
    b3 = provider.vertices[3]
    assert set(genus_three.terms) == {
        (QQ(1), (b3,)),
        (QQ(1), (b2, b1)),
        (QQ(1) / 6, (b1, b1, b1)),
    }

    checks = provider.verify(values)
    assert all(
        item["diagonal_nonzero"]
        and item["bernoulli_value"]
        and item["all_q_coefficients"]
        for item in checks.values()
    )

    # Check several coefficients beyond the equations' degree-zero input.
    genus_three_reconstructed = provider.reconstructed_stationary_series(3, values)
    polynomial, genus_three_known = provider.known_stationary_series(3)
    assert str(polynomial) == "1/6*E2^3 + 1/12*E2*E4 + 1/360*E6"
    assert genus_three_reconstructed == genus_three_known

    report = provider.report()
    assert report["vertices"][3]["genus"] == 4
    assert report["vertices"][3]["value"] == "1/9676800"
    assert report["vertices"][3]["equation_terms"] == 5


run_tests()
print("all concrete CJR equation-provider tests passed")
