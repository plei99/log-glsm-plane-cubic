"""Tests for the triangular log-GLSM infinity-vertex DP engine."""

load("log_glsm_infinity_dp.sage")


def run_tests():
    genus_one = EffectiveVertex(1, 0, (-1,), label="basic")
    genus_two_31 = EffectiveVertex(2, 0, (-3, -1), label="profile_31")
    genus_two_22 = EffectiveVertex(2, 0, (-2, -2), label="profile_22")
    positive_degree = EffectiveVertex(3, 1, (-2,), label="positive_degree")

    assert genus_one.is_balanced()
    assert genus_two_31.is_balanced()
    assert genus_two_22.is_balanced()
    assert positive_degree.is_balanced()
    assert EffectiveVertex(2, 1, (-1,), label="impossible").balance_defect() > 0

    dimension_one = EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(1,)
    )
    too_small = EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(0,)
    )
    too_large = EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(2,)
    )
    assert dimension_one.reduced_virtual_dimension == 1
    assert dimension_one.insertion_codimension == 1
    assert dimension_one.is_dimension_zero()
    assert too_small.dimension_defect() == 1
    assert too_large.dimension_defect() == -1

    contact_descendant = EffectiveVertex(
        2, 0, (-1, -2, -2), psi_min=0,
        insertions=(1, 0, 0), contact_psi=(2, 0, 0),
    )
    assert contact_descendant.contacts == (-2, -2, -1)
    assert contact_descendant.insertions == (0, 0, 1)
    assert contact_descendant.contact_psi == (0, 0, 2)
    assert contact_descendant.insertion_codimension == 3
    assert contact_descendant.is_dimension_zero()
    assert EffectiveVertex.from_record(contact_descendant.to_record()) \
        == contact_descendant
    assert EffectiveVertex(
        2, 0, (-2, -2, -1), psi_min=2,
        insertions=(0, 0, 1), contact_psi=(0, 0, 0),
    ) == EffectiveVertex(
        2, 0, (-2, -2, -1), psi_min=2, insertions=(0, 0, 1)
    )

    demo = genus_two_bo_demo(4)
    b2 = demo["genus_one_vertex"]
    b4 = demo["genus_two_vertex"]
    assert demo["values"][b2] == -QQ(1) / 24
    assert demo["values"][b4] == QQ(1) / 2880
    assert demo["dependency_graph"][b4] == (b2,)

    # Same-rank contact profiles need a matrix block rather than a scalar
    # lexicographic step: 2*x+y=5 and x+3*y=7.
    def block_provider(vertex):
        if vertex == genus_two_31:
            return LocalizationEquation(
                genus_two_31,
                known_gw=5,
                terms=((2, (genus_two_31,)), (1, (genus_two_22,))),
                probe_label="profile-31 probe",
            )
        if vertex == genus_two_22:
            return LocalizationEquation(
                genus_two_22,
                known_gw=7,
                terms=((1, (genus_two_31,)), (3, (genus_two_22,))),
                probe_label="profile-22 probe",
            )
        raise KeyError(vertex)

    block_solver = InfinityVertexDP(block_provider)
    block = block_solver.solve_block((genus_two_31, genus_two_22))
    assert block[genus_two_31] == QQ(8) / 5
    assert block[genus_two_22] == QQ(9) / 5
    checkpoint = block_solver.checkpoint()
    restored_solver = InfinityVertexDP(block_provider)
    restored = restored_solver.restore_checkpoint(checkpoint)
    assert restored == block
    assert restored_solver.checkpoint() == checkpoint

    # A single probe cannot solve a same-rank dependency outside its block.
    scalar_solver = InfinityVertexDP(block_provider)
    try:
        scalar_solver.solve(genus_two_31)
        raise AssertionError("a non-triangular scalar probe should fail")
    except NonTriangularLocalizationError:
        pass

    def singular_provider(vertex):
        if vertex == genus_two_31:
            return LocalizationEquation(
                vertex, 1,
                terms=((1, (genus_two_31,)), (1, (genus_two_22,))),
            )
        return LocalizationEquation(
            vertex, 2,
            terms=((2, (genus_two_31,)), (2, (genus_two_22,))),
        )

    try:
        InfinityVertexDP(singular_provider).solve_block(
            (genus_two_31, genus_two_22)
        )
        raise AssertionError("a singular probe block should fail")
    except SingularProbeError:
        pass


run_tests()
print("all log-GLSM infinity DP tests passed")
