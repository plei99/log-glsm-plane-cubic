"""Tests for the genus-three infinity-vertex enumeration and driver."""

load("genus_three_infinity_vertices.sage")


class StubProvider(SageObject):
    """A minimal candidate source with no graph enumeration behind it."""

    def __init__(self, relations_by_level):
        self.relations_by_level = dict(relations_by_level)
        self.calls = []

    def candidate_relations(self, genus, ambient_degree, max_markings=4,
                            t_powers=(0,), require_complete=False,
                            include_unit_relatives=False,
                            max_unit_insertions=1,
                            include_chow_relations=False,
                            max_chow_unit_insertions=1,
                            chow_primary_only=False,
                            effective_basis_only=False):
        self.calls.append((int(genus), int(ambient_degree), int(max_markings)))
        return self.relations_by_level.get(
            (int(genus), int(ambient_degree)), tuple()
        )


def test_allowed_degrees():
    # Genus three is the first genus admitting a positive infinity degree.
    assert allowed_infinity_degrees(0) == tuple()
    assert allowed_infinity_degrees(1) == (0,)
    assert allowed_infinity_degrees(2) == (0,)
    assert allowed_infinity_degrees(3) == (0, 1)
    assert allowed_infinity_degrees(4) == (0, 1, 2)


def test_profiles_generalize_the_degree_zero_function():
    for genus in range(1, 6):
        for valence in range(1, 5):
            assert infinity_contact_profiles_in_degree(genus, valence, 0) \
                == infinity_contact_profiles(genus, valence)

    assert infinity_contact_profiles_in_degree(3, 1, 0) == ((-5,),)
    assert infinity_contact_profiles_in_degree(3, 2, 0) == (
        (-5, -1), (-4, -2), (-3, -3)
    )
    # Contact excess drops from 4 to 1 in the new positive-degree sector.
    assert infinity_contact_profiles_in_degree(3, 1, 1) == ((-2,),)
    assert infinity_contact_profiles_in_degree(3, 2, 1) == ((-2, -1),)
    # Degree two is beyond the balance bound in genus three.
    assert infinity_contact_profiles_in_degree(3, 1, 2) == tuple()


def test_enumeration_is_balanced_canonical_and_deduplicated():
    vertices = enumerate_infinity_vertices(3, max_valence=3)
    assert vertices
    signatures = [vertex.signature() for vertex in vertices]
    assert len(set(signatures)) == len(signatures)

    for vertex in vertices:
        assert vertex.genus == 3
        assert vertex.is_balanced()
        assert vertex.ambient_degree in (0, 1)
        assert len(vertex.contacts) == len(vertex.insertions)
        assert 1 <= len(vertex.contacts) <= 3
        # The compiler stores sorted (contact, insertion) pairs.
        pairs = tuple(zip(vertex.contacts, vertex.insertions))
        assert tuple(sorted(pairs)) == pairs
        assert all(0 <= power <= 2 for power in vertex.insertions)
        required = PlaneCubicDimension.infinity_required_psi_min_power(
            3, vertex.ambient_degree, vertex.contacts,
            sum(vertex.insertions)
        )
        assert vertex.psi_min == required
        assert vertex.is_dimension_zero()

    degrees = set(vertex.ambient_degree for vertex in vertices)
    assert degrees == {0, 1}


def test_enumeration_truncation_controls():
    assert enumerate_infinity_vertices(3, max_valence=1, max_infinity_degree=0) \
        == enumerate_infinity_vertices(
            3, max_valence=1, max_infinity_degree=0
        )
    only_degree_zero = enumerate_infinity_vertices(
        3, max_valence=2, max_infinity_degree=0
    )
    assert all(vertex.ambient_degree == 0 for vertex in only_degree_zero)

    capped = enumerate_infinity_vertices(3, max_valence=2, max_psi_min=0)
    assert all(vertex.psi_min == 0 for vertex in capped)

    dimension_zero_only = enumerate_infinity_vertices(
        3, max_valence=1, max_psi_min=0, insertion_powers=(1,)
    )
    assert len(dimension_zero_only) == 1
    assert dimension_zero_only[0] == EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(1,)
    )

    both_degrees = enumerate_infinity_vertices(
        3, max_valence=1, max_psi_min=3, insertion_powers=(1,)
    )
    assert len(both_degrees) == 2
    assert {
        (vertex.ambient_degree, vertex.contacts, vertex.psi_min)
        for vertex in both_degrees
    } == {(0, (-5,), 0), (1, (-2,), 3)}

    # Regression: for (g,D,c)=(3,0,-5) and psi_min=0, virtual
    # dimension permits H but prunes 1 and H^2.
    one_mark_degree_zero = enumerate_infinity_vertices(
        3, max_valence=1, max_infinity_degree=0, max_psi_min=0
    )
    assert one_mark_degree_zero == (EffectiveVertex(
        3, 0, (-5,), psi_min=0, insertions=(1,)
    ),)

    # A larger valence cap only ever adds targets.
    small = set(enumerate_infinity_vertices(3, max_valence=1))
    large = set(enumerate_infinity_vertices(3, max_valence=2))
    assert small.issubset(large)


def test_probe_ambient_degree_uses_zero_hypersurface_sectors():
    assert [probe_ambient_degree(d) for d in range(7)] == list(range(7))
    for degree in range(7):
        probe = ProbeSpec.stationary(3, probe_ambient_degree(degree), (4,))
        if degree % 3:
            assert EllipticProbeValueBackend().value(probe) == 0


def test_driver_solves_a_synthetic_genus_three_block():
    root = EffectiveVertex(3, 0, (-5,), psi_min=0, insertions=(1,))
    lower = EffectiveVertex(1, 0, (-1,), psi_min=0, insertions=(1,))
    probe_lower = ProbeSpec.stationary(1, 0, (0,), label="lower")
    probe_root = ProbeSpec.stationary(3, 0, (4,), label="root")
    relations = {
        (1, 0): (CompiledLocalizationRelation(
            probe_lower, 0, QQ, -1, terms=((1, (lower,)),)
        ),),
        (3, 0): (CompiledLocalizationRelation(
            probe_root, 0, QQ, 5, terms=((2, (root,)), (3, (lower,))),
        ),),
    }
    provider = StubProvider(relations)
    computation = GenusThreeVertexComputation(
        max_valence=1,
        max_infinity_degree=0,
        max_psi_min=0,
        insertion_powers=(1,),
        max_markings=1,
        t_powers=(0,),
        include_unit_relatives=False,
        include_punctured_axioms=False,
        provider=provider,
        require_complete=False,
    )
    assert computation.roots == (root,)
    report = computation.run()
    assert report["status"] == "complete"
    assert report["solved_count"] == 1
    # 5 = 2*root + 3*(-1)  =>  root = 4.
    assert computation.orchestrator.solver.values[root] == 4
    assert computation.orchestrator.solver.values[lower] == -1
    assert computation.config.max_vertices == 20000
    assert computation.config.max_relations == 20000
    assert report["targets"][0]["value"] == "4"
    assert json.loads(json.dumps(report)) == report


def test_driver_wires_the_effective_invariant_basis():
    provider = StubProvider({})
    computation = GenusThreeVertexComputation(
        max_valence=1,
        max_infinity_degree=0,
        max_psi_min=0,
        insertion_powers=(1,),
        max_markings=1,
        t_powers=(-4,),
        include_unit_relatives=False,
        include_chow_relations=True,
        max_chow_unit_insertions=0,
        chow_primary_only=True,
        effective_basis_only=True,
        include_punctured_axioms=False,
        provider=provider,
        require_complete=False,
    )
    assert computation.config.effective_basis_only
    assert not computation.config.include_unit_relatives
    assert computation.config.include_chow_relations
    assert computation.config.chow_primary_only
    computation.run()
    assert provider.calls


def test_driver_uses_direct_zero_probe_for_positive_degree_levels():
    # max_infinity_degree is a cap, so both balance-allowed sectors are roots.
    degree_zero = EffectiveVertex(3, 0, (-5,), psi_min=0, insertions=(1,))
    degree_one = EffectiveVertex(3, 1, (-2,), psi_min=3, insertions=(1,))
    constant_probe = ProbeSpec.stationary(3, 0, (4,), label="degree zero")
    curve_probe = ProbeSpec.stationary(3, 1, (4,), label="degree one sector")
    provider = StubProvider({
        (3, 0): (CompiledLocalizationRelation(
            constant_probe, 0, QQ, 8, terms=((4, (degree_zero,)),)
        ),),
        (3, 1): (CompiledLocalizationRelation(
            curve_probe, 0, QQ, 6, terms=((3, (degree_one,)),)
        ),),
    })
    computation = GenusThreeVertexComputation(
        max_valence=1,
        max_infinity_degree=1,
        max_psi_min=3,
        insertion_powers=(1,),
        max_markings=1,
        t_powers=(0,),
        include_unit_relatives=False,
        include_punctured_axioms=False,
        provider=provider,
        require_complete=False,
    )
    assert computation.max_probe_degree == 1
    assert computation.roots == (degree_zero, degree_one)
    report = computation.run()
    assert report["status"] == "complete"
    assert computation.orchestrator.solver.values[degree_zero] == 2
    assert computation.orchestrator.solver.values[degree_one] == 2

    # The D=1 level uses its direct ambient class.  The elliptic GW side is
    # zero because one is not divisible by three, but localization still
    # gives a valid and substantially smaller equation.
    assert (3, 1, 1) in provider.calls
    assert (3, 0, 1) in provider.calls
    assert not any(call[1] == 3 for call in provider.calls)


def test_future_probe_degree_ceiling_preserves_rich_configuration():
    provider = StubProvider({})
    computation = GenusThreeVertexComputation(
        max_valence=1, max_infinity_degree=0, max_psi_min=0,
        insertion_powers=(1,), max_markings=1, t_powers=(0,),
        include_unit_relatives=False, include_punctured_axioms=False,
        additional_probe_ambient_degrees=(1, 2),
        provider=provider, require_complete=False,
    )
    assert computation.config.additional_probe_ambient_degrees == (1, 2)
    computation.restrict_future_probe_degree(0)
    computation._candidates(
        3, 0, ProbeExpansionStage(1, False), computation.config
    )
    assert provider.calls == [(3, 0, 1)]
    # The serialized mathematical configuration is unchanged, so a richer
    # checkpoint remains compatible even though future work is targeted.
    assert computation.config.additional_probe_ambient_degrees == (1, 2)


def test_driver_reports_an_unsolved_target_instead_of_guessing():
    provider = StubProvider({})
    computation = GenusThreeVertexComputation(
        max_valence=1,
        max_infinity_degree=0,
        max_psi_min=0,
        insertion_powers=(1,),
        max_markings=1,
        t_powers=(0,),
        include_unit_relatives=False,
        include_punctured_axioms=False,
        provider=provider,
        require_complete=False,
    )
    report = computation.run()
    assert report["status"] == "blocked"
    assert report["solved_count"] == 0
    assert report["stages_exhausted"]
    assert report["targets"][0]["value"] is None
    assert report["orchestration"]["frontier"] is not None
    assert json.loads(json.dumps(report)) == report


def test_time_budget_is_honoured_between_blocks():
    r"""A slow block must not run the whole stage past the budget.

    The first driver checked the clock only between probe stages, so one
    long-running stage could overrun a three-hour budget many times over.
    """
    root = EffectiveVertex(3, 0, (-5,), psi_min=0, insertions=(1,))
    lower = EffectiveVertex(1, 0, (-1,), psi_min=0, insertions=(1,))
    probe_lower = ProbeSpec.stationary(1, 0, (0,), label="lower")
    probe_root = ProbeSpec.stationary(3, 0, (4,), label="root")
    provider = StubProvider({
        (1, 0): (CompiledLocalizationRelation(
            probe_lower, 0, QQ, -1, terms=((1, (lower,)),)
        ),),
        (3, 0): (CompiledLocalizationRelation(
            probe_root, 0, QQ, 5, terms=((2, (root,)), (3, (lower,))),
        ),),
    })
    computation = GenusThreeVertexComputation(
        max_valence=1, max_infinity_degree=0, max_psi_min=0,
        insertion_powers=(1,), max_markings=1, t_powers=(0,),
        include_unit_relatives=False, include_punctured_axioms=False,
        provider=provider,
        require_complete=False,
    )

    # A block that is slow and never reaches the roots.  Without a mid-stage
    # clock check this loops max_solve_rounds times regardless of the budget,
    # which is exactly how a three-hour budget became an overnight run.
    blocks = []

    def slow_block():
        blocks.append(1)
        time.sleep(float(0.3))
        return {}

    computation.solve_next_block_with_report = slow_block
    started = time.time()
    report = computation.run(time_budget=float(0.2))
    elapsed = time.time() - started

    assert report["timed_out"]
    assert report["status"] == "blocked"
    assert report["solved_count"] == 0
    # Exactly one block ran, so the overrun is bounded by a single block
    # rather than by max_solve_rounds.
    assert len(blocks) == 1
    assert elapsed < 2, elapsed
    assert json.loads(json.dumps(report)) == report


def test_enumeration_report_is_serializable():
    report = enumeration_report(3, max_valence=2)
    assert report["genus"] == 3
    assert report["allowed_infinity_degrees"] == [0, 1]
    assert report["target_count"] == len(report["targets"])
    assert set(report["counts_by_degree"]) == {"0", "1"}
    assert json.loads(json.dumps(report)) == report


def test_driver_exposes_progressive_mixed_unit_depths():
    provider = StubProvider({})
    computation = GenusThreeVertexComputation(
        max_valence=1, max_infinity_degree=0, max_psi_min=0,
        insertion_powers=(1,), max_markings=1, t_powers=(0,),
        include_unit_relatives=True, max_unit_insertions=2,
        provider=provider, require_complete=False,
    )
    stages = computation.config.stages()
    assert [int(stage.max_unit_insertions) for stage in stages] == [0, 1, 2]
    assert json.loads(json.dumps(computation.config.to_record())) \
        == computation.config.to_record()


def run_tests():
    test_allowed_degrees()
    test_profiles_generalize_the_degree_zero_function()
    test_enumeration_is_balanced_canonical_and_deduplicated()
    test_enumeration_truncation_controls()
    test_probe_ambient_degree_uses_zero_hypersurface_sectors()
    test_driver_solves_a_synthetic_genus_three_block()
    test_driver_wires_the_effective_invariant_basis()
    test_driver_uses_direct_zero_probe_for_positive_degree_levels()
    test_driver_reports_an_unsolved_target_instead_of_guessing()
    test_time_budget_is_honoured_between_blocks()
    test_enumeration_report_is_serializable()
    test_driver_exposes_progressive_mixed_unit_depths()


run_tests()
print("all genus-three infinity-vertex tests passed")
