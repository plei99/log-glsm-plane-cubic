"""Tests for finite-closure infinity-vertex orchestration."""

load("log_glsm_infinity_orchestrator.sage")


def synthetic_fixture(max_markings=2):
    a = EffectiveVertex(1, 0, (-1,), label="a")
    x = EffectiveVertex(2, 0, (-3, -1), label="x")
    y = EffectiveVertex(2, 0, (-2, -2), label="y")

    probe_a = ProbeSpec.stationary(1, 0, (0,), label="solve a")
    probe_1 = ProbeSpec.stationary(2, 0, (2,), label="rank one")
    probe_2 = ProbeSpec.stationary(2, 0, (1, 1), label="rank two")
    relation_a = CompiledLocalizationRelation(
        probe_a, 0, QQ, 1, terms=((1, (a,)),)
    )
    relation_1 = CompiledLocalizationRelation(
        probe_1, 0, QQ, 5,
        terms=((2, (x, a)), (1, (y,))),
    )
    relation_2 = CompiledLocalizationRelation(
        probe_2, 0, QQ, 7,
        terms=((1, (x,)), (3, (y,))),
    )
    calls = []

    def candidates(genus, ambient_degree, stage, config):
        calls.append((int(genus), int(ambient_degree),
                      int(stage.max_markings)))
        if genus == 1:
            return (relation_a,)
        if genus == 2 and stage.max_markings == 1:
            return (relation_1,)
        if genus == 2:
            return relation_1, relation_2
        return tuple()

    config = InfinityOrchestrationConfig(
        max_markings=max_markings,
        t_powers=(0,),
        include_unit_relatives=False,
        max_genus=2,
        max_ambient_degree=0,
    )
    provider = object()
    orchestrator = InfinityVertexOrchestrator(
        provider,
        (x, y),
        config=config,
        coefficient_field=QQ,
        candidate_provider=candidates,
    )
    return orchestrator, (a, x, y), calls, candidates, config


def run_tests():
    orchestrator, vertices, calls, candidates, config = synthetic_fixture(2)
    a, x, y = vertices
    report = orchestrator.run()
    assert report["status"] == "complete"
    assert orchestrator.solver.values[a] == 1
    assert orchestrator.solver.values[x] == QQ(8) / 5
    assert orchestrator.solver.values[y] == QQ(9) / 5
    assert len(orchestrator.block_history) == 2
    assert orchestrator.current_stage == 1
    assert any(call == (2, 0, 1) for call in calls)
    assert any(call == (2, 0, 2) for call in calls)

    checkpoint = orchestrator.checkpoint()
    assert json.loads(json.dumps(checkpoint)) == checkpoint
    restored = InfinityVertexOrchestrator(
        object(),
        (x, y),
        config=config,
        coefficient_field=QQ,
        candidate_provider=lambda *args: (_ for _ in ()).throw(
            AssertionError("a complete checkpoint must not recompile")
        ),
    )
    restored_report = restored.restore_checkpoint(checkpoint)
    assert restored_report["status"] == "complete"
    assert restored.run()["status"] == "complete"
    assert restored.solver.values == orchestrator.solver.values
    assert restored.checkpoint() == checkpoint

    singular, singular_vertices, _, _, _ = synthetic_fixture(1)
    singular_report = singular.run()
    json.dumps(singular_report)
    assert singular_report["status"] == "blocked"
    deficiency = singular_report["frontier"]["rank_deficiencies"][0]
    assert deficiency["rank"] == 1
    assert deficiency["columns"] == 2
    assert deficiency["kernel"]

    # A richer marking schedule reuses the blocked checkpoint and compiles
    # only the newly available stage.
    singular_checkpoint = singular.checkpoint()
    extended, extended_vertices, extended_calls, _, extended_config = \
        synthetic_fixture(2)
    extended.restore_checkpoint(singular_checkpoint)
    assert extended.run()["status"] == "complete"
    assert extended.solver.values[extended_vertices[1]] == QQ(8) / 5
    assert extended.solver.values[extended_vertices[2]] == QQ(9) / 5
    assert extended_calls
    assert all(call[2] == 2 for call in extended_calls)

    # The real plane-cubic provider solves the contact-resolved genus-one
    # base vertex through the same orchestration path.
    provider = PlaneCubicFullEquationProvider(laurent_precision=8)
    genus_one = EffectiveVertex(
        1, 0, (-1,), psi_min=0, insertions=(1,)
    )
    real_config = InfinityOrchestrationConfig(
        max_markings=1,
        t_powers=(0,),
        include_unit_relatives=False,
        max_genus=1,
        max_ambient_degree=0,
    )
    real = InfinityVertexOrchestrator(
        provider, (genus_one,), config=real_config
    )
    real_report = real.run()
    assert real_report["status"] == "complete"
    relation = provider.relation(
        ProbeSpec.stationary(
            1, 0, (0,), label="stationary candidate"
        ),
        0,
    )
    reconstructed = relation.twisted_zero_level + sum(
        coefficient * prod(real.solver.values[factor] for factor in factors)
        for coefficient, factors in relation.terms
    )
    assert reconstructed == relation.known_gw

    restored_provider = PlaneCubicFullEquationProvider(laurent_precision=8)
    restored_real = InfinityVertexOrchestrator(
        restored_provider, (genus_one,), config=real_config
    )
    restored_real.restore_checkpoint(real.checkpoint())
    assert not restored_provider._compilations
    restored_provider.candidate_relations(
        1, 0, max_markings=1, t_powers=(0,),
        require_complete=True, include_unit_relatives=False,
    )
    assert not restored_provider._compilations

    # A truncated real genus-two run discovers all base contact sectors and
    # returns a certified rank frontier instead of conflating the profiles.
    genus_two_provider = PlaneCubicFullEquationProvider(laurent_precision=8)
    genus_two = EffectiveVertex(
        2, 0, (-3,), psi_min=0, insertions=(1,)
    )
    truncated_config = InfinityOrchestrationConfig(
        max_markings=1,
        t_powers=(2, 1, 0),
        include_unit_relatives=False,
        max_genus=2,
        max_ambient_degree=0,
    )
    truncated = InfinityVertexOrchestrator(
        genus_two_provider, (genus_two,), config=truncated_config
    )
    truncated_report = truncated.run()
    assert truncated_report["status"] == "blocked"
    genus_two_profiles = set(
        vertex.contacts for vertex in truncated.vertices
        if vertex.genus == 2
    )
    assert {
        (-3,), (-3, -1), (-2, -2), (-2, -2, -1)
    }.issubset(genus_two_profiles)
    assert truncated_report["frontier"]["rank_deficiencies"]


run_tests()
print("all infinity-orchestrator tests passed")
