"""Tests for finite-closure infinity-vertex orchestration."""

load("log_glsm_infinity_orchestrator.sage")


def synthetic_fixture(max_markings=2, include_kernel_basis=True):
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
        include_kernel_basis=include_kernel_basis,
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
    deep_config = InfinityOrchestrationConfig(
        max_markings=2, t_powers=(0,), include_unit_relatives=True,
        max_unit_insertions=2,
        include_chow_relations=True,
        max_chow_unit_insertions=2,
        chow_primary_only=True,
        additional_probe_ambient_degrees=(1, 2),
    )
    assert [
        (int(stage.max_markings), int(stage.max_unit_insertions))
        for stage in deep_config.stages()
    ] == [(1, 0), (2, 0), (1, 1), (2, 1), (1, 2), (2, 2)]
    assert InfinityOrchestrationConfig.from_record(
        deep_config.to_record()
    ).to_record() == deep_config.to_record()
    assert deep_config.additional_probe_ambient_degrees == (1, 2)
    assert deep_config.include_chow_relations
    assert deep_config.max_chow_unit_insertions == 2
    assert deep_config.chow_primary_only

    effective_config = InfinityOrchestrationConfig(
        max_markings=2,
        t_powers=(-4, -5),
        include_unit_relatives=False,
        max_unit_insertions=0,
        include_chow_relations=True,
        max_chow_unit_insertions=0,
        chow_primary_only=True,
        effective_basis_only=True,
    )
    assert InfinityOrchestrationConfig.from_record(
        effective_config.to_record()
    ).to_record() == effective_config.to_record()
    assert effective_config.effective_basis_only
    try:
        InfinityOrchestrationConfig(
            include_unit_relatives=False,
            include_chow_relations=False,
            chow_primary_only=True,
            effective_basis_only=True,
        )
        raise AssertionError(
            "effective-basis mode must require primary Chow probes"
        )
    except ValueError:
        pass

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

    stale_checkpoint = dict(checkpoint)
    stale_checkpoint["version"] = 1
    try:
        restored.restore_checkpoint(stale_checkpoint)
        raise AssertionError(
            "pre-dimension-filter checkpoints must be rejected"
        )
    except ValueError:
        pass

    singular, singular_vertices, _, _, _ = synthetic_fixture(1)
    singular_report = singular.run()
    json.dumps(singular_report)
    assert singular_report["status"] == "blocked"
    deficiency = singular_report["frontier"]["rank_deficiencies"][0]
    assert deficiency["rank"] == 1
    assert deficiency["columns"] == 2
    assert deficiency["kernel"]

    compact, _, _, _, _ = synthetic_fixture(
        1, include_kernel_basis=False
    )
    compact_deficiency = compact.run()["frontier"]["rank_deficiencies"][0]
    assert compact_deficiency["rank"] == 1
    assert compact_deficiency["columns"] == 2
    assert compact_deficiency["kernel_dimension"] == 1
    assert compact_deficiency["kernel"] is None

    # A singular component may still determine selected coordinates.  Here
    # x+y+z=6 and y+z=4 leave y-z free but force x=2.
    partial_x = EffectiveVertex(2, 0, (-3,), label="partial_x")
    partial_y = EffectiveVertex(2, 0, (-3,), label="partial_y")
    partial_z = EffectiveVertex(2, 0, (-3,), label="partial_z")
    partial_relations = (
        CompiledLocalizationRelation(
            ProbeSpec.stationary(2, 0, (2,), label="partial all"),
            0, QQ, 6,
            terms=((1, (partial_x,)), (1, (partial_y,)),
                   (1, (partial_z,))),
        ),
        CompiledLocalizationRelation(
            ProbeSpec.stationary(2, 0, (1, 1), label="partial pair"),
            0, QQ, 4,
            terms=((1, (partial_y,)), (1, (partial_z,))),
        ),
    )
    partial_config = InfinityOrchestrationConfig(
        max_markings=1, t_powers=(0,), include_unit_relatives=False,
        max_genus=2, max_ambient_degree=0,
    )
    partial = InfinityVertexOrchestrator(
        object(), (partial_x,), config=partial_config,
        coefficient_field=QQ,
        candidate_provider=lambda *args: partial_relations,
    )
    assert partial.run()["status"] == "complete"
    assert partial.solver.values[partial_x] == 2
    assert partial_y not in partial.solver.values
    assert partial_z not in partial.solver.values
    assert partial.block_history[-1]["partial"]

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
        include_punctured_axioms=False,
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

    obsolete = real.checkpoint()
    obsolete["version"] = 6
    try:
        restored_real.restore_checkpoint(obsolete)
        raise AssertionError(
            "checkpoints without Chow-degree filtering must be rejected"
        )
    except ValueError:
        pass
    restored_provider.candidate_relations(
        1, 0, max_markings=1, t_powers=(0,),
        require_complete=True, include_unit_relatives=False,
    )
    assert not restored_provider._compilations

    # Adding a probe ambient degree revisits completed marking stages;
    # otherwise the old stage key would hide the new dimension-zero rows.
    degree_base_config = InfinityOrchestrationConfig(
        max_markings=1, t_powers=(0,), include_unit_relatives=False,
        max_genus=1, max_ambient_degree=0,
    )
    degree_base = InfinityVertexOrchestrator(
        object(), (genus_one,), config=degree_base_config,
        coefficient_field=QQ,
        candidate_provider=lambda *args: tuple(),
        axiom_provider=None,
    )
    degree_base.current_stage = ZZ(0)
    degree_base.compiled_level_stages.add((1, 0, 1, False, 0))
    degree_extended_config = InfinityOrchestrationConfig(
        max_markings=1, t_powers=(0,), include_unit_relatives=False,
        additional_probe_ambient_degrees=(1,),
        max_genus=1, max_ambient_degree=1,
    )
    degree_extended = InfinityVertexOrchestrator(
        object(), (genus_one,), config=degree_extended_config,
        coefficient_field=QQ,
        candidate_provider=lambda *args: tuple(),
        axiom_provider=None,
    )
    degree_extended.restore_checkpoint(degree_base.checkpoint())
    assert degree_extended.current_stage == -1
    assert not degree_extended.compiled_level_stages

    # A truncated real genus-two run now includes the stabilization-corrected
    # t^0 row.  Its contact-descendant support does not falsely touch the old
    # undecorated (-3) target, so this deliberately tiny family still blocks.
    genus_two_provider = PlaneCubicFullEquationProvider(laurent_precision=8)
    genus_two = EffectiveVertex(
        2, 0, (-3,), psi_min=0, insertions=(1,)
    )
    truncated_config = InfinityOrchestrationConfig(
        max_markings=1,
        t_powers=(2, 1, 0),
        include_unit_relatives=False,
        include_punctured_axioms=False,
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
    assert genus_two_profiles == {(-3,)}
    assert truncated_report["active_relation_count"] == 0
    assert {relation.t_power for relation in truncated.relations.values()} \
        == {0}
    assert all(
        relation.probe.psi_convention == "stabilized"
        for relation in truncated.relations.values()
    )
    assert truncated_report["frontier"]["vertices_without_linear_incidence"]


def test_time_budget_stops_between_blocks():
    """An expired budget stops the run cleanly with resumable state."""
    orchestrator, vertices, calls, candidates, config = synthetic_fixture(2)
    report = orchestrator.run(time_budget=float(0))
    assert report["timed_out"]
    assert report["status"] == "blocked"
    # State is intact and checkpointable; a fresh unbudgeted run finishes.
    json.dumps(orchestrator.checkpoint())
    resumed = orchestrator.run()
    assert resumed["status"] == "complete"
    assert not resumed["timed_out"]
    a, x, y = vertices
    assert orchestrator.solver.values[x] == QQ(8) / 5


def _legacy_checkpoint_record():
    """A minimal synthetic version-8 checkpoint for the import test."""
    good_probe = ProbeSpec.stationary(1, 0, (0,), label="legacy good")
    classy_probe = ProbeSpec.stationary(1, 0, (0,), label="legacy classy")
    axiom_probe = ProbeSpec(
        1, 0, (), label="CJR 8.3 string remove=0 vertex=synthetic",
    )
    vertex = EffectiveVertex(1, 0, (-1,), label="a")
    good = CompiledLocalizationRelation(
        good_probe, 0, QQ, 1, terms=((1, (vertex,)),)
    )
    # Positive residual Chow dimension: the historical invalid row class.
    classy = CompiledLocalizationRelation(
        classy_probe, 2, QQ, 0, terms=((1, (vertex,)),)
    )
    axiom = CompiledLocalizationRelation(
        axiom_probe, 0, QQ, 0, terms=((1, (vertex,)),)
    )
    return {
        "format": InfinityVertexOrchestrator.CHECKPOINT_FORMAT,
        "version": int(8),
        "coefficient_field": str(QQ),
        "relations": [
            good.to_record(), classy.to_record(), axiom.to_record(),
        ],
    }, vertex


def test_legacy_relation_import_filters_and_activates():
    """v8 import keeps Chow-scalar rows, drops invalid and axiom rows."""
    import tempfile
    orchestrator, vertices, calls, candidates, config = synthetic_fixture(2)
    checkpoint, vertex = _legacy_checkpoint_record()
    with tempfile.NamedTemporaryFile(
            "w", suffix=".json", delete=False) as stream:
        json.dump(checkpoint, stream)
        legacy_path = stream.name

    before = len(orchestrator.relations)
    imported = orchestrator.import_legacy_relations(legacy_path)
    assert imported == 1
    assert len(orchestrator.relations) == before + 1
    labels = [relation.probe.label
              for relation in orchestrator.relations.values()]
    assert "legacy good" in labels
    assert "legacy classy" not in labels
    assert not any(" remove=" in label for label in labels)
    # Importing again is a no-op: the key already exists.
    assert orchestrator.import_legacy_relations(legacy_path) == 0

    # A modern checkpoint must be refused by the legacy path.
    modern = dict(checkpoint)
    modern["version"] = int(InfinityVertexOrchestrator.CHECKPOINT_VERSION)
    with tempfile.NamedTemporaryFile(
            "w", suffix=".json", delete=False) as stream:
        json.dump(modern, stream)
        modern_path = stream.name
    try:
        orchestrator.import_legacy_relations(modern_path)
        raise AssertionError("modern checkpoints must use restore_checkpoint")
    except ValueError:
        pass

    # The imported row participates in solving: it determines vertex a.
    report = orchestrator.run()
    assert report["status"] == "complete"
    assert orchestrator.solver.values[vertex] == 1


run_tests()
test_time_budget_stops_between_blocks()
test_legacy_relation_import_filters_and_activates()
print("all infinity-orchestrator tests passed")
