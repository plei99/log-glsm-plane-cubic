r"""Finite-closure orchestration for contact-resolved infinity vertices.

The low-level DP engine solves an equation after its probe has been chosen.
This module supplies the missing controller:

* generate progressively richer known probes;
* activate only relations incident to the requested dependency closure;
* discover every exact ``EffectiveVertex`` factor in those relations;
* substitute solved values and find linear-ready blocks;
* select exact full-rank rows and solve them with ``InfinityVertexDP``;
* report a certified kernel/frontier when the configured probes stall; and
* checkpoint extracted relations as well as solved values.

The algorithm is finite for every configured truncation.  It never resolves
a singular block by a lexicographic tie-break: lack of rank is returned as a
diagnostic requiring a richer probe family or additional mathematical input.
"""

import argparse
import json
import os
import sys

load("cjr_full_equation_provider.sage")


class OrchestrationLimitError(RuntimeError):
    """A caller-supplied finite resource bound was exceeded."""


class ProbeExpansionStage(SageObject):
    """One progressively richer family of dimension-compatible probes."""

    def __init__(self, max_markings, include_unit_relatives=False):
        self.max_markings = ZZ(max_markings)
        self.include_unit_relatives = bool(include_unit_relatives)
        if self.max_markings < 1:
            raise ValueError("a probe stage needs at least one marking")

    def to_record(self):
        return {
            "max_markings": int(self.max_markings),
            "include_unit_relatives": self.include_unit_relatives,
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["max_markings"],
            record.get("include_unit_relatives", False),
        )

    def __repr__(self):
        suffix = "+unit" if self.include_unit_relatives else ""
        return "ProbeExpansionStage(markings<=%s%s)" % (
            self.max_markings, suffix
        )


class InfinityOrchestrationConfig(SageObject):
    """Finite bounds and probe schedule for an orchestration run."""

    def __init__(self, max_markings=2, t_powers=(2, 1, 0, -1, -2),
                 include_unit_relatives=True, max_genus=None,
                 max_ambient_degree=None, max_vertices=2000,
                 max_relations=4000, max_solve_rounds=1000,
                 require_complete=True, fail_on_unsupported=False):
        self.max_markings = ZZ(max_markings)
        self.t_powers = tuple(ZZ(power) for power in t_powers)
        self.include_unit_relatives = bool(include_unit_relatives)
        self.max_genus = None if max_genus is None else ZZ(max_genus)
        self.max_ambient_degree = (
            None if max_ambient_degree is None
            else ZZ(max_ambient_degree)
        )
        self.max_vertices = ZZ(max_vertices)
        self.max_relations = ZZ(max_relations)
        self.max_solve_rounds = ZZ(max_solve_rounds)
        self.require_complete = bool(require_complete)
        self.fail_on_unsupported = bool(fail_on_unsupported)

        if self.max_markings < 1:
            raise ValueError("max_markings must be positive")
        if not self.t_powers:
            raise ValueError("at least one Laurent coefficient is required")
        if len(set(self.t_powers)) != len(self.t_powers):
            raise ValueError("t_powers cannot contain duplicates")
        if self.max_genus is not None and self.max_genus < 0:
            raise ValueError("max_genus must be nonnegative")
        if self.max_ambient_degree is not None \
                and self.max_ambient_degree < 0:
            raise ValueError("max_ambient_degree must be nonnegative")
        if self.max_vertices < 1 or self.max_relations < 1:
            raise ValueError("resource bounds must be positive")
        if self.max_solve_rounds < 1:
            raise ValueError("max_solve_rounds must be positive")

    def stages(self):
        stages = [
            ProbeExpansionStage(markings, False)
            for markings in range(1, self.max_markings + 1)
        ]
        if self.include_unit_relatives:
            stages.extend(
                ProbeExpansionStage(markings, True)
                for markings in range(1, self.max_markings + 1)
            )
        return tuple(stages)

    def to_record(self):
        return {
            "max_markings": int(self.max_markings),
            "t_powers": [int(power) for power in self.t_powers],
            "include_unit_relatives": self.include_unit_relatives,
            "max_genus": (
                None if self.max_genus is None else int(self.max_genus)
            ),
            "max_ambient_degree": (
                None if self.max_ambient_degree is None
                else int(self.max_ambient_degree)
            ),
            "max_vertices": int(self.max_vertices),
            "max_relations": int(self.max_relations),
            "max_solve_rounds": int(self.max_solve_rounds),
            "require_complete": self.require_complete,
            "fail_on_unsupported": self.fail_on_unsupported,
        }

    @classmethod
    def from_record(cls, record):
        return cls(**record)


class InfinityVertexOrchestrator(SageObject):
    r"""Discover and solve a finite contact-resolved dependency closure.

    ``provider`` normally is a ``PlaneCubicFullEquationProvider``.  Tests and
    extensions may inject

    ``candidate_provider(genus, ambient_degree, stage, config)``

    returning exact ``CompiledLocalizationRelation`` objects.
    """

    CHECKPOINT_FORMAT = "log-glsm-infinity-orchestrator"
    CHECKPOINT_VERSION = 1

    def __init__(self, provider, roots, initial_values=None, config=None,
                 coefficient_field=None, candidate_provider=None):
        self.provider = provider
        self.roots = tuple(roots)
        if not self.roots:
            raise ValueError("the orchestrator needs at least one root vertex")
        if any(not isinstance(vertex, EffectiveVertex)
               for vertex in self.roots):
            raise TypeError("every orchestration root must be an EffectiveVertex")

        self.config = config or InfinityOrchestrationConfig()
        if not isinstance(self.config, InfinityOrchestrationConfig):
            raise TypeError("config must be an InfinityOrchestrationConfig")
        if coefficient_field is None:
            rings = getattr(provider, "rings", None)
            coefficient_field = getattr(rings, "base_field", QQ)
        self.coefficient_field = coefficient_field
        self.candidate_provider = (
            candidate_provider or self._plane_cubic_candidates
        )

        initial_values = dict(initial_values or {})
        self._assignments = {}
        self.solver = InfinityVertexDP(
            self._equation_for,
            initial_values=initial_values,
            coefficient_field=self.coefficient_field,
        )
        self.vertices = set(self.roots) | set(initial_values)
        self.relations = {}
        self.active_relation_keys = set()
        self.compiled_level_stages = set()
        self.failed_level_stages = {}
        self.current_stage = -1
        self.block_history = []
        self._last_frontier = None

        inferred_vertices = tuple(self.vertices)
        self._max_genus = (
            max(vertex.genus for vertex in inferred_vertices)
            if self.config.max_genus is None
            else self.config.max_genus
        )
        self._max_ambient_degree = (
            max(vertex.ambient_degree for vertex in inferred_vertices)
            if self.config.max_ambient_degree is None
            else self.config.max_ambient_degree
        )

    def _plane_cubic_candidates(self, genus, ambient_degree, stage, config):
        return self.provider.candidate_relations(
            genus,
            ambient_degree,
            max_markings=stage.max_markings,
            t_powers=config.t_powers,
            require_complete=config.require_complete,
            include_unit_relatives=stage.include_unit_relatives,
        )

    def _equation_for(self, target):
        if target not in self._assignments:
            raise KeyError("the orchestrator has not assigned a probe to %r" % target)
        return self._assignments[target].equation_for(target)

    @staticmethod
    def _relation_key(relation):
        return relation.probe.signature(), int(relation.t_power)

    def _cache_restored_relation_in_provider(self, relation):
        """Let the concrete provider reuse a checkpointed extraction."""
        cache = getattr(self.provider, "_relations", None)
        if isinstance(cache, dict):
            cache[(relation.probe, ZZ(relation.t_power))] = relation

    @staticmethod
    def _relation_vertices(relation):
        return set(
            factor
            for coefficient, factors in relation.terms
            for factor in factors
        )

    @staticmethod
    def _level_stage_key(genus, ambient_degree, stage):
        return (
            int(genus),
            int(ambient_degree),
            int(stage.max_markings),
            bool(stage.include_unit_relatives),
        )

    def _level_allowed(self, genus, ambient_degree):
        return (
            ZZ(0) <= genus <= self._max_genus
            and ZZ(0) <= ambient_degree <= self._max_ambient_degree
        )

    def _compile_level_stage(self, level, stage_index):
        genus, ambient_degree = level
        stage = self.config.stages()[stage_index]
        key = self._level_stage_key(genus, ambient_degree, stage)
        if key in self.compiled_level_stages or key in self.failed_level_stages:
            return False
        if not self._level_allowed(genus, ambient_degree):
            self.failed_level_stages[key] = "level is outside configured bounds"
            return False

        try:
            relations = tuple(self.candidate_provider(
                genus, ambient_degree, stage, self.config
            ))
        except (UnsupportedGeometryError, TruncationError, ValueError) as error:
            if self.config.fail_on_unsupported:
                raise
            self.failed_level_stages[key] = str(error)
            return False

        new_relations = 0
        for relation in relations:
            if not isinstance(relation, CompiledLocalizationRelation):
                raise TypeError(
                    "candidate_provider must return compiled localization relations"
                )
            if relation.coefficient_field is not self.coefficient_field:
                raise TypeError(
                    "candidate relation and orchestration coefficient fields differ"
                )
            relation_key = self._relation_key(relation)
            if relation_key in self.relations:
                continue
            if len(self.relations) >= self.config.max_relations:
                raise OrchestrationLimitError(
                    "relation bound %s exceeded" % self.config.max_relations
                )
            self.relations[relation_key] = relation
            new_relations += 1
        self.compiled_level_stages.add(key)
        return bool(new_relations)

    def _activate_incident_relations(self):
        changed = False
        while True:
            local_change = False
            for key, relation in self.relations.items():
                if key in self.active_relation_keys:
                    continue
                factors = self._relation_vertices(relation)
                if not factors.intersection(self.vertices):
                    continue
                self.active_relation_keys.add(key)
                new_vertices = factors.difference(self.vertices)
                if len(self.vertices) + len(new_vertices) \
                        > self.config.max_vertices:
                    raise OrchestrationLimitError(
                        "vertex bound %s exceeded" % self.config.max_vertices
                    )
                self.vertices.update(new_vertices)
                local_change = True
                changed = True
            if not local_change:
                return changed

    def expand_closure(self):
        """Compile the current stage and close under incident relations."""
        if self.current_stage < 0:
            raise ValueError("advance to a probe stage before expanding")
        changed = True
        any_change = False
        while changed:
            changed = self._activate_incident_relations()
            any_change = any_change or changed
            levels = sorted(set(
                (vertex.genus, vertex.ambient_degree)
                for vertex in self.vertices
                if vertex not in self.solver.values
            ))
            for level in levels:
                for stage_index in range(self.current_stage + 1):
                    compiled = self._compile_level_stage(level, stage_index)
                    changed = compiled or changed
                    any_change = compiled or any_change
        return any_change

    def advance_stage(self):
        """Advance one probe-family stage and expand the active closure."""
        if self.current_stage + 1 >= len(self.config.stages()):
            return False
        self.current_stage += 1
        self.expand_closure()
        return True

    def _relation_frontier(self, relation):
        support = set()
        nonlinear_monomials = []
        for coefficient, factors in relation.terms:
            unsolved = tuple(
                factor for factor in factors
                if factor not in self.solver.values
            )
            if len(unsolved) > 1:
                nonlinear_monomials.append(unsolved)
            elif len(unsolved) == 1:
                support.add(unsolved[0])
        return support, tuple(nonlinear_monomials)

    def _relation_priority(self, relation):
        t_order = {
            power: index for index, power in enumerate(self.config.t_powers)
        }
        has_unit = any(
            insertion.kind == "unit"
            for insertion in relation.probe.insertions
        )
        return (
            t_order.get(relation.t_power, len(t_order)),
            bool(has_unit),
            int(relation.probe.marking_count),
            relation.probe.signature(),
        )

    def _ready_components(self):
        ready = []
        nonlinear = []
        active_relations = sorted(
            (self.relations[key] for key in self.active_relation_keys),
            key=self._relation_priority,
        )
        for relation in active_relations:
            support, nonlinear_monomials = self._relation_frontier(relation)
            if nonlinear_monomials:
                nonlinear.append((relation, nonlinear_monomials))
            elif support:
                ready.append((relation, support))

        adjacency = {}
        for relation, support in ready:
            for vertex in support:
                adjacency.setdefault(vertex, set()).update(
                    support.difference((vertex,))
                )
        components = []
        unseen = set(adjacency)
        while unseen:
            start = unseen.pop()
            component = {start}
            pending = [start]
            while pending:
                vertex = pending.pop()
                neighbours = adjacency.get(vertex, set()).intersection(unseen)
                unseen.difference_update(neighbours)
                component.update(neighbours)
                pending.extend(neighbours)
            component_relations = tuple(
                relation for relation, support in ready
                if support.intersection(component)
            )
            components.append((component, component_relations))
        return components, tuple(nonlinear)

    @staticmethod
    def _component_priority(component):
        return max(
            (
                vertex.genus,
                vertex.ambient_degree,
                len(vertex.contacts),
                vertex.psi_min,
                vertex.insertions,
                vertex.contacts,
            )
            for vertex in component
        )

    def _selection_for_component(self, component, relations):
        targets = tuple(sorted(component, key=lambda item: item.order_key()))
        factory = ProbeFactory(self.coefficient_field)
        selection = factory.select_relations(
            targets,
            relations,
            lower_values=self.solver.values,
        )
        return targets, selection

    def solve_next_block(self):
        """Solve one highest-priority full-rank linear-ready block."""
        components, nonlinear = self._ready_components()
        candidates = []
        deficiencies = []
        for component, relations in components:
            targets, selection = self._selection_for_component(
                component, relations
            )
            item = (component, relations, targets, selection)
            if selection.is_full_rank:
                candidates.append(item)
            else:
                deficiencies.append(item)

        self._last_frontier = self._frontier_record(
            deficiencies, nonlinear
        )
        if not candidates:
            return None
        candidates.sort(
            key=lambda item: (
                self._component_priority(item[0]),
                tuple(target.order_key() for target in item[2]),
            ),
            reverse=True,
        )
        component, relations, targets, selection = candidates[0]
        for target, relation in zip(targets, selection.relations):
            self._assignments[target] = relation
        values = self.solver.solve_block(targets)
        record = {
            "targets": [target.to_record() for target in targets],
            "probes": [relation.probe.to_record()
                       for relation in selection.relations],
            "t_powers": [int(relation.t_power)
                         for relation in selection.relations],
            "matrix": [[str(value) for value in row]
                       for row in selection.matrix.rows()],
            "values": [str(values[target]) for target in targets],
        }
        self.block_history.append(record)
        return values

    def solve_until_stalled(self):
        """Repeatedly solve every block made ready by prior substitutions."""
        solved_blocks = 0
        for _ in range(self.config.max_solve_rounds):
            values = self.solve_next_block()
            if values is None:
                return solved_blocks
            solved_blocks += 1
            if self.roots_solved():
                return solved_blocks
        raise OrchestrationLimitError(
            "solve-round bound %s exceeded" % self.config.max_solve_rounds
        )

    def roots_solved(self):
        return all(root in self.solver.values for root in self.roots)

    def run(self):
        """Expand probe stages until the roots solve or the schedule stalls."""
        while True:
            if self.roots_solved():
                break
            if self.current_stage < 0:
                if not self.advance_stage():
                    break
            self.solve_until_stalled()
            if self.roots_solved():
                break
            if not self.advance_stage():
                break
        return self.report()

    def _frontier_record(self, deficiencies, nonlinear):
        deficiency_records = []
        ready_vertices = set()
        for component, relations, targets, selection in deficiencies:
            ready_vertices.update(component)
            deficiency_records.append({
                "targets": [target.to_record() for target in targets],
                "rank": int(selection.rank),
                "columns": len(targets),
                "kernel": [[str(value) for value in vector]
                           for vector in selection.kernel_basis()],
                "candidate_rows": len(relations),
            })
        unsolved = set(self.vertices).difference(self.solver.values)
        nonlinear_vertices = set(
            factor
            for relation, monomials in nonlinear
            for monomial in monomials
            for factor in monomial
        )
        no_linear_incidence = unsolved.difference(
            ready_vertices | nonlinear_vertices
        )
        return {
            "rank_deficiencies": deficiency_records,
            "nonlinear_relation_count": len(nonlinear),
            "vertices_without_linear_incidence": [
                vertex.to_record()
                for vertex in sorted(
                    no_linear_incidence, key=lambda item: item.order_key()
                )
            ],
        }

    def frontier_report(self):
        if self._last_frontier is None:
            components, nonlinear = self._ready_components()
            deficiencies = []
            for component, relations in components:
                targets, selection = self._selection_for_component(
                    component, relations
                )
                if not selection.is_full_rank:
                    deficiencies.append(
                        (component, relations, targets, selection)
                    )
            self._last_frontier = self._frontier_record(
                deficiencies, nonlinear
            )
        return self._last_frontier

    def report(self):
        status = "complete" if self.roots_solved() else "blocked"
        return {
            "status": status,
            "stage": int(self.current_stage),
            "stage_count": len(self.config.stages()),
            "roots": [
                {
                    "vertex": root.to_record(),
                    "value": (
                        str(self.solver.values[root])
                        if root in self.solver.values else None
                    ),
                }
                for root in self.roots
            ],
            "vertex_count": len(self.vertices),
            "solved_vertex_count": len(self.solver.values),
            "compiled_relation_count": len(self.relations),
            "active_relation_count": len(self.active_relation_keys),
            "solved_blocks": list(self.block_history),
            "failed_level_stages": [
                {
                    "genus": key[0],
                    "ambient_degree": key[1],
                    "max_markings": key[2],
                    "include_unit_relatives": key[3],
                    "reason": reason,
                }
                for key, reason in sorted(self.failed_level_stages.items())
            ],
            "frontier": None if status == "complete" else self.frontier_report(),
        }

    def checkpoint(self):
        """Serialize compiled relations, closure state, and exact DP values."""
        relations = sorted(
            self.relations.values(),
            key=lambda relation: (
                relation.probe.signature(), relation.t_power
            ),
        )
        return {
            "format": self.CHECKPOINT_FORMAT,
            "version": int(self.CHECKPOINT_VERSION),
            "coefficient_field": str(self.coefficient_field),
            "config": self.config.to_record(),
            "roots": [root.to_record() for root in self.roots],
            "current_stage": int(self.current_stage),
            "current_stage_spec": (
                None if self.current_stage < 0
                else self.config.stages()[self.current_stage].to_record()
            ),
            "vertices": [
                vertex.to_record()
                for vertex in sorted(
                    self.vertices, key=lambda item: item.order_key()
                )
            ],
            "relations": [relation.to_record() for relation in relations],
            "compiled_level_stages": [
                {
                    "genus": key[0],
                    "ambient_degree": key[1],
                    "stage": ProbeExpansionStage(
                        key[2], key[3]
                    ).to_record(),
                }
                for key in sorted(self.compiled_level_stages)
            ],
            "failed_level_stages": [
                {
                    "genus": key[0],
                    "ambient_degree": key[1],
                    "stage": ProbeExpansionStage(
                        key[2], key[3]
                    ).to_record(),
                    "reason": reason,
                }
                for key, reason in sorted(self.failed_level_stages.items())
            ],
            "solver": self.solver.checkpoint(),
            "block_history": list(self.block_history),
        }

    def restore_checkpoint(self, checkpoint):
        """Restore without recompiling relations already in the checkpoint."""
        if checkpoint.get("format") != self.CHECKPOINT_FORMAT \
                or checkpoint.get("version") != self.CHECKPOINT_VERSION:
            raise ValueError("unsupported orchestration checkpoint format")
        roots = tuple(
            EffectiveVertex.from_record(record)
            for record in checkpoint.get("roots", ())
        )
        if roots != self.roots:
            raise ValueError("checkpoint roots do not match this orchestrator")
        saved_config = InfinityOrchestrationConfig.from_record(
            checkpoint["config"]
        )
        if not self._config_extends(saved_config):
            raise ValueError(
                "checkpoint configuration is not compatible with this run"
            )
        if checkpoint.get("coefficient_field") != str(self.coefficient_field):
            raise ValueError("checkpoint coefficient field does not match")

        self.relations = {}
        for record in checkpoint.get("relations", ()):
            relation = CompiledLocalizationRelation.from_record(
                record, self.coefficient_field
            )
            self.relations[self._relation_key(relation)] = relation
            self._cache_restored_relation_in_provider(relation)
        current_stage_spec = checkpoint.get("current_stage_spec")
        if current_stage_spec is None:
            self.current_stage = ZZ(-1)
        else:
            matches = [
                index for index, stage in enumerate(self.config.stages())
                if stage.to_record() == current_stage_spec
            ]
            if not matches:
                raise ValueError(
                    "checkpoint current stage is absent from this schedule"
                )
            self.current_stage = ZZ(matches[0])
        self.vertices = set(
            EffectiveVertex.from_record(record)
            for record in checkpoint.get("vertices", ())
        )
        self.vertices.update(self.roots)
        self.compiled_level_stages = set(
            self._level_stage_key(
                entry["genus"],
                entry["ambient_degree"],
                ProbeExpansionStage.from_record(entry["stage"]),
            )
            for entry in checkpoint.get("compiled_level_stages", ())
        )
        self.failed_level_stages = {
            self._level_stage_key(
                entry["genus"],
                entry["ambient_degree"],
                ProbeExpansionStage.from_record(entry["stage"]),
            ): entry["reason"]
            for entry in checkpoint.get("failed_level_stages", ())
        }
        self.block_history = list(checkpoint.get("block_history", ()))
        self._assignments = {}
        self.solver = InfinityVertexDP(
            self._equation_for,
            coefficient_field=self.coefficient_field,
        )
        self.solver.restore_checkpoint(checkpoint["solver"])
        self.active_relation_keys = set()
        self._activate_incident_relations()
        self._last_frontier = None
        return self.report()

    def _config_extends(self, saved):
        """Allow checkpoint reuse under a richer marking/unit schedule."""
        if tuple(saved.t_powers) != tuple(self.config.t_powers):
            return False
        if self.config.max_markings < saved.max_markings:
            return False
        if saved.include_unit_relatives \
                and not self.config.include_unit_relatives:
            return False
        saved_genus = (
            max(root.genus for root in self.roots)
            if saved.max_genus is None else saved.max_genus
        )
        saved_degree = (
            max(root.ambient_degree for root in self.roots)
            if saved.max_ambient_degree is None
            else saved.max_ambient_degree
        )
        if self._max_genus < saved_genus \
                or self._max_ambient_degree < saved_degree:
            return False
        if self.config.max_vertices < saved.max_vertices \
                or self.config.max_relations < saved.max_relations:
            return False
        if self.config.require_complete != saved.require_complete:
            return False
        return True

    def save_checkpoint(self, path):
        with open(path, "w") as stream:
            json.dump(self.checkpoint(), stream, indent=2, sort_keys=True)

    def load_checkpoint(self, path):
        with open(path) as stream:
            checkpoint = json.load(stream)
        return self.restore_checkpoint(checkpoint)


def orchestrate_plane_cubic_vertices(roots, initial_values=None, config=None,
                                     provider=None):
    """Convenience entry point returning the orchestrator and final report."""
    provider = provider or PlaneCubicFullEquationProvider()
    orchestrator = InfinityVertexOrchestrator(
        provider,
        roots,
        initial_values=initial_values,
        config=config,
    )
    return orchestrator, orchestrator.run()


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, required=True)
    parser.add_argument("--ambient-degree", type=int, default=0)
    parser.add_argument("--contacts", type=int, nargs="+", required=True)
    parser.add_argument("--psi-min", type=int, default=0)
    parser.add_argument("--insertions", type=int, nargs="+", required=True)
    parser.add_argument("--max-markings", type=int, default=2)
    parser.add_argument(
        "--t-powers", type=int, nargs="+", default=(2, 1, 0, -1, -2)
    )
    parser.add_argument("--laurent-precision", type=int, default=8)
    parser.add_argument("--no-unit-relatives", action="store_true")
    parser.add_argument("--checkpoint-in")
    parser.add_argument("--checkpoint-out")
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()

    if len(arguments.insertions) != len(arguments.contacts):
        parser.error(
            "provide one ambient insertion power for every contact"
        )

    target = EffectiveVertex(
        arguments.genus,
        arguments.ambient_degree,
        tuple(arguments.contacts),
        psi_min=arguments.psi_min,
        insertions=tuple(arguments.insertions),
    )
    if not target.is_balanced():
        parser.error("the requested plane-cubic infinity vertex is not balanced")
    provider = PlaneCubicFullEquationProvider(
        laurent_precision=arguments.laurent_precision
    )
    config = InfinityOrchestrationConfig(
        max_markings=arguments.max_markings,
        t_powers=tuple(arguments.t_powers),
        include_unit_relatives=not arguments.no_unit_relatives,
        max_genus=arguments.genus,
        max_ambient_degree=arguments.ambient_degree,
    )
    orchestrator = InfinityVertexOrchestrator(
        provider, (target,), config=config
    )
    if arguments.checkpoint_in:
        orchestrator.load_checkpoint(arguments.checkpoint_in)
    report = orchestrator.run()
    if arguments.checkpoint_out:
        orchestrator.save_checkpoint(arguments.checkpoint_out)

    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print("Infinity-vertex orchestration status:", report["status"])
    for root in report["roots"]:
        print("target:", root["vertex"])
        print("value:", root["value"])
    print("vertices:", report["vertex_count"])
    print("solved vertices:", report["solved_vertex_count"])
    print("compiled relations:", report["compiled_relation_count"])
    print("active relations:", report["active_relation_count"])
    if report["status"] != "complete":
        frontier = report["frontier"]
        print("rank-deficient blocks:", len(frontier["rank_deficiencies"]))
        print("nonlinear relations:", frontier["nonlinear_relation_count"])
        print("vertices without linear incidence:",
              len(frontier["vertices_without_linear_incidence"]))
    if arguments.checkpoint_out:
        print("checkpoint:", arguments.checkpoint_out)


if __name__ == "__main__" and os.path.basename(sys.argv[0]) in (
        "log_glsm_infinity_orchestrator.sage",
        "log_glsm_infinity_orchestrator.sage.py"):
    _main()
