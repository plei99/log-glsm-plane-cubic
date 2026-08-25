r"""Finite-closure orchestration for contact-resolved infinity vertices.

The low-level DP engine solves an equation after its probe has been chosen.
This module supplies the missing controller:

* generate progressively richer convention-compatible probes;
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
import time

load("cjr_full_equation_provider.sage")
load("cjr_punctured_axioms.sage")


class OrchestrationLimitError(RuntimeError):
    """A caller-supplied finite resource bound was exceeded."""


class InconsistentLocalizationRelationError(ArithmeticError):
    """A fully evaluated localization/axiom row has nonzero residual."""


class ProbeExpansionStage(SageObject):
    """One progressively richer family of dimension-compatible probes."""

    def __init__(self, max_markings, include_unit_relatives=False,
                 max_unit_insertions=None):
        self.max_markings = ZZ(max_markings)
        self.include_unit_relatives = bool(include_unit_relatives)
        if max_unit_insertions is None:
            max_unit_insertions = 1 if self.include_unit_relatives else 0
        self.max_unit_insertions = ZZ(max_unit_insertions)
        if self.max_markings < 1:
            raise ValueError("a probe stage needs at least one marking")
        if self.max_unit_insertions < 0:
            raise ValueError("max_unit_insertions must be nonnegative")
        if self.include_unit_relatives != bool(self.max_unit_insertions):
            raise ValueError(
                "unit stages need a positive max_unit_insertions depth"
            )

    def to_record(self):
        return {
            "max_markings": int(self.max_markings),
            "include_unit_relatives": self.include_unit_relatives,
            "max_unit_insertions": int(self.max_unit_insertions),
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["max_markings"],
            record.get("include_unit_relatives", False),
            record.get("max_unit_insertions"),
        )

    def __repr__(self):
        suffix = (
            "+unit<=%s" % self.max_unit_insertions
            if self.include_unit_relatives else ""
        )
        return "ProbeExpansionStage(markings<=%s%s)" % (
            self.max_markings, suffix
        )


class InfinityOrchestrationConfig(SageObject):
    """Finite bounds and probe schedule for an orchestration run."""

    def __init__(self, max_markings=2, t_powers=(0, -1, -2),
                 include_unit_relatives=True, max_unit_insertions=None,
                 include_chow_relations=False,
                 max_chow_unit_insertions=1,
                 chow_primary_only=False,
                 effective_basis_only=False,
                 additional_probe_ambient_degrees=(),
                 include_punctured_axioms=True,
                 max_genus=None,
                 max_ambient_degree=None, max_vertices=2000,
                 max_relations=4000, max_solve_rounds=1000,
                 require_complete=True, fail_on_unsupported=False,
                 include_kernel_basis=True):
        self.max_markings = ZZ(max_markings)
        self.t_powers = tuple(ZZ(power) for power in t_powers)
        self.include_unit_relatives = bool(include_unit_relatives)
        if max_unit_insertions is None:
            max_unit_insertions = 1 if self.include_unit_relatives else 0
        self.max_unit_insertions = ZZ(max_unit_insertions)
        self.include_chow_relations = bool(include_chow_relations)
        self.max_chow_unit_insertions = ZZ(max_chow_unit_insertions)
        self.chow_primary_only = bool(chow_primary_only)
        self.effective_basis_only = bool(effective_basis_only)
        self.additional_probe_ambient_degrees = tuple(sorted(set(
            ZZ(degree) for degree in additional_probe_ambient_degrees
        )))
        self.include_punctured_axioms = bool(include_punctured_axioms)
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
        self.include_kernel_basis = bool(include_kernel_basis)

        if self.max_markings < 1:
            raise ValueError("max_markings must be positive")
        if self.max_unit_insertions < 0:
            raise ValueError("max_unit_insertions must be nonnegative")
        if self.max_chow_unit_insertions < 0:
            raise ValueError(
                "max_chow_unit_insertions must be nonnegative"
            )
        if any(degree < 0 for degree in self.additional_probe_ambient_degrees):
            raise ValueError("additional probe degrees must be nonnegative")
        if self.include_unit_relatives != bool(self.max_unit_insertions):
            raise ValueError(
                "include_unit_relatives and max_unit_insertions disagree"
            )
        if self.effective_basis_only and (
                not self.include_chow_relations
                or not self.chow_primary_only
                or self.include_unit_relatives):
            raise ValueError(
                "effective-basis-only reconstruction requires primary Chow "
                "relations and excludes descendant unit relatives"
            )
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
            for unit_depth in range(1, self.max_unit_insertions + 1):
                stages.extend(
                    ProbeExpansionStage(markings, True, unit_depth)
                    for markings in range(1, self.max_markings + 1)
                )
        return tuple(stages)

    def to_record(self):
        return {
            "max_markings": int(self.max_markings),
            "t_powers": [int(power) for power in self.t_powers],
            "include_unit_relatives": self.include_unit_relatives,
            "max_unit_insertions": int(self.max_unit_insertions),
            "include_chow_relations": self.include_chow_relations,
            "max_chow_unit_insertions": int(
                self.max_chow_unit_insertions
            ),
            "chow_primary_only": self.chow_primary_only,
            "effective_basis_only": self.effective_basis_only,
            "additional_probe_ambient_degrees": [
                int(degree)
                for degree in self.additional_probe_ambient_degrees
            ],
            "include_punctured_axioms": self.include_punctured_axioms,
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
            "include_kernel_basis": self.include_kernel_basis,
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
    # Version 9 adds the CJR-I stabilization-boundary comparison and hence
    # log-domain stationary rows to effective-basis stages.  A version-8
    # checkpoint with the same configuration is missing those relations and
    # cannot be treated as a completed stage.
    CHECKPOINT_VERSION = 9

    def __init__(self, provider, roots, initial_values=None, config=None,
                 coefficient_field=None, candidate_provider=None,
                 axiom_provider=None):
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
        self.axiom_provider = axiom_provider
        if self.config.include_punctured_axioms \
                and self.axiom_provider is None:
            self.axiom_provider = PlaneCubicPuncturedAxioms(
                coefficient_field
            )
        if not self.config.include_punctured_axioms:
            self.axiom_provider = None

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
        self.compiled_axiom_vertices = set()
        self.failed_level_stages = {}
        self.current_stage = -1
        self.block_history = []
        self._last_frontier = None
        self._timed_out = False

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
        if self.config.additional_probe_ambient_degrees:
            self._max_ambient_degree = max(
                self._max_ambient_degree,
                max(self.config.additional_probe_ambient_degrees),
            )

    def _plane_cubic_candidates(self, genus, ambient_degree, stage, config):
        probe_degrees = tuple(sorted(set(
            (ZZ(ambient_degree),)
            + config.additional_probe_ambient_degrees
        )))
        relations = []
        for probe_degree in probe_degrees:
            if probe_degree < ambient_degree:
                continue
            keyword_arguments = {
                "max_markings": stage.max_markings,
                "t_powers": config.t_powers,
                "require_complete": config.require_complete,
                "include_unit_relatives": stage.include_unit_relatives,
                "max_unit_insertions": stage.max_unit_insertions,
                "include_chow_relations": (
                    config.include_chow_relations
                    and not stage.include_unit_relatives
                ),
                "max_chow_unit_insertions": (
                    config.max_chow_unit_insertions
                ),
                "chow_primary_only": config.chow_primary_only,
            }
            if config.effective_basis_only:
                keyword_arguments["effective_basis_only"] = True
            relations.extend(self.provider.candidate_relations(
                genus, probe_degree, **keyword_arguments
            ))
        return tuple(relations)

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
            int(stage.max_unit_insertions),
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

    def _compile_axiom_relations(self):
        """Add CJR unit-removal equations and closed base values once."""
        if self.axiom_provider is None:
            return False
        changed = False
        pending = sorted(
            self.vertices.difference(self.compiled_axiom_vertices),
            key=lambda item: item.order_key(),
        )
        for vertex in pending:
            known = self.axiom_provider.known_value(vertex)
            if known is not None:
                known = self.coefficient_field(known)
                if vertex in self.solver.values \
                        and self.solver.values[vertex] != known:
                    raise ArithmeticError(
                        "a solved infinity value contradicts the CJR closed "
                        "formula: %r has %s, expected %s"
                        % (vertex, self.solver.values[vertex], known)
                    )
                if vertex not in self.solver.values:
                    self.solver.values[vertex] = known
                    changed = True

            for relation in self.axiom_provider.relations_for(vertex):
                relation_key = self._relation_key(relation)
                if relation_key in self.relations:
                    continue
                if len(self.relations) >= self.config.max_relations:
                    raise OrchestrationLimitError(
                        "relation bound %s exceeded" % self.config.max_relations
                    )
                self.relations[relation_key] = relation
                changed = True
            self.compiled_axiom_vertices.add(vertex)
        return changed

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
                new_vertices = factors.difference(self.vertices)
                if len(self.vertices) + len(new_vertices) \
                        > self.config.max_vertices:
                    raise OrchestrationLimitError(
                        "vertex bound %s exceeded" % self.config.max_vertices
                    )
                # Keep activation atomic with respect to the safety bound so
                # a checkpoint taken after the exception can be resumed.
                self.active_relation_keys.add(key)
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
            changed = self._compile_axiom_relations()
            changed = self._activate_incident_relations() or changed
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
                vertex.contact_psi,
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

    def _relation_right_hand_side(self, relation):
        """Evaluate the constant side of a linear-ready relation."""
        off_diagonal = self.coefficient_field.zero()
        for coefficient, factors in relation.terms:
            unsolved = tuple(
                factor for factor in factors
                if factor not in self.solver.values
            )
            if unsolved:
                continue
            value = coefficient
            for factor in factors:
                value *= self.solver.values[factor]
            off_diagonal += value
        return (
            relation.known_gw
            - relation.twisted_zero_level
            - off_diagonal
        )

    def _check_solved_relations(self):
        """Reject a fully substituted row whose exact residual is nonzero."""
        for key in sorted(self.active_relation_keys):
            relation = self.relations[key]
            if any(
                    factor not in self.solver.values
                    for coefficient, factors in relation.terms
                    for factor in factors):
                continue
            residual = self._relation_right_hand_side(relation)
            if residual:
                raise InconsistentLocalizationRelationError(
                    "fully solved relation %r [t^%s] has residual %s"
                    % (relation.probe, relation.t_power, residual)
                )

    def _identifiable_values(self, targets, selection):
        r"""Recover coordinates fixed by a rank-deficient affine system.

        A component need not be entirely invertible.  In reduced row-echelon
        form, a pivot coordinate is nevertheless unique when its row has no
        free-variable entries.  This detects exactly the coordinates that
        vanish on the right kernel, without materialising a kernel basis.
        """
        if not selection.rows or selection.is_full_rank:
            return {}
        right_hand_side = matrix(
            self.coefficient_field, len(selection.relations), 1,
            [self._relation_right_hand_side(relation)
             for relation in selection.relations],
        )
        reduced = selection.matrix.augment(right_hand_side).rref()
        pivots = selection.matrix.pivots()
        free_columns = tuple(
            column for column in range(len(targets))
            if column not in set(pivots)
        )
        values = {}
        for row, pivot in enumerate(pivots):
            if all(not reduced[row, column] for column in free_columns):
                values[targets[pivot]] = self.coefficient_field(
                    reduced[row, len(targets)]
                )
        return values

    def solve_next_block(self):
        """Solve one highest-priority full-rank linear-ready block."""
        self._check_solved_relations()
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
            partial_candidates = []
            for component, relations, targets, selection in deficiencies:
                values = self._identifiable_values(targets, selection)
                if values:
                    partial_candidates.append((component, selection, values))
            if not partial_candidates:
                return None
            partial_candidates.sort(
                key=lambda item: self._component_priority(item[0]),
                reverse=True,
            )
            component, selection, values = partial_candidates[0]
            self.solver.values.update(values)
            self._check_solved_relations()
            self._last_frontier = None
            self.block_history.append({
                "partial": True,
                "component_size": len(component),
                "rank": int(selection.rank),
                "targets": [target.to_record() for target in sorted(
                    values, key=lambda item: item.order_key()
                )],
                "values": [str(values[target]) for target in sorted(
                    values, key=lambda item: item.order_key()
                )],
            })
            return values
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
        self._check_solved_relations()
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

    def solve_until_stalled(self, deadline=None):
        """Repeatedly solve every block made ready by prior substitutions.

        ``deadline`` is an absolute ``time.time()`` value.  It is checked
        between blocks, so a run overshoots by at most one exact block solve
        rather than by a whole stage; a single genus-three stage has been
        observed to occupy a machine for hours.
        """
        solved_blocks = 0
        for _ in range(self.config.max_solve_rounds):
            if deadline is not None and time.time() > deadline:
                self._timed_out = True
                return solved_blocks
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

    def run(self, time_budget=None):
        """Expand probe stages until the roots solve or the schedule stalls.

        ``time_budget`` is a wall-clock bound in seconds.  When it expires the
        run stops cleanly after the current block, checkpointable state is
        intact, and the report carries ``timed_out=True``.
        """
        deadline = None if time_budget is None \
            else time.time() + float(time_budget)
        self._timed_out = False
        if self.current_stage < 0 and self.axiom_provider is not None:
            changed = True
            while changed:
                changed = self._compile_axiom_relations()
                changed = self._activate_incident_relations() or changed
            self.solve_until_stalled(deadline=deadline)
        while not self._timed_out:
            if self.roots_solved():
                break
            if deadline is not None and time.time() > deadline:
                self._timed_out = True
                break
            if self.current_stage < 0:
                if not self.advance_stage():
                    break
            self.solve_until_stalled(deadline=deadline)
            if self.roots_solved():
                break
            if self._timed_out:
                break
            if not self.advance_stage():
                break
        return self.report()

    def _frontier_record(self, deficiencies, nonlinear):
        deficiency_records = []
        ready_vertices = set()
        for component, relations, targets, selection in deficiencies:
            ready_vertices.update(component)
            record = {
                "targets": [target.to_record() for target in targets],
                "rank": int(selection.rank),
                "columns": len(targets),
                "kernel_dimension": len(targets) - int(selection.rank),
                "kernel": None,
                "candidate_rows": len(relations),
            }
            if self.config.include_kernel_basis:
                record["kernel"] = [
                    [str(value) for value in vector]
                    for vector in selection.kernel_basis()
                ]
            deficiency_records.append(record)
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
            "timed_out": bool(getattr(self, "_timed_out", False)),
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
                    "max_unit_insertions": key[4],
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
                        key[2], key[3], key[4]
                    ).to_record(),
                }
                for key in sorted(self.compiled_level_stages)
            ],
            "failed_level_stages": [
                {
                    "genus": key[0],
                    "ambient_degree": key[1],
                    "stage": ProbeExpansionStage(
                        key[2], key[3], key[4]
                    ).to_record(),
                    "reason": reason,
                }
                for key, reason in sorted(self.failed_level_stages.items())
            ],
            "solver": self.solver.checkpoint(),
            "block_history": list(self.block_history),
        }

    LEGACY_IMPORT_VERSIONS = (8,)

    def import_legacy_relations(self, path):
        r"""Load the still-valid relations of an older checkpoint as rows.

        A version-8 checkpoint is rejected by ``restore_checkpoint`` because
        its *stage bookkeeping* cannot be trusted: version 9 added the
        stabilization-boundary comparison, so a v8 "completed" stage is
        missing rows.  The individual relation records, however, are exact
        statements independent of that bookkeeping, and re-deriving them
        repeats hours of graph compilation.

        This method imports only the subset that is still known to be sound:

        * records whose ``(probe, t_power)`` passes the Chow scalar gate
          (residual dimension ``<= 0``) -- positive-dimensional coefficients
          were the historical genus-one contradiction and are skipped;
        * records not already present under the current relation key; and
        * no solver values, no vertex closure, and no stage completions.

        Axiom rows (string/divisor/dilaton pseudo-probes) are also skipped:
        their pseudo-probe dimension data does not describe a localization
        extraction, and the version-9 axiom provider regenerates them
        exactly.  Returns the number of imported relations.
        """
        with open(path) as stream:
            checkpoint = json.load(stream)
        if checkpoint.get("format") != self.CHECKPOINT_FORMAT:
            raise ValueError("unrecognized checkpoint format")
        version = checkpoint.get("version")
        if version not in self.LEGACY_IMPORT_VERSIONS:
            raise ValueError(
                "legacy import supports versions %s; use restore_checkpoint "
                "for version %s"
                % (self.LEGACY_IMPORT_VERSIONS, self.CHECKPOINT_VERSION)
            )
        if checkpoint.get("coefficient_field") != str(self.coefficient_field):
            raise ValueError("legacy checkpoint coefficient field mismatch")

        chow = getattr(self.provider, "chow_backend", None) \
            or CJRInfinityChowBackend()
        imported = 0
        for record in checkpoint.get("relations", ()):
            relation = CompiledLocalizationRelation.from_record(
                record, self.coefficient_field
            )
            if not relation.is_complete:
                continue
            if " remove=" in relation.probe.label:
                # Axiom pseudo-probes (string/divisor/dilaton rows carry
                # labels like "CJR 8.3 string remove=1 vertex=...").  Their
                # pseudo-probe dimension data is not a localization
                # extraction, and version 9 regenerates them exactly.  Rows
                # labeled "CJR Chow probe" are genuine localization rows and
                # must be kept.
                continue
            degree = chow.degree(relation.probe, relation.t_power)
            if degree.is_class_valued:
                continue
            key = self._relation_key(relation)
            if key in self.relations:
                continue
            if len(self.relations) >= self.config.max_relations:
                raise OrchestrationLimitError(
                    "relation bound %s exceeded during legacy import"
                    % self.config.max_relations
                )
            self.relations[key] = relation
            self._cache_restored_relation_in_provider(relation)
            imported += 1
        if imported:
            self._activate_incident_relations()
            self._last_frontier = None
        return imported

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
        chow_schedule_richer = (
            self.config.include_chow_relations
            and (
                not saved_config.include_chow_relations
                or self.config.max_chow_unit_insertions
                > saved_config.max_chow_unit_insertions
                or (
                    saved_config.chow_primary_only
                    and not self.config.chow_primary_only
                )
            )
        )
        probe_schedule_richer = not set(
            self.config.additional_probe_ambient_degrees
        ).issubset(set(saved_config.additional_probe_ambient_degrees))
        current_stage_spec = checkpoint.get("current_stage_spec")
        if current_stage_spec is None:
            self.current_stage = ZZ(-1)
        else:
            saved_stage = ProbeExpansionStage.from_record(current_stage_spec)
            matches = [
                index for index, stage in enumerate(self.config.stages())
                if stage.to_record() == saved_stage.to_record()
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
        self.compiled_axiom_vertices = set()
        self.failed_level_stages = {
            self._level_stage_key(
                entry["genus"],
                entry["ambient_degree"],
                ProbeExpansionStage.from_record(entry["stage"]),
            ): entry["reason"]
            for entry in checkpoint.get("failed_level_stages", ())
        }
        if chow_schedule_richer or probe_schedule_richer:
            # Revisit the finite stage schedule when extra probe degrees are
            # added; restored rows are cached in the provider, so old probes
            # do not trigger fixed-locus recompilation.
            self.current_stage = ZZ(-1)
            self.compiled_level_stages = set()
            self.failed_level_stages = {}
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
        if self.config.effective_basis_only != saved.effective_basis_only:
            return False
        if self.config.max_markings < saved.max_markings:
            return False
        if saved.include_unit_relatives \
                and not self.config.include_unit_relatives:
            return False
        if saved.include_punctured_axioms \
                and not self.config.include_punctured_axioms:
            return False
        if self.config.max_unit_insertions < saved.max_unit_insertions:
            return False
        if saved.include_chow_relations \
                and not self.config.include_chow_relations:
            return False
        if self.config.max_chow_unit_insertions \
                < saved.max_chow_unit_insertions:
            return False
        if not set(saved.additional_probe_ambient_degrees).issubset(set(
                self.config.additional_probe_ambient_degrees)):
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
    parser.add_argument(
        "--contact-psi", type=int, nargs="+",
        help=("ordinary cotangent power at each contact; defaults to zero "
              "at every contact"),
    )
    parser.add_argument("--max-markings", type=int, default=2)
    parser.add_argument(
        "--t-powers", type=int, nargs="+", default=(0, -1, -2)
    )
    parser.add_argument("--laurent-precision", type=int, default=8)
    parser.add_argument("--no-unit-relatives", action="store_true")
    parser.add_argument(
        "--chow-relations", action="store_true",
        help=("add under-dimensioned probes, using only Laurent "
              "coefficients of nonpositive residual Chow dimension"),
    )
    parser.add_argument(
        "--max-chow-unit-insertions", type=int, default=1,
        help="largest number of unit markings in a Chow probe",
    )
    parser.add_argument(
        "--primary-chow-relations-only", action="store_true",
        help="use primary under-dimensioned Chow probes only",
    )
    parser.add_argument(
        "--additional-probe-ambient-degrees", type=int, nargs="*", default=(),
        help=("extra ambient probe degrees; nonmultiples of three have zero "
              "compact side (the genus-one assembly is regression-tested "
              "through ambient degree three)"),
    )
    parser.add_argument(
        "--no-punctured-axioms", action="store_true",
        help="disable the CJR contact-minus-one string/divisor/dilaton rows",
    )
    parser.add_argument(
        "--max-unit-insertions", type=int, default=1,
        help="largest number of reducible unit insertions in a mixed probe",
    )
    parser.add_argument("--checkpoint-in")
    parser.add_argument("--checkpoint-out")
    parser.add_argument(
        "--time-budget", type=float, default=None,
        help=("wall-clock seconds after which the run stops cleanly between "
              "blocks; combine with --checkpoint-out to resume later"),
    )
    parser.add_argument(
        "--import-v8-relations", metavar="PATH",
        help=("also load the Chow-scalar relations of a version-8 checkpoint "
              "as candidate rows, discarding its stage bookkeeping"),
    )
    parser.add_argument(
        "--zero-vertex-cache",
        help="persistent SQLite (recommended) or JSON O(3) zero-vertex cache",
    )
    parser.add_argument(
        "--hodge-cache",
        help="shared SQLite cache for exact Hodge integrals",
    )
    parser.add_argument(
        "--zero-vertex-strategy",
        choices=("localization", "hybrid", "givental"),
        default="localization",
        help=("Givental uses the calibrated quantum connection; hybrid keeps "
              "localization as an unsupported-geometry fallback"),
    )
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()

    if len(arguments.insertions) != len(arguments.contacts):
        parser.error(
            "provide one ambient insertion power for every contact"
        )
    if arguments.contact_psi is not None \
            and len(arguments.contact_psi) != len(arguments.contacts):
        parser.error(
            "provide one contact-psi power for every contact"
        )

    target = EffectiveVertex(
        arguments.genus,
        arguments.ambient_degree,
        tuple(arguments.contacts),
        psi_min=arguments.psi_min,
        insertions=tuple(arguments.insertions),
        contact_psi=(
            tuple(arguments.contact_psi)
            if arguments.contact_psi is not None else tuple()
        ),
    )
    if not target.is_balanced():
        parser.error("the requested plane-cubic infinity vertex is not balanced")
    if not target.is_dimension_zero():
        parser.error(
            "the requested infinity invariant violates reduced virtual dimension; "
            "expected psi_min + sum(insertions) + sum(contact_psi) = %s"
            % target.reduced_virtual_dimension
        )
    provider = PlaneCubicFullEquationProvider(
        laurent_precision=arguments.laurent_precision,
        zero_vertex_cache=arguments.zero_vertex_cache,
        autosave_zero_vertices=arguments.zero_vertex_cache is not None,
        hodge_cache=arguments.hodge_cache,
        zero_vertex_strategy=arguments.zero_vertex_strategy,
    )
    config = InfinityOrchestrationConfig(
        max_markings=arguments.max_markings,
        t_powers=tuple(arguments.t_powers),
        include_unit_relatives=not arguments.no_unit_relatives,
        include_chow_relations=arguments.chow_relations,
        max_chow_unit_insertions=arguments.max_chow_unit_insertions,
        chow_primary_only=arguments.primary_chow_relations_only,
        additional_probe_ambient_degrees=(
            arguments.additional_probe_ambient_degrees
        ),
        include_punctured_axioms=not arguments.no_punctured_axioms,
        max_unit_insertions=(
            0 if arguments.no_unit_relatives
            else arguments.max_unit_insertions
        ),
        max_genus=arguments.genus,
        max_ambient_degree=arguments.ambient_degree,
    )
    orchestrator = InfinityVertexOrchestrator(
        provider, (target,), config=config
    )
    if arguments.checkpoint_in:
        orchestrator.load_checkpoint(arguments.checkpoint_in)
    if arguments.import_v8_relations:
        imported = orchestrator.import_legacy_relations(
            arguments.import_v8_relations
        )
        print("imported %s legacy relation(s)" % imported)
    report = orchestrator.run(time_budget=arguments.time_budget)
    if arguments.checkpoint_out:
        orchestrator.save_checkpoint(arguments.checkpoint_out)

    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print("Infinity-vertex orchestration status:", report["status"])
    if report.get("timed_out"):
        print("stopped early: wall-clock budget reached")
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
