r"""Compute every genus-three log-GLSM infinity vertex in a finite truncation.

Despite the module name, the driver is genus-generic: ``--genus`` (default 3)
selects the family, and the enumeration, probe rounding, and orchestration
below make no genus-three assumption.  The `(4,1,0)` graph family already
enumerates in minutes, so genus four runs are configuration, not new code.

Genus three is the first genus in which the CJR hypersurface balancing
equation

    3*D - (2g-2) = sum_i (c_i + 1),        c_i <= -1,

admits a *positive* infinity degree.  Writing ``k_i=-c_i>=1`` the equation
becomes ``sum_i(k_i-1) = 2g-2-3D``, so a vertex exists exactly when

    0 <= D <= floor((2g-2)/3).

For ``g<=2`` this forces ``D=0``, which is the simplification used throughout
the genus-two calculation.  For ``g=3`` it allows ``D=0`` (contact excess 4)
and ``D=1`` (contact excess 1).  The ``D=1`` sector is ROADMAP risk 4 and has
no genus-two analogue.

"All genus-three vertices" is an infinite family without a truncation: a
contact ``c=-1`` contributes ``c+1=0`` to the balance, so arbitrarily many
``-1`` legs may be appended.  This module therefore enumerates the finite
family cut out by

* every balance-allowed infinity degree ``D``;
* contact valence at most ``max_valence``;
* one evaluation insertion ``alpha in {0,1,2}`` per contact; and
* the unique nonnegative ``psi_min`` power required by reduced virtual
  dimension, optionally capped by ``max_psi_min``.

The last two bounds are dimension-theoretic rather than arbitrary.  The root
enumeration chooses no ordinary contact descendants; stabilization-corrected
probe equations may add such decorated vertices to the dependency closure,
where they are retained and dimension-pruned exactly.

The enumerated vertices are handed to a single ``InfinityVertexOrchestrator``
as a common root set, so compiled relations and solved lower vertices are
shared across all of them.  Values are reported only when an exact full-rank
block determined them; otherwise the certified frontier is reported.  Nothing
is aggregated and no singular block is resolved by a tie-break.

Ordinary stationary GW descendants are localized through the stabilization
pullback.  On a contracted marked zero tail their cotangent powers become
ordinary descendants at the corresponding infinity contact.  These enlarged
vertices are kept distinct from CJR's log-domain descendants.

Run from the repository root:

    sage genus_three_infinity_vertices.sage --list-only
    sage genus_three_infinity_vertices.sage --max-valence 2 \
        --checkpoint-out genus3-infinity.json

or use the API after

    load("genus_three_infinity_vertices.sage")
"""

import argparse
import json
import os
import sys
import time
from itertools import product

load("log_glsm_infinity_orchestrator.sage")


HYPERSURFACE_DEGREE = ZZ(3)


def allowed_infinity_degrees(genus, hypersurface_degree=HYPERSURFACE_DEGREE):
    r"""Return every ambient degree an infinity vertex of ``genus`` may carry.

    Each contact satisfies ``c_i+1 <= 0``, so the balance equation forces
    ``hypersurface_degree*D <= 2g-2``.  Genus three is the first genus for
    which this permits ``D>0``.
    """
    genus = ZZ(genus)
    hypersurface_degree = ZZ(hypersurface_degree)
    if genus < 0:
        raise ValueError("genus must be nonnegative")
    if hypersurface_degree < 1:
        raise ValueError("the hypersurface degree must be positive")
    excess = 2 * genus - 2
    if excess < 0:
        return tuple()
    return tuple(range(excess // hypersurface_degree + 1))


def infinity_contact_profiles_in_degree(genus, valence, ambient_degree=0,
                                        hypersurface_degree=HYPERSURFACE_DEGREE):
    r"""List contact profiles allowed at fixed genus, valence, and degree.

    This is the positive-degree generalization of
    ``infinity_contact_profiles`` in ``log_glsm_infinity_vertices.sage``; at
    ``ambient_degree=0`` the two agree.  Entries are the negative CJR contact
    orders sorted increasingly.
    """
    genus = ZZ(genus)
    valence = ZZ(valence)
    ambient_degree = ZZ(ambient_degree)
    hypersurface_degree = ZZ(hypersurface_degree)
    if valence <= 0:
        raise ValueError("valence must be positive")
    if ambient_degree < 0:
        raise ValueError("ambient_degree must be nonnegative")

    excess = 2 * genus - 2 - hypersurface_degree * ambient_degree
    if excess < 0:
        return tuple()

    profiles = set()
    for allocation in IntegerVectors(excess, valence):
        orders = sorted((ZZ(value) + 1 for value in allocation), reverse=True)
        profiles.add(tuple(-order for order in orders))
    return tuple(sorted(profiles))


def enumerate_infinity_vertices(genus=3, max_valence=3,
                                max_infinity_degree=None,
                                insertion_powers=(0, 1, 2),
                                max_psi_min=None,
                                hypersurface_degree=HYPERSURFACE_DEGREE):
    r"""Return every balanced vertex of ``genus`` inside a finite truncation.

    A vertex is stored in the canonical form used by the graph compiler: the
    ``(contact, insertion)`` pairs are sorted, so ``contacts`` is increasing
    and ``insertions`` is aligned with it.  For each assignment, ``psi_min``
    is the unique power satisfying

        psi_min + sum(insertions) = valence + 3*ambient_degree.

    ``max_psi_min`` may discard that power but never substitutes a lower one.
    These are primary-contact roots (`contact_psi=0`); ordinary contact
    descendants introduced by stabilization enter through relation closure.
    """
    genus = ZZ(genus)
    max_valence = ZZ(max_valence)
    if max_valence < 1:
        raise ValueError("max_valence must be positive")
    insertion_powers = tuple(ZZ(power) for power in insertion_powers)
    if any(power < 0 or power > 2 for power in insertion_powers):
        raise ValueError("ambient insertions on P^2 are H^0, H^1, or H^2")
    if max_psi_min is not None:
        max_psi_min = ZZ(max_psi_min)
        if max_psi_min < 0:
            raise ValueError("max_psi_min must be nonnegative")

    degrees = allowed_infinity_degrees(genus, hypersurface_degree)
    if max_infinity_degree is not None:
        max_infinity_degree = ZZ(max_infinity_degree)
        if max_infinity_degree < 0:
            raise ValueError("max_infinity_degree must be nonnegative")
        degrees = tuple(
            degree for degree in degrees if degree <= max_infinity_degree
        )

    seen = set()
    vertices = []
    for ambient_degree in degrees:
        for valence in range(1, max_valence + 1):
            profiles = infinity_contact_profiles_in_degree(
                genus, valence, ambient_degree, hypersurface_degree
            )
            for profile in profiles:
                for assignment in product(insertion_powers, repeat=valence):
                    pairs = tuple(sorted(zip(profile, assignment)))
                    contacts = tuple(contact for contact, _ in pairs)
                    insertions = tuple(power for _, power in pairs)
                    codimension = sum(insertions)
                    psi_min = PlaneCubicDimension.infinity_required_psi_min_power(
                        genus, ambient_degree, contacts, codimension
                    )
                    if psi_min < 0 or (
                            max_psi_min is not None
                            and psi_min > max_psi_min):
                        continue
                    vertex = EffectiveVertex(
                        genus, ambient_degree, contacts,
                        psi_min=psi_min, insertions=insertions,
                    )
                    if not vertex.is_balanced(hypersurface_degree):
                        raise AssertionError(
                            "enumerated vertex %r is unbalanced" % (vertex,)
                        )
                    if not vertex.is_dimension_zero():
                        raise AssertionError(
                            "enumerated vertex %r violates virtual dimension"
                            % (vertex,)
                        )
                    signature = vertex.signature()
                    if signature in seen:
                        continue
                    seen.add(signature)
                    vertices.append(vertex)
    return tuple(sorted(vertices, key=lambda item: item.order_key()))


def probe_ambient_degree(vertex_degree, hypersurface_degree=HYPERSURFACE_DEGREE):
    r"""Return the smallest ambient probe degree seeing ``vertex_degree``.

    The ambient log-GLSM localization problem exists in every P2 degree.  If
    that degree is not divisible by three, the hypersurface stable-map side is
    empty and its stabilized invariant is exactly zero.  Such a zero probe is
    still a useful localization equation, and it is smaller than rounding up
    to the next nonempty elliptic-curve degree.  A vertex of degree ``D`` can
    already occur in total degree ``D``, so the optimal probe degree is ``D``.

    ``hypersurface_degree`` remains in the signature for API compatibility and
    to validate the convention, but does not alter the answer.
    """
    vertex_degree = ZZ(vertex_degree)
    hypersurface_degree = ZZ(hypersurface_degree)
    if vertex_degree < 0:
        raise ValueError("an infinity degree must be nonnegative")
    if hypersurface_degree < 1:
        raise ValueError("the hypersurface degree must be positive")
    return vertex_degree


class GenusThreeVertexComputation(SageObject):
    r"""Drive one orchestration over the whole enumerated genus-three family.

    Sharing a single orchestrator across every root matters: the genus-three
    vertices depend on a common closure of genus-one and genus-two data, and
    on the same compiled probe relations.  Solving them one at a time would
    recompile that closure once per target.
    """

    def __init__(self, genus=3, max_valence=3, max_infinity_degree=None,
                 insertion_powers=(0, 1, 2), max_psi_min=None,
                 max_markings=2, t_powers=(0, -1, -2),
                 include_unit_relatives=True, max_unit_insertions=None,
                 include_chow_relations=False,
                 max_chow_unit_insertions=1,
                 chow_primary_only=False,
                 effective_basis_only=False,
                 additional_probe_ambient_degrees=(),
                 include_punctured_axioms=True,
                 laurent_precision=8,
                 provider=None, require_complete=True, progress=None,
                 initial_values=None, include_kernel_basis=False,
                 zero_vertex_cache=None, hodge_cache=None,
                 zero_vertex_strategy="localization",
                 max_vertices=20000, max_relations=20000):
        self.genus = ZZ(genus)
        self.roots = enumerate_infinity_vertices(
            genus=self.genus,
            max_valence=max_valence,
            max_infinity_degree=max_infinity_degree,
            insertion_powers=insertion_powers,
            max_psi_min=max_psi_min,
        )
        if not self.roots:
            raise ValueError(
                "no balanced infinity vertex exists in the requested truncation"
            )
        self.progress = progress
        self.provider = provider or PlaneCubicFullEquationProvider(
            laurent_precision=laurent_precision,
            zero_vertex_cache=zero_vertex_cache,
            autosave_zero_vertices=zero_vertex_cache is not None,
            hodge_cache=hodge_cache,
            zero_vertex_strategy=zero_vertex_strategy,
        )

        root_degree = max(root.ambient_degree for root in self.roots)
        self.max_probe_degree = max(
            (probe_ambient_degree(root_degree),)
            + tuple(ZZ(degree)
                    for degree in additional_probe_ambient_degrees)
        )
        self.config = InfinityOrchestrationConfig(
            max_markings=max_markings,
            t_powers=tuple(t_powers),
            include_unit_relatives=include_unit_relatives,
            max_unit_insertions=max_unit_insertions,
            include_chow_relations=include_chow_relations,
            max_chow_unit_insertions=max_chow_unit_insertions,
            chow_primary_only=chow_primary_only,
            effective_basis_only=effective_basis_only,
            additional_probe_ambient_degrees=(
                additional_probe_ambient_degrees
            ),
            include_punctured_axioms=include_punctured_axioms,
            max_genus=self.genus,
            max_ambient_degree=self.max_probe_degree,
            max_vertices=max_vertices,
            max_relations=max_relations,
            require_complete=require_complete,
            include_kernel_basis=include_kernel_basis,
        )
        self._candidate_cache = {}
        # Values proved at a lower genus enter as constants rather than as
        # further unknowns.  Every such substitution removes a column from the
        # blocks below, which is what makes a staged bottom-up run different
        # from asking for genus three directly.
        self.orchestrator = InfinityVertexOrchestrator(
            self.provider,
            self.roots,
            initial_values=dict(initial_values or {}),
            config=self.config,
            candidate_provider=self._candidates,
        )

    def solved_values(self):
        """Return every exact value known so far, for use as a lower layer."""
        return dict(self.orchestrator.solver.values)

    def restrict_future_probe_degree(self, ambient_degree):
        r"""Restrict only probe relations compiled after the current point.

        Saved relations and values remain available.  This permits a richer
        checkpoint to add a targeted marking stage without recompiling all of
        its positive ambient probe degrees.
        """
        ambient_degree = ZZ(ambient_degree)
        if ambient_degree < 0:
            raise ValueError("a future probe-degree ceiling is nonnegative")
        configured_degree = self.config.max_ambient_degree
        if configured_degree is not None \
                and ambient_degree > configured_degree:
            raise ValueError(
                "a future probe-degree ceiling cannot exceed the configured bound"
            )
        self.max_probe_degree = ambient_degree

    def _report_progress(self, message):
        if self.progress is not None:
            self.progress(message)

    def _candidates(self, genus, ambient_degree, stage, config):
        r"""Compile direct ambient-degree probes, including zero GW sectors."""
        genus = ZZ(genus)
        if genus < 1:
            return tuple()
        degree = probe_ambient_degree(ambient_degree)
        if degree > self.max_probe_degree:
            return tuple()
        key = (
            int(genus), int(degree), int(stage.max_markings),
            bool(stage.include_unit_relatives),
            int(stage.max_unit_insertions),
        )
        if key in self._candidate_cache:
            return self._candidate_cache[key]
        self._report_progress(
            "compiling probes: genus %s, ambient degree %s, markings<=%s%s"
            % (genus, degree, stage.max_markings,
               ", unit depth<=%s" % stage.max_unit_insertions
               if stage.include_unit_relatives else "")
        )
        started = time.time()
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
            "max_chow_unit_insertions": config.max_chow_unit_insertions,
            "chow_primary_only": config.chow_primary_only,
        }
        if config.effective_basis_only:
            keyword_arguments["effective_basis_only"] = True
        probe_degrees = tuple(sorted(set(
            (degree,) + config.additional_probe_ambient_degrees
        )))
        relations = []
        for probe_degree in probe_degrees:
            if probe_degree < degree:
                continue
            if probe_degree > self.max_probe_degree:
                continue
            relations.extend(self.provider.candidate_relations(
                genus, probe_degree, **keyword_arguments
            ))
        relations = tuple(relations)
        self._report_progress(
            "  %s relation(s) in %.1fs" % (len(relations), time.time() - started)
        )
        self._candidate_cache[key] = relations
        return relations

    def _solve_stage(self, deadline, checkpoint_path=None):
        r"""Solve ready blocks one at a time, checking the clock between them.

        ``InfinityVertexOrchestrator.solve_until_stalled`` runs to exhaustion,
        and a single genus-three stage can occupy a machine for hours inside
        it.  Driving ``solve_next_block`` directly bounds the overrun to one
        block instead of one stage.
        """
        solved = 0
        for _ in range(self.config.max_solve_rounds):
            if deadline is not None and time.time() > deadline:
                self._report_progress("  time budget reached mid-stage")
                return solved, True
            started = time.time()
            values = self.solve_next_block_with_report()
            if values is None:
                return solved, False
            solved += 1
            if checkpoint_path:
                self.orchestrator.save_checkpoint(checkpoint_path)
            self._report_progress(
                "  block %s solved in %.1fs (%s vertex value(s))"
                % (solved, time.time() - started, len(values))
            )
            if self.orchestrator.roots_solved():
                return solved, False
        raise OrchestrationLimitError(
            "solve-round bound %s exceeded" % self.config.max_solve_rounds
        )

    def solve_next_block_with_report(self):
        """Hook for subclasses and tests; delegates to the orchestrator."""
        return self.orchestrator.solve_next_block()

    def run(self, checkpoint_path=None, time_budget=None):
        r"""Advance probe stages until the roots solve, stall, or time out.

        This mirrors ``InfinityVertexOrchestrator.run`` but checkpoints after
        every stage and honours a wall-clock budget between individual solved
        blocks, because a genus-three probe compilation and its exact block
        algebra are both far more expensive than their genus-two analogues.
        """
        started = time.time()
        deadline = None if time_budget is None else started + time_budget
        exhausted = False
        timed_out = False

        # A checkpoint can be written after relation compilation but before
        # (or midway through) exact block solving.  Finish that saved stage
        # before advancing; otherwise a time-limited resume silently skips
        # equations that are already present in the checkpoint.
        if self.orchestrator.current_stage >= 0 \
                and not self.orchestrator.roots_solved():
            solved, timed_out = self._solve_stage(
                deadline, checkpoint_path=checkpoint_path
            )
            self._report_progress(
                "  resumed stage solved %s block(s); %s/%s roots known"
                % (solved, self.solved_root_count(), len(self.roots))
            )
            if checkpoint_path:
                self.orchestrator.save_checkpoint(checkpoint_path)
        while not self.orchestrator.roots_solved():
            if timed_out:
                break
            if deadline is not None and time.time() > deadline:
                self._report_progress("time budget reached; stopping early")
                timed_out = True
                break
            try:
                advanced = self.orchestrator.advance_stage()
            except OrchestrationLimitError:
                # Candidate relations may have taken hours to compile before
                # closure activation discovers that a safety cap is too low.
                # They are valid checkpoint data even though the stage did
                # not finish activating, provided activation itself is atomic.
                if checkpoint_path:
                    self.orchestrator.save_checkpoint(checkpoint_path)
                raise
            if not advanced:
                exhausted = True
                break
            self._report_progress(
                "stage %s/%s" % (
                    self.orchestrator.current_stage + 1,
                    len(self.config.stages()),
                )
            )
            # Relation compilation can dominate a genus-three stage.  Save it
            # before exact block solving so an inconsistency or interruption
            # never forces the fixed-locus work to be repeated.
            if checkpoint_path:
                self.orchestrator.save_checkpoint(checkpoint_path)
            solved, timed_out = self._solve_stage(
                deadline, checkpoint_path=checkpoint_path
            )
            self._report_progress(
                "  solved %s block(s); %s/%s roots known"
                % (solved, self.solved_root_count(), len(self.roots))
            )
            if checkpoint_path:
                self.orchestrator.save_checkpoint(checkpoint_path)
            if timed_out:
                break
        if checkpoint_path:
            self.orchestrator.save_checkpoint(checkpoint_path)
        return self.report(
            stages_exhausted=exhausted,
            elapsed=time.time() - started,
            timed_out=timed_out,
        )

    def load_checkpoint(self, path):
        return self.orchestrator.load_checkpoint(path)

    def solved_root_count(self):
        values = self.orchestrator.solver.values
        return sum(1 for root in self.roots if root in values)

    def root_records(self):
        values = self.orchestrator.solver.values
        return tuple(
            {
                "vertex": root.to_record(),
                "description": str(root),
                "solved": root in values,
                "value": str(values[root]) if root in values else None,
            }
            for root in self.roots
        )

    def report(self, stages_exhausted=False, elapsed=None, timed_out=False):
        records = self.root_records()
        solved = sum(1 for record in records if record["solved"])
        orchestration = self.orchestrator.report()
        return {
            "genus": int(self.genus),
            "timed_out": bool(timed_out),
            "truncation": {
                "infinity_degrees": [
                    int(degree)
                    for degree in sorted(set(
                        root.ambient_degree for root in self.roots
                    ))
                ],
                "max_valence": int(max(
                    len(root.contacts) for root in self.roots
                )),
                "max_probe_ambient_degree": int(self.max_probe_degree),
                "config": self.config.to_record(),
            },
            "target_count": len(records),
            "solved_count": int(solved),
            "status": "complete" if solved == len(records) else "blocked",
            "stages_exhausted": bool(stages_exhausted),
            "elapsed_seconds": None if elapsed is None else float(elapsed),
            "targets": list(records),
            "orchestration": orchestration,
        }


def enumeration_report(genus=3, max_valence=3, max_infinity_degree=None,
                       insertion_powers=(0, 1, 2), max_psi_min=None):
    """Return the truncated target family without computing any value."""
    vertices = enumerate_infinity_vertices(
        genus=genus,
        max_valence=max_valence,
        max_infinity_degree=max_infinity_degree,
        insertion_powers=insertion_powers,
        max_psi_min=max_psi_min,
    )
    by_degree = {}
    for vertex in vertices:
        by_degree.setdefault(int(vertex.ambient_degree), []).append(vertex)
    return {
        "genus": int(genus),
        "allowed_infinity_degrees": [
            int(degree) for degree in allowed_infinity_degrees(genus)
        ],
        "max_valence": int(max_valence),
        "target_count": len(vertices),
        "counts_by_degree": {
            str(degree): len(items) for degree, items in sorted(by_degree.items())
        },
        "profiles_by_degree": {
            str(degree): [
                [int(contact) for contact in profile]
                for valence in range(1, ZZ(max_valence) + 1)
                for profile in infinity_contact_profiles_in_degree(
                    genus, valence, degree
                )
            ]
            for degree in allowed_infinity_degrees(genus)
        },
        "targets": [vertex.to_record() for vertex in vertices],
    }


def print_enumeration_report(report):
    print("Genus-%s infinity vertices in the requested truncation"
          % report["genus"])
    print("allowed infinity degrees:", report["allowed_infinity_degrees"])
    print("max valence:", report["max_valence"])
    print("targets:", report["target_count"])
    for degree, count in sorted(report["counts_by_degree"].items()):
        print("  ambient degree %s: %s target(s)" % (degree, count))
    print("contact profiles:")
    for degree, profiles in sorted(report["profiles_by_degree"].items()):
        unique = sorted(set(tuple(profile) for profile in profiles))
        print("  ambient degree %s: %s" % (
            degree, ", ".join(str(profile) for profile in unique)
        ))


def print_computation_report(report):
    print("Genus-%s infinity-vertex computation: %s"
          % (report["genus"], report["status"]))
    print("targets:", report["target_count"],
          "solved:", report["solved_count"])
    if report["elapsed_seconds"] is not None:
        print("elapsed: %.1fs" % report["elapsed_seconds"])
    print("probe stages exhausted:", report["stages_exhausted"])
    if report.get("timed_out"):
        print("stopped early: wall-clock budget reached")
    orchestration = report["orchestration"]
    print("closure vertices:", orchestration["vertex_count"],
          "solved:", orchestration["solved_vertex_count"])
    print("compiled relations:", orchestration["compiled_relation_count"],
          "active:", orchestration["active_relation_count"])
    print()
    for record in report["targets"]:
        if record["solved"]:
            print("  %s = %s" % (record["description"], record["value"]))
        else:
            print("  %s = UNDETERMINED" % record["description"])
    frontier = orchestration.get("frontier")
    if frontier:
        print()
        print("frontier:")
        print("  rank-deficient blocks:", len(frontier["rank_deficiencies"]))
        print("  nonlinear relations:",
              frontier["nonlinear_relation_count"])
        print("  vertices without linear incidence:",
              len(frontier["vertices_without_linear_incidence"]))
    failed = orchestration.get("failed_level_stages")
    if failed:
        print("  levels without a usable probe family:", len(failed))


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, default=3)
    parser.add_argument(
        "--max-valence", type=int, default=3,
        help="largest number of infinity contacts on one vertex",
    )
    parser.add_argument(
        "--max-infinity-degree", type=int, default=None,
        help="cap the infinity degree; default is every balance-allowed degree",
    )
    parser.add_argument(
        "--max-psi-min", type=int, default=None,
        help="discard targets whose dimension-required psi_min exceeds this cap",
    )
    parser.add_argument("--max-markings", type=int, default=2)
    parser.add_argument(
        "--t-powers", type=int, nargs="+", default=(0, -1, -2)
    )
    parser.add_argument("--laurent-precision", type=int, default=8)
    parser.add_argument("--no-unit-relatives", action="store_true")
    parser.add_argument(
        "--chow-relations", action="store_true",
        help=("add under-dimensioned probes at Laurent powers having "
              "nonpositive residual Chow dimension"),
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
        "--effective-basis-only", action="store_true",
        help=("use boundary-compared log stationary and primary Chow probes, "
              "requiring every infinity factor to involve only psi_min and "
              "evaluations, as in CJR Theorem 10.1"),
    )
    parser.add_argument(
        "--additional-probe-ambient-degrees", type=int, nargs="*", default=(),
        help=("extra ambient probe degrees; nonmultiples of three have zero "
              "compact side (the genus-one assembly is regression-tested "
              "through ambient degree three)"),
    )
    parser.add_argument(
        "--future-probe-ambient-degree-ceiling", type=int,
        help=("after loading a checkpoint, restrict only newly compiled "
              "probe degrees while retaining every saved relation"),
    )
    parser.add_argument(
        "--no-punctured-axioms", action="store_true",
        help="disable the CJR contact-minus-one string/divisor/dilaton rows",
    )
    parser.add_argument(
        "--max-unit-insertions", type=int, default=1,
        help="largest number of reducible unit insertions in a mixed probe",
    )
    parser.add_argument(
        "--full-kernel", action="store_true",
        help="materialize exact kernel vectors in blocked frontier reports",
    )
    parser.add_argument(
        "--max-vertices", type=int, default=20000,
        help="dependency-closure safety cap (default: 20000)",
    )
    parser.add_argument(
        "--max-relations", type=int, default=20000,
        help="compiled-relation safety cap (default: 20000)",
    )
    parser.add_argument(
        "--time-budget", type=float, default=None,
        help="wall-clock seconds after which to stop and report the frontier",
    )
    parser.add_argument("--checkpoint-in")
    parser.add_argument("--checkpoint-out")
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
    parser.add_argument(
        "--list-only", action="store_true",
        help="enumerate the truncated family without computing values",
    )
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()

    if arguments.list_only:
        report = enumeration_report(
            genus=arguments.genus,
            max_valence=arguments.max_valence,
            max_infinity_degree=arguments.max_infinity_degree,
            max_psi_min=arguments.max_psi_min,
        )
        if arguments.json:
            print(json.dumps(report, indent=2, sort_keys=True))
        else:
            print_enumeration_report(report)
        return

    def progress(message):
        print(message, flush=True)

    computation = GenusThreeVertexComputation(
        genus=arguments.genus,
        max_valence=arguments.max_valence,
        max_infinity_degree=arguments.max_infinity_degree,
        max_psi_min=arguments.max_psi_min,
        max_markings=arguments.max_markings,
        t_powers=tuple(arguments.t_powers),
        include_unit_relatives=(
            not arguments.no_unit_relatives
            and not arguments.effective_basis_only
        ),
        include_chow_relations=(
            arguments.chow_relations or arguments.effective_basis_only
        ),
        max_chow_unit_insertions=arguments.max_chow_unit_insertions,
        chow_primary_only=(
            arguments.primary_chow_relations_only
            or arguments.effective_basis_only
        ),
        effective_basis_only=arguments.effective_basis_only,
        additional_probe_ambient_degrees=(
            arguments.additional_probe_ambient_degrees
        ),
        include_punctured_axioms=not arguments.no_punctured_axioms,
        max_unit_insertions=(
            0 if (arguments.no_unit_relatives
                  or arguments.effective_basis_only)
            else arguments.max_unit_insertions
        ),
        laurent_precision=arguments.laurent_precision,
        include_kernel_basis=arguments.full_kernel,
        zero_vertex_cache=arguments.zero_vertex_cache,
        hodge_cache=arguments.hodge_cache,
        zero_vertex_strategy=arguments.zero_vertex_strategy,
        max_vertices=arguments.max_vertices,
        max_relations=arguments.max_relations,
        progress=None if arguments.quiet else progress,
    )
    if arguments.checkpoint_in:
        computation.load_checkpoint(arguments.checkpoint_in)
    if arguments.future_probe_ambient_degree_ceiling is not None:
        computation.restrict_future_probe_degree(
            arguments.future_probe_ambient_degree_ceiling
        )
    report = computation.run(
        checkpoint_path=arguments.checkpoint_out,
        time_budget=arguments.time_budget,
    )
    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print_computation_report(report)
    if arguments.checkpoint_out:
        print("checkpoint:", arguments.checkpoint_out)


if __name__ == "__main__" and os.path.basename(sys.argv[0]) in (
        "genus_three_infinity_vertices.sage",
        "genus_three_infinity_vertices.sage.py"):
    _main()
