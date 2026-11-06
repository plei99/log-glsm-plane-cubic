r"""cli / log_glsm_infinity_orchestrator."""

from sage.all import *

import argparse
import json
from log_glsm.infinity.config import (
    InfinityOrchestrationConfig,
)
from log_glsm.infinity.equations import (
    PlaneCubicFullEquationProvider,
)
from log_glsm.infinity.orchestrator import (
    InfinityVertexOrchestrator,
)
from log_glsm.infinity.solver import (
    EffectiveVertex,
)

__all__ = []


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
