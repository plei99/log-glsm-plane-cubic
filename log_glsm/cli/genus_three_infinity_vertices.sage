r"""cli / genus_three_infinity_vertices."""

from sage.all import *

import argparse
import json
from log_glsm.infinity.genus_three import (
    GenusThreeVertexComputation,
    enumeration_report,
    print_computation_report,
    print_enumeration_report,
)

__all__ = []


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
