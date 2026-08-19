r"""Compare degree-zero localization and Givental--Teleman evaluations.

Run from the repository root.  With no insertion options this profiles the
unmarked twisted genus-three constant map.  ``--h-power`` adds one insertion.
Pass ``--givental`` to include the stable-graph reconstruction; a small
genus-one request is recommended for a quick comparison.
"""

import argparse
import json
import time

load("o3_fixed_locus_graphs.sage")
load("o3_givental_teleman.sage")


def profile_degree_zero_vertex(genus=3, h_power=None, psi_power=0,
                               scale=1, include_twist=True,
                               compare_givental=False):
    insertions = () if h_power is None else (
        TwistedInsertion(h_power, psi_power, QQ(scale)),
    )
    request = TwistedZeroVertexRequest(genus, 0, insertions)

    direct = P2FixedLocusEvaluator(include_twist=include_twist)
    started = time.time()
    direct_value = direct.evaluate(request)
    direct_seconds = time.time() - started

    explicit = P2FixedLocusEvaluator(include_twist=include_twist)
    scalar, normalized = request.normalized_scales()
    started = time.time()
    planned = explicit.planned_insertions(normalized)
    graphs = explicit.fixed_graphs(normalized)
    explicit_value = scalar * sum(
        explicit.graph_contribution(graph, planned) for graph in graphs
    )
    explicit_seconds = time.time() - started

    report = {
        "genus": int(genus),
        "markings": int(request.valence),
        "twisted": bool(include_twist),
        "fixed_graphs": int(len(graphs)),
        "equal": bool(direct_value == explicit_value),
        "direct_seconds": float(direct_seconds),
        "explicit_seconds": float(explicit_seconds),
        "explicit_over_direct": (
            None if direct_seconds == 0
            else float(explicit_seconds / direct_seconds)
        ),
        "direct_cache": direct.cache_info(),
        "explicit_cache": explicit.cache_info(),
    }
    if compare_givental:
        givental = O3DegreeZeroGiventalBackend(include_twist=include_twist)
        started = time.time()
        givental_value = givental.evaluate(request)
        givental_seconds = time.time() - started
        report.update({
            "givental_equal": bool(direct_value == givental_value),
            "givental_seconds": float(givental_seconds),
            "givental_over_direct": (
                None if direct_seconds == 0
                else float(givental_seconds / direct_seconds)
            ),
        })
    return report


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, default=3)
    parser.add_argument("--h-power", type=int, choices=(0, 1, 2))
    parser.add_argument("--psi-power", type=int, default=0)
    parser.add_argument("--scale", default="1")
    parser.add_argument("--ordinary", action="store_true")
    parser.add_argument("--givental", action="store_true")
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    report = profile_degree_zero_vertex(
        genus=arguments.genus,
        h_power=arguments.h_power,
        psi_power=arguments.psi_power,
        scale=QQ(arguments.scale),
        include_twist=not arguments.ordinary,
        compare_givental=arguments.givental,
    )
    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print("degree-zero vertex: g=%s, n=%s, twisted=%s" % (
            report["genus"], report["markings"], report["twisted"]
        ))
        print("fixed graphs:", report["fixed_graphs"])
        print("exact values agree:", report["equal"])
        print("direct seconds: %.6f" % report["direct_seconds"])
        print("explicit seconds: %.6f" % report["explicit_seconds"])
        print("explicit/direct: %.3fx" % report["explicit_over_direct"])
        if arguments.givental:
            print("Givental value agrees:", report["givental_equal"])
            print("Givental seconds: %.6f" % report["givental_seconds"])
            print("Givental/direct: %.3fx" % report["givental_over_direct"])
        print("direct Hodge cache:", report["direct_cache"]["hodge_cache"])


if __name__ == "__main__":
    _main()
