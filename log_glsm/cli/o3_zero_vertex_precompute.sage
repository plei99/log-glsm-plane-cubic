r"""cli / o3_zero_vertex_precompute."""

from sage.all import *

import argparse
import json
from log_glsm.localization.precompute import (
    O3ZeroVertexPrecomputer,
    probe_family,
)

__all__ = []


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, required=True)
    parser.add_argument("--ambient-degree", type=int, default=0)
    parser.add_argument("--max-point-markings", type=int, default=1)
    parser.add_argument("--max-unit-insertions", type=int, default=0)
    parser.add_argument("--cache", required=True)
    parser.add_argument(
        "--inventory",
        help="request-manifest JSON; default is CACHE.inventory.json",
    )
    parser.add_argument(
        "--hodge-cache",
        help="shared SQLite Hodge table; default is CACHE.hodge.sqlite",
    )
    parser.add_argument("--time-budget", type=float)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--collect-only", action="store_true")
    parser.add_argument(
        "--include-requests", action="store_true",
        help="include the full request manifest in JSON output",
    )
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    probes = probe_family(
        arguments.genus,
        arguments.ambient_degree,
        arguments.max_point_markings,
        arguments.max_unit_insertions,
    )
    precomputer = O3ZeroVertexPrecomputer(
        arguments.cache,
        inventory_path=arguments.inventory,
        hodge_cache_path=arguments.hodge_cache,
        progress=None if arguments.json else lambda message: print(
            message, flush=True
        ),
    )
    report = precomputer.run(
        probes,
        time_budget=arguments.time_budget,
        collect_only=arguments.collect_only,
        include_requests=arguments.include_requests,
        shard_count=arguments.shard_count,
        shard_index=arguments.shard_index,
    )
    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print("zero-vertex requests:", report["request_count"])
        print("cached for this family:", report["cached_for_family"])
        print("complete:", report["complete"])
        print("cache:", report["cache"])
