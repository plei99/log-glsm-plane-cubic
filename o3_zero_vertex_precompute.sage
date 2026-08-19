r"""Inventory and precompute exact ``O(3)``-twisted plane zero vertices.

The CJR compiler can first be run with a collecting backend.  This discovers
the individual twisted invariants needed by a finite probe family without
paying for admcycles.  A second pass evaluates the distinct requests once and
checkpoints every completed value in a reusable SQLite (recommended) or JSON
cache.  Deterministic shards may safely share the SQLite stores.
"""

import argparse
import hashlib
import json
import os
import sys
import time
from collections import Counter

load("cjr_graph_contributions.sage")
load("cjr_probe_factory.sage")


class O3ZeroVertexPrecomputer(SageObject):
    """Discover and evaluate the zero-vertex closure of explicit probes."""

    INVENTORY_FORMAT = "log-glsm-o3-zero-vertex-inventory"
    INVENTORY_VERSION = 1

    def __init__(self, cache_path, rings=None, progress=None,
                 inventory_path=None, hodge_cache_path=None):
        self.rings = rings or PlaneCubicCoefficientRing(8)
        self.cache_path = os.path.abspath(
            specialized_zero_vertex_cache_path(
                cache_path, "nonequivariant"
            )
        )
        self.inventory_path = os.path.abspath(
            inventory_path or (self.cache_path + ".inventory.json")
        )
        self.hodge_cache_path = os.path.abspath(
            hodge_cache_path or (self.cache_path + ".hodge.sqlite")
        )
        self.progress = progress
        self.collector = CollectingTwistedZeroVertexBackend(self.rings)
        self.compiler = PlaneCubicGraphContributionCompiler(
            self.rings, self.collector
        )
        self.backend = FullTwistedZeroVertexBackend(
            self.rings,
            cache_path=self.cache_path,
            autosave=True,
            hodge_cache_path=self.hodge_cache_path,
            base_weight_specialization="nonequivariant",
        )
        self.probes = []
        self.inventory_hit = False

    def _report(self, message):
        if self.progress is not None:
            self.progress(message)

    @staticmethod
    def _probe_records(probes):
        return [probe.to_record() for probe in probes]

    def _family_digest(self, probes):
        encoded = json.dumps(
            self._probe_records(probes), sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def _inventory(self):
        if not os.path.exists(self.inventory_path):
            return {
                "format": self.INVENTORY_FORMAT,
                "version": int(self.INVENTORY_VERSION),
                "families": {},
            }
        with open(self.inventory_path) as stream:
            inventory = json.load(stream)
        if inventory.get("format") != self.INVENTORY_FORMAT \
                or inventory.get("version") != int(self.INVENTORY_VERSION):
            raise ValueError("unsupported O(3) request-inventory format")
        return inventory

    def _save_inventory(self, inventory):
        directory = os.path.dirname(self.inventory_path)
        if directory and not os.path.isdir(directory):
            os.makedirs(directory)
        temporary = self.inventory_path + ".tmp"
        with open(temporary, "w") as stream:
            json.dump(inventory, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary, self.inventory_path)

    def collect(self, probes):
        probes = tuple(probes)
        digest = self._family_digest(probes)
        inventory = self._inventory()
        saved = inventory.get("families", {}).get(digest)
        if saved is not None:
            if saved.get("probes") != self._probe_records(probes):
                raise ValueError("zero-vertex inventory digest collision")
            self.inventory_hit = True
            self.probes.extend(probes)
            requests = tuple(
                TwistedZeroVertexRequest.from_record(record)
                for record in saved.get("requests", ())
            )
            self._report(
                "loaded %s zero-vertex requests from %s"
                % (len(requests), self.inventory_path)
            )
            return requests

        self.inventory_hit = False
        self.collector = CollectingTwistedZeroVertexBackend(self.rings)
        self.compiler = PlaneCubicGraphContributionCompiler(
            self.rings, self.collector
        )
        for index, probe in enumerate(probes, start=1):
            started = time.time()
            # Stabilization changes only contracted marked zero tails, which
            # contain no stable O(3)-twisted zero request.  Compile the actual
            # probe convention anyway so the inventory path exercises the
            # same graph support as the main provider.
            self.compiler.compile_probe(probe)
            self.probes.append(probe)
            self._report(
                "collected probe %s/%s in %.1fs; %s distinct requests"
                % (index, len(probes), time.time() - started,
                   len(self.collector.requests))
            )
        requests = self.collector.requests
        inventory.setdefault("families", {})[digest] = {
            "probes": self._probe_records(probes),
            "requests": [request.to_record() for request in requests],
        }
        self._save_inventory(inventory)
        return requests

    def run(self, probes, time_budget=None, collect_only=False,
            include_requests=False, shard_count=1, shard_index=0):
        probes = tuple(probes)
        shard_count = int(shard_count)
        shard_index = int(shard_index)
        if shard_count < 1:
            raise ValueError("shard_count must be positive")
        if shard_index < 0 or shard_index >= shard_count:
            raise ValueError("shard_index must lie in [0, shard_count)")
        started = time.time()
        requests = self.collect(probes)
        shard_requests = tuple(
            request for index, request in enumerate(requests)
            if index % shard_count == shard_index
        )
        collection_seconds = time.time() - started
        deadline = (
            None if time_budget is None else started + float(time_budget)
        )
        evaluated = 0
        skipped = 0
        timed_out = False
        if not collect_only:
            for index, request in enumerate(shard_requests, start=1):
                if self.backend.is_cached(request):
                    skipped += 1
                    continue
                if deadline is not None and time.time() >= deadline:
                    timed_out = True
                    break
                request_started = time.time()
                self.backend.evaluate(request)
                evaluated += 1
                self._report(
                    "evaluated request %s/%s in %.1fs: %s"
                    % (index, len(shard_requests),
                       time.time() - request_started,
                       request)
                )
        if not os.path.exists(self.cache_path):
            self.backend.save_cache()
        cached = sum(self.backend.is_cached(request) for request in requests)
        shard_cached = sum(
            self.backend.is_cached(request) for request in shard_requests
        )
        type_counts = Counter(
            (int(request.genus), int(request.valence), int(request.degree))
            for request in requests
        )
        report = {
            "probe_count": int(len(probes)),
            "request_count": int(len(requests)),
            "shard_count": int(shard_count),
            "shard_index": int(shard_index),
            "shard_request_count": int(len(shard_requests)),
            "evaluated": int(evaluated),
            "already_cached": int(skipped),
            "cached_for_family": int(cached),
            "complete": bool(cached == len(requests)),
            "shard_complete": bool(shard_cached == len(shard_requests)),
            "timed_out": bool(timed_out),
            "collection_seconds": float(collection_seconds),
            "elapsed_seconds": float(time.time() - started),
            "cache": self.cache_path,
            "inventory": self.inventory_path,
            "hodge_cache": self.hodge_cache_path,
            "inventory_hit": bool(self.inventory_hit),
            "cache_info": self.backend.cache_info(),
            "request_types": [
                {
                    "genus": key[0],
                    "markings": key[1],
                    "degree": key[2],
                    "count": int(count),
                }
                for key, count in sorted(type_counts.items())
            ],
        }
        if include_requests:
            report["requests"] = [request.to_record() for request in requests]
        return report


def probe_family(genus, ambient_degree, max_point_markings=1,
                 max_unit_insertions=0):
    """Return the finite stationary/mixed family used by the CLI."""
    factory = ProbeFactory(QQ)
    probes = list(factory.stationary_candidates(
        genus, ambient_degree, max_markings=max_point_markings
    ))
    if max_unit_insertions:
        probes.extend(factory.mixed_unit_candidates(
            genus, ambient_degree,
            max_point_markings=max_point_markings,
            max_unit_insertions=max_unit_insertions,
        ))
    return tuple(probes)


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


if __name__ == "__main__" and os.path.basename(sys.argv[0]) in (
        "o3_zero_vertex_precompute.sage",
        "o3_zero_vertex_precompute.sage.py"):
    _main()
