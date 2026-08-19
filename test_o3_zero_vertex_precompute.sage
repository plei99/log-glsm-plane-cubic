"""Tests for zero-vertex request inventory and precomputation plumbing."""

load("o3_zero_vertex_precompute.sage")

import os
import tempfile
import json


def run_tests():
    probes = probe_family(1, 0, max_point_markings=1,
                          max_unit_insertions=1)
    assert len(probes) == 3
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "zero-cache.json")
        precomputer = O3ZeroVertexPrecomputer(path)
        report = precomputer.run(probes, collect_only=True)
        assert report["probe_count"] == 3
        assert report["request_count"] > 0
        assert not report["complete"]
        assert not report["inventory_hit"]
        assert report["cache"].endswith(
            "zero-cache.weights-nonequivariant-quadratic.json"
        )
        assert os.path.exists(report["cache"])
        assert os.path.exists(report["inventory"])
        assert os.path.exists(report["hodge_cache"])
        assert json.loads(json.dumps(report)) == report

        resumed = O3ZeroVertexPrecomputer(path)
        resumed_report = resumed.run(probes, collect_only=True)
        assert resumed_report["inventory_hit"]
        assert resumed_report["request_count"] == report["request_count"]

        left = resumed.run(
            probes, collect_only=True, shard_count=2, shard_index=0
        )
        right = resumed.run(
            probes, collect_only=True, shard_count=2, shard_index=1
        )
        assert left["shard_request_count"] \
            + right["shard_request_count"] == report["request_count"]


run_tests()
print("all O(3) zero-vertex precomputation tests passed")
