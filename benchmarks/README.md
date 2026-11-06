# November 2026 hill climb

The fixed six-workload suite improved from 33.81 s to 5.32 s (6.35x),
using the ratio of median total compute times from three alternating
before/after runs. Individual speedups range from 1.21x to 8.03x.
All six mathematical outputs agree exactly. See
[the detailed audit](../docs/PERFORMANCE_AND_DRAFT_AUDIT.md) and
[raw measurements](results-2026-11-05/summary.json).

Run from the repository root with the pinned Sage environment's Python:

```bash
python benchmarks/benchmark_hill_climb.py --output /tmp/measurement.json
python benchmarks/compare_hill_climb.py MEASUREMENT_DIRECTORY OUTPUT_DIRECTORY
```

The comparison expects `before-1.json` through `before-3.json` and
`after-1.json` through `after-3.json`. It checks every graph representative
and automorphism order, every compiled relation coefficient, and exact Sage
fraction-field equality for O(3) values. Rational-function string hashes can
differ even when the fractions are equal; the comparator checks the actual
fractions. Setup/import time and result serialization are outside the timed
region. Each case creates a fresh evaluator without disk result caches.

The unchanged-output performance cases use a primary CJR relation. The
separately corrected log-descendant relation intentionally differs from the
old, incorrect output and is checked against CJR III (8.21).

# Genus-three performance experiments

These scripts separate the combinatorial CJR graph layer from the exact
`admcycles` integration layer.

Profile a graph family:

```bash
sage benchmarks/profile_cjr_enumerator.sage --genus 3 --markings 0 --degree 3
```

Compare the direct all-genus constant-map path with the explicit auxiliary
fixed-locus sum:

```bash
sage benchmarks/profile_degree_zero_vertex.sage
sage benchmarks/profile_degree_zero_vertex.sage \
  --h-power 1 --psi-power 4 --scale 1/3 --json
```

Add the exact degree-zero Givental--Teleman stable-graph reconstruction to a
small comparison with:

```bash
sage benchmarks/profile_degree_zero_vertex.sage \
  --genus 1 --h-power 0 --givental
```

The Givental path is primarily an independent reconstruction. It need not
beat the specialized constant-map formulas, and its cost grows with the
number of stable graphs.

Profile the positive-degree calibration, including the S-matrix descendant
correction, against localization with:

```bash
sage benchmarks/profile_calibrated_vertex.sage
sage benchmarks/profile_calibrated_vertex.sage --degree 2 --json
```

Compare the same deterministic zero-vertex classification kernel in Python
and Go:

```bash
python3 benchmarks/hotloop_python.py
go run benchmarks/hotloop_go.go
```

On the development machine (2026-08-07), cProfile measured the `(3,0,3)` CJR
family at 11.20 seconds for 478 graphs.  `_zero_vertex_types` accounted for
3.02 seconds (27%); even an infinitely fast replacement of only that function
would therefore improve this graph family by at most 1.37x.  The native
microbenchmark measures the attainable speedup of that isolated kernel, not
an end-to-end speedup.  With 2,048,000 deterministic calls, Python took 1.68
seconds and Go took 0.092 seconds (18.2x faster), with identical checksums.
Applied only to the measured 27% hot spot, Amdahl's law predicts about a 1.34x
end-to-end improvement for the `(3,0,3)` enumerator.

The larger one-mark `(3,1,3)` stationary probe took 153 seconds merely to
inventory its zero vertices with all integrations disabled.  Canonicalizing
marked-insertion permutations reduced the distinct zero-vertex requests from
1,510 to 625.  This reduction is more important for complete relations: every
eliminated request avoids a full fixed-locus sum and potentially many
`admcycles` tautological products.

For ambient degree zero, the stationary plus unit-depth-two one-point family
contains 6 probes. It took 328 seconds to inventory 1,502 requests; 920 are
genus zero and 582 have higher genus. The precompute command persists this
manifest so the five-minute enumeration is paid only once.

The evaluator now has an all-genus direct constant-map path. It reduces Hodge
monomials with Mumford's relation through genus three and uses closed formulas
for pure psi, `lambda_g`, and top-Hodge-triple descendants before falling back
to admcycles. On the same machine, the unmarked twisted `(g,d)=(3,0)` request
changed from more than five minutes without completion to 0.068 seconds; the
explicit three-fixed-point sum took 0.100 seconds and returned the same exact
value. A sparse one-divisor top-descendant request remained small but favored
the explicit sum (0.218 seconds versus 0.384 seconds), illustrating that the
main gain is removal of pathological tautological products, not a uniform
speedup.

Final zero vertices and individual Hodge integrals can be stored in separate
SQLite WAL databases. `o3_zero_vertex_precompute.sage --shard-count N
--shard-index I` partitions the sorted request manifest, so multiple Sage
processes can safely accumulate one shared exact cache.

A wholesale rewrite would need to keep Sage for exact rational functions,
graph canonicalization, and `admcycles`, or reimplement those systems.  The
current evidence supports a small native extension for the graph enumerator
only after the persistent zero-vertex cache is populated; it does not support
rewriting the mathematical backend in Go.
