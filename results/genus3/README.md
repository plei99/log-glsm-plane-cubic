# Genus-three result files

Files named `stabilized-*` use the corrected stabilization pullback and contain
valid dimension-zero scalar localization relations.

Files named `g3-d0-minus5-class-*` are diagnostic artifacts from an attempted
scalar projection of under-dimensioned Chow-class probes. They are deliberately
incompatible with the current configuration schema and **must not** be used to
solve infinity vertices. Their missing terms are positive-dimensional
effective cycles. See `docs/COMPUTING_INFINITY_VERTICES.md`, “Current status of
the `(3,0,-5)` vertex”.

Files named `chow-v7-*` use the graded infinity backend, but their compiled
relations predate the corrected stable zero-side flag-descendant bound. They
are retained as diagnostic history and are rejected by checkpoint version 8.
Their Hodge-integral and zero-vertex SQLite caches remain reusable.

`chow-v7-genus2-primary-v1-d3.json` records the former positive-degree
diagnostic. Its residual `-3/2` is now understood: the compiler omitted
equivariant zero-side flag descendants. Fresh ambient-degree-one, two, and
three relations pass with residual zero.

`chow-v8-positive-d1.json` is the first certified checkpoint made with the
corrected flag-descendant bound. It contains the one-marking stage with ambient
probe degrees zero and one: 419 compiled relations, 341 active relations, 251
vertices, and 90 solved values. The `(3,0,-5;H)` target lies in a 161-column,
rank-108 component (nullity 53), so the checkpoint records it as undetermined.

`chow-v8-positive-d2-stage2-convolution.json` is the completed degree-zero,
one, and two two-marking system after incremental zero-state convolution.  It
contains 1789 compiled relations and 401 solved values.  The target component
has rank 351 on 518 columns (nullity 167).

`chow-v8-d2-plus-degree0-m3.json` adds primary three-marking, ambient-degree
zero relations.  These solve five more values but leave the target component
at rank 351 on 515 columns.

`chow-v8-d2-m3-unit1-degree0.json` additionally includes degree-zero unit
relatives through three markings.  It reduces the nonlinear frontier from
five relations to one and solves 412 values, but enlarges the target component
to 846 columns of rank 355.

`chow-v8-d3-stage1.json` adds the one-marking ambient-degree-three
relations, whose compact side is nonzero.  The exact target component has rank
357 on 886 columns (nullity 529), so `(3,0,-5;H)` is still undetermined.  This
experiment belongs to the generic stabilized-descendant branch.  Its ordinary
`contact_psi` columns are not part of the effective-invariant basis in CJR
Theorem 10.1, so it is retained as a compiler diagnostic rather than the
current hypersurface-reconstruction frontier.

Files named `effective-basis-*` are historical version-8 runs containing only primary Chow probes.  Their
configuration sets `effective_basis_only=true`, and compilation fails if an
infinity factor has a nonempty `contact_psi` tuple.  The sequence
`effective-basis-d0-m1.json` through `effective-basis-d0-m5.json` explores
ambient degree zero through five primary markings.  The `d1` and `d2`
sequences add the vanishing compact sectors in ambient degrees one and two.

`effective-basis-d2-m2.json` is the strongest of these historical primary-only
checkpoints. It has 1259 compiled relations, 1205 active relations, 692 vertices,
and 349 solved values.  The final target component has rank 278 on 343 columns
(nullity 65).  In the computed kernel basis, exactly one vector moves
`(3,0,-5;H)`; it has support on 18 genuine effective invariants and no ordinary
contact descendants.  The
target is therefore an undetermined basic effective-invariant direction at
this truncation, not evidence for a missing contact-psi backend.

Checkpoint version 9 added log-domain stationary rows. Those rows used an
incorrect specialization of CJR III (8.21): the entire weight `(3H-t)` must
be divided by contact. Version 10 fixes it. The genus-three one-point row
now has compact side `-31/967680` and target coefficient `-1/72`.
Previously reported rank gains from the faulty descendant rows are
withdrawn. These saved results are historical; rebuild descendant relations
before drawing conclusions about the corrected target block. The legacy
import keeps unaffected Chow-scalar rows and discards affected log-descendant
rows, all solver values, and all stage bookkeeping.

The corrected November 5, 2026 rerun is in
`corrected-2026-11-05/d1-extension/`. This is the authoritative completed
extension for that rerun: 1,411 relations, 408 determined values (171 in
genus three, all zero), and a remaining exact 490-by-302 linear system of
rank 248. Projection to the 35 genus-three primary candidates without
(-1) contacts leaves 30 parameters, hence five independent affine
constraints. All five one-contact coordinates remain independently free
in this finite system. `linear-certificate.json` records a particular
solution and a kernel verified against every row. See
[the detailed results](../../docs/GENUS_THREE_RESULTS_2026-11-05.md).
The parent directory also contains earlier stage snapshots; consult
`run-manifest.json` rather than treating those as the final result.

The class-valued extension is in `lower-cycles-2026-11-05/`. It computes
unmarked genus-three dual-Euler O(3) classes in ambient degrees zero and
one, and codimensions three, four, and five. Nine auxiliary pole relations
have exact Pixton certificates; six finite-part comparisons verify agreement
between two choices of torus weights. The matrices, strata ordering, rational
certificate weights, expanded classes, and all nine readable equations are
saved together. See [the derivation and results](../../docs/GENUS_THREE_LOWER_CYCLES_2026-11-05.md).
The CJR consequences identify the full infinity-containing aggregates;
they do not reconstruct individual positive-dimensional effective cycles.
