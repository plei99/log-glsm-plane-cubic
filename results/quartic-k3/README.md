# Ordinary quartic K3 calculations

The completed calculation is under `2026-11-05/`. It uses the dual Euler
twist of O(4) on P3, with ordinary K3 GW vanishing as the compact CJR input.
Reduced K3 GW invariants are a different theory.

See [the derivation and results](../../docs/QUARTIC_K3_RESULTS_2026-11-05.md)
and [the primary-source audit](../../docs/QUARTIC_K3_SOURCES.md).

The run covers unmarked genera two and three, all Chow codimensions, and
ambient degrees zero through three. `independent-relations.md` prints the
retained relation bases, with a dictionary of decorated graphs.
`cjr-aggregate-classes.md` gives all twenty CJR aggregate identities.
`coverage.json` records exact relation ranks. `basis-relations-g*-r*.json.gz`
contains the row-reduced rational strata vectors and their Pixton row coefficients;
the matching matrices and strata order are saved alongside them.
The original formal localization expressions and certificates are retained
losslessly as `final-relations-g*-r*.json.gz`; the earlier combined stage is
also compressed. The scripts read either the plain or gzip representation.

`summary.json` records the degree-one scan; `degree-two-summary.json` and
`degree-three-summary.json` record successive extensions. The degree-three
summary retains the initial eight-ray scan under `previous_scan` when the
second eight rays were run as a continuation. All final relation files are
cumulative through degree three. Earlier files remain as stage snapshots.

From the project root, in the pinned Sage environment:

```bash
python calculations/quartic_k3_relations.py --rays 64
python calculations/quartic_k3_degree_two.py --degree 2
python calculations/quartic_k3_degree_two.py --degree 3
python calculations/quartic_k3_export.py
```

The source hashes and output hashes are in `2026-11-05/run-manifest.json`.
No individual positive-dimensional effective cycle is asserted to vanish
or to have been reconstructed; the CJR identities determine the sums of
all infinity-containing contributions.
