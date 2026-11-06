# Regression fixtures

These small extracts make the suite independent of local scan directories.
They record results obtained before the module refactor; tests recompute the
mathematical quantities and compare against them.

- `p2-3-g2-classes.json.gz`: the elliptic cubic, genus two, up to three markings,
  from `results/higher-marked/2026-11-06`. Used to compare relation spans from
  the original full permutation action with the adjacent-transposition action.
- `pixton-g2-n4-r2.sobj`: the exact rational Pixton matrix for genus two, four
  markings, codimension two, from `results/marked-frontier/2026-11-06`. Used to
  verify descendant forgetting and rational-tail section identities.
- `f1-numbers.json`: the first two numerical rows and target metadata of each
  recorded Hirzebruch F1 target in the higher-marked run. Used to check that
  the convex-GKM extension preserves existing Chern numbers and GW values.

The fixtures retain the original conventions and limitations; none is a new
claim about individual infinity cycles or completeness of Pixton relations.
