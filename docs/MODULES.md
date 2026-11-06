# Module layout

The implementation now lives in `log_glsm`, with separate Python namespaces
and explicit imports. The root `.sage` files and `elliptic_cubic_gw.py` are
compatibility entry points for existing notebooks and command lines. New code
should import the objects it uses from their owning modules.

| Package | Responsibility | Typical interface |
| --- | --- | --- |
| `log_glsm.core` | Degree and sign conventions, Hodge integration, sparse tautological polynomials | `ProbeSpec`, `HodgeIntegralBackend`, `GenusFourHodgeIntegrals` |
| `log_glsm.curves` | Exact elliptic GW theory, rational power series, Bloch–Okounkov coefficients | `EllipticCurveStationaryTheory`, `connected_coefficient` |
| `log_glsm.localization` | Fixed graphs, vertex evaluation, persistent stores, Givental reconstruction | `P2FixedLocusGraphEnumerator`, `FullTwistedZeroVertexBackend` |
| `log_glsm.infinity` | CJR graphs and factors, punctured axioms, equation selection and reconstruction | `PlaneCubicFullEquationProvider`, `InfinityVertexOrchestrator` |
| `log_glsm.targets` | Projective, quartic K3, toric, Grassmannian and partial-flag presentations, target-specific localization and known ordinary inputs | `ProjectiveCompleteIntersection`, `ToricCompleteIntersection`, `convex_gkm_targets` |
| `log_glsm.cycles` | Stabilization, marked and unmarked cycle calculations, descendant reductions | `genus_four_backend`, `higher_marked_backend` |
| `log_glsm.cli` | Argument parsing and command-line reports | Installed commands listed below |

`calculations/` contains reproducible batch drivers and certificate/report tools.
They import the mathematical modules and use package-qualified imports for
other calculation drivers. `benchmarks/` contains performance experiments.
`tests/` contains the regression suite and its small, recorded fixtures.

## Using the modules

The current full suite is run with SageMath 10.9 on the ARM64 Linux `agents`
machine. Activate the Sage Python environment, then install the package:

```bash
python -m pip install --no-deps -e .
./run_tests.sh
```

The runner also works directly from a source checkout without installation: it
sets the repository import path and working directory itself. `FAST=1` skips
the explicitly listed long suites; run the complete suite before publishing.

For example, in Sage's Python environment:

```python
from log_glsm.targets.projective import ProjectiveCompleteIntersection
from log_glsm.cycles.marked import genus_four_backend

quartic = ProjectiveCompleteIntersection(3, (4,))
assert quartic.euler_characteristic() == 24
backend = genus_four_backend(quartic, genus=2)
value = backend.numerical_localization((1,))
```

The Fraction-based curve code also works in ordinary Python without importing
Sage:

```python
from log_glsm.curves.elliptic import O3TwistedPlane
assert O3TwistedPlane.nonequivariant_i0_coefficient(3) == 1680
```

Installed commands are `log-glsm-elliptic`, `log-glsm-infinity`,
`log-glsm-genus-three`, `log-glsm-equations`, and `log-glsm-precompute`.
Their arguments match the corresponding historical commands. For example:

```bash
log-glsm-infinity --help
sage log_glsm_infinity_orchestrator.sage --help
python -m calculations.marked_frontier --help
sage tests/test_marked_forgetting.sage
```

## Sage source and compatibility

Sage source remains `.sage`: numerical literals and `^` retain Sage's exact
arithmetic semantics. A package-scoped importer passes these files through
Sage's preparser. It handles only `log_glsm.*`, and Python's normal import
cache owns module identity. Repeated imports do not redefine classes or clear
their caches. The modules can be imported from outside the repository after
installation, and vendored `admcycles` is located relative to the package.

The compatibility adapter reconstructs the historical transitive exports for
`load('old_filename.sage')`. It references the same classes and functions as
the new modules. New implementation code does not use this adapter or execute
other source files with `load()`. Data loads such as Sage `.sobj` matrices are
unchanged. Command-line parsing runs only from an entry point, not as a side
effect of importing the mathematical module.

The largest files were split along their existing responsibilities:

- Fixed-locus graphs, tautological polynomials, SQLite storage and evaluation
  now have separate modules.
- Orchestration configuration and command-line parsing are separate from the
  reconstruction algorithm.
- Genus-four Hodge recursion and the target catalogue are separate from marked
  localization.
- Pure-Python power series, stationary elliptic theory, genus-two
  reconstruction and the elliptic command line are separate modules.

## Results, caches and mathematical scope

The refactor does not change formulas, admissible degrees, unknown-value
handling, or the distinction between ambient classes and individual infinity
cycles. In particular, it does not turn aggregate comparisons into
contact-resolved infinity-cycle calculations or add Pixton certificates.

Large run outputs, caches, databases, generated Sage Python files and local
visualizations remain on disk and are ignored by Git. Research reports and
source audits are versioned. The small regression fixtures are versioned under
`tests/fixtures`, so tests do not require those local scan outputs. The marking
frontier regression writes fresh output in a temporary directory.

New report hashes and the stabilization cache fingerprint use
`log_glsm.provenance` to include the implementations, not just old compatibility
filenames. Historical result provenance is left as a record of the original
run; reorganizing source files does not make those old hashes describe the new
tree. Reproducing an entire historical search may require its earlier output
artifacts as documented in the corresponding report.

Verification includes the full mathematical regression suite, legacy command
entry points, exact Hodge requests across repeated legacy loads and pickle
round trips, a spawned worker importing a Sage function, and a vendored Chow
calculation from outside the checkout.

On 2026-11-07, all 33 Sage suites and 20 Python unit tests passed. The full
suite took 223 seconds in a clean export of the staged tree, without local
scan outputs. A separately built wheel also passed vendored Hodge arithmetic
and all-target catalogue imports outside the checkout. Comparing the syntax
trees of the 199 moved functions and classes found only the intended change
to the vendored-library path resolver; mathematical function bodies were
unchanged by the split.
