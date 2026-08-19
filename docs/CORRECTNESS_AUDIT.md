# Correctness audit

This note records the mathematical scope checked against the Chen--Janda--Ruan
papers and the independent regression checks in the repository.

## Formula cross-check

The contact-resolved compiler uses the following specializations for
`X=P^2` and `E=O(3)`.

| Code component | Paper formula | Audited specialization |
|---|---|---|
| Infinity balance | CJR III, Section 9.3 | `sum(abs(c_i)) = 2g-2+n-3D` |
| Edge factor | CJR III, Section 9.3 | `c^c/(c! (t-3H)^c)` |
| Nonspecial zero vertex | CJR III, Section 9.3 | `(t-3H)^2/c`, with `c>1` |
| Marked zero vertex | CJR III, (8.21) and Section 9.3 | `(t-3H)(3H/c-t)^b` for the log-domain psi class |
| Nodal zero vertex | CJR III, Section 9.3 | `cc'/(c+c') (t-3H)` after the diagonal identification |
| Stable zero vertex | CJR III, Section 9.3 | `e((R pi_*f^*O(3))^vee)` and the flag factors `1/((t-3H)/c-psi)` |
| Stable infinity vertex | CJR III, (8.23) | `t/(-t-psi_min)` |
| Reduced comparison | CJR III, (9.7) | `(-1)^(1-g+3D)` times the cubic GW class |
| Negative Laurent coefficients | CJR III, Remark 8.5 | exact vanishing relations |
| Contact `-1` relations | *Punctured logarithmic R-maps*, Theorem 8.1 and Corollary 8.3 | string, divisor, and dilaton relations with `psi_DF=-3H` |

The auxiliary `P^2` torus is removed only after summing a complete fixed-locus
family.  The default path `(0,epsilon,epsilon^2)` avoids resonant flag weights;
specializing three distinct constants is an equivariant diagnostic, not the
nonequivariant limit.

## Descendant conventions

CJR's `psi_i` in (8.21) is the cotangent line on the coarse log-GLSM curve
before stabilization.  An ordinary elliptic-curve GW descendant instead uses
the stable-map cotangent line.  The code keeps these distinct:

- `psi_convention="log"` uses (8.21);
- `psi_convention="stabilized"` restricts the pullback of the stable-map
  class; on a contracted marked tail it becomes `contact_psi` on the adjacent
  infinity marking;
- `effective_basis_only=True` rejects `contact_psi`, so its unknowns have
  only evaluation insertions and `psi_min`, as in Theorem 10.1.

Equation (9.7) identifies virtual cycles but does not identify the two
pre-push-forward psi classes.  CJR I, equation (1.10) and Theorem 1.4,
nevertheless give a numerical virtual comparison in the compact zero sector
after the stabilization-boundary terms are included.  The code now packages
that result in `StabilizationBoundaryComparison`: Bloch--Okounkov evaluates a
separate stabilized probe and the graph compiler retains the log-domain
probe.  No class-level identification is made.

## Defects fixed by this audit

1. The graph front end admitted the globally unstable connected type
   `(g,n,D)=(1,0,0)`.  `ProbeSpec`, the known-value backend, and graph
   compilation now reject it, and the enumerator returns no loci.
2. Punctured-axiom relation metadata multiplied a positive infinity degree by
   three even though both `EffectiveVertex` and `ProbeSpec` already store the
   ambient `P^2` degree.  The duplicate conversion was removed.
3. The conventional Euler-twisted I-function (`s+3H`) is now documented as
   distinct from the dual Euler class (`t-3H`) on CJR stable zero vertices.
4. The effective-basis probe family formerly omitted the numerical
   stabilization-boundary comparison and therefore used only primary Chow
   rows.  It now also contains log-domain stationary rows with rigorously
   transferred compact values.  Orchestrator checkpoint version 9 prevents
   a version-8 stage from being mistaken for the enlarged relation family.

## Verification and limits

The regression suite checks ordinary `P^2` invariants, sparse versus standard
equivariant lifts, string and dilaton equations, constant-map Hodge formulas,
the full `O(3)` fixed-locus evaluator, Givental--Teleman values against direct
localization, CJR graph factors, punctured axioms, exact block solving, and
checkpoint behavior.

For the degree-zero genus-three target
`V(3,0,(-5); psi_min=0, H)`, the primary one-point `t^-4` relation has 85
fixed loci, 50 collected infinity monomials, and target coefficient
`-625/72`.  There is one locus with that exact `(-5)` topology but 12 loci
containing a genus-three infinity vertex.  Thus the relation is not a
one-unknown equation.  The saved effective-basis calculation remains
underdetermined.  The new boundary-compared one-point descendant row has 85
fixed loci, 53 collected monomials, target coefficient `-625/72`, compact
side `-31/967680`, and no `contact_psi`.  Added to the version-8
`effective-basis-d2-m2` system, it raises the target component rank from 278
to 279 but one kernel direction still moves the target.
rank-deficient in one direction involving this target; reporting it as
undetermined is the mathematically correct behavior, not a computational
failure that may be replaced by a guessed value.

The implementation is finite for any configured truncation, but it is not a
proof that localization alone determines every effective invariant.  Larger
probe families, additional punctured relations, or independent basic
effective-invariant input may still be necessary.
