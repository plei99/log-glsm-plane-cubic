# Plane-cubic Gromov-Witten invariants from R-localization

Exact SageMath code for logarithmic-GLSM localization of a plane cubic,
including contact-resolved CJR graph enumeration, full
`O(3)`-twisted zero vertices, tautological integration, and reconstruction
from known elliptic-curve Gromov-Witten invariants.

## Quick start

The Sage modules are tested with SageMath 10.7. The required `admcycles` 1.4
package is vendored in this repository, so no separate Python installation is
needed. Run commands from the repository root:

```bash
./run_tests.sh                # every regression suite, with timings
FAST=1 ./run_tests.sh         # skip the slow localization suites
sage cjr_full_equation_provider.sage --max-degree 8
```

Long orchestration runs should always carry their caches and a wall-clock
budget; running without them is the known-expensive path:

```bash
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 --max-psi-min 0 \
  --max-markings 2 --effective-basis-only \
  --zero-vertex-cache results/genus3/zero-vertices.sqlite \
  --hodge-cache results/genus3/hodge.sqlite \
  --time-budget 3600 --checkpoint-out genus3.json
```

The same `--time-budget` flag now exists on
`log_glsm_infinity_orchestrator.sage` itself; a budgeted run stops cleanly
between blocks with checkpointable state and reports `timed_out`.  Compiled
relations locked inside version-8 checkpoints can be recovered with
`--import-v8-relations PATH`, which keeps only Chow-scalar localization rows
and discards the untrustworthy v8 stage bookkeeping.

For a step-by-step guide to aggregate and contact-resolved infinity vertices,
see [Computing infinity vertices](docs/COMPUTING_INFINITY_VERTICES.md).
The paper-to-code formula cross-check and exact support boundary are recorded
in [Correctness audit](docs/CORRECTNESS_AUDIT.md).

### Audited scope

The contact-resolved pipeline implements CJR's compact-sector hypersurface
specialization for connected, globally stable data.  It now rejects the
exceptional unmarked constant genus-one type `(g,n,d)=(1,0,0)` before graph
enumeration or known-value lookup.  Positive degrees stored by `ProbeSpec`
and `EffectiveVertex` are always ambient `P^2` degrees; only the
Bloch--Okounkov series converts a multiple of three to intrinsic elliptic
degree.

There are two exact but different descendant workflows.  A `log` probe uses
CJR III (8.21).  A `stabilized` probe pulls an ordinary stable-map cotangent
class through `st`; contracted marked tails then produce `contact_psi`.
The numerical stabilization-boundary comparison of CJR I, equation (1.10),
now transfers known stabilized cubic descendants to separate `log` probes;
it does not identify the two line classes.  This is the descendant workflow
used by `effective_basis_only=True`.
For the effective-invariant basis of Theorem 10.1, use
`effective_basis_only=True`, which permits only evaluation insertions and
`psi_min` at infinity.  A blocked rank report means that the chosen CJR
relations do not determine the requested basic effective invariant; the code
does not assign a conjectural value.

`elliptic_cubic_gw.py` computes the connected, unmarked, positive-degree
genus-one Gromov-Witten invariants of a smooth cubic

\[
E=(3)\subset \mathbb P^2.
\]

It uses exact rational arithmetic and separates two degree conventions:

- `cover_degree = r` is the degree of the isogeny \(C\to E\);
- `ambient_degree = d=\int_C f^*H=3r` is the stable-map degree in
  \(\mathbb P^2\).

Consequently,

\[
N_{1,d}(E)=
\begin{cases}
\sigma_1(r)/r,&d=3r,\\
0,&3\nmid d.
\end{cases}
\]

## R-equivariant twisted plane theory

Let \(\lambda_0,\lambda_1,\lambda_2\) be the fixed-point restrictions of
\(H\), let \(s\) be the character scaling the fiber of \(\mathcal O(3)\),
and let \(z\) be the loop parameter. `O3TwistedPlane` evaluates

\[
I_{i,d}(z;s)=
\frac{\prod_{m=1}^{3d}(s+3\lambda_i+mz)}
{\prod_{j=0}^2\prod_{m=1}^{d}(\lambda_i-\lambda_j+mz)}.
\]

This is the fixed-point restriction of the
\((\mathbb C^*)^3\times\mathbb C_R^*\)-equivariant
\(\mathcal O(3)\)-twisted \(I\)-function. The same class on a degree-\(d\)
invariant edge from \(p_i\) to \(p_j\) has section weights

\[
s+\frac{3d-k}{d}\lambda_i+\frac{k}{d}\lambda_j,
\qquad 0\le k\le 3d.
\]

For example:

```bash
python3 elliptic_cubic_gw.py twisted-i \
  --max-degree 3 \
  --fixed-point 0 \
  --base-weights 0 1 3 \
  --fiber-weight 5 \
  --z 2
```

The fiber parameter is called `s` here to avoid confusing it with the
logarithmic mirror coordinate used below.

## Reduced log-GLSM localization

For \(X=\mathbb P^2\) and \(\mathcal E=\mathcal O(3)\), equation (9.7) of
Chen-Janda-Ruan III specializes the reduced log-GLSM class at R-weight zero:

\[
\left.\operatorname{st}_*[\mathcal U_{g}(\mathfrak P,d)]^{\rm red}
\right|_{s=0}
=(-1)^{1-g+3d}\,
i_*[\overline M_g(E,d)]^{\rm vir}.
\]

The stable zero-level vertices in their localization formula carry the
\(\mathcal O(3)\)-twisted class, while the infinity vertices provide the
effective/reduced part. The word "reduced" here refers to the log-GLSM
obstruction theory. Zinger's reduced genus-one GW class is a different
construction; for a curve target with no descendants, its invariant equals
the standard genus-one invariant. In this case, the resummed fixed-locus
series is

\[
F_1(q)=\frac18(T-\log q)
-\frac1{24}\log(1-27q)
-\frac12\log I_0(q),
\]

where

\[
I_0(q)=\sum_{d\ge0}\frac{(3d)!}{(d!)^3}q^d,
\qquad
T-\log q=\frac{I_1(q)}{I_0(q)},
\]

and

\[
[q^d]I_1(q)=3\frac{(3d)!}{(d!)^3}(H_{3d}-H_d).
\]

After the mirror change of variables
\(Q=e^T=q\exp(I_1/I_0)\), exact series reversion gives

\[
F_1(Q)=-\sum_{r\ge1}\log(1-Q^{3r}).
\]

The program reports both the reduced log-GLSM coefficient (including the
sign above) and the recovered GW invariant.

### Descendant convention

There are two different psi classes in this calculation. CJR III inserts the
cotangent line on the coarse log-GLSM universal curve before the stabilization
map `st`. Bloch--Okounkov computes the ordinary stable-map descendant after
stabilization. These classes differ when `st` contracts a marked zero tail.
Equation (9.7) identifies the pushed-forward virtual cycles, but it does not
identify those two pre-push-forward descendant insertions.

For compact zero-sector probes, CJR I's virtual comparison supplies the
missing boundary correction at the level of numerical invariants.  The code
implements this in `StabilizationBoundaryComparison`: it evaluates a
stabilized copy with Bloch--Okounkov, while compiling the original probe with
`psi_convention="log"`.  Hence CJR (8.21), rather than `contact_psi`, appears
on the graph side.  Reports retain both probe records so this transfer is
auditable and cannot be confused with equality of cotangent classes.

For example, the effective-basis equation associated to the known
genus-three one-point descendant is obtained by

```sage
load("cjr_full_equation_provider.sage")

provider = PlaneCubicFullEquationProvider()
stabilized = ProbeSpec.stationary(3, 0, (4,))
relation = provider.boundary_compared_relation(stabilized)
```

Here `relation.probe.psi_convention == "log"`, its known side is
`-31/967680`, and every infinity factor has `contact_psi=()`.

The code records this distinction as `ProbeSpec.psi_convention` and implements
both cases. On a stable zero component the two restrictions agree. On an
unstable marked zero tail of contact `c`, the log-domain class contributes
`3H/c-t`, while the pullback of the stable-map class becomes the ordinary
cotangent line at the adjacent infinity contact. `EffectiveVertex.contact_psi`
records those powers. Primary probes are unchanged.

## Genus two from log-GLSM localization

For an infinity vertex of genus $h\leq2$ and ambient plane degree $D$,
the hypersurface balancing equation in CJR I, equation (7.2), becomes

\[
3D-(2h-2)=\sum_E(c(E_\infty)+1).
\]

Every infinity contact satisfies $c(E_\infty)\leq-1$, so the right hand
side is nonpositive. For $D>0$ and $h\leq2$, the left hand side is
positive. Thus **every infinity vertex appearing in the genus-two
calculation has degree zero**.

`GenusTwoLogGLSMLocalization` separates the graph sum into the positive-degree
zero-level blocks

\[
A_2(q)=\sum_{d>0}\sigma_1(d)q^d,
\qquad
A_4(q)=\frac1{12}\sum_{d>0}\sigma_3(d)q^d,
\]

and the degree-zero infinity Hodge weights

\[
b_2=-\frac1{24},\qquad b_4=\frac1{2880}.
\]

The latter are calculated from

\[
\log\!\left(\frac{z/2}{\sinh(z/2)}\right)
=\sum_{h\geq1}\frac{-B_{2h}}{2h(2h)!}z^{2h}.
\]

The four specialized CJR graph types give

\[
A_4+\frac12A_2^2+b_2A_2+\left(b_4+\frac12b_2^2\right).
\]

In particular, the degree-zero term is

\[
b_4+\frac12b_2^2=\frac7{5760}.
\]

Run this calculation with:

```bash
python3 elliptic_cubic_gw.py log-glsm-g2 --max-degree 8
python3 elliptic_cubic_gw.py log-glsm-g2 --max-degree 8 --json
```

This is a specialized numerical reconstruction of the CJR graph sum. It
does not attempt to construct the punctured R-map moduli spaces themselves.
The program checks its result against the independent stationary formula
described next.

### Reconstructing the infinity vertices from known GW theory

The Sage module `log_glsm_infinity_vertices.sage` performs the inverse
calculation requested by the virtual localization formula. It loads the
genus-one and genus-two series from `bo_coefficient.sage`, obtaining

\[
C_{[0]}=E_2,
\qquad
C_{[2]}=\frac12E_2^2+\frac1{12}E_4,
\]

and first solves the genus-one localization equation

\[
C_{[0]}=A_2+b_2.
\]

It then computes the positive-degree zero-level O(3)-twisted contribution

\[
A_4+\frac12A_2^2,
\]

and solves

\[
C_{[2]}=A_4+\frac12A_2^2+b_2A_2+
\left(b_4+\frac12b_2^2\right)
\]

The genus-one residual is supported only in degree zero and gives
$b_2=-1/24$. In genus two, after removing $b_2A_2$, the residual again has
only degree zero and gives $b_4=1/2880$.

```bash
sage log_glsm_infinity_vertices.sage --max-degree 8
sage log_glsm_infinity_vertices.sage --max-degree 8 --json
```

The balance equation also enumerates the allowed degree-zero profiles. The
basic genus-one profile is `(-1)`. At genus two, the primitive profiles begin
with `(-3)` and `(-2,-2)` (with possible additional `-1` contacts). The
one-point resummed equation determines their aggregate primitive contribution
$b_4$. The full R-equivariant compiler keeps these profiles separate.
Both the constant Bloch--Okounkov rows and the nonconstant Laurent rows now
use the stabilization-corrected contact-descendant state space; see
[`docs/COMPUTING_INFINITY_VERTICES.md`](docs/COMPUTING_INFINITY_VERTICES.md).

### Dynamic programming for general infinity vertices

`log_glsm_infinity_dp.sage` implements the triangular algebra needed for a
general reconstruction. An unknown is indexed by

\[
(g,D,\mathbf c,k,\boldsymbol\alpha,\mathbf b),
\]

where $D$ is ambient degree, $\mathbf c$ is the negative contact profile,
$k$ is the power of $\psi_{\min}$, $\boldsymbol\alpha$ records evaluation
insertions, and $\mathbf b$ records ordinary cotangent powers at the contact
markings. Keeping only $(g,n,D)$ is insufficient because different contact
profiles, cohomology insertions, and descendants can mix in one equation.

For each requested vertex, an equation provider must construct a probe
correlator and enumerate the finite CJR graph sum in the form

\[
\mathrm{GW}(\text{probe})=
\mathrm{Tw}_{\mathcal O(3)}(\text{probe})+
\sum_m a_m\prod_{V\in m}\mathrm{Eff}(V).
\]

The DP engine then:

1. recursively solves and substitutes every strictly lower vertex;
2. collects the coefficient of the requested vertex;
3. solves the resulting linear equation and memoizes it;
4. solves a finite matrix block when several same-rank contact profiles mix.

The block step is important: a scalar lexicographic recursion can depend on
an arbitrary tie-breaking order and fail even when the full family of probes
is invertible. The implementation raises an error for non-triangular,
nonlinear, or singular probe systems instead of silently assigning a value.

The generic engine and its genus-two Bloch--Okounkov base-case demonstration
can be loaded with

```sage
load("log_glsm_infinity_dp.sage")
demo = genus_two_bo_demo(8)
print(demo["values"])
```

The full O(3)-twisted descendant classes and gluing pairings are attached to
the enumerated graphs. `log_glsm_infinity_orchestrator.sage` constructs a
finite dependency closure, expands stabilization-corrected probe families,
selects exact full-rank blocks, back-substitutes solved values, and checkpoints
the extracted relations. The former determinant-9 witness is still withdrawn:
although its psi convention is now fixed, its old matrix omitted the new
contact-descendant columns. Infinite computing power addresses fixed-graph
growth, but not a genuinely rank-deficient probe matrix; blocked runs report
the exact kernel.

There is a second, mathematically different use of string and dilaton. The
module `cjr_punctured_axioms.sage` implements the unit-removal axioms of
[Chen--Janda--Ruan, *Punctured logarithmic R-maps*](https://arxiv.org/abs/2208.04519)
for the *punctured infinity vertices themselves*. In the plane-cubic theory
the unit sector has contact `-1`, and

\[
\psi_{\mathrm{DF}}=c_1(\mathcal O(-\infty))=-3H.
\]

Writing \(V^{(k)}_{\mathbf c}(\alpha_1,\ldots,\alpha_n)\) for an effective
vertex with `psi_min=k`, the implemented scalar forms are

\[
\begin{aligned}
V^{(k)}_{\mathbf c\cup(-1)}(\boldsymbol\alpha,1)
 &=\sum_j |c_j|\sum_{a=0}^{k-1}
 V^{(k-1-a)}_{\mathbf c}
 (\alpha_1,\ldots,\alpha_j(-3H)^a,\ldots,\alpha_n),\\
V^{(k)}_{\mathbf c\cup(-1)}(\boldsymbol\alpha,H)
 &=D V^{(k)}_{\mathbf c}(\boldsymbol\alpha)
 +\sum_j |c_j|\sum_{a=0}^{k-1}
 V^{(k-1-a)}_{\mathbf c}
 (\alpha_1,\ldots,\alpha_jH(-3H)^a,\ldots,\alpha_n),\\
V^{(k+1)}_{\mathbf c\cup(-1)}(\boldsymbol\alpha,1)
 &+3V^{(k)}_{\mathbf c\cup(-1)}(\boldsymbol\alpha,H)
 =(2g-2+n)V^{(k)}_{\mathbf c}(\boldsymbol\alpha).
\end{aligned}
\]

The first two are CJR equations (8.3) and (8.4); the third is their dilaton
combination, Corollary 8.3. The empty sum at `k=0` gives the fundamental-class
vanishing for a unit insertion, while the divisor equation becomes
`V(...,H)=D*V(...)`. Terms containing \(H^3\) are discarded. The provider also
loads CJR's genus-zero vanishing, genus-one emptiness away from degree zero
with all contacts `-1`, and the two base values

\[
V^{(0)}_{1,0,(-1)}(H)=0,\qquad
V^{(1)}_{1,0,(-1)}(1)=\frac18.
\]

These punctured axioms are enabled by default and are applied before expensive
compact probes. Use `--no-punctured-axioms` only for a legacy localization-row
comparison. Product relations for disconnected infinity contributions are
already represented by graph monomials. No splitting axiom is invented: CJR
explicitly leaves such a formula to future work.

For example, the contact-resolved genus-one base vertex is now returned
directly from the CJR genus-one formula by

```bash
sage log_glsm_infinity_orchestrator.sage \
  --genus 1 --contacts=-1 --insertions 1 \
  --max-markings 1 --t-powers 0 --no-unit-relatives
```

This prints `0`. Replacing `--psi-min 0 --insertions 1` by
`--psi-min 1 --insertions 0` prints `1/8`.

### All genus-three infinity vertices

`genus_three_infinity_vertices.sage` enumerates the complete genus-three
target family and drives one shared orchestration over it. Genus three is the
first genus in which the balance equation

\[
3D-(2g-2)=\sum_i(c_i+1),\qquad c_i\leq-1,
\]

permits a positive infinity degree: it forces \(3D\leq 2g-2\), so
\(D=0\) for \(g\leq2\) but \(D\in\{0,1\}\) for \(g=3\). The degree-zero sector
has contact excess \(4\) and the new degree-one sector has contact excess
\(1\):

```text
D=0:  (-5), (-5,-1), (-4,-2), (-3,-3), (-5,-1,-1), (-4,-2,-1), ...
D=1:  (-2), (-2,-1), (-2,-1,-1), ...
```

Because a contact \(c=-1\) contributes \(c+1=0\), the family is infinite
without a truncation; `--max-valence` bounds it. Within a fixed valence the
enumeration is exactly the set of vertices the graph compiler can emit:
insertions run over \(\{1,H,H^2\}\), and reduced virtual dimension selects
the unique nonnegative power

\[
\psi_{\min}=\operatorname{valence}+3D-\sum_i\operatorname{codim}(\alpha_i).
\]

Assignments with a negative required power vanish. `--max-psi-min` discards
targets above the requested cap; it does not replace the required power by a
smaller one.

```bash
sage genus_three_infinity_vertices.sage --list-only --max-valence 2
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 --max-psi-min 0 \
  --max-markings 1 --effective-basis-only \
  --max-chow-unit-insertions 0 \
  --t-powers 0 -1 -2 -3 -4 -5 -6 \
  --zero-vertex-cache genus3-o3-zero-vertices.sqlite \
  --hodge-cache genus3-o3-hodge.sqlite \
  --time-budget 600 \
  --checkpoint-out genus3-degree-zero.json
```

The provider tracks the ordinary Chow dimension of every Laurent
coefficient. If a probe has residual dimension `delta`, its `t^k`
coefficient has dimension `delta+k`: it is a zero-cycle when this is zero,
is a forced-vanishing numerical row when it is negative, and remains a Chow
class requiring an additional test class when it is positive. Thus extending
`--t-powers` adds scalar rows only after this grading check.
`--max-unit-insertions 2` adds canonical mixed unit/point probes
whose compact invariants are reduced to stationary theory by repeated string
and dilaton equations. This compact probe depth is independent of the
punctured contact-`-1` axioms, which are applied automatically at every stage.
Unit depth is staged (`0`, then `1`, then `2`) so a smaller family is always
attempted first.

For the hypersurface reconstruction of Section 9.3 and Theorem 10.1, use
`--effective-basis-only`.  This mode combines primary Chow probes with
stationary log-domain probes whose known sides are supplied by the CJR-I
stabilization-boundary comparison.  It rejects a compiled row if any infinity
factor contains an ordinary cotangent class at a contact.  Its infinity
unknowns are therefore exactly evaluation insertions and powers of
`psi_min`, as in (10.1).  The broader direct stabilized-descendant mode is
retained as a diagnostic compiler.

The ambient log-GLSM problem exists in every plane degree, and the compact
hypersurface side is zero when that degree is not divisible by three. Extra
probe degrees can be requested with `--additional-probe-ambient-degrees`.
The genus-one rows in ambient degrees one, two, and three are exact regression
checks: after substituting the CJR genus-one infinity theory, their complete
sums are `0`, `0`, and `-1`. The exact solved-relation check remains enabled
for every larger run.

All targets share a single `InfinityVertexOrchestrator`, so the common
genus-one and genus-two closure and the compiled probe relations are built
once rather than once per target. Unsolved targets are reported as
`UNDETERMINED` alongside the exact frontier; no value is ever assigned by a
tie-break.

This is a driver, not a new mathematical shortcut. Historical profiling on
2026-08-07 preceded the exact infinity virtual-dimension filter and therefore
overcounted the closure by treating dimension-incompatible insertions as
independent unknowns. Those old rank and nullity figures are not mathematical
statements about the corrected system. In particular, for the one-contact
degree-zero vertex `(g,D,c)=(3,0,-5)` with `psi_min=0`, only the `H` insertion
survives; the `1` and `H^2` variants are pruned before orchestration.

A corrected Chow-graded retry gives the target linear incidence. The
primary point probe has defect four, so `t^-4` is its first numerical
coefficient and contains
`(-625/72) V(3,0,(-5); psi_min=0,H)`. The complete row has 50 infinity terms;
retaining all of them removes the former genus-one contradiction.

The clean effective-basis run through primary two-marking probes in ambient
degrees zero, one, and two is saved as
`results/genus3/effective-basis-d2-m2.json`.  It compiles 1259 relations,
activates 1205, and proves 349 values.  The final target component has rank
278 on 343 columns.  In the computed kernel basis, exactly one of its 65
directions moves
`V(3,0,(-5);psi_min=0,H)`; that direction has support on 18 effective
invariants, all with empty `contact_psi`.  Thus the current obstruction is a
basic effective-invariant direction, not a missing ordinary-contact-psi
backend.  This also reflects the direction of Theorem 10.1: it computes
hypersurface GW invariants *in terms of* effective invariants; it does not
claim that localization alone determines every effective invariant.  Negative
Laurent coefficients provide additional relations, as noted in Remark 10.7,
but a rank computation must decide whether they determine a chosen basic
invariant.

Coefficients with `delta+k>0` are still genuine positive-dimensional Chow
classes. The numerical backend rejects them and reports the missing test-class
codimension; arbitrary tautological pairings remain a future extension.

Positive infinity degree remains harder, but positive ambient probe degrees
are now usable. The former genus-one contradiction was a truncation bug: an
equivariant O(3) Euler class may supply powers of `t`, so stable zero-side
flag descendants must be bounded by the ordinary ambient stable-map virtual
dimension, not the dimension after subtracting the Euler rank. Work degree
zero first and checkpoint it, then add ambient degrees one through three. A
time budget is checked between stages and solve blocks; it cannot interrupt
one graph compilation or one admcycles integral.
Omit `--full-kernel` for large runs: the default compact frontier records exact
rank and nullity without spending minutes printing huge kernel vectors.

The expensive individual zero vertices now have a second, independent
checkpoint. Precompute them once and reuse them in every orchestration:

```bash
sage o3_zero_vertex_precompute.sage \
  --genus 3 --ambient-degree 1 --max-point-markings 1 \
  --cache genus3-o3-zero-vertices.sqlite \
  --hodge-cache genus3-o3-hodge.sqlite --time-budget 3600
```

Every completed exact value is committed separately. SQLite is recommended:
it uses WAL mode, permits multiple deterministic shards to share the same
zero-vertex and Hodge-integral stores, and avoids rewriting the whole cache.
The older JSON format remains readable when the cache name ends in `.json`.
The equation provider evaluates fixed loci along the nonresonant path
`(0,epsilon,epsilon^2)` and sets `epsilon=0` only after the complete
fixed-locus sum. A linear ray such as `epsilon*(0,1,3)` makes two flag weights
cancel at a degree-three nodal locus. Setting the three weights to fixed
distinct numbers is also not the nonequivariant limit while the R-weight `t`
remains. To keep cache conventions distinct, passing `foo.sqlite` as
`--zero-vertex-cache` uses
`foo.weights-nonequivariant-quadratic.sqlite`; the precompute command applies
the same derivation, so pass the same base cache name to both commands.
Insertion scales and provenance labels are factored out of persistent keys,
marking permutations are canonicalized, and psi profiles are sorted before
admcycles; for the one-mark `(3,1,3)` probe this reduced the request inventory
from 1,510 to 625 distinct zero vertices. Request discovery is separately
stored in `CACHE.inventory.json`, so a resumed precompute does not repeat
graph enumeration.

Degree-zero requests use the direct constant-map model
`Mbar_(g,n)(P2,0) = Mbar_(g,n) x P2` in every genus. Before admcycles, the
Hodge backend applies Mumford's relations through genus three, the lambda-g
theorem, and the closed top-Hodge-triple formula (with string/dilaton
descendants). On the development machine, the unmarked twisted `(g,d)=(3,0)`
request fell from more than five minutes in an admcycles multiplication to
about 0.07 seconds; the direct formula agreed with the three fixed-locus sum.

A Go prototype makes the isolated zero-type classifier about 18.2x faster,
but cProfile assigns that function only 27% of the `(3,0,3)` enumeration.
That predicts roughly a 1.34x end-to-end graph-enumeration improvement, while
it does not accelerate admcycles. The benchmark and reproducible profiler are
in [`benchmarks/`](benchmarks/README.md); a wholesale rewrite is therefore not
the current bottleneck-removal strategy.

A staged implementation plan, including component interfaces and completion
criteria, is in [`docs/ROADMAP.md`](docs/ROADMAP.md).

The roadmap's computational architecture is implemented in the following
modules:

- `log_glsm_conventions.sage`: `ProbeSpec`, elliptic-to-ambient insertion
  conversion, exact equivariant rings, dimensions, signs, and truncations;
- `cjr_graph_factors.sage`: edge, unstable-zero, diagonal, stable-flag, and
  infinity-descendant factors with `H_infinity=3H`;
- `hodge_integrals.sage`: all-genus DVV psi intersections, the lambda-g and
  top-Hodge-triple formulas, low-genus Mumford reductions, a concurrent exact
  SQLite cache, and arbitrary remaining lambda/psi monomials through vendored
  `admcycles`;
- `o3_twisted_plane_vertices.sage`: request objects, the lightweight partial
  backend, twisted I-function restrictions, and resummed stationary blocks;
- `o3_fixed_locus_graphs.sage`: full decorated stable-map fixed graphs on
  `P^2`, internal automorphisms and deck groups, tautological expansions, and
  the complete `FullTwistedZeroVertexBackend` with a persistent exact cache;
- `givental_teleman.sage`: truncated R-matrices and the unit-preserving
  Givental stable-graph action with exact DVV vertex integrals;
- `o3_givental_teleman.sage`: the degree-zero QRR initial condition, calibrated
  positive-q quantum connection, R- and S-matrices, and hybrid backend;
- `o3_zero_vertex_precompute.sage`: request inventory, resumable batch
  evaluation, and cache prepopulation for expensive twisted zero vertices;
- `cjr_graph_contributions.sage`: auditable per-graph compilation into
  contact-resolved `EffectiveVertex` polynomials;
- `cjr_infinity_chow.sage`: residual Chow-degree tracking and the numerical
  projection guard for zero-cycle and forced-vanishing Laurent rows;
- `elliptic_probe_values.sage`: stabilized Bloch--Okounkov values with string
  and dilaton reduction, plus guards against use as CJR log-psi values;
- `cjr_probe_factory.sage`: incremental exact row-echelon screening, diagonal
  matrices, pivot selection, and kernel reporting;
- `cjr_full_equation_provider.sage`: coefficient extraction, field-valued DP,
  convention-safe reports, and the aggregate genus-two end-to-end command;
- `log_glsm_infinity_orchestrator.sage`: finite dependency closure, staged
  probe expansion, exact block solving, frontier reports, and checkpoints;
- `genus_three_infinity_vertices.sage`: positive-degree balance enumeration,
  the complete genus-three target family, and a shared orchestration driver.

Run the aggregate genus-two computation with

```bash
sage cjr_full_equation_provider.sage --max-degree 8
sage cjr_full_equation_provider.sage --max-degree 8 --json
```

The compiler never silently discards geometry. Its default backend now
computes arbitrary finite higher-genus and positive-degree individual
twisted zero vertices by base-torus localization. The former partial backend
can still be injected explicitly when testing failure semantics.

### Full O(3)-twisted zero vertices

For an internal fixed vertex of colour `i`, genus `g`, and internal edge
valence `r`, put

\[
w_{ij}=\lambda_i-\lambda_j,\qquad W_i=t-3\lambda_i.
\]

The stable vertex integrand implemented in
`o3_fixed_locus_graphs.sage` is

\[
\prod_{j\ne i}\Lambda_g^\vee(w_{ij})w_{ij}^{r-1}
\prod_F\frac1{\omega_F-\psi_F}\cdot
\frac{W_i^{1-r}}{\Lambda_g^+(W_i)}.
\]

The final factor is the normalization of
`e((R pi_* f^* O(3))^vee)`. In genus one and degree zero it is precisely
`(t-3H)/(t-3H+lambda_1)`, as in CJR III, Example 9.2. Each invariant-line
edge carries the standard stable-map virtual-normal factor and all `3d+1`
section weights of `O(3d)`. The global fixed-locus denominator is
`|Aut(Gamma)| * product_e d_e`. Nodal unstable vertices receive the extra
factor `W_i^-1`, and a psi class on an unstable marked endpoint restricts to
`-omega`.

All tautological rational functions are expanded only through the complex
dimension of the relevant `Mbar_(g,n)`. Equivariant coefficients stay in
`QQ(lambda_0,lambda_1,lambda_2,t)` while `admcycles` sees only individual
lambda/psi monomials. This separation avoids symbolic-coefficient ambiguity
inside the tautological ring.

The full evaluator uses sparse equivariant lifts by default. For a divisor it
cycles through

\[
H-\lambda_i,
\]

which vanishes at `p_i`. For a point class it cycles through

\[
\delta_i=\prod_{j\ne i}(H-\lambda_j),
\]

which is supported only at `p_i`. These specialize respectively to `H` and
`H^2` when the auxiliary base weights vanish. Marking placements at which a
chosen lift restricts to zero are rejected before fixed graphs are built.
Thus the degree-two five-point calculation uses supports `(0,1,2,0,1)` and
has only three fixed graphs, one for each possible middle colour, instead of
1557. Both calculations give `1`:

```sage
load("o3_fixed_locus_graphs.sage")
rings = PlaneCubicCoefficientRing()
request = TwistedZeroVertexRequest(
    0, 2, tuple(TwistedInsertion(2) for _ in range(5))
)

fast = P2FixedLocusEvaluator(rings, include_twist=False)
assert fast.graph_count(request) == 3
assert fast.evaluate(request) == 1

reference = P2FixedLocusEvaluator(
    rings, include_twist=False, lift_strategy="standard"
)
assert reference.graph_count(request) == 1557
assert reference.evaluate(request) == 1
```

Explicit choices are available as
`TwistedInsertion.H_vanishing_at(i)` and
`TwistedInsertion.point_supported_at(i)`. The planner preserves explicit
lifts. Use `lift_strategy="standard"` to disable automatic sparsification.

### Givental--Teleman backend

The second backend is implemented without replacing localization. The generic
engine applies

\[
R^{-1}(\psi),\qquad
\frac{\eta^{-1}-R^{-1}(\psi')\eta^{-1}R^{-1}(\psi'')^T}
     {\psi'+\psi''},\qquad
T(z)=z(1-R^{-1}(z)1)
\]

to legs, edges, and translation markings of every stable graph. Its vertices
are semisimple TFT vertices, so their integrals are pure psi intersections
rather than stable-map or admcycles calculations.

At degree zero, the three auxiliary fixed points give an exact semisimple
algebra with

\[
\eta_i=\frac{t-3\lambda_i}
{\prod_{j\ne i}(\lambda_i-\lambda_j)}.
\]

Mumford GRR determines its diagonal R-matrix. The implementation agrees with
the localization backend in genera zero, one, and two. The raw Picard--Fuchs
principal symbol is retained as a diagnostic by Hensel lifting the roots of

\[
\prod_i(p-\lambda_i)=q(t-3p)^3.
\]

Its residue pairing varies with (q), so the code does not mistake this raw
symbol for flat Frobenius data. The production calibration instead computes
the exact small quantum product in the flat basis (1,H,H^2) from genus-zero
three-point twisted invariants. After q-adic diagonalization, it solves

\[
[P,R_{k+1}]=D R_k+A R_k-R_k B,
\qquad D=q\frac{d}{dq},
\]

where (P) is multiplication by (H), (A=E^{-1}DE), and
(B=\frac12\operatorname{diag}(D\log\eta_i)). The diagonal integration
constants are fixed by the degree-zero quantum-Riemann--Roch matrix.

Teleman reconstructs ancestors. Stable-map descendants are recovered from

\[
D S_{m+1}=C_H S_m-S_m C_H(0),\qquad S_0=1,
\quad
\tau_k(\gamma)=\sum_{m=0}^k\bar\tau_{k-m}(S_m\gamma).
\]

Primary and descendant positive-degree coefficients are checked exactly
against virtual localization in degrees one and two. Use the calibrated path
with:

```bash
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 \
  --zero-vertex-strategy hybrid --time-budget 600
```

`--zero-vertex-strategy givental` uses the calibrated reconstruction directly.
`hybrid` uses it where supported and keeps localization as a fallback. The
default remains `localization` while larger positive-degree workloads are
profiled.

Run the focused checks with

```bash
DOT_SAGE=/tmp/codex-sage sage test_givental_teleman.sage
DOT_SAGE=/tmp/codex-sage sage test_hodge_integrals.sage
DOT_SAGE=/tmp/codex-sage sage test_o3_fixed_locus_graphs.sage
DOT_SAGE=/tmp/codex-sage sage test_cjr_graph_contributions.sage
DOT_SAGE=/tmp/codex-sage sage test_cjr_full_equation_provider.sage
```

### Decorated bipartite graph enumerator

`cjr_bipartite_graphs.sage` completes the finite combinatorial layer for the
smooth plane-cubic, trivial-sector specialization with compact-type ordinary
markings. A graph records

- zero- and infinity-level vertex genera and ambient degrees;
- the zero vertex carrying each labeled ordinary marking;
- every edge endpoint and positive contact order;
- stable, nonspecial, marked, and nodal zero-vertex types;
- the graph genus, infinity balancing equations, and automorphism order.

The search enumerates connected bipartite multigraphs, distributes genus and
degree, solves the positive contact compositions independently at each
infinity vertex, and removes isomorphic copies by a colored incidence-graph
canonical label. Edge instances are incidence vertices, so permutations of
parallel equal-contact edges are included automatically in
`automorphism_order()` and `localization_weight()` returns
`1/|Aut(Gamma)|`.

For mixed-level graphs the stability and balance constraints give the finite
conservative bound

\[
|V(\Gamma)|\leq 4g+2n+d-1.
\]

That bound is very loose — genus-two graphs never exceed four zero vertices
while it permits ten — so the search is organized to avoid paying for it.

*Incidence matrices are generated up to relabeling of the zero vertices.*
Row `z` of the matrix records how many edges join zero vertex `z` to each
infinity vertex, and `_canonical_incidence_matrices` emits rows in
non-increasing `(sum, row)` order, one representative per orbit. This loses
no isomorphism class: the permutation sorting an arbitrary matrix's rows can
be applied to that graph's genus, degree, and marking decorations too, and
every decoration of every emitted matrix is enumerated. Each row sum is
positive, since a zero vertex with no incident edge would disconnect a
mixed-level graph.

*The contact budget is checked before contacts are assigned.* Contacts at an
infinity vertex sum to `2g-2+val-3d`; every edge contributes at least one and
an edge meeting a nonspecial zero vertex contributes at least two. Testing
that necessary condition in `_contact_budget_available` discards
configurations whose budget can never balance, instead of building each of
their contact assignments first.

*Inner-loop decorations are reused.* Zero-vertex types are computed from one
pass over edge and marking valences instead of constructing a provisional
graph for every decoration. Degree compositions are cached by vertex count,
and the graph compiler shares one enumerated family across every probe with
the same `(g,n,d)`. The O(3)-twisted backend likewise shares internal `P^2`
fixed-locus graphs when only insertion coefficients change.

Both reductions are exact. `test_cjr_bipartite_graphs.sage` cross-checks the
canonical generator against the brute-force labeled search, and the resulting
isomorphism classes, automorphism orders, and infinity dependencies are
unchanged. Measured:

| `(g,n,d)` | before | after |
|---|---|---|
| (2,1,0) | 2.5 s | 0.02 s |
| (2,0,3) | 32.7 s | 0.31 s |
| (2,2,0) | 129.5 s | 0.24 s |
| (2,3,0) | ~1.4 h (projected) | 5.5 s |
| (3,1,0) | ~1.2 d (projected) | 0.34 s |
| (3,2,0) | ~95 d (projected) | 3.31 s |
| (3,0,3) | not measured | 6.17 s |
| (3,1,3) | >5 min | 115.58 s |

Genus three is therefore reachable: `(3,0,0)` has 32 classes, `(3,1,0)` has
85, and `(3,2,0)` has 300.

One-level exceptional graphs are checked separately. Run, for example,

```bash
sage cjr_bipartite_graphs.sage 2 0 0
sage cjr_bipartite_graphs.sage 2 1 0
```

The first command reproduces the six graph types in CJR I, Figure 4, with
contact profiles

\[
\varnothing,(1),(1,1),(1,1),(2,2),(3).
\]

The second finds eleven types after adding one labeled compact marking.
Programmatically:

```sage
load("cjr_bipartite_graphs.sage")
graphs = PlaneCubicGraphEnumerator(2, 1, 0).graphs()
for graph in graphs:
    print(graph.signature(), graph.localization_weight())
    print(graph.all_infinity_vertex_data())
```

`all_infinity_vertex_data()` uses the negative CJR convention and returns
triples `(genus, ambient_degree, contact_profile)`. The one-point equation
provider exposes their union as
`combinatorial_infinity_dependencies(genus, ambient_degree)`; descendant
powers and evaluation insertions are added only when a particular graph
contribution is expanded.

This is a complete graph enumerator within the stated compact plane-cubic
scope. It does not yet evaluate a graph's vertex integrals. Ramified target
sectors, nontrivial inertia labels, noncompact ordinary markings, and a
general bundle in place of `O(3)` require extending the decoration and
balance data.

### First concrete equation provider

`cjr_plane_cubic_equation_provider.sage` supplies a genuine provider for the
one-point stationary sector in arbitrary genus. It uses

\[
\sum_{g\geq0}\left\langle\tau_{2g-2}(\mathrm{pt})\right\rangle_gz^{2g}
=\exp\left(\sum_{m\geq1}\frac{2E_{2m}(q)}{(2m)!}z^{2m}\right)
\]

and splits each primitive block as

\[
\frac{2E_{2m}(q)}{(2m)!}=A_{2m}(q)+b_{2m}.
\]

Here $A_{2m}$ is the positive-degree resummed O(3)-twisted block and
$b_{2m}$ is degree-zero infinity data. Integer partitions of $g$ enumerate
the graph monomials in the resulting triangular equation. The provider
solves these equations with `InfinityVertexDP` and checks every requested
$q$-coefficient against `bo_coefficient.sage`.

```bash
sage cjr_plane_cubic_equation_provider.sage \
  --max-genus 4 --max-degree 6
```

The first primitive cumulants are

\[
-\frac1{24},\qquad \frac1{2880},\qquad
-\frac1{181440},\qquad \frac1{9676800}.
\]

The same module implements the contact-resolved diagonal for CJR basic
vertices with only contact `-2`. For

\[
n_2=2g-2-3D,
\]

the $n_2$ edge/unstable-leaf pairs cancel and the dominant star graph has
coefficient $1/n_2!$. For other contact orders the uncancelled moving-weight
factor is retained explicitly; their full equation still requires the
general bipartite graph enumerator and twisted descendant integrations.

## Independent stationary-theory check

The one-point stationary formula for an elliptic curve is

\[
\sum_{g\ge0}\left\langle\tau_{2g-2}(\mathrm{pt})\right\rangle_g z^{2g}
=\exp\left(\sum_{k\ge1}C_{2k}(q)z^{2k}\right).
\]

The virtual dimension is \(3\), exactly the codimension of
\(\mathrm{pt}\,\psi^2\). Although the correlation-function theorem is most
naturally stated in normalized disconnected theory, for one marked
insertion the normalization cancels all unmarked components, so this is the
connected invariant requested here.

It follows that

\[
\left\langle\mathrm{pt}\,\psi^2\right\rangle_{2,1}
=C_4+\frac12C_2^2,
\]

where

\[
C_2=-\frac1{24}+\sum_{d\ge1}\sigma_1(d)q^d,
\qquad
C_4=\frac1{2880}+\frac1{12}\sum_{d\ge1}\sigma_3(d)q^d.
\]

Thus, for \(d>0\),

\[
\left\langle\mathrm{pt}\,\psi^2\right\rangle_{2,1,d}
=\frac{\sigma_3(d)}{12}-\frac{\sigma_1(d)}{24}
+\frac12\sum_{a=1}^{d-1}\sigma_1(a)\sigma_1(d-a).
\]

The degree-zero value is \(7/5760\). The first positive-degree values are

\[
\frac1{24},\quad \frac98,\quad \frac{31}{6},\quad
\frac{343}{24},\quad \frac{117}{4},\quad \frac{111}{2}.
\]

Here \(d\) is intrinsic degree on \(E\). Its degree as a stable map to
\(\mathbb P^2\) is \(3d\). Also
\([\mathrm{pt}]=\frac13i^*H\), so a log-GLSM calculation with ambient
insertions should use \(H/3\). The `stationary-g2` command reports the CJR
reduced-class sign separately.

## Usage

```bash
python3 elliptic_cubic_gw.py invariants --max-cover-degree 8
python3 elliptic_cubic_gw.py invariants --max-cover-degree 8 --json
python3 elliptic_cubic_gw.py log-glsm-g2 --max-degree 8
python3 elliptic_cubic_gw.py stationary-g2 --max-degree 8
```

Python API:

```python
from elliptic_cubic_gw import (
    EllipticCubicLocalization,
    EllipticCurveStationaryTheory,
    GenusTwoLogGLSMLocalization,
    O3TwistedPlane,
)

calculation = EllipticCubicLocalization(max_cover_degree=8)
print(calculation.gw_invariants_by_cover_degree())
print(calculation.localization_terms_flat())
print(calculation.verify())

twisted = O3TwistedPlane(base_weights=(0, 1, 3), fiber_weight=5)
print(twisted.fixed_point_i_series(fixed_point=0, max_degree=3, z=2))
print(twisted.edge_euler_class(0, 1, degree=2))

stationary = EllipticCurveStationaryTheory(max_degree=8)
print(stationary.genus_two_tau2_point_invariants())
print(stationary.verify())

genus_two_glsm = GenusTwoLogGLSMLocalization(max_degree=8)
print(genus_two_glsm.localization_graph_terms())
print(genus_two_glsm.gw_invariants())
print(genus_two_glsm.verify())
```

Run the tests with:

```bash
python3 -m unittest -v test_elliptic_cubic_gw.py
sage test_log_glsm_infinity_vertices.sage
sage test_log_glsm_infinity_dp.sage
sage test_cjr_plane_cubic_equation_provider.sage
```

## Scope

This module computes the numerical positive-degree genus-one primary theory
and the one-point genus-two invariant by the specialized log-GLSM graph
reconstruction above. It does not yet implement arbitrary descendants or
the full all-genus effective-vertex formalism.

## References

- Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged
  Linear Sigma Models I: The Tropical Decomposition Formula*, especially
  Theorem 1.4 and equation (1.10).
- Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged
  Linear Sigma Models III: The Localization Formula*, especially equations
  (9.7) and the hypersurface contributions in Section 9.3.
- A. Zinger, *The Reduced Genus-One Gromov-Witten Invariants of Calabi-Yau
  Hypersurfaces*, equations (0.18)-(0.20), for the closed genus-one
  resummation used by the implementation.
- A. Pixton, *The Gromov-Witten Theory of an Elliptic Curve and Quasimodular
  Forms*, Theorem 3.2.2 and equation (4.3), for the stationary formula and
  normalized Eisenstein series.
