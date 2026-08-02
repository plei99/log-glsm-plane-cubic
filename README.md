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
sage test_o3_fixed_locus_graphs.sage
sage cjr_full_equation_provider.sage --max-degree 8
```

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
one-point equation determines their aggregate primitive contribution
$b_4$; it does not separate the two punctured-vertex integrals individually.

### Dynamic programming for general infinity vertices

`log_glsm_infinity_dp.sage` implements the triangular algebra needed for a
general reconstruction. An unknown is indexed by

\[
(g,D,\mathbf c,k,\boldsymbol\alpha),
\]

where $D$ is ambient degree, $\mathbf c$ is the negative contact profile,
$k$ is the power of $\psi_{\min}$, and $\boldsymbol\alpha$ records evaluation
insertions. Keeping only $(g,n,D)$ is insufficient because different contact
profiles and cohomology insertions can mix in one localization equation.

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

The full O(3)-twisted descendant classes and gluing pairings are now attached
to the enumerated graphs. The remaining reconstruction problem is to select
enough known probes to make each contact-resolved diagonal block nonsingular.
Infinite computing power addresses fixed-graph growth, but not a zero or
rank-deficient probe matrix.

A staged implementation plan, including component interfaces and completion
criteria, is in [`docs/ROADMAP.md`](docs/ROADMAP.md).

The roadmap's computational architecture is implemented in the following
modules:

- `log_glsm_conventions.sage`: `ProbeSpec`, elliptic-to-ambient insertion
  conversion, exact equivariant rings, dimensions, signs, and truncations;
- `cjr_graph_factors.sage`: edge, unstable-zero, diagonal, stable-flag, and
  infinity-descendant factors with `H_infinity=3H`;
- `hodge_integrals.sage`: all-genus DVV psi intersections, the lambda-g
  theorem, and arbitrary lambda/psi monomials through vendored `admcycles`;
- `o3_twisted_plane_vertices.sage`: request objects, the lightweight partial
  backend, twisted I-function restrictions, and resummed stationary blocks;
- `o3_fixed_locus_graphs.sage`: full decorated stable-map fixed graphs on
  `P^2`, internal automorphisms and deck groups, tautological expansions, and
  the complete `FullTwistedZeroVertexBackend`;
- `cjr_graph_contributions.sage`: auditable per-graph compilation into
  contact-resolved `EffectiveVertex` polynomials;
- `elliptic_probe_values.sage`: Bloch--Okounkov values with string and dilaton
  reduction;
- `cjr_probe_factory.sage`: exact diagonal matrices, pivot selection, and
  kernel reporting;
- `cjr_full_equation_provider.sage`: coefficient extraction, field-valued DP,
  reports, and the genus-two end-to-end command.

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

Run the focused checks with

```bash
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
