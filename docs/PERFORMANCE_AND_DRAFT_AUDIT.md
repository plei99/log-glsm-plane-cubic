The fixed six-workload benchmark is **6.35 times faster**: median total
compute time decreased from **33.81 seconds to 5.32 seconds**. All six
mathematical outputs agree exactly with the preserved starting code. A
separate correctness fix changes log-domain descendants: CJR III (8.21)
requires the entire cotangent weight to be divided by contact. Previously,
only its hyperplane term was divided.

This audit used the existing Linux checkout in the OrbStack `agents`
machine, SageMath 10.9, Python 3.14.10, and the vendored admcycles package.
The starting submodule commit was
`e5e591a3d5d087482c1d2d9ee4394f7267176c7a`, including the four existing
uncommitted files listed in [the provenance manifest](../benchmarks/results-2026-11-05/provenance.json).
Those existing changes were preserved. A complete starting snapshot remains
at `/home/patrick/agent-work/benchmarks/log-glsm-20261105/baseline`.
No commit or push was made.

The measured workloads were selected before optimizing. Each case creates
fresh evaluator objects and uses no persistent result cache. Imports and
output serialization are outside the timed region. There were three
alternating before/after process pairs on the same machine. The table gives
per-case medians; the total is the median of each run's summed compute time.

| Workload | Before, seconds | After, seconds | Speedup |
|---|---:|---:|---:|
| CJR graph family `(g,n,D)=(3,0,3)`, 478 classes | 2.3402 | 0.6478 | 3.61x |
| CJR graph family `(3,1,3)`, 2,088 classes | 29.5747 | 3.6822 | 8.03x |
| Dual Euler O(3), genus 3, degree 0, unmarked | 0.3202 | 0.2050 | 1.56x |
| Dual Euler O(3), genus 3, degree 0, `H psi^4` | 0.2497 | 0.1709 | 1.46x |
| Dual Euler O(3), genus 2, degree 1, `H psi^3` | 0.0580 | 0.0281 | 2.06x |
| Complete genus-3 primary relation, `t^-4` coefficient | 0.7428 | 0.6124 | 1.21x |
| Total fixed workload | 33.8057 | 5.3224 | 6.35x |

The O(3) microbenchmarks retain the auxiliary base torus. The complete CJR
relation takes the nonequivariant base limit. The aggregate gain is weighted
heavily toward graph enumeration; it is not a claim that every call, the
whole test suite, or every reconstruction is twice as fast. Raw timing
ranges and every exact result are saved beside
[summary.json](../benchmarks/results-2026-11-05/summary.json).
The comparator checks graph canonical keys and automorphism orders,
relation coefficients, and exact Sage fraction-field equality. Two equal
rational functions can print differently, so their string hashes alone are
not used to establish equality.

The actual iCloud drafts are CJR I and CJR III dated **October 17, 2025**.
The supplementary punctured-R-map source is dated **April 25, 2023**.
The provenance manifest records their filenames and SHA-256 hashes. I read
the relevant text and inspected the rendered CJR III pages 82 and 90,
including the parentheses in (8.21). No private draft was uploaded.

For the plane cubic, write `D` for ambient degree in P2 and set
\(w=t-3H\). The formulas checked against CJR III, Section 9.3, are

\[
\sum_{e\ni v}c_e=2g_v-2+\operatorname{val}(v)-3D_v,
\qquad
\operatorname{Cont}_e=\frac{c_e^{c_e}}{c_e!w^{c_e}},
\]

\[
\operatorname{Cont}_{0,\mathrm{ns}}=\frac{w^2}{c},\qquad
\operatorname{Cont}_{0,\mathrm{nd}}=\frac{cc'}{c+c'}w,
\qquad
\operatorname{Cont}_{\infty}=\frac{t}{-t-\psi_{\min}}.
\]

A nonspecial zero tail requires \(c>1\). Stable zero vertices carry

\[
e\!\left((R\pi_*f^*\mathcal O(3))^\vee\right)
\prod_e\frac{1}{w_e/c_e-\psi_e}.
\]

The reduced comparison has sign \((-1)^{1-g+3D}\), as in (9.7).
An elliptic point insertion is \(H/3\); intrinsic elliptic degree is
\(D/3\) when \(3\mid D\). The existing code respects these conventions.
The punctured string and divisor formulas retain \(|c_j|\),
\(\psi_{\rm DF}=-3H\), and the ambient divisor degree \(D\), as prescribed
by Theorem 8.1 and Corollary 8.3 of the punctured-R-map paper.

The substantive formula correction is in
[cjr_graph_factors.sage](../cjr_graph_factors.sage). Equation (8.21), with
all orbifold factors equal to one, gives

\[
\operatorname{Cont}_{0,\mathrm{mk}}
=\operatorname{ev}^*(\alpha)\,(t-3H)
\left(\frac{3H-t}{c}\right)^b
=\operatorname{ev}^*(\alpha)\,
\frac{(-1)^b}{c^b}w^{b+1}.
\]

The old expression was \(w(3H/c-t)^b\). The two agree for primary
insertions or contact one, which explains why many earlier tests passed.
For \(c=2\), the corrected first two descendants are

\[
-\frac{w^2}{2}
=-\frac{t^2}{2}+3tH-\frac92H^2,
\qquad
\frac{w^3}{4}
=\frac{t^3}{4}-\frac94t^2H+\frac{27}{4}tH^2
\pmod{H^3}.
\]

The new expected values were tested against the old implementation first
and failed. For the genus-three one-point `tau_4(pt)` probe, the marked
contact-five star contributes

\[
-\frac{H}{3}\,
\frac{5^5}{5!w^5}\,\frac{w^5}{5^4}
=-\frac{H}{72}.
\]

Consequently the coefficient of
\(V(3,0,(-5);\psi_{\min}=0,H)\) is **\(-1/72\)**, replacing
\(-625/72\). The compact side remains \(-31/967680\).
The primary `t^-4` relation still has coefficient \(-625/72\); it is a
different equation and is unchanged.

CJR I, equation (1.10) and Theorem 1.4, supplies the numerical comparison
with hypersurface descendants. This does not identify log-domain and
stabilized cotangent classes on the source moduli space. The code continues
to distinguish those conventions. Only the log-domain marked-tail factor
was corrected.

Checkpoint version **10** rejects older solver states. Explicit legacy
import accepts versions 8 and 9 but drops log-domain descendant rows,
non-scalar Chow rows, and axiom pseudo-probes. It imports neither solved
values nor completed-stage bookkeeping. Tests cover both rejection and
salvage, including preservation of an unaffected stabilized descendant.
The O(3) and Hodge caches remain mathematically valid. Historical rank
claims that used the faulty log-descendant rows have been withdrawn from
the documentation. The large corrected genus-three system has not been
rerun, so this audit supplies no new value or rank claim for its unresolved
effective invariant.

The performance work proceeded in three measured rounds. The initial
profile counted 6,080,529 zero-vertex classifications and 7,341,796 calls to
the marking-assignment generator on `(3,1,3)`. Most of those decorations
could never satisfy contact balance.

First, suppose \(L_i\) is the number of currently unmarked, genus-zero,
degree-zero leaves adjacent to infinity vertex \(i\). Its excess contact
budget is

\[
B_i=2g_i-2-3D_i.
\]

Each leaf left unmarked costs at least one excess unit because it needs
contact at least two. Thus it requires at least
\(\max(0,L_i-B_i)\) markings. Since these leaves each meet only one
infinity vertex, the entire decoration is impossible when

\[
\sum_i\max(0,L_i-B_i)>n.
\]

The enumerator now makes this test before generating marking assignments.
The first experimental `(3,1,3)` run fell from 37.62 to 8.79 seconds, with
the exact same graph orbits. These exploratory timings are retained
separately from the final repeated measurements.

Second, a tighter vertex bound removes whole impossible graph families.
Let \(Z,I,L,b\) be the numbers of zero vertices, infinity vertices, zero
leaves, and independent cycles. Let \(g_0,D_0,D_\infty\) be the zero-genus,
zero-degree, and infinity-degree sums. Summing contact balance bounds
nonspecial leaves by \(2(g-b-g_0)-2I-3D_\infty\). Every other leaf
consumes a marking, positive genus, or positive degree. Therefore

\[
L\le 2g-2b-2I+n+D-g_0-4D_\infty
\le 2g-2b-2I+n+D.
\]

Connectedness gives \(E=Z+I-1+b\), while summing zero-vertex valences gives
\(E\ge L+2(Z-L)\). Combining these yields

\[
\boxed{|V|+b\le 2g+n+D-1}.
\]

This replaces the looser bound \(|V|\le4g+2n+D-1\). One-level graphs
are treated separately. The argument preserves every allowed graph and its
symmetry factor; it changes the search space, not localization weights.
An independent unpruned enumeration under the old bound checks small
families, and the benchmark compares all 478 and 2,088 graph orbits in the
larger families.

Third, [o3_fixed_locus_graphs.sage](../o3_fixed_locus_graphs.sage) now uses
Mumford's identity instead of expanding an inverse Hodge polynomial as a
long geometric series. Define

\[
\Lambda_g^\pm(w)=\sum_{j=0}^g(\pm1)^j\lambda_jw^{g-j}.
\]

Since \(c(\mathbb E)c(\mathbb E^\vee)=1\),

\[
\Lambda_g^+(w)\Lambda_g^-(w)=w^{2g},
\qquad
\boxed{\frac1{\Lambda_g^+(w)}=
\frac{\Lambda_g^-(w)}{w^{2g}}}.
\]

This is an identity in the Chow ring, not a numerical approximation. The
new representative has at most \(g+1\) terms. Its justification also
follows from the vanishing even Chern characters recalled in Coates and
Givental, Section 3, Remark (2), p. 7.

Ordinary psi insertions are now accounted for before expanding the
remaining characteristic classes. A stable vertex with \(m\) flags and
markings has dimension \(3g-3+m\). If prescribed ordinary descendants
consume \(B=\sum_j b_j\), only degree at most

\[
3g-3+m-B
\]

can contribute from the remaining Hodge and flag factors. Negative residual
dimension returns zero. The final polynomial restores the prescribed psi
powers before integration. This avoids terms already forced to vanish by
dimension.

The twist convention was checked explicitly. CJR uses the dual Euler
class, with fiber weight \(t-3H\). The ordinary Euler twist with character
\(-t\) has fiber weight \(3H-t\). Riemann--Roch gives virtual rank
\(\chi(f^*\mathcal O(3))=3D+1-g\), hence

\[
e\!\left((R\pi_*f^*\mathcal O(3))^\vee\otimes\mathbb C_t\right)
=(-1)^{3D+1-g}
e\!\left(R\pi_*f^*\mathcal O(3)\otimes\mathbb C_{-t}\right).
\]

The existing conversion and ancestor/descendant separation were preserved.
The nonequivariant base limit is taken after the complete fixed-locus sum;
three distinct numerical base weights compute an equivariant theory.

Validation includes all **22 Sage suites** and the elliptic-cubic Python
suite, plus **8 standalone t-twisted tests**. Independent new O(3) checks use
\(\int_{\overline M_{g,1}}\psi^{3g-2}=1/(24^g g!)\) for genera 1--4,
where dimension forces all positive Hodge terms to vanish. Genus-one
mapping-to-a-point calculations separately give the nonequivariant dual
values \(-1/(8t)\), \(-1/8\), and zero for insertions \(1,H,H^2\).
The genus-zero three-point insertion `(H,1,1)` gives \(-3\), fixing the
Euler-twist sign.

The extended independent comparison also checks nine ordinary t-twisted
ancestor coefficients in genera 2--5 and degrees 1--2 against the fully
equivariant calibrated backend, applying the rank sign and then removing
base weights. Its exact results and log are stored with the benchmarks.
For example, genus 5, degree 2, insertion \(H\bar\psi^8\) agrees at
\(186937/3406233600\), and genus 4, degree 1, insertions
\((H\bar\psi^6,H)\) agree at \(-1/50400\).

To reproduce the benchmark, use the pinned Sage environment's Python from
each repository root:

```bash
python benchmarks/benchmark_hill_climb.py --output /tmp/before-1.json
# Repeat for before-2/3 and after-1/2/3 in the respective snapshots.
python benchmarks/compare_hill_climb.py /tmp benchmarks/reproduced-results
bash run_tests.sh
```

For the separate t-twisted suite, run
`python -m unittest -v test_t_twisted_cubic` from the parent `test` directory.
The benchmark code and mathematical tests use exact arithmetic throughout.
The audit covers the implemented compact plane-cubic sector and the stated
finite checks; it does not certify unsupported orbifold sectors or prove
that localization determines every effective invariant.

References:

1. Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models I: The Tropical Decomposition Formula*, October 17, 2025 iCloud draft, Theorem 1.4 and (1.10), pp. 7--8. [Local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_1_the_Tropical_Decomposition_Formula.pdf>).
2. Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, October 17, 2025 iCloud draft, (8.21), p. 82; (8.22)--(8.23), p. 83; (9.7) and Section 9.3, pp. 89--90. [Local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>).
3. Q. Chen, F. Janda, Y. Ruan, *Punctured logarithmic R-maps*, April 25, 2023 local version, Theorem 8.1, (8.3)--(8.4), and Corollary 8.3. [Local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Punctured_Logarithmic_R-Maps.pdf>); [arXiv](https://arxiv.org/abs/2208.04519).
4. T. Coates and A. Givental, *Quantum Riemann--Roch, Lefschetz and Serre*, Section 3, Remark (2), p. 7; Section 4 for equivariant Euler conventions; Appendix 2 for ancestors and descendants. [Author's PDF](https://math.berkeley.edu/~giventh/papers/rls.pdf), [arXiv](https://arxiv.org/abs/math/0110142).
