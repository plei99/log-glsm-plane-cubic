The new [complete-intersection module](../projective_complete_intersections.sage) treats all nine nonlinear Fano hypersurface types of dimensions two through four and all five smooth CY3 complete-intersection types in ordinary projective space. The calculation covers connected unmarked genera two and three, ambient curve degrees \(D=0,1,2,3\), and every codimension on \(\overline{\mathcal M}_g\). Two nonresonant auxiliary torus-weight choices are used for each target. The complete numerical tables, independent relation bases, and exact Pixton certificates are linked below.

Here \(d_j\) always denotes a defining polynomial degree, while \(D\) denotes a curve's degree against the hyperplane class. The program is parameterized by arbitrary positive defining degrees and ambient dimension; this particular scan is restricted to the fourteen requested types. Linear defining equations can be eliminated. A degree-one Fano hypersurface is projective space itself, which the parameterized module also accepts.

The scan retained 2463 relations across individual targets. Removing repetitions across all fourteen targets leaves 463 independent relations: 23 in genus two and 440 in genus three. Every retained pure tautological relation has a rational Pixton membership certificate. This is a statement about the computed finite search. It does not assert that every possible log-GLSM relation lies in Pixton's system, or that Pixton's system is complete in general. The CJR identities involving the aggregate of infinity contributions have a separate geometric meaning, explained below.

Write
\[
Y=Y_{d_1,\ldots,d_c}\subset\mathbb P^N,\qquad
E=\bigoplus_{j=1}^c\mathcal O(d_j),\qquad
m=N-c,\quad s=\sum_jd_j,\quad I=N+1-s.
\]
The Fano condition for a hypersurface of dimension \(m\) is \(d\le m+1\). Thus the surface cases are degrees 2 and 3; the threefold cases are degrees 2, 3, and 4; the fourfold cases are degrees 2, 3, 4, and 5. A CY3 complete intersection satisfies \(N=c+3\) and \(\sum(d_j-1)=4\). The partitions of 4 give exactly
\[
(5)\subset\mathbb P^4,\quad (2,4),(3,3)\subset\mathbb P^5,\quad
(2,2,3)\subset\mathbb P^6,\quad (2,2,2,2)\subset\mathbb P^7.
\]
This list concerns a single ordinary projective space; products and weighted projective spaces have additional types.

The CJR conventions were checked against the supplied October 17, 2025 draft of *Structural Formulae in Logarithmic Gauged Linear Sigma Models III*: Theorem 8.4 and Remark 8.5 give localization and regularity, §9.3 treats hypersurfaces, and §9.4 explicitly treats higher-rank bundles. The [source audit](PROJECTIVE_COMPLETE_INTERSECTION_SOURCES.md) records the page numbers and links to the private PDF. The public arbitrary-rank compact comparison is also in [Chen–Janda–Ruan, Proposition 4.1 and Corollary 1.10](https://arxiv.org/abs/1906.04345v2); the fields comparison sign is in [Chen–Janda–Webb, Corollary 1.2.1](https://arxiv.org/abs/1911.09825v2). Use of the private draft means use of its stated formula, not an independent proof of that manuscript.

For a compact target, the ordinary virtual dimension and resulting codimension after stabilization are
\[
\nu=(1-g)(m-3)+ID,
\qquad
\alpha=3g-3-\nu=m(g-1)-ID.
\]
Define the ambient dual-Euler pushforward
\[
Z_D(t)=\operatorname{st}_*\left(
 e_t\big((R\pi_*f^*E)^\vee\big)
 [\overline{\mathcal M}_g(\mathbb P^N,D)]^{\rm vir}\right)
=\sum_{r=0}^{3g-3} t^{\alpha-r}Z_{D,r}.
\]
The Euler class is a ratio for the derived complex, including its \(R^1\) term. Every summand of \(E\) has the same parameter \(t\). In particular, the higher-genus calculation does not substitute the genus-zero convex quantum Lefschetz rule.

Let \(B_D(t)\) be the sum of all CJR localization contributions containing infinity data, including mixed graphs and all splittings. The full reduced pushforward is \(Z_D+B_D\). Regularity and the compact comparison give
\[
\boxed{B_{D,r}=-Z_{D,r}\quad(r>\alpha)},
\qquad
\boxed{B_{D,\alpha}=(-1)^{c(1-g)+sD}G_D-Z_{D,\alpha}},
\]
where the second equality applies when \(0\le\alpha\le3g-3\), and
\[
G_D=\sum_{\beta:\,i_*\beta=D[\ell]}
\operatorname{st}_*[\overline{\mathcal M}_g(Y,\beta)]^{\rm vir}
\]
is the ordinary compact-target class. Positive powers of \(t\) impose neither of these vanishing identities. The code keeps an unknown positive-degree compact class as unknown. For surfaces, the sum over \(\beta\) matters: a quadric has two rulings, and a cubic surface has many curve classes. The calculations sum all classes of a fixed ambient degree.

For higher rank, the infinity target is \(\mathbb P_{\mathbb P^N}(E^\vee)\), as in CJR III §9.4. The implementation multiplies the zero-vertex twisting factors. It does not replace the infinity theory by a product of hypersurface infinity theories. The displayed formulas determine the stabilized aggregate \(B_{D,r}\) in the stated ranges; they do not reconstruct each punctured effective vertex individually.

The fixed-graph calculation uses [Graber–Pandharipande, §4](https://arxiv.org/abs/alg-geom/9708001v2). At a contracted vertex of genus \(h\), colour \(i\), and valence \(v\), put
\[
w_{ik}=a_i-a_k,\qquad q_{ij}=t-d_ja_i,
\qquad \Lambda_h^\vee(x)=\sum_{\ell=0}^h(-1)^\ell\lambda_\ell x^{h-\ell}.
\]
Before flag denominators, the vertex factor is
\[
\left(\prod_{k\ne i}w_{ik}\right)^{v-1}
\prod_{k\ne i}\Lambda_h^\vee(w_{ik})
\prod_{j=1}^c\left[q_{ij}^{1-v-2h}\Lambda_h^\vee(q_{ij})\right].
\]
Each flag contributes \((w_{ik}/e-\psi)^{-1}\), where \(e\) is its edge degree. The expression uses Mumford's identity
\(1/\Lambda_h^+(q)=\Lambda_h^\vee(q)/q^{2h}\) in rational Chow. A degree-\(e\) edge from \(i\) to \(k\) contributes
\[
\frac{(-1)^e e^{2e}}{(e!)^2(a_i-a_k)^{2e}}
\frac{\displaystyle\prod_{j=1}^c\prod_{u=0}^{d_je}
 \left(t-\frac{d_je-u}{e}a_i-\frac{u}{e}a_k\right)}
{\displaystyle\prod_{\ell\ne i,k}\prod_{u=0}^{e}
 \left(\frac{u}{e}a_i+\frac{e-u}{e}a_k-a_\ell\right)}.
\]
Deck factors \(1/e\) and graph automorphisms are supplied separately by the graph enumerator. Univalent and bivalent unstable vertices receive their usual unstable factors; the bivalent twisting denominator is \(\prod_jq_{ij}\). The source audit derives these conventions and the independent checks.

To obtain classes on the correct moduli space, the program retains stable contracted components, turns rational tails into temporary markings, glues the remaining flags, and forgets those temporary markings. Components that become unstable are then contracted. This includes the correction terms from forgetting several tails and the degree-three case in which a rational contracted component destabilizes. No integral is substituted for a class before stabilization.

Set \(a_i=\epsilon b_i\), and use homogeneity to set \(t=1\) during the exact rational calculation. Proper equivariant pushforward is regular at \(\epsilon=0\). Hence every negative coefficient
\[
[\epsilon^{-k}]Z_{D,r}(1;\epsilon b)=0\quad(k>0)
\]
is a tautological relation. Two different choices of \(b\) have the same finite limit in Chow, so their finite-part difference is another relation. The weight rays use the binary digits of \(i\), interpreted in bases 4 and 5; for \(\mathbb P^7\) these are
\[
(0,1,4,5,16,17,20,21),\qquad
(0,1,5,6,25,26,30,31).
\]
They avoid edge and bivalent-vertex resonances through curve degree three. The negative \(\epsilon\) powers here are distinct from the negative CJR \(t\) powers above.

For each genus and codimension, let \(M\) be the matrix of generated Pixton 3-spin relations in the saved strata basis. Each retained localization relation \(v\) is certified by exact rational coefficients \(u\) satisfying
\[
uM=v.
\]
The export script reloads and replays every certificate, checks the claimed ranks, and removes repetitions across targets. Membership uses the Chow validity proved by [Janda, Theorem 1.1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf), for the system constructed by [Pandharipande–Pixton–Zvonkine](https://arxiv.org/abs/1303.1043). No completeness conjecture is needed. “Independent relations” means independent vectors in the strata algebra before taking the Chow quotient; their images in Chow are all zero.

The exact tables follow. A count of negative-\(t\) coefficient slots includes identities whose two sides are zero; it must not be added to the ranks of pure tautological relations.

All results use unmarked genera two and three, ambient degrees D=0,1,2,3, and two nonresonant weight rays. Each rank counts independent strata expressions before quotienting by tautological relations.

| Target | Genus 2 ranks r=1,2,3 | Genus 3 ranks r=2,3,4,5,6 | Negative-t coefficient slots (g2,g3) |
|---|---|---|---|
| quadric surface | 1, 6, 16 | 6, 15, 21, 21, 27 | 12, 19 |
| cubic surface | 1, 5, 16 | 6, 15, 21, 21, 27 | 10, 14 |
| quadric threefold | 1, 6, 16 | 6, 26, 39, 45, 51 | 11, 16 |
| cubic threefold | 1, 6, 16 | 6, 26, 39, 45, 51 | 10, 12 |
| quartic threefold | 1, 6, 16 | 6, 30, 39, 45, 51 | 6, 6 |
| quadric fourfold | 1, 6, 16 | 6, 25, 39, 45, 51 | 11, 15 |
| cubic fourfold | 1, 6, 16 | 6, 28, 39, 45, 51 | 10, 12 |
| quartic fourfold | 1, 6, 16 | 6, 30, 39, 45, 51 | 8, 6 |
| quintic fourfold | 1, 6, 16 | 6, 30, 39, 45, 51 | 3, 1 |
| CY3 (5,) | 1, 6, 16 | 6, 30, 39, 45, 51 | 0, 0 |
| CY3 (2, 4) | 1, 6, 16 | 6, 30, 39, 45, 51 | 0, 0 |
| CY3 (3, 3) | 1, 6, 16 | 6, 30, 39, 45, 51 | 0, 0 |
| CY3 (2, 2, 3) | 1, 6, 16 | 6, 30, 40, 45, 51 | 0, 0 |
| CY3 (2, 2, 2, 2) | 1, 6, 16 | 6, 18, 22, 22, 27 | 0, 0 |

Across all targets, the independent ranks are:

| Genus | Codimension | Found | Generated Pixton rank |
|---|---|---:|---:|
| 2 | 1 | 1 | 1 |
| 2 | 2 | 6 | 6 |
| 2 | 3 | 16 | 16 |
| 3 | 2 | 6 | 6 |
| 3 | 3 | 30 | 30 |
| 3 | 4 | 86 | 114 |
| 3 | 5 | 150 | 259 |
| 3 | 6 | 168 | 377 |

Exact ambient twisted top integrals, genus 2:

| Target | D=0 | D=1 | D=2 | D=3 |
|---|---|---|---|---|
| quadric surface | -1/288 | 17/120 | 0 | 0 |
| cubic surface | -1/288 | 1/72 | 0 | 0 |
| quadric threefold | 0 | 53/192 | 0 | 0 |
| cubic threefold | 0 | 25/1152 | 0 | 0 |
| quartic threefold | 0 | 65/1152 | -25/16 | 0 |
| quadric fourfold | 0 | 343/720 | 0 | 0 |
| cubic fourfold | 0 | 1/32 | 0 | 0 |
| quartic fourfold | 0 | 29/96 | 0 | 0 |
| quintic fourfold | 0 | 3509/576 | 9361/4 | 65250 |
| CY3 (5,) | 0 | 1645/1152 | -351685/144 | 1652182775/384 |
| CY3 (2, 4) | 0 | 467/480 | 78323/144 | 1689847/5 |
| CY3 (3, 3) | 0 | 49/64 | 44011/160 | 4040769/32 |
| CY3 (2, 2, 3) | 0 | 469/480 | -42593/192 | 5479821/80 |
| CY3 (2, 2, 2, 2) | 0 | 181/180 | 16933/90 | 559072/15 |

Exact ambient twisted top integrals, genus 3:

| Target | D=0 | D=1 | D=2 | D=3 |
|---|---|---|---|---|
| quadric surface | 1/72576 | -179/6048 | -491/3780 | 0 |
| cubic surface | 1/72576 | -21293/362880 | 0 | 0 |
| quadric threefold | 0 | -22249/290304 | -1541/3024 | 0 |
| cubic threefold | 0 | -1087/10752 | 0 | 0 |
| quartic threefold | 0 | -21017/290304 | -1525/13824 | 125/96 |
| quadric fourfold | 0 | -949/6912 | -10009/6480 | 0 |
| cubic fourfold | 0 | -737/5760 | 0 | 0 |
| quartic fourfold | 0 | -121/1728 | -89/192 | 0 |
| quintic fourfold | 0 | -397/6912 | -136723687/82944 | -49574937/64 |
| CY3 (5,) | 0 | -3265/290304 | -4368535/72576 | 53253995/256 |
| CY3 (2, 4) | 0 | -12091/69120 | 948043/362880 | -1388524609/90720 |
| CY3 (3, 3) | 0 | -1637/5760 | 323153/53760 | -22057333/4480 |
| CY3 (2, 2, 3) | 0 | -609/1280 | -1152557/138240 | 23035193/6048 |
| CY3 (2, 2, 2, 2) | 0 | -221/864 | -53243/22680 | -245269/135 |

Ordinary CY3 GW numbers (distinct from the ambient twisted values above):

| CY3 | Lines | N₂,₀ | N₃,₀ | N₂,₁ | N₃,₁ |
|---|---:|---|---|---|---|
| CY3 (5,) | 2875 | -5/144 | 5/36288 | 575/48 | 2875/6048 |
| CY3 (2, 4) | 1280 | -11/360 | 11/90720 | 16/3 | 40/189 |
| CY3 (3, 3) | 1053 | -1/40 | 1/10080 | 351/80 | 39/224 |
| CY3 (2, 2, 3) | 720 | -1/40 | 1/10080 | 3 | 5/42 |
| CY3 (2, 2, 2, 2) | 512 | -1/45 | 1/11340 | 32/15 | 16/189 |

Integrated CY3 infinity aggregates at t⁰, obtained from the signed compact comparison:

| CY3 | g2,D0 | g3,D0 | g2,D1 | g3,D1 |
|---|---|---|---|---|
| CY3 (5,) | 5/144 | 5/36288 | 12155/1152 | -134735/290304 |
| CY3 (2, 4) | -11/360 | 11/90720 | 2093/480 | 187037/483840 |
| CY3 (3, 3) | -1/40 | 1/10080 | 1159/320 | 18479/40320 |
| CY3 (2, 2, 3) | 1/40 | 1/10080 | 971/480 | 9589/26880 |
| CY3 (2, 2, 2, 2) | -1/45 | 1/11340 | 203/180 | 2059/6048 |


Two short expressions in the computed relation span are
\[
\kappa_1=2\lambda_1+\frac12\xi_{1,1,*}(1\boxtimes1)
\quad\text{on }\overline{\mathcal M}_2,
\]
and
\[
\boxed{\kappa_2-\lambda_1\kappa_1+4\lambda_2
=\xi_{1,2,*}\bigl(1\boxtimes(\psi_1-\lambda_1)\bigr)}
\quad\text{on }\overline{\mathcal M}_3.
\]
Here \(\xi_{1,2}\) attaches the marked point of a genus-one curve to that of a genus-two curve. The right side is a class supported on the separating boundary, with the cotangent/Hodge decoration on the genus-two factor. The genus-two formula has a factor \(1/2\) because the ordered genus-one gluing map has degree two onto its boundary divisor. The export includes exact coordinates expressing these short relations in the computed union basis, in `short-geometric-relations.json`. They are familiar tautological identities recovered by the calculation; no novelty claim is intended.

The surface examples have \(\alpha=2-ID\) in genus two and \(\alpha=4-ID\) in genus three. The quadric surface has \(I=2\), and the cubic surface has \(I=1\). Their ordinary constant-map classes, derived directly from the obstruction bundle, are respectively
\[
G_{2,0}=8\lambda_2,\quad G_{3,0}=8\lambda_1\lambda_3,
\qquad
G_{2,0}=3\lambda_2,\quad G_{3,0}=3\lambda_1\lambda_3.
\]
Their ambient constant-map classes agree because both use \(\mathbb P^3\):
\[
Z_0(t)=-20\lambda_1\lambda_2\,t^{-1}\quad(g=2),
\qquad
Z_0(t)=20\lambda_1\lambda_2\lambda_3\,t^{-2}\quad(g=3).
\]
Thus, in genus three, both have the explicit infinity-aggregate identities
\[
B_{0,5}=0,\qquad B_{0,6}=-20\lambda_1\lambda_2\lambda_3,
\]
and their \(t^0\) comparisons give \(B_{0,4}=8\lambda_1\lambda_3\) or \(3\lambda_1\lambda_3\). In genus two the compact comparison has sign minus, so \(B_{0,2}=-8\lambda_2\) or \(-3\lambda_2\), while \(B_{0,3}=20\lambda_1\lambda_2\). These are class identities, not just integrals.

The three Fano threefolds have indices 3, 2, and 1. Their genus-three codimensions are \(6-3D\), \(6-2D\), and \(6-D\). For a threefold, put \(C_Y=\int_Y(c_3-c_1c_2)\). The quadric, cubic, and quartic give \(C_Y=-20,-30,-80\), respectively, and
\[
G_{2,0}=C_Y\lambda_1\lambda_2,
\qquad G_{3,0}=-C_Y\lambda_1\lambda_2\lambda_3.
\]
Their ordinary genus-three degree-zero numbers are \(1/72576\), \(1/48384\), and \(1/18144\). The ambient degree-zero class is zero because \(\mathbb P^4\) has negative degree-zero virtual dimension for \(g\ge2\). Consequently the infinity aggregate supplies the entire signed compact constant-map class. A zero in the ambient twisted table does not imply a zero ordinary GW invariant.

The four Fano fourfolds have indices 4, 3, 2, and 1, giving genus-three codimensions \(8-4D\), \(8-3D\), \(8-2D\), and \(8-D\). The ordinary constant-map virtual dimension is negative, so \(G_{g,0}=0\). At degree one, the quadric has negative-\(t\) identities in codimensions 5 and 6; the cubic has one in codimension 6; the quartic has its compact comparison in codimension 6; the quintic has \(\alpha=7\), outside the Chow range. At degree two the quintic reaches its compact comparison in codimension 6, and at degree three it gives \(B_{3,6}=-Z_{3,6}\). These distinctions are why the grading is recomputed for each target and degree.

All five CY3s have \(I=0\) and \(\alpha=3g-3\), independently of \(D\). There are no negative-\(t\) coefficient slots on unmarked \(\overline{\mathcal M}_g\). Their pure tautological relations still come from the auxiliary \(\epsilon\)-pole and finite-limit equations. The top \(t^0\) comparison retains ordinary GW theory:
\[
\int B_{D,3g-3}=(-1)^{c(1-g)+sD}N^Y_{g,D}
-\int Z_{D,3g-3}.
\]
The five Euler characteristics are \(-200,-176,-144,-144,-128\). They follow independently from
\[
c(TY)=\frac{(1+H)^{N+1}}{\prod_j(1+d_jH)},\qquad
\int_YH^3=\prod_jd_j.
\]
By the constant-map obstruction bundle and [Faber–Pandharipande, Theorem 4](https://arxiv.org/abs/math/9810173),
\[
N^Y_{2,0}=\frac{\chi(Y)}{5760},\qquad
N^Y_{3,0}=-\frac{\chi(Y)}{1451520}.
\]
The quintic has 2875 lines; the \((2,4)\) intersection has 1280; the bicubic has 1053; the \((2,2,3)\) intersection has 720; four quadrics have 512. The code obtains each count both by fixed-line localization and by the independent Grassmannian formula
\[
L_Y=[x^Ny^{N-1}]\left((x-y)
\prod_{j=1}^c\prod_{u=0}^{d_j}(ux+(d_j-u)y)\right).
\]
The four non-quintic counts appear in [Libgober–Teitelbaum, Table 1](https://homepages.math.uic.edu/~libgober/otherpapers/export/1993linesoncy.pdf); all five also appear in the [Macaulay2 EnumerationCurves documentation](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EnumerationCurves/html/_rational__Curve.html).

On a general member these lines are isolated and transverse. Their normal bundles are \(\mathcal O(-1)^{\oplus2}\), since the normal bundle has degree \(-2\) and no sections. [Faber–Pandharipande, Theorem 3](https://arxiv.org/abs/math/9810173) gives the contribution of degree-\(k\) covers of one such line:
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3}.
\]
Every degree-one map has a line as its image, so
\[
N^Y_{2,1}=L_Y/240,\qquad N^Y_{3,1}=L_Y/6048.
\]
These include contracted higher-genus components. The formulas extend by deformation invariance to smooth members. Together with the computed twisted integrals they give the exact degree-one integrated infinity aggregates in the table. The raw class scan conservatively records positive-degree compact classes as unspecified; this numerical refinement is stored separately in `ordinary-cy3-gw.json`. Ordinary degree-two and degree-three invariants are not inferred from the ambient twisted table.

The implementation uses Sage 10.9 in the native Linux OrbStack machine `agents`. The shared stabilization cache stores curve geometry independently of the target, and worker processes evaluate different target Euler factors. The new module recovers the quartic K3 backend through degree one in every codimension, and reproduces its genus-two degree-two and degree-three integrals, \(-67/36\) and \(133/2\). Regression tests also check ordinary/reduced distinctions, Chern constants, all five line counts by two methods and two torus choices, compact-comparison signs, and rejection of unsupported ordinary GW numbers. The original quartic test suite passes after the small shared bivalent-vertex hook was introduced.

The batch computation integrates top classes using the exact functional annihilating the Pixton rows and normalized by
\[
\int_{\overline{\mathcal M}_g}\kappa_{3g-3}
=\langle\tau_{3g-2}\rangle_g=\frac{1}{24^g g!}.
\]
The export verification independently checks degree-one integrals in both genera for every ambient dimension \(N=3,4,5,6,7\), using direct vertex Hodge integration rather than that functional.

Reproduce the calculation from the project directory with Sage's Python:

```bash
python -c 'from sage.all import *; load("tests/test_projective_complete_intersections.sage"); load("tests/test_quartic_k3.sage")'
python calculations/projective_ci_scan.py --workers 3 --max-degree 3 --bases 4 5
python calculations/projective_ci_export.py
```

From macOS, run those commands inside `orb -m agents -w /home/patrick/agent-work/repos/research/test/log-glsm-plane-cubic /bin/bash -lc`, with `/home/patrick/agent-work/tools/mamba/envs/sage/bin` on `PATH`. No floating-point arithmetic enters a reported coefficient or certificate.

For an individual example, the module exposes class-valued and numerical calls:

```python
load('projective_complete_intersections.sage')
Y = ProjectiveCompleteIntersection(5, (5,))  # quintic fourfold
C = CompleteIntersectionCycles(Y, genus=3)
Z = C.cycle(3, 6)                           # Z_{D=3,r=6}
comparison = C.aggregate_relation(3, 6)     # known infinity aggregate = -Z
quintic = ProjectiveCompleteIntersection(4, (5,))
cy3_ordinary_gw_number(quintic, 3, 1)        # 2875/6048
```

A comparison with an unknown compact class returns `compact_known=False` and `infinity_aggregate=None`; it never silently replaces that class by zero.

The [result directory](../results/complete-intersections/2026-11-06) contains one folder per target. `finite-classes-g*.json.gz` stores the formal localization terms, stabilized strata vectors, grading, and known degree-zero compact terms. `relations-g*-r*.json.gz` stores independent relation vectors with their ambient-degree/weight origins and sparse Pixton certificates. The common `strata-g*-r*.json` dictionaries and `pixton-g*-r*.sobj` matrices fix the certificate conventions. `union-relations-g*-r*.json.gz` removes repetitions across targets; [`independent-relations.txt`](../results/complete-intersections/2026-11-06/independent-relations.txt) gives the complete readable union basis and its strata dictionary (also supplied as gzip). `verification.json` records the replayed certificate counts and direct Hodge checks.

This finite search does not exhaust higher degrees, additional weight rays, markings, insertions, or genera. It computes the specified stabilized zero-vertex classes and aggregate comparisons. Individual infinity effective cycles and ordinary positive-degree Fano GW classes remain separate reconstruction problems.
