# Projective complete intersections: CJR comparison, ambient localization, and independent constants

Research date: 2026-11-06. Work over \(\mathbb C\), with rational Chow groups and connected, unmarked genus \(g\ge2\), principally \(g=2,3\). This note checks geometry and derives constants independently; it does not audit implementation code. The optional genus-zero line counts are separate validation targets.

The requested extension is mathematically supported. In particular, the supplied private CJR III draft **explicitly treats arbitrary-rank complete intersections in §9.4, p. 92**. Its compact comparison has sign \((-1)^{c(1-g)+sD}\), and its zero-vertex twist is the dual Euler ratio for the entire rank-\(c\) bundle. Thus the higher-rank aggregate interpretation is directly supported by this draft, not merely inferred from its hypersurface §9.3.

Ordinary GW theory does **not** vanish for these targets in general. Negative-\(t\) cancellation and the \(t^0\) ordinary-GW comparison are different statements. For the five CY3 targets, \(\alpha=3g-3\): there are **no negative-\(t\) identities after the existing unmarked pushforward to \(\overline{\mathcal M}_g\)**. There is a top-codimension \(t^0\) comparison, whose right side must retain ordinary GW theory.

## 1. Primary sources and exact locators

**[CJR-III]** Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, supplied private draft dated **October 17, 2025**. Read through the supplied text extraction /tmp/log-glsm-audit/cjr-iii.txt; [supplied PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>). Page numbers below are printed draft pages.

- §4.2.2–§4.2.3: canonical and reduced obstruction theories; Assumption 4.2 is the compact-type marking condition.
- Theorem 8.4, equation (8.28), and Remark 8.5, p. 84: general localization and regularity in \(t\). The theorem is not rank-one-only.
- §8.7, p. 85 onward: the **additional simplification** when the bundle has rank one.
- §9.3, pp. 89–90, equation (9.7): hypersurface compact comparison and zero-vertex twist.
- **§9.4, p. 92:** rank-\(R\) bundle \(E=E_1\); the displayed compact comparison immediately after its setup; stable zero-vertex factor; infinity target \(\mathbb P_X(E^\vee)\).
- §10.1, Theorem 10.1, p. 93: ambient GW invariants use the sum over target curve classes with a specified ambient pushforward.

This note uses the localization theorem as stated in that private version. It does not certify the draft's proofs or turn its unpublished status into a public theorem citation.

**[CJR-GLSM]** Q. Chen, F. Janda, Y. Ruan, *The logarithmic gauged linear sigma model*, Invent. Math. **225** (2021), 1077–1154; [arXiv:1906.04345v2](https://arxiv.org/abs/1906.04345v2), [PDF](https://arxiv.org/pdf/1906.04345v2). Corollary 1.10, preprint p. 7, and Proposition 4.1 with its proof, pp. 34–35, give the arbitrary-rank reduced-log-GLSM/ordinary-target comparison. In our scheme setting there are no orbifold-age corrections.

**[CJW]** Q. Chen, F. Janda, R. Webb, *Virtual cycles of stable (quasi-)maps with fields*, Adv. Math. **385** (2021), 107781; [arXiv:1911.09825v2](https://arxiv.org/abs/1911.09825v2), [PDF](https://arxiv.org/pdf/1911.09825v2). Theorem 1.1.1 and Corollary 1.2.1, pp. 3–4: the cosection-localized fields class equals the target virtual class with sign \((-1)^{\chi(f^*E)}\). Theorem 1.2.3, pp. 4–5: convex quantum Lefschetz. The introduction explicitly distinguishes the genus-zero Euler-class principle from its failure in positive genus.

**[GP]** T. Graber, R. Pandharipande, *Localization of virtual classes*, Invent. Math. **135** (1999), 487–518; [arXiv:alg-geom/9708001v2](https://arxiv.org/abs/alg-geom/9708001v2), [PDF](https://arxiv.org/pdf/alg-geom/9708001v2). Equation (1) and §4, preprint pp. 13–18: projective-space fixed graphs, inverse normal Euler factors, and graph/deck automorphisms. Section 6 below applies this framework and independently normalizes the extra bundle twist; its combined formulas are not quoted verbatim from GP.

**[FP]** C. Faber, R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. **139** (2000), 173–199; [arXiv:math/9810173](https://arxiv.org/abs/math/9810173), [PDF](https://arxiv.org/pdf/math/9810173). Introduction equation (1) and the next paragraph, p. 1: the constant-map obstruction bundle. Theorem 4, p. 4, and §4.3, pp. 20–21: the Bernoulli integral for \(\lambda_{g-1}^3\) and its relation to \(2\lambda_g\lambda_{g-1}\lambda_{g-2}\).

**[M]** D. Mumford, *Towards an Enumerative Geometry of the Moduli Space of Curves*, in *Arithmetic and Geometry II*, Progress in Math. **36** (1983), 271–328; [author-hosted scan](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf). §5, equations (5.2), (5.4), Corollary (5.3), printed pp. 305–307: algebraic GRR gives \(\operatorname{ch}_{2k}(\mathbb E)=0\), hence \(c(\mathbb E)c(\mathbb E^\vee)=1\) in rational Chow, including the compactification.

**[LT]** A. Libgober, J. Teitelbaum, *Lines on Calabi–Yau complete intersections, mirror symmetry, and Picard–Fuchs equations*, IMRN **1993**, 29–39; [arXiv:alg-geom/9301001](https://arxiv.org/abs/alg-geom/9301001), [author PDF](https://homepages.math.uic.edu/~libgober/otherpapers/export/1993linesoncy.pdf). Preprint p. 2 computes the bicubic Euler characteristic; Table 1, p. 9, gives the four non-quintic degree-one counts.

**[M2]** Official Macaulay2 EnumerationCurves documentation, [rationalCurve](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EnumerationCurves/html/_rational__Curve.html), “For lines,” examples i1–i3: all five counts in the order used here.

**[Stacks]** [Chow Homology and Chern Classes, Lemma 42.44.1, tag 0FA8](https://stacks.math.columbia.edu/tag/0FA8): a regular section realizes the top Chern class. This supplies the ordinary intersection-theoretic interpretation of the independent Grassmannian calculation in §7.

## 2. Scope, completeness, and Chern data

Write
\[
Y=Y_{d_1,\ldots,d_c}\subset X=\mathbb P^N,\qquad
E=\bigoplus_{j=1}^c\mathcal O_X(d_j),\qquad
s=\sum_jd_j,\quad m=N-c,\quad I=N+1-s.
\]
Let \(H=c_1(\mathcal O_Y(1))\), \(h=\int_YH^m=\prod_jd_j\), and \(D=\int_\beta H\). An ambient-degree invariant sums over **all** effective target classes with that ambient pushforward. This matters for surfaces: a quadric has two ruling classes, and a cubic surface has many curve classes. One must not replace target curve-class information by a presumed rank-one \(H_2(Y)\). This is the convention in [CJR-III, §10.1].

### Completeness within the requested families

The Euler and normal exact sequences give, independently,
\[
c(TY)=\frac{(1+H)^{N+1}}{\prod_j(1+d_jH)},\qquad
K_Y=\mathcal O_Y(s-N-1).
\]
For hypersurfaces \(N=m+1\), the Fano condition is \(d<m+2\). Removing degree-one equations leaves
\[
(m,d)=(2,2),(2,3),(3,2),(3,3),(3,4),
 (4,2),(4,3),(4,4),(4,5).
\]
There are \(2+3+4=9\) types. This is completeness among nonlinear Fano **hypersurfaces** in the stated dimensions, not a classification of all Fano complete intersections or all Fano varieties.

For a CY3 complete intersection in one ordinary projective space, \(N=c+3\) and \(s=N+1=c+4\). After eliminating linear equations, \(d_j\ge2\), and
\[
\sum_j(d_j-1)=4.
\]
The five partitions of \(4\) give exactly
\[
(5)\mathbb P^4,\quad (2,4)\mathbb P^5,\quad
(3,3)\mathbb P^5,\quad (2,2,3)\mathbb P^6,\quad
(2,2,2,2)\mathbb P^7.
\]
Eliminating a linear equation lowers both \(N\) and \(c\) by one, so it preserves the variety and the condition above. Smooth examples exist by choosing general transverse sections. The familiar additional CY conditions follow from the complete-intersection Koszul resolution. This argument is restricted to ordinary projective space; weighted projective spaces and products have further types.

### Independent Chern calculation

The following tables were obtained by expanding the displayed total Chern class, with exact integer arithmetic. They do not depend on localization output. For a reproducible alternative, put \(p_k=N+1-\sum_jd_j^k\), so \(p_1=I\), and write \(c_k(TY)=b_kH^k\). Newton identities give
\[
\begin{aligned}
b_1&=I,\\
b_2&=(I^2-p_2)/2,\\
b_3&=(I^3-3Ip_2+2p_3)/6,\\
b_4&=(I^4-6I^2p_2+3p_2^2+8Ip_3-6p_4)/24.
\end{aligned}
\]
Every Chern number is \(h\) times the corresponding product of \(b_k\)'s; the topological Euler characteristic is \(hb_m\).

| Fano target | \(N\) | \(I\) | \((b_1,\ldots,b_m)\) | \(\chi_{\rm top}(Y)\) |
|---|---:|---:|---|---:|
| Quadric surface | 3 | 2 | \((2,2)\) | 4 |
| Cubic surface | 3 | 1 | \((1,3)\) | 9 |
| Quadric threefold | 4 | 3 | \((3,4,2)\) | 4 |
| Cubic threefold | 4 | 2 | \((2,4,-2)\) | -6 |
| Quartic threefold | 4 | 1 | \((1,6,-14)\) | -56 |
| Quadric fourfold | 5 | 4 | \((4,7,6,3)\) | 6 |
| Cubic fourfold | 5 | 3 | \((3,6,2,9)\) | 27 |
| Quartic fourfold | 5 | 2 | \((2,7,-8,47)\) | 188 |
| Quintic fourfold | 5 | 1 | \((1,10,-30,165)\) | 825 |

Here are all Chern numbers in each of these dimensions; integration over \(Y\) is understood in every numeric column.

| Surface | \(c_1^2\) | \(c_2\) |
|---|---:|---:|
| Quadric | 8 | 4 |
| Cubic | 3 | 9 |

| Threefold | \(c_1^3\) | \(c_1c_2\) | \(c_3\) | \(c_3-c_1c_2\) |
|---|---:|---:|---:|---:|
| Quadric | 54 | 24 | 4 | -20 |
| Cubic | 24 | 24 | -6 | -30 |
| Quartic | 4 | 24 | -56 | -80 |

| Fourfold | \(c_1^4\) | \(c_1^2c_2\) | \(c_2^2\) | \(c_1c_3\) | \(c_4\) |
|---|---:|---:|---:|---:|---:|
| Quadric | 512 | 224 | 98 | 48 | 6 |
| Cubic | 243 | 162 | 108 | 18 | 27 |
| Quartic | 64 | 112 | 196 | -64 | 188 |
| Quintic | 5 | 50 | 500 | -150 | 825 |

For the CY3 targets, \(c_1=0\), so \(c_1^3=c_1c_2=0\). The useful remaining data are:

| CY3 type | \(N\) | \(c\) | \(s\) | \(h=\int H^3\) | \(b_2\) | \(\int c_2H\) | \(b_3\) | \(\chi_{\rm top}=\int c_3\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| \((5)\) | 4 | 1 | 5 | 5 | 10 | 50 | -40 | -200 |
| \((2,4)\) | 5 | 2 | 6 | 8 | 7 | 56 | -22 | -176 |
| \((3,3)\) | 5 | 2 | 6 | 9 | 6 | 54 | -16 | -144 |
| \((2,2,3)\) | 6 | 3 | 7 | 12 | 5 | 60 | -12 | -144 |
| \((2,2,2,2)\) | 7 | 4 | 8 | 16 | 4 | 64 | -8 | -128 |

As an external primary-source check, [LT, p. 2] derives \(c(TY)(1+3H)^2=(1+H)^6\) for the bicubic and obtains \(\chi_{\rm top}=-144\).

## 3. Grading and the two different localization conclusions

Let \(\mathrm{st}_X:\overline{\mathcal M}_g(X,D)\to\overline{\mathcal M}_g\) stabilize the domain, and let \(\pi,f\) be the universal curve and map. Define
\[
Z_D(t)=\mathrm{st}_{X,*}\left[
 \prod_{j=1}^c e_t\!\left((R\pi_*f^*\mathcal O(d_j))^\vee\right)
 \cap[\overline{\mathcal M}_g(X,D)]^{\rm vir}
\right].
\]
The Euler symbol is the multiplicative Euler ratio for a derived complex:
\[
e_t(V^0-V^1)=e_t(V^0)/e_t(V^1).
\]
In particular, the \(R^1\) term cannot be discarded. All summands use the **same** \(t\), corresponding to CJR's \(E=E_1\), not different R-charges indexed by the polynomial degrees \(d_j\). The degrees belong to line bundles on \(X\). This is exactly the stable zero-vertex specialization in [CJR-III, §9.4].

The compact target's ordinary virtual dimension is
\[
\boxed{\nu=(1-g)(m-3)+ID.}
\]
Indeed, deformation theory gives the source contribution \(3g-3\) and the map index
\(\chi(f^*TY)=m(1-g)+c_1(TY)\cdot\beta\). Their sum is \(\nu\).
Equivalently, the ambient virtual dimension and twisting rank are
\[
v_X=(1-g)(N-3)+(N+1)D,\qquad
\chi(f^*E)=c(1-g)+sD,
\]
and \(v_X-\chi(f^*E)=\nu\).

Since \(\dim\overline{\mathcal M}_g=3g-3\), the total twisted pushforward codimension is
\[
\boxed{\alpha=3g-3-\nu=m(g-1)-ID.}
\]
After forgetting the auxiliary projective-space torus, homogeneity therefore gives
\[
Z_D(t)=\sum_{r=0}^{3g-3}t^{\alpha-r}Z_{D,r},
\qquad Z_{D,r}\in A^r(\overline{\mathcal M}_g).
\]
The auxiliary weights are removed first; this is a statement in \(A^*(\overline{\mathcal M}_g)[t,t^{-1}]\).

Let \(B_D(t)\) denote the pushforward of the sum of **all** CJR localization terms containing infinity data. It includes mixed zero/infinity graphs, their gluing maps, curve-class splittings and multiplicities. It is not an individual punctured effective cycle. The all-zero term is \(Z_D(t)\), so write
\[
F_D^{\rm red}(t)=Z_D(t)+B_D(t).
\]
By [CJR-III, Theorem 8.4 and Remark 8.5], \(F_D^{\rm red}(t)\) has no negative powers of \(t\). With \(B_{D,r}=[t^{\alpha-r}]B_D(t)\), this implies only
\[
\boxed{B_{D,r}=-Z_{D,r}\quad\text{when }0\le r\le3g-3,\ r>\alpha.}
\]

Define the ordinary compact pushforward
\[
G_D=\sum_{\beta:\,i_*\beta=D[\ell]}
 \mathrm{st}_{Y,*}[\overline{\mathcal M}_g(Y,\beta)]^{\rm vir}
 \in A^\alpha(\overline{\mathcal M}_g).
\]
The compact comparison instead gives
\[
\boxed{[t^0](Z_D+B_D)=
 (-1)^{c(1-g)+sD}G_D.}
\]
When \(0\le\alpha\le3g-3\), this means
\[
B_{D,\alpha}=(-1)^{c(1-g)+sD}G_D-Z_{D,\alpha}.
\]
It does not give \(B_{D,\alpha}=-Z_{D,\alpha}\) without a separate vanishing theorem. If \(\alpha\) lies outside the Chow range, the stabilized \(G_D\) is zero by degree; this alone is not a general vanishing of target GW theory with arbitrary markings and insertions.

For orientation, the genus-specific values are:

| Target dimension | \(g=2\): \(\alpha\) | \(g=3\): \(\alpha\) |
|---|---|---|
| Fano surface | \(2-ID\) | \(4-ID\) |
| Fano threefold | \(3-ID\) | \(6-ID\) |
| Fano fourfold | \(4-ID\) | \(8-ID\) |
| CY3 | \(3\) | \(6\) |

For any of the five CY3s, \(I=0\) and \(\alpha=3g-3\), so \(\alpha-r\ge0\) throughout the available range. Thus this unmarked setup supplies no negative-\(t\) aggregate relations at any \(D\). The comparison is at \(r=3g-3\); lower-codimension terms have positive \(t\)-powers and need not cancel. This does not rule out relations from other markings, insertions, pushforwards, or auxiliary-torus poles.

## 4. What the higher-rank CJR statement establishes

The rank-\(c\) choice in [CJR-III, §9.4] is
\[
P=\mathbb P_{X\times B\mathbb C_\omega^*}
 (E^\vee\otimes L_\omega\oplus\mathcal O),
\qquad r=a=\widetilde r=1,\qquad L=\mathcal O_X,
\]
where \(E\) is pulled back from \(X\). Its infinity fiber over a point of \(B\mathbb C_\omega^*\) is \(\mathbb P_X(E^\vee)\). The compact-type marking condition is vacuous for the present unmarked data. Its displayed nonequivariant comparison reads, in the present notation,
\[
\left.\mathrm{st}_*[U_{g,\emptyset}(P,D)]^{\rm red}\right|_{t=0}
 =(-1)^{c(1-g)+sD}\,
 i_*[\overline{\mathcal M}_g(Y,D)]^{\rm vir},
\]
where the right side includes the ambient-degree union of target curve classes. The sign is also independently explained by Riemann–Roch:
\[
\chi(f^*E)=c(1-g)+sD.
\]
The public [CJR-GLSM, Proposition 4.1/Corollary 1.10] and [CJW, Corollary 1.2.1] support this arbitrary-rank comparison independently of the private localization formula.

The canonical and reduced **log-GLSM** theories are distinct. The canonical theory is constructed from the hybrid target; the reduced theory uses the superpotential and its cosection modification. The stated ordinary-target comparison uses the reduced log-GLSM class. The word “reduced” here does not mean reduced K3 GW theory, nor a reduced higher-genus virtual class of \(Y\). The class of \(Y\) on the right is its **ordinary** stable-map virtual class.

For \(c>1\), \(\mathbb P_X(E^\vee)\) has positive-dimensional fibers. The infinity theory retains its curve classes and tautological \(\mathcal O(1)\), and the rank-one simplifications of §8.7/§9.3 cannot simply be multiplied together. In contrast, the zero-vertex twist factors over the direct sum \(E=\bigoplus_j\mathcal O(d_j)\). This distinction permits the same ambient GP backend for \(Z_D\), but does not identify each infinity invariant from the aggregate \(B_D\).

Consequently, the answer to the source-scope question is **yes for this supplied private draft**: its §9.4 expressly establishes the needed higher-rank specialization, and Theorem 8.4 supplies the aggregate localization identity. No “hypersurface-only conjectural extension” caveat is appropriate here. The caveat that remains is the private-draft provenance of the localization statement.

## 5. Degree zero: compact targets versus ambient zero-vertex theory

### Independent compact obstruction calculation

For \(g\ge2\), the constant-map space and virtual class are
\[
\overline{\mathcal M}_g(Y,0)=\overline{\mathcal M}_g\times Y,\qquad
[\overline{\mathcal M}_g(Y,0)]^{\rm vir}
=c_{mg}(\mathbb E^\vee\boxtimes TY)
 \cap[\overline{\mathcal M}_g\times Y].
\]
This is the constant-map input of [FP, introduction]. To derive the requested pushforward independently, let \(x_1,\ldots,x_m\) be formal Chern roots of \(TY\), and set
\[
P_g(x)=\sum_{k=0}^g(-1)^k\lambda_kx^{g-k}.
\]
Then \(G_0=\int_Y\prod_{i=1}^mP_g(x_i)\), meaning the coefficient of target degree \(m\). Mumford's Chow relation [M] gives, in the genera used here,
\[
\lambda_1^2=2\lambda_2,\qquad
\lambda_2^2=2\lambda_1\lambda_3,\qquad
\lambda_3^2=0,
\]
with classes beyond the Hodge-bundle rank set to zero.

For a surface in genus two, the target-degree-two part of
\(\prod_{i=1}^2(x_i^2-\lambda_1x_i+\lambda_2)\) is
\[
\lambda_2(x_1^2+x_2^2)+\lambda_1^2x_1x_2
=\lambda_2c_1(Y)^2.
\]
In genus three, the degree-two part of
\(\prod_{i=1}^2(x_i^3-\lambda_1x_i^2+\lambda_2x_i-\lambda_3)\)
is
\[
\lambda_1\lambda_3(x_1^2+x_2^2)+\lambda_2^2x_1x_2
=\lambda_1\lambda_3c_1(Y)^2.
\]
Therefore
\[
\boxed{
G_0^{g=2}=\left(\int_Yc_1^2\right)\lambda_2,\qquad
G_0^{g=3}=\left(\int_Yc_1^2\right)\lambda_1\lambda_3
\quad(m=2).}
\]
The coefficients are \(8\) for the quadric surface and \(3\) for the cubic surface. These have codimensions \(2\) and \(4\), respectively; they are not unmarked scalar invariants until paired with complementary classes.

For a threefold in genus two, put
\(m_{21}=\sum_{i\ne j}x_i^2x_j=c_1c_2-3c_3\). The target-degree-three part is
\[
-\lambda_1\lambda_2m_{21}-\lambda_1^3c_3
=\lambda_1\lambda_2(c_3-c_1c_2).
\]
In genus three, the terms of target degree three are
\[
\lambda_3^2\sum_i x_i^3+
\lambda_1\lambda_2\lambda_3m_{21}+\lambda_2^3c_3
=\lambda_1\lambda_2\lambda_3(c_1c_2-c_3).
\]
Thus, for \(g=2,3\),
\[
\boxed{
G_0=(-1)^g\,C_Y\,\lambda_g\lambda_{g-1}\lambda_{g-2},
\qquad C_Y=\int_Y(c_3-c_1c_2),\quad \lambda_0=1.}
\]
This directly verifies the sign and the absence of an extra factor \(1/2\) in this basis. In the alternative basis \(\lambda_{g-1}^3\), there **is** a factor \(1/2\), since
\(\lambda_{g-1}^3=2\lambda_g\lambda_{g-1}\lambda_{g-2}\) in these genera.

| Threefold | \(C_Y\) |
|---|---:|
| Fano quadric | -20 |
| Fano cubic | -30 |
| Fano quartic | -80 |
| CY \((5)\) | -200 |
| CY \((2,4)\) | -176 |
| CY \((3,3)\) | -144 |
| CY \((2,2,3)\) | -144 |
| CY \((2,2,2,2)\) | -128 |

For CY3s, \(c_1=0\), so \(C_Y=\chi_{\rm top}(Y)\). As an optional scalar normalization check, [FP, Theorem 4] yields
\[
\int_{\overline{\mathcal M}_2}\lambda_2\lambda_1=\frac1{5760},
\qquad
\int_{\overline{\mathcal M}_3}\lambda_3\lambda_2\lambda_1
=\frac1{1451520}.
\]
Hence the **ordinary** degree-zero threefold GW numbers are \(C_Y/5760\) in genus two and \(-C_Y/1451520\) in genus three. For example, the quintic gives \(-5/144\) and \(5/36288\). These nonzero values already preclude a blanket ordinary-GW vanishing assumption.

For a fourfold, \(\nu=1-g<0\) at \(D=0\). Its ordinary virtual class, and therefore \(G_0\), vanishes for every \(g\ge2\). Equivalently, the predicted codimension \(4g-4\) exceeds \(3g-3\).

### Ambient degree-zero vanishing

For \(X=\mathbb P^N\), \(D=0\), the ambient virtual dimension is
\[
v_X=(1-g)(N-3).
\]
If \(N\ge4\) and \(g\ge2\), it is negative; thus the **nonequivariant ambient virtual class is zero**. Keeping the GLSM parameter \(t\) inverted, the dual Euler ratio is a well-defined operational class, so
\[
\boxed{Z_0(t)=0\quad(N\ge4,\ g\ge2)}
\]
after the auxiliary projective-space torus is forgotten. Individual auxiliary fixed graphs can be nonzero and singular in their weights; they must be summed before taking that limit.

This argument concerns the ambient zero-vertex contribution, not the compact target. For all threefolds here, compact \(G_0\) is nonzero even though ambient \(Z_0=0\). Their entire compact constant term comes from the infinity aggregate:
\[
[t^0]B_0(t)=(-1)^{c(1-g)}G_0.
\]
For threefolds in genus two this is \((-1)^c C_Y\lambda_1\lambda_2\); in genus three it is \(-C_Y\lambda_1\lambda_2\lambda_3\). Thus, for example, the quintic infinity constant term is \(+200\lambda_1\lambda_2\) in genus two and \(+200\lambda_1\lambda_2\lambda_3\) in genus three, before numerical integration. These signs include the CJR comparison sign.

The negative-ambient-dimension argument does not apply to the two surface examples, whose ambient space is \(\mathbb P^3\). In fact, for \(g=2,3\) one can additionally derive
\[
Z_0(t)=
-20(-1)^g\,t^{1-g}\lambda_g\lambda_{g-1}\lambda_{g-2}
\qquad(X=\mathbb P^3,\ c=1).
\]
Here the ambient constant-map virtual class has dimension zero; all positive-degree terms in the Euler-ratio expansion annihilate it, leaving its leading term \(t^{1-g}\). The remaining pushforward is the threefold formula for \(\mathbb P^3\), with \(\int(c_3-c_1c_2)=4-24=-20\). In particular, its \(t^0\) term is zero, while its negative-power term is not. The surface compact constant term is again supplied by infinity, with sign \((-1)^{1-g}\).

## 6. General projective-space GP graph factors for the dual twist

These are **auxiliary torus fixed graphs of maps to \(\mathbb P^N\)** used to calculate \(Z_D\). They are distinct from the CJR zero/infinity graphs defining \(B_D\). The following formulas are obtained by combining the GP deformation/normalization calculation with the normalization sequence for each line-bundle summand.

Choose \(H|_{p_i}=a_i\) and coordinate weights \(-a_i\). Define
\[
w_{ik}=a_i-a_k,\qquad T_i=\prod_{k\ne i}w_{ik},
\qquad q_{ij}=t-d_ja_i.
\]
The tangent weight at \(p_i\) along the line to \(p_k\) is \(w_{ik}\). For a flag \(F\) at \(i\), on an edge of degree \(e\) ending at \(k\), let \(w_F=(a_i-a_k)/e\). Write
\[
\Lambda_h^\vee(z)=\sum_{\ell=0}^h(-1)^\ell\lambda_\ell z^{h-\ell},
\qquad
\Lambda_h^+(z)=\sum_{\ell=0}^h\lambda_\ell z^{h-\ell}.
\]

### Stable contracted vertex

For a vertex at \(p_i\) of genus \(h\) and unmarked valence \(n\), with \(2h-2+n>0\), its factor on \(\overline{\mathcal M}_{h,n}\) is
\[
\boxed{
V_{i,h,n}=
\left[\prod_{k\ne i}\Lambda_h^\vee(w_{ik})\right]
T_i^{\,n-1}
\left[\prod_{j=1}^c
 \frac{q_{ij}^{\,1-n}}{\Lambda_h^+(q_{ij})}\right]
\prod_{F\ni v}\frac1{w_F-\psi_F}.}
\]
The tangent normalization contributes \(T_i^{n-1}\) and the obstruction Hodge factors \(\Lambda_h^\vee\). For one twisting summand, a contracted component has derived dual
\((1-\mathbb E)\otimes\mathcal O(-d_j)|_{p_i}\), with Euler ratio \(q_{ij}/\Lambda_h^+(q_{ij})\). Normalizing each attached node removes one copy of its fiber, hence \(q_{ij}^{-n}\). This proves the displayed \(q_{ij}^{1-n}\) and explains why the Hodge polynomial in the twist has **plus** signs.

### Edge

For a degree-\(e\ge1\) edge between \(p_i,p_k\), keep the deck factor separate. Its remaining contribution is
\[
\boxed{
A_{ik,e}=
\frac{(-1)^e e^{2e}}
 {(e!)^2(a_i-a_k)^{2e}}
\left[
 \prod_{\ell\ne i,k}\prod_{u=0}^{e}
 \left(\frac{u}{e}a_i+\frac{e-u}{e}a_k-a_\ell\right)
\right]^{-1}
\prod_{j=1}^c\prod_{u=0}^{d_je}
\left(t-\frac{d_je-u}{e}a_i-\frac{u}{e}a_k\right).}
\]
The first two factors are inverse moving tangent weights. The final product is independently read from the \(d_je+1\) monomial weights of \(H^0(\mathbb P^1,\mathcal O(d_je))^\vee\). Endpoint weights are included. Both the tangent product and the twist are invariant under interchanging the endpoints and reversing the monomial index.

### Unstable vertices and multiplicities

With the same edge normalization, an unmarked genus-zero valence-one vertex contributes \(w_F\). A genus-zero valence-two vertex, representing the node between two edges rather than a contracted component, contributes
\[
\boxed{
V_{i,0,2}=
\frac{T_i}{(w_{F_1}+w_{F_2})\prod_{j=1}^c q_{ij}}.}
\]
This can also be obtained from the stable expression using the formal conventions
\[
\int_{\overline{\mathcal M}_{0,1}}\frac1{w-\psi}=w,\qquad
\int_{\overline{\mathcal M}_{0,2}}
 \frac1{(w_1-\psi_1)(w_2-\psi_2)}
=\frac1{w_1+w_2}.
\]
These are conventions encoding unstable deformation factors, not integrations over actual stable-curve stacks. A degree-zero graph in the present genus range has one stable genus-\(g\) vertex of valence zero, so no unstable-valence-zero convention is needed.

For each colored, degree-decorated graph \(\Gamma\), multiply the edge and vertex factors and include
\[
\frac1{|\operatorname{Aut}(\Gamma)|}\prod_{e\in E(\Gamma)}\frac1{e}.
\]
The separate \(1/e\) accounts for deck transformations; it must not be duplicated inside \(A_{ik,e}\). Adjacent vertex colors differ, \(D=\sum_e e\), and
\(g=\sum_v g(v)+b_1(\Gamma)\). To obtain a class on \(\overline{\mathcal M}_g\), use the actual stabilization/gluing pushforward from the fixed-locus product, including contraction of unstable rational components.

The auxiliary nonequivariant limit is taken only after the complete graph sum, while \(t\) remains a formal invertible parameter. These formulas compute the dual-twisted ambient contribution, not the full ordinary positive-genus CI theory. The latter additionally requires the infinity aggregate described above.

## 7. Optional genus-zero, degree-one check

For a general smooth CY3 in the five families, the degree-one ordinary GW invariants are:

| Type | Degree-one invariant |
|---|---:|
| \((5)\mathbb P^4\) | 2875 |
| \((2,4)\mathbb P^5\) | 1280 |
| \((3,3)\mathbb P^5\) | 1053 |
| \((2,2,3)\mathbb P^6\) | 720 |
| \((2,2,2,2)\mathbb P^7\) | 512 |

The four non-quintic counts appear in [LT, Table 1]; all five are also explicitly returned by [M2, examples i1–i3]. Literal enumeration of isolated lines is a general-member statement; on special smooth members the same numbers retain their virtual/deformation-invariant meaning.

For an independent exact check, lines are parameterized by \(\mathrm{Gr}(2,N+1)\), with tautological rank-two bundle \(S\). The defining equations give a section of
\(\bigoplus_j\mathrm{Sym}^{d_j}S^\vee\). Its rank \(\sum_j(d_j+1)=2N-2\) equals the Grassmannian dimension. If \(x,y\) are roots of \(S^\vee\), the top Chern number can be evaluated as
\[
\int_{\mathrm{Gr}(2,N+1)}
 c_{\rm top}\!\left(\bigoplus_j\mathrm{Sym}^{d_j}S^\vee\right)
=
[x^Ny^{N-1}]\,
(x-y)\prod_j\prod_{u=0}^{d_j}\bigl(ux+(d_j-u)y\bigr).
\]
Expanding this polynomial independently gives \(2875,1280,1053,720,512\), in the displayed order. The section/top-Chern interpretation uses [Stacks, Lemma 42.44.1].

This tests the **ordinary, nondual** genus-zero Euler twist \(e(R^0\pi_*f^*E)\), as allowed by convexity and [CJW, Theorem 1.2.3]. At \(g=0,D=1\), the dual twist differs at \(t=0\) by \((-1)^{c+s}\); for these five CY3 types, \(c+s=2c+4\) is even, so the numerical check coincides. This coincidence must not be used to erase the general dual-twist sign. The existing genus-\(2,3\) backend is not asserted to implement this optional genus-zero test.

## 8. Consequences for the requested generalization

- Use the complete lists in §2 and the single rank-\(c\) bundle \(E\), with one common GLSM parameter \(t\).
- Grade the unmarked output by \(t^{m(g-1)-ID-r}\). Only \(r>\alpha\) gives automatic negative-power cancellation.
- At \(r=\alpha\), retain \((-1)^{c(1-g)+sD}G_D\). Ordinary K3 vanishing is not a property of the present Fano or CY3 families.
- For unmarked CY3s, report the top compact comparison; there is no negative-\(t\) range.
- At \(D=0\), use \(Z_0=0\) for \(N\ge4\), while preserving the nonzero compact threefold constants and their CJR signs in the infinity aggregate.
- The product formulas in §6 justify the ambient localization extension. Higher-rank infinity factors still require the rank-\(c\) geometry of CJR §9.4; an aggregate identity does not determine each effective cycle separately.

All numerical Chern data, the genus-two/three obstruction expansions, and the five Grassmannian line counts were independently derived in this note. No implementation audit, additional file, agent, or commit is part of this research task.

## 9. Ordinary CY3 degree-one numbers

Faber–Pandharipande [FP, Theorem 3, printed p. 4] proves that a fixed rational curve with normal bundle \(\mathcal O(-1)\oplus\mathcal O(-1)\) contributes
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3}
\]
to the connected ordinary genus-\(g\), degree-\(k\) invariant. In particular, \(C(2,1)=1/240\) and \(C(3,1)=1/6048\).

For the five complete intersections, the bundle \(\bigoplus_j\operatorname{Sym}^{d_j}S^\vee\) on \(\operatorname{Gr}(2,N+1)\) is globally generated and has rank equal to the Grassmannian dimension. General defining equations give a transverse zero-dimensional scheme of lines. For each such line, \(H^0(N_{\ell/Y})=0\); adjunction gives \(\deg N_{\ell/Y}=-2\), so the normal bundle is \(\mathcal O(-1)^{\oplus2}\). Every degree-one map has a line as its image. Therefore, if \(L_Y\) is the line count in §7,
\[
N^Y_{2,1}=\frac{L_Y}{240},\qquad N^Y_{3,1}=\frac{L_Y}{6048}.
\]
These equalities first hold for a general smooth member and extend to smooth deformations by GW deformation invariance. They include contracted higher-genus components and do not count embedded higher-genus degree-one curves. This is an application of FP's local-curve theorem and the line geometry, independent of the ambient twisted computation. No degree-two or degree-three ordinary GW numbers are inferred from the line counts.
