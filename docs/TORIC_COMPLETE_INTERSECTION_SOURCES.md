# Toric complete intersections: bounded source audit

Research date: 2026-11-06. Work over \(\mathbb C\), in rational Chow, with connected unmarked stable maps of genus \(g=2,3\). This audit covers exactly the twelve proposed targets, ambient Mori tuples with nonnegative coefficients of sum at most two, and two generic auxiliary-torus specializations. It does not classify toric Fanos, audit the implementation, or compute the proposed localization/Pixton scan.

All twelve families admit general smooth members of the stated Fano, K3, or CY3 type. The supplied CJR III draft expressly supports their arbitrary-rank ambient comparison. The proposed toric edge and vertex formulas are correct with the conventions below. Independent intersection calculations confirm the CY3 Euler characteristics \(-168,-162,-144\) and the requested primitive curve counts \(64,189,189,9\).

The distinction between ordinary target GW theory, the dual-twisted ambient term, and the full CJR infinity aggregate remains essential. Unknown ordinary positive-degree classes must remain unknown. In particular, the three CY3s give no negative-\(t\) range after this unmarked stabilization to \(\overline{\mathcal M}_g\).

## 1. Sources and provenance

The mathematical statements below are labeled as source inputs or deductions. The numeric Chern tables, weight example, curve counts, and applications of the local-curve formula are independent deductions in this audit. The existing [projective audit](PROJECTIVE_COMPLETE_INTERSECTION_SOURCES.md) and [quartic K3 audit](QUARTIC_K3_SOURCES.md) were read for conventions; the relevant formulas were checked against the sources below. The vertex and edge conventions in the existing projective and quartic fixed-graph files were also inspected, without changing or testing that code.

- **[CJR-III]** Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, supplied private draft dated October 17, 2025. [Supplied PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>), read through the supplied extraction /tmp/log-glsm-audit/cjr-iii.txt. Printed-page locators: §2.1.2, p. 7 (smooth proper base and graded vector bundle); Assumption 4.2 (compact-type markings); Theorem 8.4, equation (8.28), and Remark 8.5, p. 84 (localization and regularity in \(t\)); §9.4, p. 92 (arbitrary-rank complete intersections, compact comparison, stable zero-vertex twist); §10.1, p. 93 (ambient-class sums and Theorem 10.1). Corollary 10.2, p. 93, adds ampleness for a finiteness assertion. Claims involving this localization theorem retain the provenance of this private version.
- **[CJR]** Chen–Janda–Ruan, *The logarithmic gauged linear sigma model*, Invent. Math. 225 (2021), 1077–1154; [arXiv:1906.04345v2](https://arxiv.org/pdf/1906.04345v2). Corollary 1.10, p. 7, and Proposition 4.1, pp. 34–35, give the reduced-log-GLSM comparison with the ordinary complete-intersection virtual class.
- **[CJW]** Q. Chen, F. Janda, R. Webb, *Virtual cycles of stable (quasi-)maps with fields*, Adv. Math. 385 (2021), 107781; [arXiv:1911.09825v2](https://arxiv.org/pdf/1911.09825v2). Theorem 1.1.1, p. 3, and Corollary 1.2.1, p. 4, give the arbitrary-rank fields comparison with sign \((-1)^{\chi(f^*E)}\), including the sum over target classes of fixed ambient pushforward. Theorem 1.2.3, pp. 4–5, gives convex quantum Lefschetz.
- **[LS]** C.-C. M. Liu, A. Sheshmani, *Equivariant Gromov–Witten Invariants of Algebraic GKM Manifolds*, SIGMA 13 (2017), 048; [publisher PDF](https://sigma-journal.com/2017/048/sigma17-048.pdf), [arXiv:1407.1370](https://arxiv.org/abs/1407.1370). §2.2 describes edge weights and normal degrees. Lemma 4.5, equations (4.7)–(4.8), pp. 16–17, gives the edge factor, including negative normal degrees. Proposition 4.4, equation (4.6), gives the contracted-vertex Hodge factor. §4.3.4 and Theorem 4.7, equation (4.11), p. 19, specify graph/deck factors and unstable conventions. These are toric/GKM formulas, so no projective-space-only extrapolation is needed.
- **[CLS]** D. Cox, J. Little, H. Schenck, *Toric Varieties*; [author book page](https://dacox.people.amherst.edu/toric.html), [consulted book PDF](https://pi.math.cornell.edu/~tsh/cornell-only/cox-little-schenck-toric.pdf). Theorem 6.3.12, printed p. 291, identifies nef, globally generated, and nonnegative invariant-curve degree for complete toric varieties; Theorem 6.3.22, p. 295, characterizes ampleness by the interior of the nef cone. Theorem 8.1.6, p. 363, is the toric Euler sequence. Example 10.4.6, p. 498, computes Hirzebruch-surface intersections.
- **[MPT]** D. Maulik, R. Pandharipande, R. P. Thomas, with appendix by A. Pixton, *Curves on K3 surfaces and modular forms*; [author PDF](https://people.math.ethz.ch/~rahul/mpt.pdf), [arXiv:1001.2719](https://arxiv.org/abs/1001.2719). Introduction §0.1, p. 2: the ordinary virtual class vanishes for every nonzero effective K3 class because the obstruction theory has a trivial quotient; reduced K3 theory removes that quotient.
- **[M]** D. Mumford, *Towards an Enumerative Geometry of the Moduli Space of Curves*; [author scan](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf). §5, equations (5.2), (5.4), Corollary (5.3), printed pp. 305–307: the algebraic GRR relation implies \(c(\mathbb E)c(\mathbb E^\vee)=1\) in rational Chow on the compactification.
- **[FP]** C. Faber, R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. 139 (2000), 173–199; [arXiv:math/9810173](https://arxiv.org/pdf/math/9810173). Introduction equation (1), p. 1: the constant-map obstruction bundle. Equation (7), p. 3, and Theorem 3, p. 4: the all-genus local \(\mathcal O(-1)\oplus\mathcal O(-1)\) contribution. Theorem 4, p. 4, and §4.3, pp. 20–21: the needed Hodge integrals.
- **[Stacks]** [Chow Homology and Chern Classes, Lemma 42.44.1, tag 0FA8](https://stacks.math.columbia.edu/tag/0FA8): a regular section represents the top Chern class. This is the intersection-theoretic input for the curve counts in §7.

## 2. Geometry of the twelve targets

Write \(A,B\) for the hyperplane divisors on a product of projective spaces. On \(F_1\times\mathbb P^n\), write \(S,F,H\) for the negative section, fiber, and projective hyperplane. The same symbols \(S,F\) denote the corresponding curve classes on the surface when used in a Mori tuple. Thus
\[
S^2=-SF,\qquad F^2=0,\qquad \int_{F_1}SF=1,
\qquad -K_{F_1}=2S+3F.
\]

For \(Y=Z(s)\subset X\), \(s\) a transverse section of \(E=\bigoplus_jL_j\), adjunction gives
\[
m=\dim Y=\dim X-\operatorname{rk}E,\qquad
c_1(TY)=\left(c_1(TX)-\sum_jL_j\right)|_Y.
\]

The following table is the complete sample, not a completeness claim about ambient varieties or complete intersections.

| ID | Ambient \(X\) | Defining bundle \(E\) | \(m\) | \(c_1(TY)\), restricted from \(X\) | Type |
|---|---|---|---:|---|---|
| T1 | \(\mathbb P^1\times\mathbb P^2\) | \(\mathcal O(1,1)\) | 2 | \(A+2B\) | Fano surface |
| T2 | \(\mathbb P^1\times\mathbb P^2\) | \(\mathcal O(2,3)\) | 2 | \(0\) | K3 |
| T3 | \(\mathbb P^1\times\mathbb P^3\) | \(\mathcal O(1,2)\) | 3 | \(A+2B\) | Fano threefold |
| T4 | \(\mathbb P^1\times\mathbb P^3\) | \(\mathcal O(2,4)\) | 3 | \(0\) | CY3 |
| T5 | \(\mathbb P^2\times\mathbb P^2\) | \(\mathcal O(1,2)\) | 3 | \(2A+B\) | Fano threefold |
| T6 | \(\mathbb P^2\times\mathbb P^2\) | \(\mathcal O(3,3)\) | 3 | \(0\) | CY3 |
| T7 | \(\mathbb P^2\times\mathbb P^2\) | \(\mathcal O(1,1)^{\oplus2}\) | 2 | \(A+B\) | Fano surface |
| T8 | \(\mathbb P^2\times\mathbb P^2\) | \(\mathcal O(1,1)\oplus\mathcal O(2,2)\) | 2 | \(0\) | K3 |
| T9 | \(F_1\times\mathbb P^1\) | \(\mathcal O(S+F+H)\) | 2 | \(S+2F+H\) | Fano surface |
| T10 | \(F_1\times\mathbb P^1\) | \(\mathcal O(2S+3F+2H)\) | 2 | \(0\) | K3 |
| T11 | \(F_1\times\mathbb P^2\) | \(\mathcal O(S+F+H)\) | 3 | \(S+2F+2H\) | Fano threefold |
| T12 | \(F_1\times\mathbb P^2\) | \(\mathcal O(2S+3F+3H)\) | 3 | \(0\) | CY3 |

All ambients are smooth toric Fanos of dimension three or four. On \(F_1\), a divisor \(aS+bF\) has degrees \(b-a\) on \(S\) and \(a\) on \(F\); it is nef when \(b\ge a\ge0\), and ample when \(b>a>0\). This follows from the two generators of its Mori cone and [CLS, Theorems 6.3.12 and 6.3.22]. In particular \(2S+3F\) and \(S+2F\) are ample. Every nonzero \(c_1(TY)\) in the table is the restriction of an ample divisor, proving the Fano assertions.

All summands of \(E\) are globally generated. For a globally generated bundle, the incidence variety of pairs \((x,s)\) with \(s(x)=0\) is the kernel bundle of the evaluation map over smooth \(X\). Generic smoothness in characteristic zero gives transverse general sections. The displayed positive top intersections ensure nonempty zero loci. This argument also covers rank two and the two nonample bundles.

For T9 and T11, \(S+F+H\) has degree zero on \(S\times\{p\}\), so it is nef but not ample. More concretely, \(S+F\) is the pullback of a line under \(F_1\to\mathbb P^2\). A general \((1,1)\) divisor in \(\mathbb P^2\times\mathbb P^n\) meets \(\{p\}\times\mathbb P^n\) transversely in \(\mathbb P^{n-1}\); its inverse image is its smooth blowup along that intersection. This also verifies connectedness for the nonample cases.

For the anticanonical hypersurfaces, the sequence \(0\to K_X\to\mathcal O_X\to\mathcal O_Y\to0\), Serre duality, and \(H^i(X,\mathcal O_X)=0\) for \(i>0\) give connectedness and \(H^i(Y,\mathcal O_Y)=0\) for \(0<i<m\). Here the ambient vanishing can be computed directly from projective spaces and the \(\mathbb P^1\)-bundle \(F_1\to\mathbb P^1\). Adjunction gives \(K_Y\simeq\mathcal O_Y\). Thus the surfaces are K3 and the threefolds have the usual CY3 cohomology vanishings. For T8, the Koszul complex with terms \(\mathcal O(-3,-3)\), \(\mathcal O(-1,-1)\oplus\mathcal O(-2,-2)\), and \(\mathcal O_X\) gives \(h^0(\mathcal O_Y)=h^2(\mathcal O_Y)=1,\ h^1(\mathcal O_Y)=0\); adjunction again gives a K3. These are deductions for the listed families.

## 3. CJR scope, degree tuples, and grading

[CJR-III, §9.4] takes a transverse section of an arbitrary-rank bundle \(E=E_1\), with \(L=\mathcal O_X\) in the hybrid-target data and \(r=a=\widetilde r=1\). This applies to the smooth projective toric ambients above, including T7 and T8 and the nef, nonample T9 and T11. The compact-type marking assumption is vacuous here. All line-bundle summands use the same GLSM parameter \(t\). Their divisor coefficients are not distinct R-charges.

This applicability does not use Corollary 10.2: that later finiteness statement assumes \(E\) ample and \(H^1(\mathbb P_X(E^\vee))=0\), and its ampleness hypothesis fails for T9 and T11. The finite auxiliary graph search here instead follows from the explicit Mori bound. For rank two, the infinity target is \(\mathbb P_X(E^\vee)\), with positive-dimensional fibers; its theory cannot be reconstructed by multiplying rank-one infinity formulas. The zero-vertex twist does factor over the summands.

Let \(\beta\) denote an ambient class. For products of projective spaces use \(\beta=(d_1,d_2)\) in the line bases, so \((aA+bB)\cdot\beta=ad_1+bd_2\). For \(F_1\times\mathbb P^n\), use \(\beta=(d_S,d_F,d_H)\), so
\[
(aS+bF+cH)\cdot\beta=(b-a)d_S+a d_F+c d_H.
\]
In particular
\[
c_1(TX)\cdot\beta=d_S+2d_F+(n+1)d_H,\qquad
(S+F+H)\cdot\beta=d_F+d_H.
\]
For the nonample Fano examples, \(c_1(TY)\cdot\beta=d_S+d_F+n d_H\). For every CY or K3 example it is zero. On T7, \(c_1(E)\cdot\beta=2(d_1+d_2)\); on T8 it is \(3(d_1+d_2)\).

Set \(c=\operatorname{rk}E\), \(m=\dim Y\), and define
\[
Z_\beta(t)=\mathrm{st}_{X,*}
\left(e_t((R\pi_*f^*E)^\vee)\cap
[\overline{\mathcal M}_g(X,\beta)]^{\mathrm{vir}}\right).
\]
The Euler class here is an Euler ratio of the derived complex; \(R^1\) is retained. Let \(B_\beta(t)\) be the pushforward of the full sum of CJR terms involving infinity vertices, with all splittings and gluing. Let
\[
G_{g,\beta}=
\sum_{\beta':\,i_*\beta'=\beta}
\mathrm{st}_{Y,*}[\overline{\mathcal M}_g(Y,\beta')]^{\mathrm{vir}}.
\]
Ambient tuples need not distinguish target classes, especially for surfaces. This sum is part of the definition, as in [CJR-III, §10.1] and [CJW, Corollary 1.2.1].

Riemann–Roch and adjunction give the following grading independently:
\[
\nu=(1-g)(m-3)+c_1(TY)\cdot\beta,\qquad
\boxed{\alpha=3g-3-\nu=m(g-1)-c_1(TY)\cdot\beta},
\]
\[
\chi(f^*E)=c(1-g)+c_1(E)\cdot\beta,\qquad
\boxed{\sigma_{g,\beta}=(-1)^{c(1-g)+c_1(E)\cdot\beta}}.
\]
Here \(\alpha\) is the codimension of \(G_{g,\beta}\). The sign agrees with the compact comparison in [CJR-III, §9.4] and the public [CJR, Proposition 4.1].

After removing the auxiliary torus,
\[
Z_\beta(t)=\sum_{r=0}^{3g-3}t^{\alpha-r}Z_{\beta,r}.
\]
Regularity of \(Z_\beta+B_\beta\) in \(t\), from [CJR-III, Remark 8.5], yields
\[
B_{\beta,r}=-Z_{\beta,r}\quad(r>\alpha),\qquad
[t^0](Z_\beta+B_\beta)=\sigma_{g,\beta}G_{g,\beta}.
\]
The \(t^0\) identity is an ordinary-GW comparison, not an automatic vanishing. If \(\alpha\) is outside \([0,3g-3]\), the stabilized compact class is zero by grading; this does not imply vanishing for other markings or insertions.

For K3 targets T2, T8, T10, [MPT, §0.1] gives ordinary virtual vanishing in every nonzero effective class. The degree-zero calculation in §6 gives vanishing in genera two and three as well. Thus \(G_{g,\beta}=0\), \(\alpha=2g-2\), and the aggregate equality extends to \(r\ge2g-2\). Reduced K3 invariants are different invariants. CJR's reduced log-GLSM class compares to ordinary K3 theory.

For CY3 targets T4, T6, T12, \(\alpha=3g-3\) at every degree. There is no \(r>\alpha\) in this pushforward. Their top-codimension comparison must retain ordinary GW theory, including the nonzero constants and the specific positive-degree checks below.

## 4. Auxiliary toric graph factors

These graphs are torus-fixed stable maps to \(X\), used to compute \(Z_\beta\). They are distinct from the CJR zero/infinity graphs. The untwisted factors agree with [LS, Lemma 4.5 and Theorem 4.7]. The twist and the equivalent Hodge-polynomial expression below follow directly from normalization and [M].

### Edges, including negative normal degree

Orient an invariant curve \(C\simeq\mathbb P^1\) from \(i\) to \(j\). Let \(w\) be its tangent weight at \(i\); the weight at \(j\) is \(-w\). For a normal summand \(\mathcal O_C(a)\) with fiber weight \(\ell\) at \(i\), compatible linearization gives
\[
\ell_j=\ell-aw.
\]
For a cover of degree \(e\), the source tangent weight is \(u=w/e\). The monomial weights of \(H^0(\mathcal O(ae))\), for \(ae\ge0\), are \(\ell-ru\), \(0\le r\le ae\). For \(ae<0\), a Čech basis of \(H^1\) gives weights with \(ae+1\le r\le-1\). Hence the inverse Euler ratio is
\[
\boxed{\mathcal N(a,\ell,w,e)=
\begin{cases}
\displaystyle\prod_{r=0}^{ae}(\ell-rw/e)^{-1},&ae\ge0,\\[4pt]
\displaystyle\prod_{r=ae+1}^{-1}(\ell-rw/e),&ae<0.
\end{cases}}
\]
Empty products equal one. There is no extra sign in the \(H^1\) numerator.

For the tangent \(\mathcal O(2e)\), omit its single fixed weight. The remaining weights are \(\{\pm kw/e:1\le k\le e\}\), so their inverse product is
\[
\boxed{\mathcal T(w,e)=
\frac{(-1)^e e^{2e}}{(e!)^2w^{2e}}}.
\]
The fixed tangent deformation cancels the fixed source automorphism in this convention. Deck transformations contribute a separate \(1/e\).

For a nef twisting line bundle \(L\), let \(d=L\cdot C\ge0\), and let \(L_i,L_j\) be its fiber weights. Then \(L_j=L_i-dw\). Dualizing the \(ed+1\) section weights gives
\[
\boxed{\mathcal Q_L(i,C,e)=
\prod_{r=0}^{ed}(t-L_i+rw/e)}.
\]
Endpoint factors are included. When \(d=0\), the product is the single factor \(t-L_i\). The edge factor before deck symmetry is
\(\mathcal T\prod_{\text{normal summands}}\mathcal N\prod_{L\subset E}\mathcal Q_L\).
Reversing orientation replaces \((w,\ell,L_i)\) by \((-w,\ell-aw,L_i-dw)\); reversing the product indices proves invariance.

### Vertices

Let a contracted stable vertex at \(i\) have genus \(h\) and valence \(n\). Write \(\tau_{i,1},\ldots,\tau_{i,N}\) for the ambient tangent weights, \(T_i=\prod_k\tau_{i,k}\), and \(q_{ij}=t-L_{j,i}\). Define
\[
\Lambda_h^\vee(z)=\sum_{k=0}^h(-1)^k\lambda_kz^{h-k},
\qquad
\Lambda_h^+(z)=\sum_{k=0}^h\lambda_kz^{h-k}.
\]
With \(u_F=w_F/e_F\) the source tangent weight for each flag,
\[
\boxed{V_{i,h,n}=
T_i^{n-1}\prod_{k=1}^N\Lambda_h^\vee(\tau_{i,k})
\prod_{j=1}^c\left(q_{ij}^{1-n-2h}\Lambda_h^\vee(q_{ij})\right)
\prod_F(u_F-\psi_F)^{-1}.}
\]

To check the twist, a contracted component contributes \(q/\Lambda_h^+(q)\) from \((R\Gamma(\mathcal O)\otimes L_i)^\vee\). Normalizing its \(n\) attached nodes supplies \(q^{-n}\). Mumford's relation gives
\[
\Lambda_h^+(q)\Lambda_h^\vee(q)=q^{2h},
\qquad
\frac{q^{1-n}}{\Lambda_h^+(q)}
=q^{1-n-2h}\Lambda_h^\vee(q).
\]
Thus the proposed dual-polynomial form agrees with the plus-polynomial denominator in the existing source audit. It requires this relation in the Hodge ring; replacing the denominator by a dual polynomial without its power of \(q\) would be wrong.

With the edge convention above, an unmarked genus-zero univalent vertex contributes \(u_F\), and a bivalent node contributes
\[
\boxed{V_{i,0,2}=
\frac{T_i}{(u_{F_1}+u_{F_2})\prod_jq_{ij}}.}
\]
The weights in the denominator are cover-adjusted flag weights. Multiply the graph by
\[
\frac{1}{|\operatorname{Aut}\Gamma|}\prod_{E\in E(\Gamma)}\frac1{e_E}.
\]
Then apply the actual stabilization/gluing pushforward; unstable components cannot be replaced by extra stable factors.

### Explicit \(F_1\) check and Mori bookkeeping

Let \(S_+=S+F\), and let \(F_0,F_\infty\) be the invariant fibers. Label the four fixed points
\[
A=S\cap F_0,\quad B=S\cap F_\infty,\quad
C=S_+\cap F_0,\quad D=S_+\cap F_\infty.
\]
Choose tangent weights \(u\) along \(S\) and \(v\) along \(F_0\) at \(A\). A compatible table is

| Point | Two tangent weights |
|---|---|
| \(A\) | \(u,v\) |
| \(B\) | \(-u,v+u\) |
| \(C\) | \(u,-v\) |
| \(D\) | \(-u,-v-u\) |

The normal degrees in \(F_1\) are \(-1\) on \(S\), \(+1\) on \(S_+\), and \(0\) on each fiber. For instance along \(A\to B\), the normal weight changes from \(v\) to \(v+u=v-(-1)u\). Along \(C\to D\), it changes from \(-v\) to \(-v-u\). Each extra product-factor tangent direction is a degree-zero normal summand.

For \(S\), a degree-one cover has \(\mathcal N(-1,v,u,1)=1\), since \(H^0(\mathcal O(-1))=H^1(\mathcal O(-1))=0\). A double cover has
\[
\mathcal N(-1,v,u,2)=v+u/2.
\]
This belongs in the numerator. The positive section instead has the two degree-one denominator factors \((-v)^{-1}(-v-u)^{-1}\). Fibers have one factor \(\ell^{-1}\) for each degree-zero normal line. For an edge in the \(\mathbb P^n\) factor, the normal degrees are \(1\) repeated \(n-1\) times and \(0\) repeated twice.

As another check, normalize the weight of \(\mathcal O(aS+bF)\) at \(A\) to zero. The weights at \((A,B,C,D)\) are
\[
(0,\;-(b-a)u,\;-av,\;-bu-av).
\]
They give degrees \(b-a,a,b\) on \(S,F,S_+\), respectively. Add \(cH\)'s projective weight to obtain a product linearization.

The curve classes of invariant edges in \(F_1\times\mathbb P^n\) are
\[
[S]=(1,0,0),\quad [F_0]=[F_\infty]=(0,1,0),\quad
[S_+]=(1,1,0),\quad [\text{projective line}]=(0,0,1).
\]
For any auxiliary graph the required class is
\[
\boxed{\beta(\Gamma)=\sum_E e_E[C_E],\qquad
g=\sum_vg(v)+b_1(\Gamma).}
\]
In particular a degree-one edge on \(S_+\) and a two-edge chain on \(S\cup F_0\) both have class \((1,1,0)\), although their sums of cover degrees are one and two. The cut \(\sum d_i\le2\) admits both; it excludes a double cover of \(S_+\), whose class is \((2,2,0)\).

There are six ambient tuples including zero for each rank-two Mori lattice and ten for each rank-three lattice. The sample has eight product targets and four \(F_1\)-product targets: 88 target/class pairs, or 352 target/class/genus/ray evaluations if all are run. Every primitive invariant edge has positive coefficient sum, so the bound implies at most two nonconstant edges. Enumeration must still retain all genus decorations and both possible \(S+F\) realizations.

The two requested torus rays are two generic one-parameter specializations of the full auxiliary character lattice, not two selected invariant curves. Check every tangent, normal cohomology, and bivalent flag-sum denominator for resonance. Sum all graphs before taking the auxiliary nonequivariant limit, with \(t\) still formal. Agreement of two rays is a useful check, not a proof of equality of unspecialized expressions or a certificate of a new tautological relation.

## 5. Independent Chern numbers

The calculations in this section use ordinary product Chow rings and the normal exact sequence, independently of localization:
\[
A^*(\mathbb P^p\times\mathbb P^q)
=\mathbb Q[A,B]/(A^{p+1},B^{q+1}),\qquad
\int A^pB^q=1,
\]
\[
A^*(F_1\times\mathbb P^n)
=\mathbb Q[S,F,H]/(S^2+SF,F^2,H^{n+1}),\qquad
\int SFH^n=1.
\]
The Euler sequence [CLS, Theorem 8.1.6] and the four invariant divisors \(S,S+F,F,F\) give
\[
c(T(\mathbb P^p\times\mathbb P^q))=(1+A)^{p+1}(1+B)^{q+1},
\]
\[
c(T(F_1\times\mathbb P^n))
=(1+S)(1+S+F)(1+F)^2(1+H)^{n+1}
=(1+2S+3F+4SF)(1+H)^{n+1}.
\]
Then
\[
c(TY)=\frac{c(TX)}{\prod_j(1+L_j)}\bigg|_Y,\qquad
\int_Y P(c(TY))=\int_X\left(\prod_jL_j\right)P(c(TY)).
\]
These equations are also a reproducible prescription: expand each \((1+L_j)^{-1}\) through ambient dimension, reduce by the displayed ideal, and extract the indicated top monomial. Exact Sage arithmetic in the native Linux guest verified every row.

| Surface | \(\int c_1^2\) | \(\int c_2=\chi_{\mathrm{top}}\) |
|---|---:|---:|
| T1: \((1,1)\subset\mathbb P^1\times\mathbb P^2\) | 8 | 4 |
| T2: \((2,3)\subset\mathbb P^1\times\mathbb P^2\) | 0 | 24 |
| T7: \([(1,1),(1,1)]\subset\mathbb P^2\times\mathbb P^2\) | 6 | 6 |
| T8: \([(1,1),(2,2)]\subset\mathbb P^2\times\mathbb P^2\) | 0 | 24 |
| T9: \(S+F+H\subset F_1\times\mathbb P^1\) | 7 | 5 |
| T10: \(2S+3F+2H\subset F_1\times\mathbb P^1\) | 0 | 24 |

| Threefold | \(\int c_1^3\) | \(\int c_1c_2\) | \(\int c_3=\chi_{\mathrm{top}}\) | \(C_Y=\int(c_3-c_1c_2)\) |
|---|---:|---:|---:|---:|
| T3: \((1,2)\subset\mathbb P^1\times\mathbb P^3\) | 32 | 24 | 4 | -20 |
| T4: \((2,4)\subset\mathbb P^1\times\mathbb P^3\) | 0 | 0 | -168 | -168 |
| T5: \((1,2)\subset\mathbb P^2\times\mathbb P^2\) | 30 | 24 | 6 | -18 |
| T6: \((3,3)\subset\mathbb P^2\times\mathbb P^2\) | 0 | 0 | -162 | -162 |
| T11: \(S+F+H\subset F_1\times\mathbb P^2\) | 42 | 24 | 8 | -16 |
| T12: \(2S+3F+3H\subset F_1\times\mathbb P^2\) | 0 | 0 | -144 | -144 |

These are all Chern numbers in dimensions two and three. For explicit verification of the three CY3 Euler values, the ambient representatives of their Chern classes are
\[
\begin{array}{c|l|l}
&c_2(TY)&c_3(TY)\\ \hline
T4&8AB+6B^2&-32AB^2-20B^3\\
T6&3A^2+9AB+3B^2&-27A^2B-27AB^2\\
T12&4SF+6SH+9FH+3H^2&-24SFH-18SH^2-27FH^2.
\end{array}
\]
Multiplying the last column by \(2A+4B\), \(3A+3B\), and \(2S+3F+3H\), respectively, gives top coefficients \(-168,-162,-144\). Additional checks are
\[
(\int c_2A,\int c_2B)_{T4}=(24,44),\qquad
(\int c_2A,\int c_2B)_{T6}=(36,36),
\]
\[
(\int c_2S,\int c_2F,\int c_2H)_{T12}=(12,24,36).
\]
For the two rank-two surfaces the ambient \(c_2\) representatives are \(3AB\) for T7 and \(A^2+5AB+B^2\) for T8, integrated against \((A+B)^2\) and \(2(A+B)^2\).

## 6. Constant-map comparisons

The source input is [FP, introduction (1)]:
\[
[\overline{\mathcal M}_g(Y,0)]^{\mathrm{vir}}
=c_{mg}(\mathbb E^\vee\boxtimes TY)
\cap[\overline{\mathcal M}_g\times Y].
\]
Here is an independent expansion for the present genera. If \(x_i\) are the target Chern roots, take the target-degree-\(m\) part of
\(\prod_i\Lambda_g^\vee(x_i)\), and integrate over \(Y\). Use [M]'s identities
\(\lambda_1^2=2\lambda_2\), \(\lambda_2^2=2\lambda_1\lambda_3\), and the rank truncations. One obtains
\[
\boxed{
\begin{aligned}
\dim Y=2:\quad&
G_{2,0}=(\int_Yc_1^2)\lambda_2,\qquad
G_{3,0}=(\int_Yc_1^2)\lambda_1\lambda_3;\\
\dim Y=3:\quad&
G_{2,0}=C_Y\lambda_1\lambda_2,\qquad
G_{3,0}=-C_Y\lambda_1\lambda_2\lambda_3.
\end{aligned}}
\]
For example, in genus two on a surface the target-degree-two polynomial is
\(\lambda_2(x_1^2+x_2^2)+\lambda_1^2x_1x_2=\lambda_2c_1^2\).
For a threefold it is
\(-\lambda_1\lambda_2(c_1c_2-3c_3)-\lambda_1^3c_3
=\lambda_1\lambda_2(c_3-c_1c_2)\).
The genus-three expansion gives the signs displayed above.

Thus all ordinary compact constants are fixed by §5. In particular the K3 constants vanish. In fact the full constant-map virtual class vanishes for these K3s in genera two and three: target-degree-zero terms contain \(\lambda_g^2=0\), target-degree-one terms contain \(c_1=0\), and the target-degree-two term above vanishes. Combined with [MPT], this covers every ambient degree.

For threefolds, [FP, Theorem 4 and §4.3] gives
\[
\int_{\overline{\mathcal M}_2}\lambda_1\lambda_2=\frac1{5760},
\qquad
\int_{\overline{\mathcal M}_3}\lambda_1\lambda_2\lambda_3
=\frac1{1451520}.
\]
Hence the ordinary scalar constants are \(C_Y/5760\) and \(-C_Y/1451520\). These signs precede the additional CJR factor \(\sigma_{g,0}=(-1)^{c(1-g)}\).

For an ambient fourfold, \(\mathrm{vdim}\,\overline{\mathcal M}_g(X,0)=1-g<0\). Consequently \(Z_0(t)=0\) after forgetting the auxiliary torus. This applies to T3–T8 and T11–T12, including the two surface complete intersections. It does not imply that \(G_{g,0}\) vanishes: the signed compact constant must then come entirely from \([t^0]B_0\).

For the ambient threefolds of T1–T2 and T9–T10, one can also check the nonzero ambient constant term. Both ambients have \(\int c_1c_2=24\); their Euler characteristics are six and eight. Set \(C_X=\int_X(c_3-c_1c_2)\), so \(C_X=-18\) for \(\mathbb P^1\times\mathbb P^2\) and \(-16\) for \(F_1\times\mathbb P^1\). Their ambient constant virtual cycle has dimension zero, so only the leading term \(t^{1-g}\) of the rank-one dual twist acts nontrivially. Therefore
\[
Z_0(t)=t^{1-g}G^X_{g,0},\qquad
G^X_{2,0}=C_X\lambda_1\lambda_2,\quad
G^X_{3,0}=-C_X\lambda_1\lambda_2\lambda_3.
\]
This is a further independently derived ambient check, valid after the full auxiliary graph sum.

## 7. Independent primitive rational-curve counts

The following are deductions for general smooth members, using parameter spaces of curves and transverse sections. The top-Chern interpretation is [Stacks, Lemma 42.44.1]. They were not imported from mirror-symmetry tables.

### T4: \(\mathbb P^1\)-fiber class \((1,0)\)

Write the equation of a \((2,4)\) hypersurface as
\[
x_0^2f_0(y)+x_0x_1f_1(y)+x_1^2f_2(y)=0,
\qquad f_i\in H^0(\mathbb P^3,\mathcal O(4)).
\]
A fiber \(\mathbb P^1\times\{y\}\) lies in \(Y\) exactly when all three quartics vanish. Three general quartics meet transversely in
\[
\int_{\mathbb P^3}c_3(\mathcal O(4)^{\oplus3})=4^3=64
\]
points. These are all images of maps of ambient class \((1,0)\), since projection to \(\mathbb P^3\) has degree zero.

For such a fiber \(C\), \(N_{C/X}=\mathcal O^{\oplus3}\) and \(L|_C=\mathcal O(2)\). Transversality makes the map on global sections in the normal sequence an isomorphism \(\mathbb C^3\to H^0(\mathcal O(2))\). Hence \(H^0(N_{C/Y})=0\). Adjunction gives \(\deg N_{C/Y}=-2\), so splitting on \(\mathbb P^1\) forces \(N_{C/Y}=\mathcal O(-1)^{\oplus2}\).

### T6: bicubic classes \((1,0)\) and \((0,1)\)

Curves of class \((1,0)\) have the form \(\ell\times\{p\}\). Their parameter space is
\(\mathrm{Gr}(2,3)\times\mathbb P^2\), of dimension four. If \(U\) is the tautological rank-two subbundle and \(h\) is the hyperplane class on the second factor, restriction of the equation is a section of
\[
V=\mathrm{Sym}^3U^\vee\boxtimes\mathcal O(3),
\qquad\operatorname{rk}V=4.
\]
The space of bicubic equations generates \(V\), so a general section has transverse isolated zeros. Let \(x,y\) be the roots of \(U^\vee\). Its symmetric-cube roots are \(3x,2x+y,x+2y,3y\). If \(p=c_1(U^\vee)\), then on \(\mathrm{Gr}(2,3)\simeq\mathbb P^{2\vee}\),
\[
x+y=p,\quad xy=p^2,\quad \int p^2=1,
\]
\[
c_2(\mathrm{Sym}^3U^\vee)
=11(x^2+y^2)+32xy=21p^2.
\]
Only the \(h^2\) term contributes to \(\int c_4(V)\); it is \(9h^2c_2(\mathrm{Sym}^3U^\vee)\). Thus the count is
\[
\boxed{189}.
\]
Equivalently, the coefficient of \(x^2yh^2\) in
\((x-y)\prod_{r=0}^3(3h+(3-r)x+ry)\) is 189; exact expansion independently checked this value. Interchanging the two factors gives 189 for \((0,1)\).

Here \(N_{C/X}=\mathcal O(1)\oplus\mathcal O^{\oplus2}\) and \(L|_C=\mathcal O(3)\). Transversality gives \(H^0(N_{C/Y})=0\), and the degree is again \(-2\), proving the \(\mathcal O(-1)^{\oplus2}\) splitting. Degree one in one factor and zero in the other exhausts the possible images.

### T12: negative-section class \((1,0,0)\)

Restrict the anticanonical equation to \(S\times\mathbb P^2\). Its bidegree there is \((1,3)\), since \((2S+3F)\cdot S=1\). It is therefore
\[
x_0g_0(y)+x_1g_1(y)=0
\]
with two cubic forms on \(\mathbb P^2\). Their coefficients can be general: under \(F_1=\operatorname{Bl}_p\mathbb P^2\), anticanonical sections are cubics through \(p\), whose first-order terms at \(p\) span \(H^0(S,\mathcal O_S(1))\). Thus restriction onto the displayed bidegree space is surjective. The curves \(S\times\{y\}\subset Y\) occur at the transverse common zeros of two cubics:
\[
\int_{\mathbb P^2}c_2(\mathcal O(3)^{\oplus2})=9.
\]
The negative section is the unique effective curve in its primitive class on \(F_1\); hence these exhaust the images of ambient class \((1,0,0)\).

Now \(N_{C/X}=\mathcal O(-1)\oplus\mathcal O^{\oplus2}\), while \(L|_C=\mathcal O(1)\). The map on global sections is the transverse derivative of the two cubic coefficients, an isomorphism \(\mathbb C^2\to H^0(\mathcal O(1))\). The same normal-sequence argument proves \(N_{C/Y}=\mathcal O(-1)^{\oplus2}\).

### Ordinary genus-zero and higher-genus implications

Nef line bundles have nonnegative degree on every component of a genus-zero tree, so \(H^1(C,f^*E)=0\). By [CJW, Theorem 1.2.3], the ordinary genus-zero invariant is computed by the nondual Euler twist. The four counts above therefore provide independent genus-zero localization checks.

At \(g=0\), the dual twist at \(t=0\) has sign \((-1)^{c+c_1(E)\cdot\beta}\). For T4 in \((1,0)\) it is \((-1)^3=-1\), for the two T6 classes it is \((-1)^4=+1\), and for T12 in \((1,0,0)\) it is \((-1)^2=+1\). A dual-twisted genus-zero check should consequently return \(-64,+189,+189,+9\). The ordinary counts are all positive.

[FP, equation (7) and Theorem 3] gives the local contribution of a fixed \((-1,-1)\) curve:
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3},
\qquad C(2,1)=\frac1{240},\quad C(3,1)=\frac1{6048}.
\]
For each primitive class above, nonnegative Mori coefficients rule out a decomposition into two positive ambient classes, and the curve geometry exhausts all possible images. Stable maps may have contracted higher-genus components; those are included in \(C(g,1)\). Therefore the connected ordinary numbers, summed over any target classes with that ambient pushforward, are

| Target and ambient class | \(N_{0,\beta}=L\) | \(N_{2,\beta}=L/240\) | \(N_{3,\beta}=L/6048\) |
|---|---:|---:|---:|
| T4, \((1,0)\) | 64 | \(4/15\) | \(2/189\) |
| T6, \((1,0)\) | 189 | \(63/80\) | \(1/32\) |
| T6, \((0,1)\) | 189 | \(63/80\) | \(1/32\) |
| T12, \((1,0,0)\) | 9 | \(3/80\) | \(1/672\) |

These are ordinary invariants, before the CJR comparison sign. Literal isolated-curve geometry is asserted for general smooth equations; the invariant values extend to smooth members by deformation invariance.

There is one additional bounded consequence with an explicit exhaustion argument. Every map of pure class \(k(1,0)\) in T4 has image one of the same 64 fibers, because its projection to \(\mathbb P^3\) is constant. Every map of pure class \(k(1,0,0)\) in T12 has image one of the same nine negative sections: zero projective degree fixes \(y\), and intersection with \(S\) forces every positive component of a pure multiple of \(S\) to be supported on \(S\). The fibers/sections are disjoint, so a connected map cannot use two of them. Thus within the requested bound \(k=2\), the full ordinary numbers are \(L\,C(g,2)\):

| Target and class | \(N_{2,\beta}\) | \(N_{3,\beta}\) |
|---|---:|---:|
| T4, \((2,0)\) | \(8/15\) | \(16/189\) |
| T12, \((2,0,0)\) | \(3/40\) | \(1/84\) |

The same conclusion has not been proved here for doubles of a bicubic ruling: plane conics are additional possible images, and the line count alone does not determine their contributions. This audit leaves T4's \((0,1),(0,2),(1,1)\), T6's \((2,0),(0,2),(1,1)\), and T12's other positive tuples unknown. In particular the \(F_1\) fiber class, the projective-line class, and \(S+F\) cannot inherit the negative-section count. No positive-degree Fano ordinary values are inferred from these CY calculations.

## 8. What is established for the bounded search

The source inputs justify applying the CJR comparison and toric ambient localization to this sample, with the rank-two infinity geometry retained. The independent calculations establish all twelve Chern-number rows, all compact constants, the four primitive CY curve counts, and the two pure double-cover checks just proved. Exact guest-side Sage arithmetic also checked reversal invariance of the normal factor for normal degrees \(-2,-1,0,1,2\) and cover degrees one and two.

A future implementation still needs to enumerate the vector-valued curve classes correctly, apply stabilization, and compare the two auxiliary specializations after summing. No scan result or Pixton reduction has been certified by this source audit. Any remaining \(t^0\) compact term must be retained symbolically unless it is covered by the K3 vanishing, grading, constant-map formulas, or explicit ordinary CY calculations above.
