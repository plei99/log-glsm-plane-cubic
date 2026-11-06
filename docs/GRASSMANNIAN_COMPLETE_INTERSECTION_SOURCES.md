# Grassmannian complete intersections: bounded source audit

Research date: 2026-11-06. Work over \(\mathbb C\), in rational Chow, with connected unmarked stable maps of genus \(g=2,3\) and Plücker degree \(D=0,1,2\). This audit covers the fourteen specified complete intersections in \(\operatorname{Gr}(2,4)\) and \(\operatorname{Gr}(2,5)\). It does not run the localization/Pixton relation search or classify other homogeneous complete intersections.

The extension is supported: the supplied CJR III draft permits smooth projective nontoric ambients and arbitrary-rank defining bundles; Grassmannians are explicit examples in the algebraic GKM localization source. All fourteen proposed families have the stated Fano, K3, or CY3 geometry. The proposed normal splitting and weight assignment are correct. Independent intersection calculations below verify every Chern-number row and give CY3 line counts \(1280,540,400\). Primary sources also supply the genus-zero conic numbers \(92288,12555,5540\), with the multiple-cover convention specified below.

The ordinary target GW class, dual-twisted ambient class, and full CJR infinity aggregate remain separate objects. For these CY3s there is **no negative-\(t\) range after unmarked stabilization to \(\overline{\mathcal M}_g\)**. The quadric realization of \(\operatorname{Gr}(2,4)\) gives valid geometric and genus-zero checks, but does not identify higher-genus ambient Euler twists.

## 1. Sources, locators, and provenance

The [toric audit](TORIC_COMPLETE_INTERSECTION_SOURCES.md) and [projective audit](PROJECTIVE_COMPLETE_INTERSECTION_SOURCES.md) supply the notation convention. Their claims are not treated as substitutes for the primary sources. Calculations explicitly labeled deductions below were performed independently for this audit.

- **[CJR-III]** Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, supplied private draft dated October 17, 2025. [Supplied PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>), consulted through the existing extraction `/tmp/log-glsm-audit/cjr-iii.txt`. Printed locators: §2.1.2, p. 7 (smooth proper DM base with projective coarse space); Assumption 4.2 (compact-type markings); Theorem 8.4, equation (8.28), and Remark 8.5, p. 84; §9.4, p. 92 (arbitrary-rank complete intersections); §10.1 and Theorem 10.1, p. 93 (ambient-class sums). This audit uses the theorem as stated in this private version; it does not certify its proof or cite it as a published result.
- **[CJW]** Q. Chen, F. Janda, R. Webb, *Virtual cycles of stable (quasi-)maps with fields*, Adv. Math. 385 (2021), 107781; [arXiv:1911.09825v2](https://arxiv.org/pdf/1911.09825v2). Theorem 1.1.1, p. 3, and Corollary 1.2.1, p. 4: arbitrary-rank fields comparison, its sign, and target-class sums. Theorem 1.2.3, p. 4: convex quantum Lefschetz. Introduction p. 2 explicitly explains the failure of naive quantum Lefschetz in positive genus.
- **[LS]** C.-C. M. Liu, A. Sheshmani, *Equivariant Gromov–Witten Invariants of Algebraic GKM Manifolds*, SIGMA 13 (2017), 048; [publisher PDF](https://sigma-journal.com/2017/048/sigma17-048.pdf). Example 2.3, p. 3, treats Grassmannians, including their fixed points and invariant lines. §2.2 describes normal degrees and endpoint weights. **Lemma** 4.4, equation (4.6), p. 15, gives vertex Hodge factors; Lemma 4.5, equations (4.7)–(4.8), pp. 16–17, gives edge factors. Theorem 4.7, equation (4.11), p. 19, includes graph/deck factors and unstable conventions. The prior toric audit calls 4.4 a proposition; the consulted PDF labels it a lemma.
- **[GP]** T. Graber, R. Pandharipande, *Localization of virtual classes*, Invent. Math. 135 (1999), 487–518; [arXiv:alg-geom/9708001v2](https://arxiv.org/pdf/alg-geom/9708001v2). Introduction equation (1): virtual localization for equivariant perfect obstruction theories; §4: stable-map graph calculations. This is not a theorem restricted to toric targets.
- **[BCFKvS]** V. Batyrev, I. Ciocan-Fontanine, B. Kim, D. van Straten, *Conifold Transitions and Mirror Symmetry for Calabi–Yau Complete Intersections in Grassmannians*, Nuclear Phys. B 514 (1998), 640–666; [arXiv:alg-geom/9710022v1](https://arxiv.org/pdf/alg-geom/9710022v1), [DOI](https://doi.org/10.1016/S0550-3213(98)00020-0). Preprint §2.1, p. 6 onward: quartics in \(G(2,4)\); §2.2: \(G(2,5)[1,1,3]\). The table on p. 29 records CY3 Euler characteristics and Hodge numbers. Theorem 6.2.1, p. 30, and tables §7.1–§7.2, p. 31, give virtual rational-curve numbers. The Yukawa expansion immediately before the theorem fixes the instanton convention.
- **[Kim]** B. Kim, *Quantum hyperplane section theorem for homogeneous spaces*, Acta Math. 183 (1999), 71–99; [article PDF](https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/6582-11511_2006_Article_BF02392947.pdf). Theorem 1, p. 73, and Remark (4), p. 74: genus-zero hyperplane theorem and verification of the Grassmannian rational-curve predictions in [BCFKvS]. §7 treats Grassmannians. This gives a mathematical genus-zero justification, rather than only an unaudited mirror prediction.
- **[M2]** Official Macaulay2 `EnumerationCurves` documentation, [rationalCurve](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/EnumerationCurves/html/_rational__Curve.html). Examples i2–i3 fix the ordering of the projective CY complete intersections and their line numbers; i5 gives degree-two GW values; i7 subtracts one eighth of the line number to obtain conic numbers. This is a first-party computational reference, not a new geometric proof in this audit.
- **[MPT]** D. Maulik, R. Pandharipande, R. P. Thomas, with appendix by A. Pixton, *Curves on K3 surfaces and modular forms*; [arXiv:1001.2719v4](https://arxiv.org/pdf/1001.2719v4). Introduction §0.1, p. 2: the ordinary virtual class vanishes in every nonzero effective K3 class; reduced K3 theory removes a trivial obstruction quotient.
- **[FP]** C. Faber, R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. 139 (2000), 173–199; [arXiv:math/9810173](https://arxiv.org/pdf/math/9810173). Introduction equation (1), p. 1: constant-map virtual class. Equation (7), p. 3, and Theorem 3, p. 4: local \((-1,-1)\)-curve contributions. Theorem 4, p. 4, and §4.3, pp. 20–21: Hodge integrals used below.
- **[M]** D. Mumford, *Towards an Enumerative Geometry of the Moduli Space of Curves*; [author-hosted scan](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf). §5, equations (5.2), (5.4), Corollary (5.3), printed pp. 305–307: \(c(\mathbb E)c(\mathbb E^\vee)=1\) in rational Chow on the compactification. The same GRR input is used in [FP, §4.3].
- **[Stacks]** [Chow Homology and Chern Classes, §42.44, Lemma 42.44.1](https://stacks.math.columbia.edu/tag/0FA8): a regular section represents the top Chern class, the intersection-theoretic input for the line calculation.

## 2. Geometry and all Chern numbers of the sample

Put
\[
X=\operatorname{Gr}(2,n),\quad N=2(n-2),\quad
H=c_1(\det S^\vee),\quad E=\bigoplus_{j=1}^{c}\mathcal O_X(d_j),
\]
\[
s=\sum_jd_j,\qquad m=N-c,\qquad \kappa=n-s.
\]
Here \(S\) is the tautological rank-two subbundle and \(\mathcal O_X(1)=\det S^\vee\) is the Plücker bundle. The tautological sequence identifies \(TX=S^\vee\otimes Q\), giving \(c_1(TX)=nH\). The normal sequence and adjunction give
\[
K_Y=\mathcal O_Y(s-n),\qquad
c(TY)=\frac{c(TX)}{\prod_j(1+d_jH)}\bigg|_Y,
\qquad [Y]=\left(\prod_jd_j\right)H^c.
\]

All \(d_j\) are positive, so all summands are very ample and general equations meet transversely. Bertini gives smooth members. For completeness, connectedness and the CY cohomology condition can also be checked with the Koszul complex: for \(0<A<n\), the bundle \(\mathcal O_X(-A)=K_X\otimes\mathcal O_X(n-A)\) has no higher cohomology by Kodaira vanishing and no global sections by negativity. The final term when \(s=n\) is \(K_X\), with its sole nonzero cohomology in degree \(N\). Thus the Koszul resolution gives \(h^0(\mathcal O_Y)=h^m(\mathcal O_Y)=1\) and \(h^i(\mathcal O_Y)=0\) for \(0<i<m\). The \(m=2\) cases are K3 surfaces, and the \(m=3\) cases are CY3s, not merely varieties with numerically trivial first Chern class. For \(s<n\), \(-K_Y=\kappa H\) is ample, proving the Fano statements. These are deductions for the specified families.

### Independent calculation prescription

Let \(x,y\) be the Chern roots of \(S^\vee\), so \(H=x+y\), \(P=xy=c_2(S^\vee)\). In K-theory,
\[
[TX]=n[S^\vee]-[S^\vee\otimes S],\qquad
\operatorname{ch}(TX)=n(e^x+e^y)-(2+e^{x-y}+e^{y-x}).
\]
Expanding and applying Newton identities gives
\[
\begin{aligned}
c_1(TX)&=nH,\\
c_2(TX)&=\frac{n^2-n+2}{2}H^2+(n-4)P,\\
c_3(TX)&=\frac{n(n^2-3n+8)}6H^3+n(n-5)HP.
\end{aligned}
\]
For example, \(c_2(T\operatorname{Gr}(2,5))=11H^2+P\): retaining only powers of \(H\) would lose information needed for the target Chern numbers.

For a symmetric polynomial of top degree,
\[
\boxed{\int_{\operatorname{Gr}(2,n)}f(x,y)
=[x^{n-1}y^{n-2}]\,(x-y)f(x,y).}
\]
This is the two-variable Schur alternant integration rule. It gives
\[
\int_{G(2,4)}H^4=2,\quad \int_{G(2,4)}H^2P=1;
\qquad
\int_{G(2,5)}H^6=5,\quad \int_{G(2,5)}H^4P=2.
\]
To reproduce the table, expand the normal-sequence quotient through degree \(m\), multiply by \((\prod d_j)H^c\), and apply this coefficient rule. Every row was checked with exact guest-side SymPy arithmetic and independently by fixed-point integration using the tangent weights of §4; this does not use the repository's stable-map localization implementation.

| ID | Ambient | Degrees | \(m\) | \(\kappa\) | \(\int_YH^m\) | Type |
|---|---|---|---:|---:|---:|---|
| G1 | \(G(2,4)\) | \((1)\) | 3 | 3 | 2 | Fano threefold, quadric \(Q^3\) |
| G2 | \(G(2,4)\) | \((2)\) | 3 | 2 | 4 | Fano threefold, \((2,2)\subset\mathbb P^5\) |
| G3 | \(G(2,4)\) | \((3)\) | 3 | 1 | 6 | Fano threefold, \((2,3)\subset\mathbb P^5\) |
| G4 | \(G(2,4)\) | \((4)\) | 3 | 0 | 8 | CY3, \((2,4)\subset\mathbb P^5\) |
| G5 | \(G(2,4)\) | \((1,1)\) | 2 | 2 | 2 | Fano surface, quadric \(Q^2\) |
| G6 | \(G(2,4)\) | \((1,2)\) | 2 | 1 | 4 | Fano surface, degree-four del Pezzo |
| G7 | \(G(2,4)\) | \((1,3)\) | 2 | 0 | 6 | K3, polarization genus 4 |
| G8 | \(G(2,4)\) | \((2,2)\) | 2 | 0 | 8 | K3, polarization genus 5 |
| G9 | \(G(2,5)\) | \((1,1,1)\) | 3 | 2 | 5 | Fano threefold \(V_5\) |
| G10 | \(G(2,5)\) | \((1,1,2)\) | 3 | 1 | 10 | Fano threefold of degree 10 |
| G11 | \(G(2,5)\) | \((1,1,3)\) | 3 | 0 | 15 | CY3 |
| G12 | \(G(2,5)\) | \((1,2,2)\) | 3 | 0 | 20 | CY3 |
| G13 | \(G(2,5)\) | \((1,1,1,1)\) | 2 | 1 | 5 | Fano surface, degree-five del Pezzo |
| G14 | \(G(2,5)\) | \((1,1,1,2)\) | 2 | 0 | 10 | K3, polarization genus 6 |

The polarization genus of a K3 is \(1+H^2/2\); it is unrelated to the stable-map genus \(g\).

| Surface | \(\int c_1^2\) | \(\int c_2=\chi_{\rm top}\) |
|---|---:|---:|
| G5 | 8 | 4 |
| G6 | 4 | 8 |
| G7 | 0 | 24 |
| G8 | 0 | 24 |
| G13 | 5 | 7 |
| G14 | 0 | 24 |

| Threefold | \(\int c_1^3\) | \(\int c_1c_2\) | \(\int c_3=\chi_{\rm top}\) | \(\int Hc_2\) | \(C_Y=\int(c_3-c_1c_2)\) |
|---|---:|---:|---:|---:|---:|
| G1 | 54 | 24 | 4 | 8 | -20 |
| G2 | 32 | 24 | 0 | 12 | -24 |
| G3 | 6 | 24 | -36 | 24 | -60 |
| G4 | 0 | 0 | -176 | 56 | -176 |
| G9 | 40 | 24 | 4 | 12 | -20 |
| G10 | 10 | 24 | -16 | 24 | -40 |
| G11 | 0 | 0 | -150 | 66 | -150 |
| G12 | 0 | 0 | -120 | 68 | -120 |

These include all Chern numbers in dimensions two and three. The extra \(Hc_2\) column provides a polarization check. For G11, the ambient representatives are \(c_2=8H^2+P\), \(c_3=-14H^3-5HP\); multiplying the latter by \(3H^3\) gives \(-150\). For G12 they are \(c_2=7H^2+P\), \(c_3=-4H^3-5HP\); multiplying by \(4H^3\) gives \(-120\). For G4, the projective model gives \(c_2=7H^2\), \(c_3=-22H^3\), hence \(-176\).

The three CY3 Euler values agree with [BCFKvS, p. 29], which also gives \(h^{2,1}=89,76,61\), respectively, with \(h^{1,1}=1\).

## 3. Why CJR applies, and what its comparison means

The base condition in [CJR-III, §2.1.2] is smooth and proper with projective coarse space, not toric. Section 9.4 takes a transverse section of a rank-\(c\) bundle on a projective variety. Thus take \(X=G(2,n)\), \(E=E_1=\bigoplus_j\mathcal O_X(d_j)\), \(L=\mathcal O_X\), and \(r=a=\widetilde r=1\). All summands have the same GLSM grading; the integers \(d_j\) are Plücker degrees, not different R-charges. The marking condition is vacuous in the unmarked problem.

Theorem 8.4 and §9.4 therefore apply directly, with no toric degeneration required. The auxiliary torus used to compute stable maps to the ambient Grassmannian is separate from CJR's parameter \(t\). General defining equations need not be invariant under that auxiliary torus. For \(c>1\), the CJR infinity target is \(\mathbb P_X(E^\vee)\); its fiber geometry must be retained. Multiplying independent rank-one infinity formulas does not reconstruct this target. Only the Euler twist factors immediately over the summands.

Use \(D=\int_\beta H\) for the ambient curve class, and define
\[
G_{g,D}=\sum_{\beta':\,i_*\beta'=D[\mathrm{line}]}
\mathrm{st}_{Y,*}[\overline{\mathcal M}_g(Y,\beta')]^{\rm vir}.
\]
This sum is necessary for surfaces; Plücker degree alone need not identify a target curve class. It matches [CJR-III, §10.1] and [CJW, Corollary 1.2.1]. Define separately
\[
Z_D(t)=\mathrm{st}_{X,*}\left(
e_t\big((R\pi_*f^*E)^\vee\big)
\cap[\overline{\mathcal M}_g(X,D)]^{\rm vir}\right),
\]
and let \(B_D(t)\) be the pushforward of the **full** CJR infinity aggregate, including its splittings and gluing. The derived Euler ratio retains \(R^1\); positivity of \(E\) does not make it convex for genus two or three.

Riemann–Roch and adjunction give
\[
\nu_Y=(1-g)(m-3)+\kappa D,\qquad
\boxed{\alpha=3g-3-\nu_Y=m(g-1)-\kappa D},
\]
\[
\boxed{\sigma_{g,D}=(-1)^{\chi(f^*E)}
=(-1)^{c(1-g)+sD}}.
\]
After removing the auxiliary torus, write
\[
Z_D(t)=\sum_{r=0}^{3g-3}t^{\alpha-r}Z_{D,r},
\quad Z_{D,r}\in A^r(\overline{\mathcal M}_g).
\]
The grading and [CJR-III, Remark 8.5] yield
\[
\boxed{B_{D,r}=-Z_{D,r}\quad(r>\alpha)},\qquad
\boxed{[t^0](Z_D+B_D)=\sigma_{g,D}G_{g,D}}.
\]
Here \(B_{D,r}\) means the coefficient of \(t^{\alpha-r}\). Cancellation of negative powers and comparison at exponent zero are distinct assertions. Neither says that the infinity aggregate or the ordinary target GW class vanishes by itself. No result of a tautological/Pixton reduction is being asserted here.

The bounded graph search does not depend on Corollary 10.2. Its hypotheses are compatible with these positive bundles, but the supplied draft's proof still contains an unresolved reference; Theorem 8.4, §9.4, and the explicit degree bound suffice for the audit's application.

## 4. Grassmannian fixed points, invariant lines, and local factors

Let \(T=(\mathbb C^*)^n\) act diagonally on \(\mathbb C^n\), with vector weights \(\lambda_1,\ldots,\lambda_n\). The diagonal subtorus acts trivially on the base. Its action on the chosen bundle linearization is harmless; this audit fixes the natural linearization \(\mathcal O(1)=\det S^\vee\).

In [LS, Example 2.3], fixed points are two-subsets \(I\subset\{1,\ldots,n\}\), and an invariant line joins \(I\) to \(J\) exactly when one element is exchanged. Write
\[
I=\{a,c\},\qquad J=\{b,c\},\qquad
K=\{1,\ldots,n\}\setminus\{a,b,c\}.
\]
Then the tangent and Plücker weights, derived from \(TX=S^\vee\otimes Q\), are
\[
\boxed{T_IX:\ \lambda_j-\lambda_i\quad(i\in I,\ j\notin I)},
\qquad w_{I\to J}=\lambda_b-\lambda_a,
\]
\[
\boxed{H_I=-\sum_{i\in I}\lambda_i},\qquad
H_J=H_I-w,\qquad L_{d,I}=dH_I.
\]
Thus the proposed signs agree: a degree-\(d\) bundle's endpoint weights differ by \(-dw\).

### Derivation of the normal splitting, with each source weight

Along this line, the varying two-plane is \(\langle e_c\rangle\oplus\ell\), where \(\ell\) varies in \(\langle e_a,e_b\rangle\). Consequently
\[
S|_C=\mathcal O\oplus\mathcal O(-1),\qquad
Q|_C=\mathcal O(1)\oplus\mathcal O^{n-3}.
\]
Tensoring gives
\[
TX|_C=\mathcal O(2)\oplus\mathcal O(1)^{n-2}
\oplus\mathcal O^{n-3},
\]
\[
\boxed{N_{C/X}=\mathcal O(1)^{n-2}\oplus\mathcal O^{n-3}.}
\]
The full equivariant assignment is as follows; it is derived here, not assumed from the nonequivariant splitting.

| Summand | Weight at \(I\) | Weight at \(J\) | Multiplicity |
|---|---|---|---:|
| Tangent \(\mathcal O(2)\) | \(w=\lambda_b-\lambda_a\) | \(-w\) | 1 |
| Normal \(\mathcal O(1)\) | \(\lambda_k-\lambda_a\), \(k\in K\) | \(\lambda_k-\lambda_b\) | \(n-3\) |
| Normal \(\mathcal O(1)\) | \(\lambda_b-\lambda_c\) | \(\lambda_a-\lambda_c\) | 1 |
| Normal \(\mathcal O(0)\) | \(\lambda_k-\lambda_c\), \(k\in K\) | same | \(n-3\) |

Each normal degree \(a_N\) satisfies \(\ell_J=\ell_I-a_Nw\). In particular, the source tangent direction \((i=c,j=b)\) is a degree-one normal, not degree zero. This is easy to misassign if only the unordered tangent-weight set is stored. The rank is \(2n-5=N-1\), and the total normal degree is \(n-2=c_1(TX)\cdot C-2\).

### Explicit edge and vertex expressions

For an edge cover of degree \(e\), the domain tangent weight is \(u=w/e\). An \(\mathcal O(a_N)\) normal with source weight \(\ell\) contributes
\[
\mathcal N(a_N,\ell,w,e)=
\prod_{k=0}^{a_Ne}(\ell-kw/e)^{-1},\qquad a_N\in\{0,1\}.
\]
The moving tangent \(\mathcal O(2e)\) factor, omitting its fixed section, is
\[
\mathcal T(w,e)=\frac{(-1)^e e^{2e}}{(e!)^2w^{2e}}.
\]
These specialize [LS, Lemma 4.5]. Directly dualizing the section weights of each defining line bundle gives the additional factor
\[
\mathcal Q_d(I,J,e)=\prod_{k=0}^{de}(t-dH_I+kw/e).
\]
Thus the edge contribution before deck symmetry is
\[
\mathcal T(w,e)
\prod_{\ell\in\{\lambda_k-\lambda_a:k\in K\}\cup\{\lambda_b-\lambda_c\}}
\prod_{k=0}^{e}(\ell-kw/e)^{-1}
\prod_{k\in K}(\lambda_k-\lambda_c)^{-1}
\prod_{j=1}^c\mathcal Q_{d_j}(I,J,e).
\]
The separate deck factor is \(1/e\). Reversing orientation sends \((w,\ell,dH_I)\) to \((-w,\ell-a_Nw,dH_I-dw)\); reversing product indices proves equality. In particular the degree-zero normal factors occur only once, irrespective of cover degree.

For a stable contracted vertex at \(I\) of genus \(h\) and valence \(v\), let
\[
T_I=\prod_{i\in I,j\notin I}(\lambda_j-\lambda_i),\quad
q_{I,j}=t-d_jH_I,\quad
\Lambda_h^\vee(z)=\sum_{k=0}^h(-1)^k\lambda_kz^{h-k}.
\]
Here the \(\lambda_k\) inside \(\Lambda\) are Hodge classes, not the auxiliary torus characters; implementations should use distinct variable names. Normalization of \(R\Gamma(f^*TX)\) and \(R\Gamma(f^*E)\), with [LS, Lemma 4.4] and Mumford's relation, gives
\[
\boxed{V_{I,h,v}=T_I^{v-1}
\prod_{i\in I,j\notin I}\Lambda_h^\vee(\lambda_j-\lambda_i)
\prod_{j=1}^c q_{I,j}^{1-v-2h}\Lambda_h^\vee(q_{I,j})
\prod_F(u_F-\psi_F)^{-1}.}
\]
Indeed the dual twist before applying Mumford is \(q^{1-v}/\Lambda_h^+(q)\), where \(\Lambda_h^+(q)=\sum_k\lambda_kq^{h-k}\), and \(\Lambda_h^+(q)\Lambda_h^\vee(q)=q^{2h}\). This explains both the dual polynomial and the exponent of \(q\).

With these endpoint-inclusive edge conventions, an unmarked genus-zero univalent vertex contributes \(u_F\), and a bivalent vertex contributes
\[
\frac{T_I}{(u_{F_1}+u_{F_2})\prod_jq_{I,j}}.
\]
Multiply a graph by \(|\operatorname{Aut}\Gamma|^{-1}\prod_Ee_E^{-1}\), then apply the actual stabilization/gluing pushforward. These auxiliary fixed-map graphs are distinct from the CJR zero/infinity graphs.

### Explicit numerical check and bounded enumeration

For \(n=4\), take weights \((0,1,4,10)\), \(I=\{1,2\}\), \(J=\{2,3\}\), so \(a=1,c=2,b=3\). Then \(w=4\), source tangent weights are \(4,10,3,9\), the two degree-one normal source weights are \(10,3\), and the degree-zero weight is \(9\). At \(J\), the corresponding weights are \(-4,6,-1,9\). Also \(H_I=-1\), \(H_J=-5\). A degree-two cover has source tangent weight 2; its degree-one normal sections have weights \((10,8,6)\) and \((3,1,-1)\). This gives a direct sign and interpolation check.

The fixed-point/invariant-line counts are \((6,12)\) for \(G(2,4)\) and \((10,30)\) for \(G(2,5)\). Every invariant line has Plücker degree one, hence
\[
D=\sum_E e_E,\qquad g=\sum_vg(v)+b_1(\Gamma).
\]
For \(D\le2\), retain degree-zero graphs, one edge of degree one or two, and two degree-one edges. The last case includes a two-edge chain with repeated endpoint labels when allowed, and two parallel edges forming a graph cycle; in the cyclic case the vertex genera sum to \(g-1\). A tree-only enumeration would miss legitimate fixed loci. All unstable conventions still apply.

Two usable specializations are
\[
(0,1,4,10,22)|_{1,\ldots,n},\qquad
(0,2,7,17,40)|_{1,\ldots,n}.
\]
Exact arithmetic checked nonvanishing of tangent and normal-cohomology denominators for \(e=1,2\), and bivalent flag sums for the permitted two-edge degree-one graphs. Use these as two auxiliary torus rays, keeping \(t\) formal. Sum before taking \(\lambda_i=\epsilon a_i\), \(\epsilon\to0\). Agreement of two rays is a diagnostic, not proof of an unspecialized identity. The full grid has \(14\cdot3\cdot2=84\) target/degree/genus cases, or 168 evaluations with two rays; no such scan was run in this audit.

## 5. Constants and the genus-two/three relation ranges

From [FP, equation (1)], the degree-zero virtual class is
\[
c_{mg}(\mathbb E^\vee\boxtimes TY)
\cap[\overline{\mathcal M}_g\times Y].
\]
Expanding its target-degree-\(m\) component and applying Mumford's relation gives the following deductions:
\[
\begin{aligned}
m=2:\quad&G_{2,0}=(\int_Yc_1^2)\lambda_2,
&G_{3,0}&=(\int_Yc_1^2)\lambda_1\lambda_3;\\
m=3:\quad&G_{2,0}=C_Y\lambda_1\lambda_2,
&G_{3,0}&=-C_Y\lambda_1\lambda_2\lambda_3.
\end{aligned}
\]
For example, the genus-two surface expression before simplification is \(\lambda_2(c_1^2-2c_2)+\lambda_1^2c_2\), and \(\lambda_1^2=2\lambda_2\). Thus §2 fixes every compact constant, including vanishing for the three K3 families. For threefolds, the scalar integrals are
\[
N_{2,0}=C_Y/5760,\qquad N_{3,0}=-C_Y/1451520
\]
by [FP, Theorem 4 and §4.3]. Apply the additional sign \(\sigma_{g,0}\) only when forming the CJR comparison.

For \(D>0\), ordinary K3 virtual classes vanish by [MPT, §0.1]. Combining this with the constants proves \(G_{g,D}=0\) for all three K3 families in the requested grid. Reduced K3 invariants are different invariants; CJR's reduced log-GLSM comparison here is with **ordinary** target GW theory.

The following table gives \((\alpha_{D=0},\alpha_{D=1},\alpha_{D=2})\). For each entry, negative-\(t\) cancellation concerns integer codimensions \(\max(0,\alpha+1)\le r\le3g-3\).

| Targets | \((m,\kappa)\) | Genus 2: \(\alpha\) | Genus 3: \(\alpha\) |
|---|---|---|---|
| G1 | \((3,3)\) | \((3,0,-3)\) | \((6,3,0)\) |
| G2, G9 | \((3,2)\) | \((3,1,-1)\) | \((6,4,2)\) |
| G3, G10 | \((3,1)\) | \((3,2,1)\) | \((6,5,4)\) |
| G4, G11, G12 | \((3,0)\) | \((3,3,3)\) | \((6,6,6)\) |
| G5 | \((2,2)\) | \((2,0,-2)\) | \((4,2,0)\) |
| G6, G13 | \((2,1)\) | \((2,1,0)\) | \((4,3,2)\) |
| G7, G8, G14 | \((2,0)\) | \((2,2,2)\) | \((4,4,4)\) |

Consequences for this bounded calculation:

1. For a K3, \(B_{D,r}=-Z_{D,r}\) also at \(r=\alpha\): the applicable range is \(r=2,3\) for genus two, and \(r=4,5,6\) for genus three.
2. For a CY3, \(\alpha=3g-3\), so no negative power survives the grading. Its \(t^0\) equation has a compact GW term, including the nonzero constants and the degree-one values in §6. It cannot be converted into a vanishing relation by deleting that term.
3. If \(\alpha<0\), the compact stabilization vanishes by grading. This occurs at \((g,D)=(2,2)\) for G1, G2, G5, G9. It is a statement about this unmarked pushforward, not every invariant with markings or insertions.
4. For the ambient Grassmannian,
   \[
   \nu_X=(1-g)(2n-7)+nD.
   \]
   It is negative at \(D=0\) for both ambients and both genera, so the nonequivariant ambient virtual class, and hence \(Z_0(t)\), vanish. It is also negative for \(n=5,g=3,D=1\), giving \(Z_1(t)=0\) in that case. A derived Euler ratio with formal invertible \(t\) does not restore a zero nonequivariant virtual class. Nonzero compact constants or degree-one CY3 values then come from \(B\).
5. Other ordinary Fano positive-degree classes remain unknown here. A computed \(Z\), or a putative relation in its decorated-graph expression, is not permission to replace an unknown \(G\) by zero.

The prescribed bounded computation should report the negative-power range, the independently justified \(t^0\) terms, and any residual unknown compact class separately. To claim a relation, it must also specify the full CJR expression used and its reduction in the tautological ring. Auxiliary localization alone computes \(Z\), not the unknown effective theory at infinity.

## 6. Independent CY3 line and conic checks

### Lines from the flag variety

A Plücker line in \(G(2,n)\) consists of planes \(A\subset V_2\subset B\), with \(\dim A=1\), \(\dim B=3\). Its parameter space is
\[
F=F(1,3;n),\qquad \dim F=3n-7.
\]
Write \(U_1\subset U_3\) for its tautological bundles and \(W=U_3/U_1\). On the universal line \(\mathbb P(W)\), the Plücker bundle is \(U_1^\vee\otimes\mathcal O_{\mathbb P(W)}(1)\). Therefore the equation of degree \(d\) restricts to a section of
\[
V_d=\operatorname{Sym}^d(W^\vee)\otimes(U_1^\vee)^{\otimes d},
\qquad \operatorname{rk}V_d=d+1.
\]
This confirms the proposed bundle including the tensor power on \(U_1^\vee\). For a CY3, \(c=2n-7\) and \(s=n\), so \(\sum_j(d_j+1)=3n-7=\dim F\). Restriction of Plücker polynomials generates each fiber of \(V_d\), since the Plücker embedding maps the line to a projective line and homogeneous polynomials restrict surjectively to it. Generic smoothness of the incidence correspondence gives transverse zeros for general equations. By [Stacks], their count is
\[
\boxed{L=\int_{F(1,3;n)}\prod_jc_{d_j+1}(V_{d_j}).}
\]
In root notation \(x=c_1(U_1^\vee)\), \(y,z\) the roots of \(W^\vee\), the integrand is exactly
\[
\boxed{\prod_j\prod_{k=0}^{d_j}
\big(d_jx+(d_j-k)y+kz\big).}
\]

Here is a finite rational formula that evaluates it without a flag-ring implementation. Choose distinct vector weights \(\mu_1,\ldots,\mu_n\). Sum over \(a\) and unordered \(\{b,c\}\subset\{1,\ldots,n\}\setminus\{a\}\), with \(U=\{a,b,c\}\):
\[
\boxed{
L=\sum_{a,\{b,c\}}
\frac{\displaystyle\prod_j\prod_{k=0}^{d_j}
[-d_j\mu_a-(d_j-k)\mu_b-k\mu_c]}
{\displaystyle(\mu_b-\mu_a)(\mu_c-\mu_a)
\prod_{i\in U}\prod_{r\notin U}(\mu_r-\mu_i)}.}
\]
The denominator is the flag tangent Euler class, from \(\operatorname{Hom}(U_1,U_3/U_1)\) and \(\operatorname{Hom}(U_3,\mathbb C^n/U_3)\). There is no extra factor of two: \(\{b,c\}\) is unordered. Exact rational evaluation with both weight lists in §4 gave
\[
L_{G4}=1280,\qquad L_{G11}=540,\qquad L_{G12}=400.
\]
These are independent intersection calculations, not outputs copied from the mirror tables. The latter two agree with [BCFKvS, §7.1–§7.2]; the first agrees with [M2, i3].

### Published degree-two values and their normalization

Let \(n_d^{(0)}\) denote genus-zero instanton numbers, and \(N_{0,D}\) the connected unmarked ordinary GW invariant. Their convention is
\[
N_{0,D}=\sum_{k\mid D}\frac{n_{D/k}^{(0)}}{k^3},\qquad
N_{0,1}=n_1^{(0)},\quad N_{0,2}=n_2^{(0)}+\frac{n_1^{(0)}}8.
\]
The normalization follows from the Yukawa expansion in [BCFKvS, p. 30], or directly from [M2, i5 versus i7].

| CY3 | \(n_1^{(0)}\), lines | \(n_2^{(0)}\), conic/instanton number | Ordinary \(N_{0,2}\) | Primary locator |
|---|---:|---:|---:|---|
| G4: \(G(2,4)[4]\) | 1280 | 92288 | 92448 | [M2], i3, i5, i7, second entry |
| G11: \(G(2,5)[1,1,3]\) | 540 | 12555 | \(25245/2\) | [BCFKvS], Theorem 6.2.1 and §7.1 |
| G12: \(G(2,5)[1,2,2]\) | 400 | 5540 | 5590 | [BCFKvS], Theorem 6.2.1 and §7.2 |

For G11 and G12 the source states virtual rational-curve numbers; [Kim, Theorem 1 and Remark (4)] supplies their genus-zero mathematical justification. This audit has not independently proved that every degree-two image for a general member is an isolated embedded conic with the required normal bundle. Accordingly the table is a verified genus-zero invariant check, not an assertion of that stronger enumerative geometry.

Genus-zero convexity does hold: each summand has nonnegative degree on every component of a genus-zero tree, so \(H^1(C,f^*E)=0\). The nondual Euler twist therefore computes these ordinary invariants by [CJW, Theorem 1.2.3]. At \(t=0\), the dual twist instead has sign \((-1)^{c+sD}\). Its degree-one/degree-two numerical checks should be
\[
\begin{array}{c|rr}
&D=1&D=2\\\hline
G4&-1280&-92448\\
G11&540&-25245/2\\
G12&400&-5590.
\end{array}
\]

### What the line calculation proves in genus two and three

For a transverse isolated line \(C\subset Y\), the normal sequence has \(H^0(N_{C/Y})=0\): its kernel is the tangent space to the zero scheme on \(F(1,3;n)\). Since \(\deg N_{C/Y}=-2\), splitting on \(\mathbb P^1\) forces \(N_{C/Y}=\mathcal O(-1)\oplus\mathcal O(-1)\). Every degree-one image in the Plücker embedding is a line. A degree-one nonconstant component cannot have positive geometric genus, and the remaining higher-genus components are contracted. These observations exhaust the degree-one images for general members.

Applying the local contribution [FP, Theorem 3], including contracted tails,
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3},
\qquad C(2,1)=1/240,\quad C(3,1)=1/6048,
\]
gives independently justified ordinary scalar checks:

| CY3 | \(N_{2,1}\) | \(N_{3,1}\) |
|---|---:|---:|
| G4 | \(16/3\) | \(40/189\) |
| G11 | \(9/4\) | \(5/56\) |
| G12 | \(5/3\) | \(25/378\) |

These are ordinary numbers before multiplying by \(\sigma_{g,1}\). Scalar integration alone is not a presentation of a class in \(A^{3g-3}(\overline{\mathcal M}_g)\); any use as a class-level comparison must justify the tautological/top-degree identification it employs.

In degree two, the double-cover-of-lines contribution is known: \(N_{2,2}^{\rm lines}=2L/240\), \(N_{3,2}^{\rm lines}=8L/6048\). This is **only a contribution**. The primary genus-zero conic table alone does not prove the full genus-two/three invariants. If one separately proves that all other images are disjoint isolated \((-1,-1)\) conics and excludes connected reducible images, then the total would be
\[
N_{2,2}=\frac{n_2^{(0)}+2n_1^{(0)}}{240},\qquad
N_{3,2}=\frac{n_2^{(0)}+8n_1^{(0)}}{6048}.
\]
That additional geometric assertion is **not verified here**, so these last total formulas are conditional and must not be installed as certified degree-two ordinary inputs. Keep the full positive-degree-two CY3 compact terms symbolic unless a further proof or independent higher-genus source is provided.

## 7. The \(G(2,4)=Q^4\subset\mathbb P^5\) regression

The Plücker coordinates satisfy one nondegenerate quadratic equation
\[
p_{12}p_{34}-p_{13}p_{24}+p_{14}p_{23}=0.
\]
Thus \(G(2,4)\) is the smooth quadric fourfold, with Plücker \(H\) equal to the projective hyperplane class. A target \(G(2,4)[d_1,\ldots,d_c]\) has the projective presentation \(\mathbb P^5[2,d_1,\ldots,d_c]\), allowing linear equations to be eliminated. The explicit equation proves the isomorphism; [BCFKvS, §2.1] uses the same quartic/quadric model.

This supplies three kinds of valid regression:

- Ordinary Chern numbers and ordinary target GW invariants agree for the same variety and polarization.
- In genus zero, convex quantum Lefschetz relates the quadric presentation and the projective presentation, with the nondual Euler twists or the correctly converted dual signs.
- Complete CJR expressions compare to the same ordinary target class, after their respective comparison signs.

It does **not** justify equating
\[
Z^{G(2,4),\,\oplus\mathcal O(d_j)}_D(t)
\quad\text{and}\quad
Z^{\mathbb P^5,\,\mathcal O(2)\oplus\bigoplus\mathcal O(d_j)}_D(t)
\]
in genus two or three. For a contracted genus-\(h>0\) component, the pullback of \(\mathcal O(2)\) is trivial and has \(H^1\) of dimension \(h\). Thus the required convexity fails on the actual compactified moduli space, even for total degree one or two. Replacing the missing higher-genus Lefschetz theorem with a formal Euler ratio would assume the very comparison that needs proof. The infinity targets and aggregates also differ between the two presentations.

There is a sign difference even before addressing those aggregates. Adding the quadric equation changes
\[
(c,s)\longmapsto(c+1,s+2),\qquad
\frac{\sigma^{\mathbb P^5}_{g,D}}{\sigma^{G(2,4)}_{g,D}}
=(-1)^{1-g+2D}=(-1)^{1-g}.
\]
Hence the sign reverses for genus two and agrees for genus three. A regression of full CJR \(t^0\) comparisons must account for this; it should compare the recovered ordinary target classes, not raw signed sums or individual zero-vertex contributions.

## 8. Established inputs and remaining work

This audit establishes the fourteen geometric types and Chern-number rows, the Grassmannian fixed-point/edge data including equivariant normal splitting, the private-draft CJR applicability and signs, all ordinary constant-map comparisons, ordinary K3 vanishing, three independently derived CY3 line counts, verified genus-zero degree-two values with multiple-cover normalization, and the ordinary genus-two/three degree-one CY3 scalar values. Exact arithmetic was run in the persistent Linux guest without modifying implementation files.

The bounded relation search still needs actual auxiliary graph enumeration, stabilization, the two-ray comparison, and any claimed tautological reduction. Full ordinary degree-two CY3 higher-genus terms and other positive-degree Fano terms remain uncomputed here. No new tautological relation, completed scan, or higher-genus quantum Lefschetz identity is certified by this source note.
