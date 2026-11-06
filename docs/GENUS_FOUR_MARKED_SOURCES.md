# Genus four with one or two markings: source and formula audit

Research date: 2026-11-06. This audit concerns connected ordinary stable maps, rational Chow classes where specified, genus \(g=4\), \(n=1,2\), and the existing ambient curve bound two. Each labeled marking carries precisely one of
\[
U=\tau_0(1),\qquad P=\tau_1(1)=\psi^{\mathrm{map}},\qquad
H=\tau_0(H_X).
\]
The letter \(P\) in an insertion profile means the true stable-map cotangent class, not an ancestor and not the log-domain cotangent class.

The grading, full-twist marking identities, and compact-type marking extension are supported. The principal qualification is that two different stabilizations occur: from log maps to stable maps, and from stable maps to stable curves. Neither allows unconditional replacement of the corresponding cotangent lines. The numerical identities below do not determine the new marked Chow classes or certify their Pixton membership.

This note audits formulas and existing geometric arguments. It does not run a genus-four localization calculation, independently rerun earlier curve counts, or certify an implementation.

## 1. Scope and provenance of the 43 presentations

The requested sample is \(1+1+14+12+14+1=43\) ambient presentations. Isomorphic compact targets in different ambients remain separate presentations because their ambient twists and CJR infinity theories differ.

| Group | Presentations retained | Curve bound and existing record |
|---|---|---|
| Elliptic cubic | Smooth cubic in \(\mathbb P^2\), defining bundle \(\mathcal O(3)\) | Ambient \(D=0,1,2\); [existing descendant audit](CORRECTNESS_AUDIT.md) |
| Quartic K3 | Smooth quartic in \(\mathbb P^3\), defining bundle \(\mathcal O(4)\) | \(D=0,1,2\); [K3 source audit](QUARTIC_K3_SOURCES.md) |
| Projective complete intersections | Nine nonlinear Fano hypersurfaces: dimensions \(m=2,3,4\), degrees respectively \(2\!:\!3,2\!:\!4,2\!:\!5\); five CY3 types \((5),(2,4),(3,3),(2,2,3),(2,2,2,2)\) | \(D=0,1,2\); [projective source audit](PROJECTIVE_COMPLETE_INTERSECTION_SOURCES.md) |
| Toric complete intersections | T1–T12, retaining the products, \(\mathbb F_1\) factors, and rank-two presentations | Nonnegative Mori tuples with coefficient sum at most two; [toric source audit](TORIC_COMPLETE_INTERSECTION_SOURCES.md) and [later geometric additions](TORIC_COMPLETE_INTERSECTION_RESULTS_2026-11-06.md) |
| Grassmannian complete intersections | G1–G14 in \(\operatorname{Gr}(2,4)\) and \(\operatorname{Gr}(2,5)\) | Plücker \(D=0,1,2\); [Grassmannian source audit](GRASSMANNIAN_COMPLETE_INTERSECTION_SOURCES.md) |
| Additional Grassmannian CY3 | \(\operatorname{Gr}(2,7)[1^7]\) | Plücker \(D=0,1,2\); [degree-two geometric audit](GRASSMANNIAN_27_SOURCES.md) |

Thus all \(3^n\) labeled profiles are included. Symmetry permits storing six unordered two-mark profiles, but recovering the labeled Chow class requires the corresponding marking permutation. There is no \(1/n!\) in the fixed-profile formulas.

On a projective space or Grassmannian, \(H_X\) is the standard hyperplane or Plücker class. A toric product has no unique divisor called “the hyperplane”: record the chosen polarization and calculate
\[
d_H=\int_\beta H_X.
\]
For \(H_X=aA+bB\), \(d_H=ad_1+bd_2\). For
\(H_X=aS+bF+cH_{\mathbb P}\) on \(\mathbb F_1\times\mathbb P^q\),
\[
d_H=(b-a)d_S+ad_F+cd_{\mathbb P}.
\]
The usual Segre choice on a product is \(A+B\); an ample scroll/product choice on the latter ambient is \(S+2F+H_{\mathbb P}\). The formulas hold for either explicitly recorded choice, or any other specified divisor. The Mori coefficient bound is not automatically the degree in that polarization. In particular \(S+F\) has tuple \((1,1,0)\), whereas \(2(S+F)\) is outside the bound.

## 2. Primary references and what they establish

**[CJR-III]** Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, supplied private draft, October 17, 2025. [Local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>). Rechecked the existing extraction at \(/tmp/log-glsm-audit/cjr-iii.txt\). Printed locators: §2.5.2–§2.5.3, pp. 14–15 (compact-type and unit sectors); Assumption 4.2; opening of §8 and (8.2) (stabilization and coarse log-domain cotangent lines); (8.21)–(8.23), pp. 82–83 (marked localization factors); Lemma 8.3, Theorem 8.4/(8.28), Remark 8.5, p. 84; §9.3/(9.6)–(9.7), pp. 89–90; §9.4, p. 92 (arbitrary rank); §10.1/Theorem 10.1, p. 93 (ambient descendants); Proposition 10.9, p. 95 (conditional effective-cycle tautologicality). These are private-draft statements, not asserted here to be published theorems.

**[CJR-I]** Same authors, *Structural Formulae in Logarithmic Gauged Linear Sigma Models I: The Tropical Decomposition Formula*, supplied October 17, 2025 draft. [Local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_1_the_Tropical_Decomposition_Formula.pdf>), existing extraction \(/tmp/log-glsm-audit/cjr-i.txt\). Theorem 1.4 and (1.10), pp. 7–8, explicitly discuss ambient descendant numerical comparison. They must not be read as equality of cotangent line classes before stabilization. For arbitrary rank, use the sign explicitly stated in CJR III §9.4 and [CJW]; the displayed hypersurface-style sign in this draft of CJR I is not a replacement for it.

**[CJR-GLSM]** Same authors, *The logarithmic gauged linear sigma model*, Invent. Math. 225 (2021), 1077–1154, [arXiv:1906.04345v2](https://arxiv.org/pdf/1906.04345v2). Corollary 1.10 and Proposition 4.1, preprint pp. 7 and 34–35, give the reduced-log-GLSM/ordinary-complete-intersection comparison.

**[CJW]** Q. Chen, F. Janda, R. Webb, *Virtual cycles of stable (quasi-)maps with fields*, Adv. Math. 385 (2021), 107781, [arXiv:1911.09825v2](https://arxiv.org/pdf/1911.09825v2). Corollary 1.2.1, p. 4, gives the Euler-characteristic sign and sum over target curve classes. Theorem 1.2.3 is convex quantum Lefschetz; it justifies the earlier genus-zero inputs, not positive-genus convexity.

**[CG]** T. Coates, A. Givental, *Quantum Riemann–Roch, Lefschetz and Serre*, Ann. of Math. 165 (2007), 15–53, [author PDF](https://math.berkeley.edu/~giventh/papers/rls.pdf), [arXiv:math/0110142](https://arxiv.org/abs/math/0110142). In the 32-page author PDF, Theorems 1 and \(1'\), pp. 6–7, concern the full twisted descendant potential. Appendix 1, formula (23) and properties (i)–(vi), pp. 23–24, supplies cotangent comparison, virtual-class naturality, forgetting compatibility of \(R\pi_*f^*E\), and normalization at nodes. Appendix 2, pp. 28–30, gives the descendant/ancestor comparison, including the virtual rational-tail divisor. Quantum Riemann–Roch itself does not identify ancestors with descendants.

**[P]** R. Pandharipande, *Rational curves on hypersurfaces [after A. Givental]*, Séminaire Bourbaki, Astérisque 252 (1998), exposé 848, [author's exposition at Numdam](https://numdam.org/item/SB_1997-1998__40__307_0.pdf). §1.2 states string, dilaton, and divisor equations, and their higher-genus form. The application to the present twist is the deduction in §5 below using [CG]'s forgetting compatibility.

**[FP]** C. Faber, R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. 139 (2000), 173–199, [arXiv:math/9810173](https://arxiv.org/pdf/math/9810173). Introduction (1) identifies the constant-map obstruction bundle. Equation (7) and Theorem 3 define and evaluate the full local \((-1,-1)\)-curve contribution. Theorem 4 and its following paragraph give the CY3 constant-map number; §4.3 relates \(\lambda_{g-1}^3\) to \(2\lambda_g\lambda_{g-1}\lambda_{g-2}\). All genus-four substitutions below are deductions, not newly sourced enumerative calculations.

**[M]** D. Mumford, *Towards an Enumerative Geometry of the Moduli Space of Curves*, Arithmetic and Geometry II (1983), 271–328, [author-hosted scan](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf). §5, (5.2), Corollary (5.3), (5.4), printed pp. 305–307: algebraic GRR and \(c(\mathbb E)c(\mathbb E^\vee)=1\) in rational Chow.

**[MPT]** D. Maulik, R. Pandharipande, R. P. Thomas, *Curves on K3 surfaces and modular forms*, [arXiv:1001.2719](https://arxiv.org/pdf/1001.2719). §0.1, p. 2: the ordinary virtual class vanishes for nonzero K3 curve classes; the reduced theory is different.

**[CJR-P]** Q. Chen, F. Janda, Y. Ruan, *Punctured logarithmic R-maps*, [arXiv:2208.04519v2](https://arxiv.org/pdf/2208.04519v2). §8.1.1 and Theorem 8.1, p. 78, define \(\psi_{\min}\) and punctured equations; Theorem 8.4/(8.8), p. 80, includes stabilized boundary corrections. Question 9.36 asks about effective-cycle tautologicality. Ordinary GW marking equations cannot simply replace these punctured equations.

**[J]** F. Janda, *Relations on \(\overline{\mathcal M}_{g,n}\) via equivariant Gromov–Witten theory of \(\mathbb P^1\)*, Algebraic Geometry 4 (2017), 311–336, [journal PDF](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf). Theorem 1.1 proves the 3-spin/Pixton relations in Chow. The discussion on pp. 312–313 separates this from completeness of their tautological closure. The original [PPZ paper](https://arxiv.org/abs/1303.1043) proves its relations in cohomology.

**[GP]** T. Graber, R. Pandharipande, *Localization of virtual classes*, [arXiv:alg-geom/9708001](https://arxiv.org/pdf/alg-geom/9708001), equation (1) and §4: virtual localization and stable-map graph factors. The earlier geometric counts also use the [regular-section/top-Chern formula](https://stacks.math.columbia.edu/tag/0FA8). Those sources provide the methods; the target-specific incidence arguments are recorded in the linked repository audits.

## 3. Marked grading and the comparison sign

Write \(Y=Z(s)\subset X\), \(E\) of rank \(c\), \(m=\dim X-c\), and use the ambient class \(\beta\). In expressions involving \(c_1(TY)\cdot\beta\), this means the ambient pairing
\[
(c_1(TX)-c_1(E))\cdot\beta.
\]
For each effective class of \(Y\) pushing forward to \(\beta\), it has the same value. Preserve the sum over such target classes, especially for surfaces.

Let
\[
s_n:\mathcal M_{n,\beta}:=\overline{\mathcal M}_{4,n}(X,\beta)
 \longrightarrow\overline{\mathcal M}_{4,n}
\]
forget the map and stabilize. Set
\[
\Theta_{n,\beta}(t)=e_t((R\pi_*f^*E)^\vee),\qquad
\mathcal Z_{\beta,\mathbf I}(t)
=s_{n,*}\left(\Theta_{n,\beta}(t)\prod_i I_i
 \cap[\mathcal M_{n,\beta}]^{\rm vir}\right).
\]
This is the full derived dual Euler ratio, with one common CJR weight \(t\). It retains \(R^1\), and is expanded with \(t\) invertible. It is not an \(R^0\)-only twist, a product of separate infinity theories, or the compact GW theory.

Let \(a=\sum_i\deg I_i\), using complex codimension:
\(\deg U=0,\ \deg P=\deg H=1\). To prevent counting \(a\) twice, distinguish
\[
\alpha_0=m(g-1)-c_1(TY)\cdot\beta,\qquad
\alpha_{\mathbf I}=\alpha_0+a.
\]
The first is the uninserted shift; the second is the total marked equivariant codimension. Both have no separate \(n\) term. Indeed
\[
\begin{aligned}
\operatorname{vdim}\overline{\mathcal M}_{g,n}(X,\beta)
 &=(1-g)(\dim X-3)+c_1(TX)\cdot\beta+n,\\
\operatorname{rk}(R\pi_*f^*E)&=c(1-g)+c_1(E)\cdot\beta,
\end{aligned}
\]
so the twisted, inserted dimension is
\[
(1-g)(m-3)+c_1(TY)\cdot\beta+n-a.
\]
Subtracting it from \(3g-3+n\) gives \(\alpha_0+a\).

At genus four, therefore,
\[
\boxed{\mathcal Z_{\beta,\mathbf I}(t)
 =\sum_{r=0}^{9+n}t^{\alpha_0+a-r}
       \mathcal Z_{\beta,\mathbf I;r},\qquad
 \alpha_0=3m-c_1(TY)\cdot\beta.}
\]
The coefficient which can be integrated on \(\overline{\mathcal M}_{4,n}\) has power
\[
\boxed{t^{\alpha_0+a-(9+n)}.}
\]
If “alpha” already includes insertion degrees, this is \(t^{\alpha_{\mathbf I}-(9+n)}\), not \(t^{\alpha_{\mathbf I}+a-(9+n)}\).

The compact comparison sign remains
\[
\boxed{\sigma_{4,\beta}=(-1)^{c(1-4)+c_1(E)\cdot\beta}
=(-1)^{-3c+c_1(E)\cdot\beta}.}
\]
There is no marking or descendant contribution to this sign. Its origin is the Euler characteristic in [CJW], and the arbitrary-rank version is explicit in [CJR-III, §9.4]. The public [CJR-GLSM, Corollary 1.10] includes orbifold age terms; these vanish for every scheme target and unit sector in the present sample.

For orientation, \(\alpha_0=3\) for the elliptic cubic, \(6\) for every K3, and \(9\) for every CY3. For a rank-one-polarized Fano of index \(I\), it is \(3m-ID\). For toric Fanos use the actual divisor pairing, not the coefficient sum of the Mori tuple.

## 4. Compact-type markings, global insertion, and the CJR aggregate

The marking assumption is no longer vacuous, but it is satisfied. Choose the discrete data \(1^n\) of [CJR-III, (9.6)]: zero contact, unit zero sector, evaluation in \(0_P\simeq X\). This is precisely compact type by §2.5.2–§2.5.3 and satisfies Assumption 4.2. “Unit sector” is discrete geometric data; it does not require the evaluation insertion to be the cohomological unit. Ambient \(H_X\), as well as \(1\), is permitted. Cotangent powers do not change the sector.

Section 9.4 explicitly reuses (9.6) for arbitrary rank. Thus the marking extension covers all 43 presentations, including the toric nef bundles. No ampleness-based finiteness corollary is required for this bounded audit. For \(c>1\), infinity still has target \(\mathbb P_X(E^\vee)\) and its own fiber classes.

To define the requested true stable-map insertion unambiguously, let
\[
\rho:\mathcal U_{4,1^n}(P,\beta)\longrightarrow
 \mathcal M_{n,\beta}
\]
be CJR's stabilization (their \(st\)), and insert
\(\rho^*(\prod I_i)\) globally before applying localization. Equivalently, cap the pushed-forward CJR identity on \(\mathcal M_{n,\beta}\) with the product of stable-map cotangent and evaluation classes, then apply \(s_{n,*}\). These are operational classes, so the projection formula applies even though the stable-map stack need not be smooth.

This construction gives
\[
\mathcal L_{\beta,\mathbf I}(t)
 =\mathcal Z_{\beta,\mathbf I}(t)+\mathcal B_{\beta,\mathbf I}(t).
\]
Here \(\mathcal B\) is the full sum of all infinity-containing localization terms with the restrictions of this global insertion. In particular, it includes stabilization corrections on contracted vertical tails. On the all-zero locus the insertion is the desired stable-map insertion, giving exactly \(\mathcal Z\) as defined above.

The global left side has no negative \(t\)-powers. Define
\[
\mathcal G_{\beta,\mathbf I}
=\sum_{i_*\beta'=\beta}s_{Y,*}\left(
 \prod_i I_i^Y\cap[\overline{\mathcal M}_{4,n}(Y,\beta')]^{\rm vir}\right).
\]
The stable-map cotangent lines agree under the closed inclusion of maps to \(Y\) into maps to \(X\): neither inclusion contracts any component. Consequently [CJR-III, (9.7), §9.4] and the projection formula yield
\[
\boxed{\mathcal B_{\beta,\mathbf I;r}
=-\mathcal Z_{\beta,\mathbf I;r}\quad(r>\alpha_0+a),}
\]
and at the constant power,
\[
\boxed{\mathcal B_{\beta,\mathbf I;\alpha_0+a}
=\sigma_{4,\beta}\mathcal G_{\beta,\mathbf I}
 -\mathcal Z_{\beta,\mathbf I;\alpha_0+a}.}
\]
Only codimensions \(0\le r\le9+n\) occur. Negative or out-of-range compact codimension means the corresponding stabilized class vanishes by grading. Positive \(t\)-coefficients are equivariant corrections and cannot be deleted using the compact comparison.

For a CY3 profile with \(a=n\), the threshold is \(9+n\), the dimension of the receiving moduli space, so there is no negative-\(t\) range. With units, \(a<n\) and a negative range opens, although all its top numerical values in these profiles vanish by §5. For K3 the negative range is \(r>6+a\), and ordinary cycle vanishing also makes the aggregate identity hold at \(r=6+a\).

This is an equality for the aggregate. It does not give individual punctured effective cycles or assert their support on the stable-curve boundary. Contracting tails can push infinity contributions onto the interior.

## 5. Numerical marking identities for the entire ambient twist

Let
\[
Z_\beta(t)=\int_{[\overline{\mathcal M}_{4,0}(X,\beta)]^{\rm vir}}
 \Theta_{0,\beta}(t).
\]
This is the integrated unmarked ambient twist, including its \(t\)-power, after removal of auxiliary target-torus parameters. It is not \(N^Y_{4,\beta}\), and is not just a selected lower-codimension coefficient of \(\mathcal Z\).

For the map \(q\) forgetting one marking,
\[
q^*(R\pi_*f^*E)=R\pi'_*f'^*E
\]
in \(K\)-theory [CG, Appendix 1, property (iv)]. Dualization and the equivariant Euler ratio preserve this identity. Virtual-class compatibility then makes the ordinary string, dilaton, and divisor arguments apply coefficientwise to this entire twist. Write brackets for these integrated twisted invariants. The three rules are
\[
\begin{aligned}
\langle U,\prod_j\tau_{k_j}(\gamma_j)\rangle
 &=\sum_{j:k_j>0}
 \langle\tau_{k_j-1}(\gamma_j)\prod_{\ell\ne j}\tau_{k_\ell}(\gamma_\ell)\rangle,\\
\langle P,\prod_{j=1}^{n}\tau_{k_j}(\gamma_j)\rangle
 &=(2g-2+n)\langle\prod_{j=1}^{n}\tau_{k_j}(\gamma_j)\rangle,\\
\langle H,\prod_j\tau_{k_j}(\gamma_j)\rangle
 &=d_H\langle\prod_j\tau_{k_j}(\gamma_j)\rangle\\
 &\quad+\sum_{j:k_j>0}
 \langle\tau_{k_j-1}(H_X\gamma_j)
       \prod_{\ell\ne j}\tau_{k_\ell}(\gamma_\ell)\rangle.
\end{aligned}
\]
These displayed formulas are the application of [P] and [CG] to the present twist.

At \(g=4\), the resulting factors are:

| Marked profile | Integrated ambient twist | Insertion degree \(a\) |
|---|---:|---:|
| \(U\) | \(0\) | 0 |
| \(P\) | \(6Z_\beta(t)\) | 1 |
| \(H\) | \(d_HZ_\beta(t)\) | 1 |
| \(U,U\) | \(0\) | 0 |
| \(U,P\) | \(0\) | 1 |
| \(U,H\) | \(0\) | 1 |
| \(P,P\) | \(42Z_\beta(t)\) | 2 |
| \(P,H\) | \(7d_HZ_\beta(t)\) | 2 |
| \(H,H\) | \(d_H^2Z_\beta(t)\) | 2 |

For clarity, \(\langle U,P\rangle=\langle U\rangle=0\); it is not \(Z_\beta\).
The two dilatons give \(7\cdot6=42\), not \(6^2\).
The divisor rule gives
\(\langle H,P\rangle=d_H(6Z_\beta)+\langle H\rangle
=7d_HZ_\beta\). The extra term is essential.
For distinct divisors at the two markings the last entry is \(d_{H_1}d_{H_2}Z_\beta\).

These are identities of integrated \(t\)-expressions. If
\[
Z_\beta(t)=z_\beta\,t^{\alpha_0-9},
\]
each potentially nonzero profile has \(a=n\), hence exactly the same power
\(t^{\alpha_0-9}\). Profiles with a unit have \(a<n\), and their candidate top power is shifted downward by \(n-a\), but its coefficient is zero. The statement “only \(a=n\) can be nonzero in this profile list” follows from the equations, not from homogeneity alone.

There are no global unstable exceptions in this task: \(g=4\) remains stable after forgetting all markings, including at \(\beta=0\). At degree zero \(d_H=0\), so all numerical entries involving \(H\) vanish; \(P\) and \(P,P\) retain factors 6 and 42. The genus-zero classical three-point exceptions and the unstable unmarked genus-one degree-zero case are outside this range.

Unstable *localization vertices* still occur. Their marked-edge cotangent weights and node conventions must be retained; global stability is not permission to omit them. Neither the table nor string/dilaton/divisor applies to individual CJR infinity vertices, a truncated selection of ambient graphs, or the \(\psi_{\min}\) theory.

The same marking factors hold for integrated ordinary compact GW theory with the twist replaced by 1. For a CY3, replace \(Z_\beta(t)\) by the ordinary number \(N^Y_{4,\beta}\). For other targets, dimension may force that number to vanish while leaving a nonzero stabilized Chow class.

## 6. The three cotangent conventions and their corrections

There are three different classes:
\[
\psi_i^{\log},\qquad \rho^*\psi_i^{\rm map},\qquad
\rho^*s_n^*\psi_i^{\overline{\mathcal M}_{4,n}}.
\]
There is also the distinct punctured class
\(\psi_{\min}=c_1(\mathcal O(-e_{\min}))\).
None of these names is interchangeable.

### Stable-map descendants versus stable-curve ancestors

Set \(\bar\psi_i=s_n^*\psi_i^{\overline{\mathcal M}_{4,n}}\). The cotangent-line comparison is
\[
\psi_i^{\rm map}=\bar\psi_i+T_i,
\]
where \(T_i\) is the contraction divisor: a rational tail containing mark \(i\), no other retained marking, and positive map degree is contracted when the map is forgotten. One can interpret the formula by capping with the virtual class and summing the corresponding gluing loci, as in [CG, Appendix 2]. The stable-map factor on such a tail is \(\overline{\mathcal M}_{0,2}(X,\beta_{\rm tail})\) with \(\beta_{\rm tail}>0\); it is not an unstable degree-zero two-point space.

Thus a one-descendant pushed class has the form
\[
\mathcal Z_{\beta,P_i,\ldots}
=\psi_i^{\overline{\mathcal M}_{4,n}}
 \mathcal Z_{\beta,U_i,\ldots}
+s_{n,*}\left(T_i\,\Theta_{n,\beta}
 \prod_{j\ne i}I_j\cap[\mathcal M_{n,\beta}]^{\rm vir}\right).
\]
For two descendants, expand the product, including the intersection term \(T_1T_2\). If the contractions are evaluated by gluing twisted invariants, use the normalization correction at the node:
\[
R\Gamma(C,f^*E)=R\Gamma(C_1,f^*E)+R\Gamma(C_2,f^*E)-E|_{\rm node}.
\]
The dual Euler twists therefore glue with the inverse node factor
\(e_t(E^\vee|_{\rm node})^{-1}\), or equivalently the twisted pairing. Ordinary untwisted diagonal contractions would give the wrong answer.

At \(\beta=0\), the domain of a stable map is already a stable curve, so \(T_i=0\). Forgetting a *marking* is a separate operation: its familiar comparison divisor is the degree-zero bubble carrying the forgotten and retained markings. That correction supplies the string/divisor identities; it is not the positive-degree \(T_i\) above.

### Log-domain descendants versus stable-map descendants

CJR III explicitly defines the \(\psi_i\) at the start of §8 on its coarse log-domain universal curve. Under \(\rho\), vertical components may contract even when their maps to the log target are nonconstant. The cotangent comparison therefore has its own contraction correction. Keeping \(\rho^*\psi_i^{\rm map}\) as a global operational insertion, as in §4, avoids silently identifying the two classes.

For example, in rank one the unstable marked-zero factor (8.21), specialized to \(\mathcal O(3)\), uses
\[
(t-3H)\left(\frac{3H-t}{c(E)}\right)^b
\]
for a log-domain cotangent power \(b\). It is not the restriction of a stable-curve ancestor. If that marked vertical tail is contracted in forming the stable \(X\)-map, the surviving marking is on the adjacent component. Its stable-map cotangent insertion must be restricted there, with further contractions accounted for. In repository terminology this can require a contact cotangent insertion at infinity, not merely \(\psi_{\min}\).

The compact numerical comparison in [CJR-I, Theorem 1.4/(1.10)] and the algorithmic statement [CJR-III, Theorem 10.1] support ambient descendants. They do not prove a termwise equality of raw log-domain and stable-map cotangent factors. Reusing a compiled “log” formula for “stabilized” data requires the appropriate comparison, as already noted in [the correctness audit](CORRECTNESS_AUDIT.md).

There is no single GRR theorem that removes these corrections. The relevant division of labor is:

- [CG, Theorems 1 and \(1'\)]: compute the full twisted descendant theory from ordinary theory.
- [CG, Appendix 1]: forgetting and node compatibility of the twisting complex.
- [CG, Appendix 2]: stable-map descendant/ancestor conversion.
- [CJR-III, Theorem 8.4 and §9.4] plus the projection formula: the globally inserted, stabilized CJR identity.
- [CJR-P, §8]: punctured identities with their own descendants and boundary terms.

## 7. Genus-four Hodge algebra and constant maps

Let \(\lambda_i=c_i(\mathbb E)\), with \(\mathbb E\) of rank four. Mumford's relation gives
\[
\boxed{\lambda_1^2=2\lambda_2,\quad
\lambda_2^2=2\lambda_1\lambda_3-2\lambda_4,\quad
\lambda_3^2=2\lambda_2\lambda_4,\quad
\lambda_4^2=0.}
\]
In particular, reusing the genus-three identity
\(\lambda_2^2=2\lambda_1\lambda_3\) would drop a necessary term.

The formal symmetric Hodge algebra
\[
A_{\rm Hodge}=\mathbb Q[\lambda_1,\lambda_2,\lambda_3,\lambda_4]/
(\lambda_1^2-2\lambda_2,\lambda_2^2-2\lambda_1\lambda_3+2\lambda_4,
 \lambda_3^2-2\lambda_2\lambda_4,\lambda_4^2)
\]
has dimension 16. Here is an algebraic justification, independent of a numerical reduction: its equations say that
\(C(z)C(-z)=1\), so their only common zero over an algebraic closure is the origin. Four homogeneous equations in four variables thus form a zero-dimensional complete intersection. With \(\deg\lambda_i=i\), its Hilbert series is
\[
\frac{(1-q^2)(1-q^4)(1-q^6)(1-q^8)}
 {(1-q)(1-q^2)(1-q^3)(1-q^4)}
=\prod_{i=1}^4(1+q^i),
\]
whose value at \(q=1\) is 16. This describes the symmetric quotient, not an unrestricted algebra in four individual Chern roots, and not the whole tautological ring. On \(\overline{\mathcal M}_4\) additionally impose dimension nine. Hodge classes on \(\overline{\mathcal M}_{4,n}\) pull back from \(\overline{\mathcal M}_4\); adding markings does not resurrect its degree-ten Hodge monomial.

For any smooth compact target \(Y\), [FP, (1)] gives
\[
[\overline{\mathcal M}_{4,n}(Y,0)]^{\rm vir}
=c_{4m}(\mathbb E^\vee\boxtimes TY)
 \cap[\overline{\mathcal M}_{4,n}\times Y].
\]
This supplies all degree-zero marked classes without pretending that their numerical vanishing is cycle vanishing.

For the ambient twist, write \(x_i\) for Hodge roots, \(y_k\) for roots of \(TX\), and \(\ell_j\) for roots of \(E\). At degree zero a reproducible symbolic formula is
\[
\mathcal Z_{0,\mathbf I}(t)=
\left(\prod_{i:I_i=P}\psi_i\right)
\int_X H_X^{\,\#\{i:I_i=H\}}
\prod_{k=1}^{\dim X}\prod_{i=1}^{4}(y_k-x_i)
\prod_{j=1}^{c}\frac{t-\ell_j}{\prod_{i=1}^{4}(t-\ell_j+x_i)}.
\]
This uses \((R\pi_*f^*E)^\vee=(\mathcal O-\mathbb E)\boxtimes E^\vee\).
Expand in ordinary Chow degree, use the Hodge relations, and integrate the target factors. Existing Chern-number tables can be used as input; this note performs no such full expansion.

If \(\dim X\ge4\), the ambient virtual dimension at degree zero is
\(-3(\dim X-3)+n<0\) for \(n\le2\). Thus the nonequivariant ambient virtual cycle, and every requested degree-zero ambient twisted class, is zero. The auxiliary target-torus limit must precede this conclusion. Compact CY3 constants may nevertheless be nonzero and then come entirely from the CJR infinity aggregate.

Useful independent unmarked compact consequences of the same obstruction polynomial are
\[
\begin{array}{ll}
\dim Y=2:&
\mathcal G_{4,0}=\left(\int_Yc_1(TY)^2\right)\lambda_2\lambda_4,\\[2pt]
\dim Y=3:&
\mathcal G_{4,0}=\left(\int_Y(c_3(TY)-c_1(TY)c_2(TY))\right)
\lambda_2\lambda_3\lambda_4.
\end{array}
\]
These are deductions from the root formula and Mumford's relations. For example, the surface target-degree-two term is
\(\lambda_2\lambda_4(c_1^2-2c_2)+\lambda_3^2c_2
=\lambda_2\lambda_4c_1^2\).
For fourfolds, the unmarked and the present one-/two-mark degree-zero virtual cycles vanish by negative dimension. Marked evaluation classes should be handled by the full obstruction formula, not by multiplying an already target-integrated number.

### Ordinary K3 vanishing

For \(\beta'\ne0\), [MPT] gives ordinary virtual-cycle vanishing; adding markings preserves it. At degree zero the separate obstruction calculation is
\[
c_8(\mathbb E^\vee\boxtimes TS)
=\lambda_4^2+c_2(TS)(\lambda_3^2-2\lambda_2\lambda_4)=0.
\]
Both summands vanish by the genus-four Hodge relations. Merely citing \(\lambda_4^2=0\) without checking the \(c_2(TS)\) term is incomplete. This proves cycle-level ordinary vanishing for every K3 presentation and every requested insertion, including degree zero. It does not vanish the ambient twist or individual infinity terms. CJR's reduced log theory compares here to ordinary K3 theory, not reduced K3 GW theory.

### Elliptic cubic: numerical and class distinctions

For a smooth plane cubic \(C\), \(T_C\) is trivial. At degree zero the virtual class on
\(\overline{\mathcal M}_{4,n}\times C\) is \(\lambda_4\) times the fundamental class. Without a target divisor its pushforward to \(\overline{\mathcal M}_{4,n}\) is zero. With one \(H\), \(\int_C H=3\), so, for example,
\[
\mathcal G_{0,(H)}=3\lambda_4,\qquad
\mathcal G_{0,(U,H)}=3\lambda_4,\qquad
\mathcal G_{0,(P,H)}=3\lambda_4\psi_1,
\]
where \(P\) is at marking 1. With two \(H\)'s it is zero because \(H^2|_C=0\).
All the displayed nonzero-codimension classes have zero numerical integral in the requested profile range. Thus a vanishing scalar cannot erase these compact classes.

The ambient degrees \(D=1,2\) have no maps into \(C\): every nonconstant component maps onto \(C\), and its ambient degree is three times its positive covering degree. This establishes compact emptiness in the present bound, not ambient twisted vanishing. These elementary deductions do not assert vanishing of elliptic GW cycles in other degrees. Any unsupported compact class outside the stated vanishing or constant-map arguments remains symbolic.

## 8. Compact CY3 numerical checks in genus four

Keep these ordinary numbers separate from \(Z_\beta(t)\), from the signed compact term \(\sigma_{4,\beta}N^Y_{4,\beta}\), and from a stabilized Chow-cycle presentation.

### Constant maps

The Bernoulli convention is \(B_6=1/42\), \(B_8=-1/30\); the FP formulas use absolute values. Substitution gives
\[
\int_{\overline{\mathcal M}_4}\lambda_3^3=\frac1{43545600},
\qquad
\int_{\overline{\mathcal M}_4}\lambda_4\lambda_3\lambda_2
=\frac1{87091200},
\]
and hence
\[
\boxed{N^Y_{4,0}=\frac{\chi(Y)}{87091200}.}
\]
For a general smooth threefold, the corresponding constant is
\(\int_Y(c_3-c_1c_2)/87091200\).
These follow from [FP, Theorem 4, §4.3] and the preceding obstruction calculation.

### Isolated rational curves

The local integral in [FP, Theorem 3] includes all connected stable domains and contracted higher-genus components:
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3}.
\]
Thus
\[
\boxed{C(4,k)=\frac{k^5}{172800},\qquad
C(4,1)=\frac1{172800},\qquad C(4,2)=\frac1{5400}.}
\]
Whenever the established image-exhaustion argument gives \(L\) isolated \((-1,-1)\) curves in a primitive class and no other possible image, the full ordinary number is \(L/172800\). Applying this local formula is valid in genus four because the incidence, normal-bundle, and image-classification arguments do not depend on domain genus.

The existing CY3 characteristic numbers and degree-one counts give the following exact inputs; the two final columns state the formulas without hiding their denominators.

| Presentation | \(\chi(Y)\) | Primitive line count \(L\) | \(N_{4,0}\) | \(N_{4,1}\) |
|---|---:|---:|---:|---:|
| \(\mathbb P^4[5]\) | \(-200\) | 2875 | \(-200/87091200\) | \(2875/172800\) |
| \(\mathbb P^5[2,4]\) | \(-176\) | 1280 | \(-176/87091200\) | \(1280/172800\) |
| \(\mathbb P^5[3,3]\) | \(-144\) | 1053 | \(-144/87091200\) | \(1053/172800\) |
| \(\mathbb P^6[2,2,3]\) | \(-144\) | 720 | \(-144/87091200\) | \(720/172800\) |
| \(\mathbb P^7[2,2,2,2]\) | \(-128\) | 512 | \(-128/87091200\) | \(512/172800\) |
| G4: \(\operatorname{Gr}(2,4)[4]\) | \(-176\) | 1280 | \(-176/87091200\) | \(1280/172800\) |
| G11: \(\operatorname{Gr}(2,5)[1,1,3]\) | \(-150\) | 540 | \(-150/87091200\) | \(540/172800\) |
| G12: \(\operatorname{Gr}(2,5)[1,2,2]\) | \(-120\) | 400 | \(-120/87091200\) | \(400/172800\) |
| \(\operatorname{Gr}(2,7)[1^7]\) | \(-98\) | 196 | \(-98/87091200\) | \(196/172800\) |

The underlying counts and their transverse line-parameter-space arguments are in the projective and Grassmannian audits linked in §1. Their use here is a genus substitution, not an independent count. Agreement of the compact G4 and projective \([2,4]\) rows does not identify their ambient Euler twists.

### Seven toric primitive classes and three pure doubles

The later toric results document extends the original toric source audit's four primitive counts to all seven. Its additional parameter-space restrictions are:
\(\mathcal O(2)\boxtimes\operatorname{Sym}^4U^\vee\) on
\(\mathbb P^1\times\operatorname{Gr}(2,4)\), giving 640;
\(\mathcal O(1,3)\oplus\mathcal O(2,3)\oplus\mathcal O(3,3)\) on
\(\mathbb P^1\times\mathbb P^2\), giving 54; and the line bundle restriction calculation
\(21\int_{\mathbb F_1}(-K)^2=168\).
Generated restriction bundles and transversality imply \(H^0(N_{C/Y})=0\); CY adjunction then gives \(N_{C/Y}=\mathcal O(-1)^{\oplus2}\). Primitive extremal classes exhaust the possible images. These geometric arguments remain valid at genus four.

| Toric target and ambient class | \(L\) | Full ordinary \(N_{4,\beta}\) |
|---|---:|---:|
| T4: \((2,4)\subset\mathbb P^1\times\mathbb P^3\), \((1,0)\) | 64 | \(64/172800\) |
| T4, \((0,1)\) | 640 | \(640/172800\) |
| T6: bicubic in \(\mathbb P^2\times\mathbb P^2\), \((1,0)\) | 189 | \(189/172800\) |
| T6, \((0,1)\) | 189 | \(189/172800\) |
| T12: anticanonical in \(\mathbb F_1\times\mathbb P^2\), \((1,0,0)\) | 9 | \(9/172800\) |
| T12, \((0,1,0)\) | 54 | \(54/172800\) |
| T12, \((0,0,1)\) | 168 | \(168/172800\) |

Their constant numbers are respectively
\(-168/87091200\), \(-162/87091200\), and \(-144/87091200\).

For T4 class \((2,0)\), and T12 classes \((2,0,0)\), \((0,2,0)\), the zero degrees in the other projections and the section/fiber geometry force the image onto one of the same isolated rational curves. Connectedness prohibits using two disjoint such fibers or sections. Therefore the full double-cover values are
\[
\boxed{
N^{T4}_{4,(2,0)}=\frac{64}{5400},\quad
N^{T12}_{4,(2,0,0)}=\frac9{5400},\quad
N^{T12}_{4,(0,2,0)}=\frac{54}{5400}.}
\]
The multiplier relative to the primitive number is \(2^5=32\), not the genus-three multiplier 8.

This argument does not cover doubles of the projective-line classes, where embedded conics may contribute, or mixed Mori classes. Their genus-four ordinary compact terms remain unknown here.

### The \(\operatorname{Gr}(2,7)\) degree-two case

The prior audit establishes, for a general smooth member, that every degree-at-most-two map is supported on an isolated \((-1,-1)\) line or smooth conic, and excludes distinct meeting lines. This is independent of domain genus. It takes the recorded genus-zero ambient computation
\(N_{0,2}=2499/2\) as input; subtracting \(196/8\) gives 1225 conics.

With that same explicitly retained input dependency,
\[
\boxed{
N_{4,2}=\frac{1225+2^5\cdot196}{172800}
=\frac{7497}{172800}
=\frac{2499}{57600}
=\frac{833}{19200}.}
\]
The proposed intermediate fraction is correct but not in lowest terms.
The primitive value is \(49/43200\), and the constant is \(-7/6220800\).

No ordinary degree-two genus-four value for the other projective or Grassmannian CY3s follows merely from a genus-zero conic table or a primitive line count. A comparable exhaustion/obstruction argument, or another all-genus theorem, is required. This audit keeps those compact terms symbolic.

For every established CY3 number above, the marked ordinary values are obtained by §5: \(6N\), \(d_HN\), \(42N\), \(7d_HN\), \(d_H^2N\), with all unit-containing numerical profiles zero. For \(D=2\) in the last example, \(d_H=2\). Apply the CJR sign only when placing these ordinary numbers into the compact comparison.

## 9. Unit pullbacks and the precise limit of Pixton inference

Let \(p:\overline{\mathcal M}_{4,n}\to\overline{\mathcal M}_4\) forget all markings. The unit axiom for the full twisted GW class gives
\[
\boxed{\mathcal Z_{\beta,(U,\ldots,U)}(t)
=p^*\mathcal Z_{\beta,\varnothing}(t).}
\]
It is compatible with the forgetting identity for the twisting complex in §5. It does not say that the marked class itself is zero. Its numerical integral vanishes because it is pulled back from a space of dimension nine.

Consequently, if an unmarked auxiliary-pole relation already has a certified expression in the chosen Pixton relation ideal, its unit pullback is certified by applying the tautological pullback operation to that expression. Multiplication by stable-curve \(\psi\)-classes or other tautological classes, and the appropriate gluing/forgetting operations, similarly preserve membership in the explicitly defined functorially closed relation system.

This inference has two hypotheses: the original relation must actually be certified, and the relation system used must include the required tautological operations. Equality of unmarked scalar integrals is not such a certificate.

For an ancestor insertion, the corresponding statement follows by the projection formula. For a true \(\psi^{\rm map}\) insertion, §6 adds rational-tail contributions with positive-degree genus-zero factors and normalization terms. For \(H\), the evaluation insertion is not generally a pullback from \(\overline{\mathcal M}_{4,n}\). Neither case is automatically just the pullback of an old relation.

Thus the following are separate assertions:

1. A negative auxiliary torus coefficient vanishes geometrically after summing all loci.
2. Its graph expression is tautological.
3. That expression belongs to the particular generated Pixton ideal.
4. It yields a relation independent of previously generated ones.

Neither the numerical table nor the quantum Riemann–Roch theorem proves all four. [J] proves the Chow validity of the named Pixton generators, not that every new vanishing graph sum is already in an arbitrary implementation's finite relation span.

A sufficient computational certificate remains an exact rational identity
\[
uM=v,
\]
where \(v\) is the candidate relation in documented strata coordinates and the rows of \(M\) are valid generators and their allowed tautological consequences. No such new genus-four membership calculation is supplied by this audit. Two auxiliary weight rays are useful implementation checks, but do not replace it.

At the CJR level, a known numerical compact term determines the numerical aggregate at \(t^0\). It does not determine its Chow class. The K3 cycle vanishing, unit pullbacks, the explicit constant-map polynomials, and the stated elliptic degree restriction are stronger inputs and can be used at class level. Other compact Fano/elliptic/CY classes must be retained until separately justified. Individual effective cycles and their general tautologicality remain outside the conclusion.

## 10. Requirements for a subsequent genus-four calculation

The formulas above support the requested marked extension provided its record distinguishes:

- the 43 ambient presentations, the actual ambient curve class, the chosen divisor \(H_X\), and ordered marking profile;
- the full derived twist, \(\alpha_0\), insertion degree \(a\), top coefficient power, and the rank-dependent compact sign;
- log-domain, stable-map, stable-curve ancestor, contact, and minimal-degeneracy cotangent classes;
- integrated marking identities from class-valued pullbacks and rational-tail corrections;
- proved compact vanishing or numerical values from compact classes still unknown;
- geometric auxiliary-pole cancellation from explicit Pixton membership certificates;
- the total CJR infinity aggregate from its individual effective constituents.

In particular, the factors \(6,d_H,42,7d_H,d_H^2\) and zero are strong numerical consistency checks, but copying an unmarked class and multiplying it by those scalars is not a construction of the requested marked classes.
