# Gr(2,7)[1^7]: degree-two proof and bounded source audit

Research date: 2026-11-06. Work over \(\mathbb C\), with connected, unmarked, ordinary stable-map invariants
\[
N_{g,D}=\int_{[\overline{\mathcal M}_{g,0}(Y,D)]^{\mathrm{vir}}}1,
\qquad Y=\operatorname{Gr}(2,7)[1^7].
\]
Here \(D\) is Plücker degree. The requested computational bound is \(g=2,3\), \(0\leq D\leq2\).

**Verdict.** The proposed geometric proof is sound after making the finite-type incidence argument and the stable-map obstruction comparison explicit. For a general smooth \(Y\), every positive-degree stable map of total degree at most two is supported on an isolated \((-1,-1)\) line or smooth conic. Consequently, using the parent's exact genus-zero calculation \(N_{0,2}=2499/2\), the full ordinary invariants are
\[
\boxed{N_{2,2}=539/80,\qquad N_{3,2}=133/288.}
\]
These are not inferred from higher-genus mirror predictions. The genus-zero degree-two computation is a supplied computational input: this audit verifies its applicable comparison theorem and arithmetic, but does not independently rerun its stable-map graph sum. The line count and characteristic numbers were independently recomputed here.

This is a separate extension of [the fourteen-target source audit](GRASSMANNIAN_COMPLETE_INTERSECTION_SOURCES.md). Its designated output directory is **results/grassmannian-27/2026-11-06**. This source audit writes no results there and makes no changes to the prior fourteen-target results.

## 1. Primary references and their precise roles

- **[FP]** C. Faber and R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. 139 (2000), 173–199; [arXiv:math/9810173](https://arxiv.org/pdf/math/9810173). Introduction equation (7), printed p. 3, defines the local normal-bundle integral; Theorem 3, p. 4, evaluates it for \(g\geq2\). The adjacent genus-zero formula is \(C(0,k)=k^{-3}\). Theorem 4 and its following paragraph give the constant-map formula.
- **[CJW]** Q. Chen, F. Janda and R. Webb, *Virtual cycles of stable (quasi-)maps with fields*, Adv. Math. 385 (2021), 107781; [arXiv:1911.09825v2](https://arxiv.org/pdf/1911.09825v2). Theorem 1.2.3, printed p. 4, gives convex quantum Lefschetz. Corollary 1.2.1 gives the fields comparison, including the Euler-characteristic sign and the sum over target classes. Remark 1.2.2 discusses the Lefschetz identification of classes.
- **[FuP]** W. Fulton and R. Pandharipande, *Notes on stable maps and quantum cohomology*; [arXiv:alg-geom/9608011v2](https://arxiv.org/pdf/alg-geom/9608011v2). Theorem 2, printed p. 12, gives smoothness and dimension of the automorphism-free genus-zero stable-map locus for a convex target. Sections 5–6 provide the deformation and boundary constructions. This also supplies a finite-type realization of the embedded-conic parameter space.
- **[GP]** T. Graber and R. Pandharipande, *Localization of virtual classes*, Invent. Math. 135 (1999), 487–518; [arXiv:alg-geom/9708001v2](https://arxiv.org/pdf/alg-geom/9708001v2). Introduction equation (1) supplies virtual localization. Section 5, printed p. 19, describes a balanced rational curve's stable-map component and its normal obstruction bundle. The all-genus numerical evaluation used here is [FP].
- **[BF]** K. Behrend and B. Fantechi, *The intrinsic normal cone*, Invent. Math. 128 (1997), 45–88; [arXiv:alg-geom/9601010](https://arxiv.org/pdf/alg-geom/9601010). Sections 5 and 7 construct virtual classes from compatible obstruction theories. The normal obstruction comparison below is an application of this construction, not a claim that [BF] counts this particular Grassmannian's conics.
- **[HK]** B. Haghighat and A. Klemm, *Topological Strings on Grassmannian Calabi-Yau manifolds*, JHEP 01 (2009), 029; [arXiv:0802.2908](https://arxiv.org/pdf/0802.2908). Section 5.4 identifies this model. Table B.5, printed p. 31 (PDF page 32), has \(n_0(1)=196\), \(n_0(2)=1225\), and zeros in degrees 1 and 2 for the displayed genera 1 through 5. These are GV/BPS numbers, not ordinary \(N_{g,D}\). The abstract explicitly describes use of a mirror conjecture. The table is a consistency check, not the proof of the ordinary invariants.
- **[CJR-III]** Q. Chen, F. Janda and Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, supplied private draft dated October 17, 2025. [Supplied PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>). Rechecked through the existing extraction /tmp/log-glsm-audit/cjr-iii.txt: Theorem 8.4, equation (8.28), Remark 8.5, printed p. 84; §9.4, p. 92; §10.1, p. 93. Applicability and comparison below refer to that version, not to an independently certified published theorem.

The incidence and image-classification arguments in §§3–5 are deductions in this audit. No cited paper is being credited with stating this exact combined proof.

## 2. Geometry and effective curve classes

Put \(X=\operatorname{Gr}(2,7)\), \(H=c_1(S^\vee)\), and \(P=c_2(S^\vee)\). Then
\[
\dim X=10,\quad TX=S^\vee\otimes Q,\quad c_1(TX)=7H,
\quad E=\mathcal O_X(1)^{\oplus7},\quad [Y]=H^7.
\]
Bertini and adjunction give a general smooth threefold with \(K_Y=\mathcal O_Y\). The Koszul resolution, using Kodaira vanishing for \(\mathcal O_X(-j)=K_X\otimes\mathcal O_X(7-j)\), \(1\leq j\leq6\), gives \(h^1(\mathcal O_Y)=h^2(\mathcal O_Y)=0\).

Expansion of \(c(TX)/(1+H)^7\) gives
\[
c_2(TY)=(H^2+3P)|_Y,\qquad c_3(TY)=-7HP|_Y.
\]
For Chern roots \(x,y\) of \(S^\vee\), integration is
\[
\int_X f(x,y)=[x^6y^5](x-y)f(x,y).
\]
In particular,
\[
\int_XH^{10}=42,\quad \int_XH^8P=14,
\qquad \boxed{H_Y^3=42,\quad c_2(Y)H_Y=84,\quad\chi(Y)=-98.}
\]
These are direct intersection calculations, independent of [HK].

Successive smooth hyperplane sections down to complex dimension three preserve integral \(H_2\) by weak Lefschetz. Thus
\[
H_2(Y,\mathbb Z)\xrightarrow{\;\sim\;}H_2(X,\mathbb Z)=\mathbb Z[\mathrm{line}].
\]
Plücker degree therefore specifies the homology class used here; there is no hidden sum over distinct target classes. Equivalently, one may retain the ambient-class-sum convention of [CJW, Corollary 1.2.1]; here it has one class. Since \(H\) is ample, every nonconstant component has strictly positive integral degree. Nonconstant degree-zero components cannot occur.

## 3. One general choice of equations works for all three incidence problems

Let \(V\) be the 21-dimensional space of ambient Plücker linear forms, and let \(A=V^{\oplus7}\) be the affine space of ordered defining equations, of dimension 147. All genericity assertions below are open conditions in this same irreducible parameter space. Intersect their finitely many dense open subsets with the open locus of independent equations and smooth complete intersections.

### 3.1 Lines

The line parameter space is the smooth projective flag variety
\[
F=F(1,3;7),\qquad \dim F=14.
\]
For its tautological bundles \(U_1\subset U_3\), restriction of one equation gives
\[
V_1=U_1^\vee\otimes(U_3/U_1)^\vee,\qquad \operatorname{rk}V_1=2.
\]
The map \(V\otimes\mathcal O_F\to V_1\) is surjective: a Plücker line is an actual projective line. Hence seven equations give a generated rank-14 bundle. A general section has a reduced zero-dimensional zero scheme.

The count is \(\int_Fc_2(V_1)^7\). An independent exact evaluation in this audit gave **196** for each weight list \(\mu_i=4^{i-1}\) and \(\mu_i=5^{i-1}\), \(1\leq i\leq7\). The reproducible fixed-point sum is
\[
\sum_{a,\{b,c\}\subset\{1,\ldots,7\}\setminus\{a\}}
\frac{[-(\mu_a+\mu_b)]^7[-(\mu_a+\mu_c)]^7}
{(\mu_b-\mu_a)(\mu_c-\mu_a)
 \prod_{i\in\{a,b,c\}}\prod_{r\notin\{a,b,c\}}(\mu_r-\mu_i)}
=196.
\]
The pair \(\{b,c\}\) is unordered. This repeats the parent's flag-variety check without using its stable-map calculation.

### 3.2 Smooth embedded conics: nonproperness is harmless

Let \(Q\) be the open parameter space of smooth embedded Plücker conics in \(X\). It is an open subscheme of the Hilbert scheme with polynomial \(2m+1\), hence finite type and quasi-projective. Equivalently it is the smooth-domain embedded locus of \(\overline{\mathcal M}_{0,0}(X,2)\), with no map automorphisms.

The homogeneous tangent bundle \(TX\) is globally generated. For \(C\simeq\mathbb P^1\) embedded in \(X\), its quotient \(N_{C/X}\) is globally generated, so \(H^1(N_{C/X})=0\). Its rank is 9 and degree is \(14-2=12\). Therefore
\[
\dim Q=h^0(N_{C/X})=9+12=21
\]
at every conic. This also agrees with [FuP, Theorem 2].

Every embedded degree-two integral curve in projective space is a smooth plane conic. Linear forms on its plane restrict isomorphically to \(H^0(C,\mathcal O_C(1))=H^0(\mathbb P^1,\mathcal O(2))\). Restriction from the full Plücker projective space onto that plane is surjective. Thus for the universal conic \(p:\mathcal C\to Q\),
\[
B=p_*\mathcal O_{\mathcal C}(1),\qquad \operatorname{rk}B=3,
\qquad V\otimes\mathcal O_Q\twoheadrightarrow B.
\]
Base change applies because \(H^1(\mathbb P^1,\mathcal O(2))=0\). The seven restrictions form a section of \(B^{\oplus7}\), of rank 21.

More precisely, the incidence
\[
I_Q=\{(C,s)\in Q\times A:s|_C=0\}
\]
is the total space of the kernel of a surjective bundle map from the trivial rank-147 bundle onto \(B^{\oplus7}\). It is smooth, finite type, and of dimension 147. Characteristic-zero generic smoothness, applied to the dominant components of \(I_Q\to A\), gives smooth zero-dimensional general fibers. Nondominant components can be removed by excluding the closures of their images. There are only finitely many components.

A zero-dimensional finite-type scheme over \(\mathbb C\) is finite. Properness of \(Q\) is not required for this assertion. The general zero locus is thus a finite reduced set, possibly empty at this stage. Existence and its count are obtained in §6. This argument does not integrate a top Chern class over the nonproper space \(Q\).

### 3.3 Two distinct meeting lines: exclude the entire locus

Let \(U=\overline{\mathcal M}_{0,1}(X,1)\), the universal line over \(F\). It is smooth of dimension 15. Evaluation \(U\to X\) is smooth of relative dimension 5: lines through a fixed two-plane \(W\) are specified by a one-plane in \(W\) and a three-plane containing \(W\), a fiber \(\mathbb P^1\times\mathbb P^4\).

Consequently \(U\times_X U\) has dimension \(15+15-10=20\). Restrict to the open locus \(R\) of **distinct** lines. Their union \(C=L_1\cup L_2\) is a nodal plane conic, with
\[
h^0(\mathcal O_C(1))=2+2-1=3,\qquad H^1(\mathcal O_C(1))=0.
\]
Linear forms restrict surjectively to this three-dimensional space. Seven equations therefore impose rank 21 over \(R\). The incidence of such pairs and equations has dimension
\[
20+(147-21)=146<\dim A.
\]
The closure of its image in \(A\) is proper, so a general \(Y\) has no two distinct meeting lines.

The exclusion of the diagonal matters: for the same line twice, restriction has rank 2 per equation, not 3. Such stable maps are not excluded and are handled by local multiple-cover theory.

## 4. Normal bundles and double-line degenerations

For every line or smooth conic in the general \(Y\) selected above, transversality identifies the kernel of the derivative of the equation-restriction section with \(H^0(N_{C/Y})\). The kernel is zero. Adjunction on the smooth CY threefold gives
\[
\operatorname{rk}N_{C/Y}=2,\qquad \deg N_{C/Y}=-2.
\]
Splitting on \(\mathbb P^1\) then forces
\[
\boxed{N_{C/Y}=\mathcal O_{\mathbb P^1}(-1)^{\oplus2}.}
\]
For a conic, this \(\mathcal O(-1)\) is intrinsic to \(\mathbb P^1\); the restricted Plücker line bundle is \(\mathcal O(2)\).

There is no need to regard a Hilbert-scheme double line as a separate stable-map image. A map from a reduced nodal curve has reduced scheme-theoretic image. A degree-two map with set-theoretic image a line factors through that reduced line, possibly with two nonconstant components and arbitrary contracted subcurves.

One can also exclude the Hilbert boundary relevant to smooth plane conics directly. A flat limit of plane conics remains a plane conic: planes and quadratic equations have a projective parameter space. The singular possibilities are a pair of distinct meeting lines, already excluded, or a planar double line \(2L\). The latter would yield a conormal quotient
\[
N_{L/Y}^{\vee}\twoheadrightarrow\mathcal O_L(-1),
\]
because the square-zero ideal of \(L\) in a planar double line is \(\mathcal O_L(-1)\). But \(N_{L/Y}^{\vee}=\mathcal O(1)^{\oplus2}\) has no nonzero map to \(\mathcal O(-1)\). Thus no such double conic lies in this general \(Y\). This concerns embedded double lines, not stable-map double covers.

## 5. Exhaustion of stable maps and identification of their virtual contributions

The following argument is independent of domain genus. For a connected stable map \(f:D_{\mathrm{dom}}\to Y\), write the positive-degree image cycle as
\[
f_*[D_{\mathrm{dom}}]=\sum_i k_i[C_i],
\qquad \sum_i k_i\deg_H(C_i)\leq2.
\]
An integral projective curve of degree one is a line; an integral degree-two curve spans a plane and is a smooth conic over \(\mathbb C\). The latter classification follows from the bound \(\deg C\geq\dim\langle C\rangle\); an integral plane quadratic is nonsingular in characteristic zero.

If there are two distinct nonconstant images, both are lines. Connectedness forces them to meet: a connected chain of contracted components has one image point and cannot connect two disjoint curves. Section 3.3 excludes this case. Contracted components also cannot add an image point outside the positive-degree support.

Thus the entire map factors through one of the finitely many lines or conics. At total degree two its degree onto that support is respectively 2 or 1. A degree-one map onto a smooth conic has rational normalization on its nonconstant component; all extra genus lies in the rest of the stable domain. Degree-two maps to a line can also have noncontracted positive-genus components. Both kinds belong to the full stable-map space to \(\mathbb P^1\).

**Why the local model computes the ordinary virtual class.** Fix a support \(C\simeq\mathbb P^1\) and a positive cover degree \(k\). For every connected stable map \(f:D_{\mathrm{dom}}\to C\),
\[
H^0(D_{\mathrm{dom}},f^*N_{C/Y})=0.
\]
Indeed, each summand has negative degree on every nonconstant irreducible component, so its sections vanish there. On contracted components it is trivial; a section is constant and must vanish at a node meeting an already vanishing component. Connectedness propagates the vanishing through every contracted subcurve, including positive-genus components.

This also proves formal rigidity as a map into \(Y\), not just rigidity of the embedded curve. Induct over small extensions of Artinian bases. If a deformation factors through \(C\) modulo the square-zero kernel, its possible failure to factor through \(C\) is a homomorphism from the pulled-back conormal sheaf into that kernel, hence a section of \(f^*N_{C/Y}\) tensored with the kernel. It vanishes. This argument allows the domain itself to deform.

The closed immersion
\[
\overline{\mathcal M}_{g,0}(C,k)\longrightarrow
\overline{\mathcal M}_{g,0}(Y,k[C])
\]
is therefore also formally étale, hence open. It identifies a proper open-and-closed substack. The finitely many support substacks exhaust all geometric points in degrees 1 and 2. They do not intersect, even if two support curves happen to intersect as curves: a positive-degree map cannot factor through their zero-dimensional intersection. No assertion that all conics are pairwise disjoint is needed.

Pulling back
\[
0\longrightarrow TC\longrightarrow TY|_C\longrightarrow N_{C/Y}\longrightarrow0
\]
gives the compatible ordinary stable-map obstruction theories. The normal complex has \(R^0=0\), and cohomology and base change give the obstruction bundle
\[
\mathcal V=R^1\pi_*f^*N_{C/Y},\qquad
\operatorname{rk}\mathcal V=2(k+g-1).
\]
The virtual dimension of \(\overline{\mathcal M}_{g,0}(C,k)\) is the same rank. The excess-obstruction formula from the virtual-class construction identifies the restriction of the \(Y\) virtual class with
\[
c_{\mathrm{top}}(\mathcal V)\cap
[\overline{\mathcal M}_{g,0}(C,k)]^{\mathrm{vir}}.
\]
This is precisely the integral defined in [FP, equation (7)], consistent with [GP, §5]; [BF, §§5,7] supplies the compatible-obstruction framework. These are the ordinary theories, including all contracted tails, reducible domains, and automorphism factors. No positive-genus convex Lefschetz assertion is used.

## 6. Conic count and full ordinary invariants

For genus zero, \(\mathcal O_X(1)^{\oplus7}\) is convex on every genus-zero stable-map tree: each summand has nonnegative component degrees, and normalization gives \(H^1=0\). Thus [CJW, Theorem 1.2.3] identifies the ordinary invariant with the **nondual** Euler-twisted ambient integral
\[
N_{0,D}=
\int_{[\overline{\mathcal M}_{0,0}(X,D)]^{\mathrm{vir}}}
e(R^0\pi_*f^*E).
\]
The parent reports exact GKM evaluation at two rays:
\[
N_{0,1}=196,\qquad N_{0,2}=2499/2.
\]
Agreement between two rays is a computational consistency check, not itself a proof of the localization implementation. This audit takes the stated exact calculation as input; [GP] and [CJW] establish the mathematical formula it must evaluate. If the stored genus-zero result is from the dual Euler twist at \(t=0\), first multiply by \((-1)^{7(D+1)}\): the raw checks are \(196\) at \(D=1\) and \(-2499/2\) at \(D=2\).

By §§3–5, a smooth conic contributes 1 to genus-zero degree two, whereas each line contributes \(1/8\). Therefore the actual number of smooth conics is
\[
Q_Y=\frac{2499}{2}-\frac{196}{8}=\boxed{1225}.
\]
This subtraction is justified by the geometric exhaustion and normal-bundle proof, not merely by labeling a mirror instanton number a curve count.

[FP, Theorem 3] gives, for \(g\geq2\),
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}\,k^{2g-3},
\qquad C(2,k)=\frac{k}{240},\quad C(3,k)=\frac{k^3}{6048}.
\]
Consequently
\[
N_{2,2}=1225C(2,1)+196C(2,2)
=\frac{1225+392}{240}=\frac{539}{80},
\]
\[
N_{3,2}=1225C(3,1)+196C(3,2)
=\frac{1225+1568}{6048}=\frac{133}{288}.
\]
The constant-map formula in [FP, following Theorem 4],
\[
N_{g,0}=\frac{(-1)^g\chi(Y)}2
\frac{|B_{2g}|}{2g}\frac{|B_{2g-2}|}{2g-2}
\frac1{(2g-2)!},
\]
completes the requested scalar bounds:

| Genus | \(N_{g,0}\) | \(N_{g,1}\) | \(N_{g,2}\) |
|---|---:|---:|---:|
| 2 | \(-49/2880\) | \(49/60\) | \(539/80\) |
| 3 | \(7/103680\) | \(7/216\) | \(133/288\) |

The argument proves the full connected ordinary numerical invariants in this bound, subject to the explicitly supplied genus-zero numerical input. It does not certify an implementation's presentation of their stabilized Chow classes merely from equality of integrals.

## 7. CJR grading, signs, ambient dimension, and scope

In the notation of the existing audit, set
\[
X=\operatorname{Gr}(2,7),\quad E=E_1=\mathcal O_X(1)^{\oplus7},
\quad L=\mathcal O_X,\quad r=a=\widetilde r=1.
\]
All seven summands have the same GLSM grading. The infinity target is \(\mathbb P_X(E^\vee)\); its full aggregate cannot be replaced by a product of seven rank-one infinity answers. Section 9.4 of the supplied draft allows this projective base and rank-seven bundle.

The target virtual dimension is zero. Riemann–Roch and stabilization to \(\overline{\mathcal M}_g\) give
\[
\boxed{\alpha=3g-3,\qquad
\sigma_{g,D}=(-1)^{7(1-g+D)},\qquad
\nu_X=(1-g)(10-3)+7D=7(D+1-g).}
\]

| \(g\) | \(D\) | \(\alpha\) | \(\sigma_{g,D}\) | Ambient \(\nu_X\) | Signed ordinary number \(\sigma_{g,D}N_{g,D}\) |
|---|---:|---:|---:|---:|---:|
| 2 | 0 | 3 | \(-1\) | \(-7\) | \(49/2880\) |
| 2 | 1 | 3 | \(+1\) | \(0\) | \(49/60\) |
| 2 | 2 | 3 | \(-1\) | \(7\) | \(-539/80\) |
| 3 | 0 | 6 | \(+1\) | \(-14\) | \(7/103680\) |
| 3 | 1 | 6 | \(-1\) | \(-7\) | \(-7/216\) |
| 3 | 2 | 6 | \(+1\) | \(0\) | \(133/288\) |

Retain the separate objects
\[
G_{g,D}=\mathrm{st}_{Y,*}[\overline{\mathcal M}_{g,0}(Y,D)]^{\mathrm{vir}},
\]
\[
Z_D(t)=\mathrm{st}_{X,*}
\left(e_t((R\pi_*f^*E)^\vee)
\cap[\overline{\mathcal M}_{g,0}(X,D)]^{\mathrm{vir}}\right),
\]
and the **full** CJR infinity aggregate \(B_D(t)\). In the conventions of [CJR-III, §9.4],
\[
[t^0](Z_D+B_D)=\sigma_{g,D}G_{g,D}.
\]
This is a class identity; the last table gives its integrated right side.

After the auxiliary torus is removed,
\[
Z_D(t)=\sum_{r=0}^{3g-3}t^{\,3g-3-r}Z_{D,r}.
\]
Thus there is **no negative-\(t\) range** in the requested stabilized unmarked CY3 problem. Remark 8.5's negative-power cancellation does not by itself produce a relation here. Negative ambient virtual dimension forces \(Z_D=0\) for \((g,D)=(2,0),(3,0),(3,1)\), but the corresponding ordinary target invariants above are nonzero and must be supplied by the full comparison.

No scan, Pixton reduction, or equality of individual infinity terms is certified by this note. The prior fourteen-target audit and results remain separate; this proof only resolves the additional Gr(2,7)[1^7] low-degree ordinary numerical inputs.
