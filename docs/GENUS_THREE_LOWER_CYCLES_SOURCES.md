# Genus-three lower-codimension CJR/O(3) cycles: sources and justification

Research date: 2026-11-05. Scope: unmarked genus three, ambient degrees \(D=0,1\), and codimensions \(r=3,4,5\), with the top class \(r=6\) as a consistency check. This note supplies mathematical justification and an independent degree-zero derivation. It does not certify the running degree-one computation or determine individual punctured effective cycles.

The proposed deliverable is valid when formulated as an **aggregate identity in rational Chow**:
\[
 B_{D,r}=-Z_{D,r},\qquad r=3,4,5,
\]
where \(B_{D,r}\) means the sum of *all* infinity-containing CJR graph contributions after pushforward to \(\overline{\mathcal M}_3\), and \(Z_{D,r}\) is the indicated coefficient of the single all-zero-vertex contribution. A separate auxiliary torus calculation expresses \(Z_{D,r}\) tautologically and supplies pure tautological relations from its auxiliary poles.

## 1. Primary sources and precise scope

**[CJR-III]** Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, private draft dated **October 17, 2025**. Consulted the supplied PDF/extraction, not a purported public substitute:
[local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>).

- Theorem 8.4, equation (8.28), and Remark 8.5, printed p. 84: localization is an equality in localized equivariant **Chow**, and every negative \(t\)-coefficient on its right vanishes.
- Section 9.3, pp. 89–90: hypersurface specialization, stable zero-vertex dual Euler twist, edge and unstable-vertex factors, and reduced infinity factor.
- Equation (9.7), pp. 89–90: the \(t^0\) specialization gives the hypersurface virtual class with sign \((-1)^{1-g+\beta\cdot c_1(E)}\).
- Example 9.2, pp. 90–92: an explicit mapping-to-a-point check and cancellation of a \(t^{-1}\) coefficient.
- Remark 10.7, p. 94: negative powers give relations between twisted GW and effective invariants.
- Proposition 10.9, p. 95, assumes tautologicality of effective cycles; it does not prove that assumption in general.

**[CJR-P]** Q. Chen, F. Janda, Y. Ruan, *Punctured logarithmic R-maps*, [arXiv:2208.04519v2](https://arxiv.org/abs/2208.04519v2), [public PDF](https://arxiv.org/pdf/2208.04519v2). The supplied extraction identifies arXiv v2, **21 April 2023**; the user's April 25 date should not replace that printed version date. Relevant items: definition of \(\psi_{\min}\), §8.1.1, p. 78; Theorem 8.1, equations (8.3)–(8.4), p. 78; Corollary 8.2, p. 79; Theorem 8.4, equation (8.8), p. 80; §9.8, Question 9.36 (tautologicality of effective cycles remains a question in the stated generality).

**[GP]** T. Graber, R. Pandharipande, *Localization of virtual classes*, Invent. Math. **135** (1999), 487–518, [arXiv:alg-geom/9708001](https://arxiv.org/abs/alg-geom/9708001), [PDF](https://arxiv.org/pdf/alg-geom/9708001). Equation (1) gives virtual localization in equivariant Chow. Section 4 (preprint pp. 13–18) gives fixed graphs for projective space, tangent weights, normal factors, and the extension
\[
1\longrightarrow\prod_e\mu_{d_e}\longrightarrow A_\Gamma
\longrightarrow\operatorname{Aut}(\Gamma)\longrightarrow1.
\]
This is the source for graph/deck factors and Hodge–cotangent factors.

**[M]** D. Mumford, *Towards an Enumerative Geometry of the Moduli Space of Curves*, in *Arithmetic and Geometry II*, Progress in Math. **36** (1983), 271–328, [author-hosted scan](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf). Section 5, formula (5.2), Corollary (5.3), and (5.4), printed pp. 305–307: GRR, \(\operatorname{ch}_{2k}(\mathbb E)=0\) for \(k>0\), and \(c(\mathbb E)c(\mathbb E^\vee)=1\). Use the algebraic GRR statement for Chow, rather than relying only on its subsequent cohomological explanation.

**[FP]** C. Faber, R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. **139** (2000), 173–199, [arXiv:math/9810173](https://arxiv.org/abs/math/9810173), [PDF](https://arxiv.org/pdf/math/9810173). Introduction, equation (1) and the following paragraph: constant maps and their obstruction bundle. Theorem 4 and §4.3: the Bernoulli formula for \(\int\lambda_{g-1}^3\), and its relation to \(2\int\lambda_g\lambda_{g-1}\lambda_{g-2}\).

**[PPZ]** R. Pandharipande, A. Pixton, D. Zvonkine, *Relations on \(\overline{\mathcal M}_{g,n}\) via 3-spin structures*, J. Amer. Math. Soc. **28** (2015), 279–309, [arXiv:1303.1043](https://arxiv.org/abs/1303.1043), [PDF](https://arxiv.org/pdf/1303.1043). Definitions 0.1–0.2, Theorem 1 and Corollary 2: 3-spin/Pixton relations in **cohomology**.

**[J]** F. Janda, *Relations on \(\overline{\mathcal M}_{g,n}\) via equivariant Gromov–Witten theory of \(\mathbb P^1\)*, Algebraic Geometry **4** (2017), 311–336, [journal PDF](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf), [arXiv:1509.08421](https://arxiv.org/abs/1509.08421). Theorem 1.1, printed p. 313, proves the relations in the **tautological Chow ring**. It distinguishes this assertion from the conjecture that their closure generates all tautological relations.

## 2. Definitions, grading, and the aggregate CJR relation

Let
\[
\mathcal M_D=\overline{\mathcal M}_{3,0}(\mathbb P^2,D),\qquad
s:\mathcal M_D\longrightarrow\overline{\mathcal M}_3,
\]
and let \(\pi,f\) be the universal curve and map. All Chow groups have rational coefficients. Equip \((f^*\mathcal O(3))^\vee\) with the CJR character of weight \(+t\); its ordinary fibre Chern class is \(-3H\). Define
\[
Z_D(t)=s_*\left(
 e_{\mathbb R}\bigl((R\pi_*f^*\mathcal O(3))^\vee\bigr)
 \cap[\mathcal M_D]^{\mathrm{vir}}\right),\qquad
Z_{D,r}=[t^{\,2-r}]Z_D(t).
\]
The Euler class here is the Euler ratio of a derived complex, **not** the Euler class of \(R^0\) alone. Its interpretation and sign agree with the Serre-duality form of the stable zero-vertex factor in [CJR-III, §9.3].

Here is the grading calculation:
\[
\operatorname{vdim}\mathcal M_D=(1-3)(2-3)+3D=2+3D,\qquad
\chi(f^*\mathcal O(3))=1-3+3D=3D-2.
\]
Capping with the Euler ratio therefore gives total equivariant dimension \(4\). Since \(\dim\overline{\mathcal M}_3=6\), its pushforward has total codimension \(2\):
\[
Z_D(t)=\sum_{r=0}^{6}t^{2-r}Z_{D,r},\qquad Z_{D,r}\in A^r(\overline{\mathcal M}_3).
\]
Tautologicality for these \(Z_{D,r}\) follows from the concrete auxiliary fixed-locus expressions below.

Write the reduced CJR localization identity, composed with curve stabilization, as
\[
L_D(t)=Z_D(t)+I_D(t),\qquad
I_D(t)=
\sum_{\Gamma:\,V_\infty(\Gamma)\ne\varnothing}
\frac{\operatorname{st}_{\Gamma,*}(\text{CJR graph class})}
     {|\operatorname{Aut}\Gamma|}.
\]
The summand without an infinity vertex is the single zero vertex carrying genus three and degree \(D\). Proper pushforward of [CJR-III, Theorem 8.4] preserves the identity, and its unlocalized left side has no negative powers of \(t\), by Remark 8.5. Thus, defining \(B_{D,r}=[t^{2-r}]I_D(t)\), one obtains
\[
\boxed{B_{D,r}+Z_{D,r}=0\quad(3\le r\le6).}
\]
The lower-codimension deliverable uses \(t^{-1},t^{-2},t^{-3}\), hence codimensions \(3,4,5\). The existing \(t^{-4}\) relation is codimension \(6\). Codimension \(2\) is the constant term, not a negative-pole relation.

This identity proves that the **aggregate** \(B_{D,r}\) equals a tautological class. It supplies neither separate values nor separate tautological presentations for its infinity-vertex constituents. Also, “infinity-containing” does not mean supported on the boundary of \(\overline{\mathcal M}_3\): stabilization may contract rational tails and yield classes on the whole moduli space.

The sign in (9.7) is \((-1)^{3D-2}\), but concerns the constant term. It must not be inserted as an additional sign into the negative-pole identity. These conclusions use [CJR-III, §9.3, Theorem 8.4 and Remarks 8.5, 10.7].

## 3. Independent degree-zero derivation

Put \(h=c_1(\mathcal O_{\mathbb P^2}(1))\), with \(h^3=0\) and \(\int_{\mathbb P^2}h^2=1\), and write \(\lambda_i=c_i(\mathbb E)\). Constant maps give
\[
\mathcal M_0=\overline{\mathcal M}_3\times\mathbb P^2,\qquad
[\mathcal M_0]^{\mathrm{vir}}
=c_6(\mathbb E^\vee\boxtimes T_{\mathbb P^2})
 \cap[\overline{\mathcal M}_3\times\mathbb P^2].
\]
This is [FP, introduction], also used in [CJR-III, Example 9.2].

If \(x_1,x_2,x_3\) are the Hodge roots, then \(c_1(T_{\mathbb P^2})=3h\) and \(c_2(T_{\mathbb P^2})=3h^2\) give the following direct calculation:
\[
\begin{aligned}
c_6(\mathbb E^\vee\boxtimes T_{\mathbb P^2})
&=\prod_{i=1}^3(x_i^2-3hx_i+3h^2)\\
&=\lambda_3^2-3h\lambda_2\lambda_3
 +3h^2(\lambda_2^2+\lambda_1\lambda_3).
\end{aligned}
\]
Mumford's relation specializes to
\[
\lambda_1^2=2\lambda_2,\qquad
\lambda_2^2=2\lambda_1\lambda_3,\qquad
\lambda_3^2=0.
\]
Consequently the obstruction class is
\[
-3h\lambda_2\lambda_3+9h^2\lambda_1\lambda_3.
\]

For constant maps, in equivariant \(K\)-theory,
\[
(R\pi_*f^*\mathcal O(3))^\vee
=(\mathcal O-\mathbb E)\boxtimes\mathcal O(-3).
\]
Define
\[
F(u)=\frac{u}{u^3+\lambda_1u^2+\lambda_2u+\lambda_3}.
\]
The dual Euler ratio is \(F(t-3h)\). Extracting \(h^2\) therefore yields the particularly useful formula
\[
\boxed{Z_0(t)=9\lambda_1\lambda_3F(t)
                    +9\lambda_2\lambda_3F'(t).}
\]
Since
\[
F(t)=t^{-2}-\lambda_1t^{-3}+\lambda_2t^{-4}+O(t^{-5}),
\qquad
F'(t)=-2t^{-3}+3\lambda_1t^{-4}+O(t^{-5}),
\]
and all terms of Chow degree above six vanish, this becomes
\[
\boxed{
Z_0(t)=9\lambda_1\lambda_3\,t^{-2}
       -36\lambda_2\lambda_3\,t^{-3}
       +36\lambda_1\lambda_2\lambda_3\,t^{-4}.}
\]
In particular,
\[
Z_{0,3}=0,\qquad
Z_{0,4}=9\lambda_1\lambda_3,\qquad
Z_{0,5}=-36\lambda_2\lambda_3,\qquad
Z_{0,6}=36\lambda_1\lambda_2\lambda_3.
\]
Thus the degree-zero lower-codimension CJR aggregates are
\[
B_{0,3}=0,\qquad
B_{0,4}=-9\lambda_1\lambda_3,\qquad
B_{0,5}=36\lambda_2\lambda_3.
\]
These are derived here, independently of the fixed-locus backend, from mapping to a point and [M, §5].

For the scalar check, [FP, Theorem 4 and §4.3] gives
\[
\int_{\overline{\mathcal M}_g}\lambda_g\lambda_{g-1}\lambda_{g-2}
=\frac{|B_{2g}|\,|B_{2g-2}|}
 {2(2g)(2g-2)(2g-2)!}.
\]
At \(g=3\), this is \(1/1451520\). Hence
\[
\int_{\overline{\mathcal M}_3}Z_{0,6}
=\frac{36}{1451520}=\frac1{40320}.
\]
The rational arithmetic was checked separately using Python's exact fractions inside OrbStack.

## 4. Degree one: an explicit class-valued fixed-graph formula

This section gives a direct specialization of [GP, §4] with the dual Euler twist, rather than a claim that the cited papers print these particular genus-three coefficients.

Write the auxiliary hyperplane restrictions as \(a_0,a_1,a_2\), to distinguish them from Hodge classes. Use the convention
\[
H|_{p_i}=a_i,\qquad
w_{ij}=a_i-a_j,\qquad q_i=t-3a_i.
\]
The tangent weights are \(a_i-a_j\), as in [GP]. For \(\{i,j,k\}=\{0,1,2\}\), set
\[
\Lambda_g^\vee(u)=\sum_{\ell=0}^g(-1)^\ell\lambda_\ell u^{g-\ell},
\qquad
\Lambda_g^+(u)=\sum_{\ell=0}^g\lambda_\ell u^{g-\ell},
\]
and define on \(\overline{\mathcal M}_{g,1}\)
\[
V_{g,i|j}=
\frac{\Lambda_g^\vee(w_{ij})\Lambda_g^\vee(w_{ik})}
     {\Lambda_g^+(q_i)(w_{ij}-\psi)}.
\]
The degree-one edge factor is
\[
E_{ij}=
-\frac{\displaystyle\prod_{m=0}^{3}
       \bigl(t-(3-m)a_i-ma_j\bigr)}
       {(a_i-a_j)^2(a_i-a_k)(a_j-a_k)}
=E_{ji}.
\]
The numerator is the Euler class of the four edge sections of the dual-twisted \(\mathcal O(3)\). Normalization at a valence-one stable vertex cancels its \(H^0\) fibre factor; this explains the absence of an additional \(q_i\) in \(V_{g,i|j}\).

There is one degree-one edge, so genus is entirely at its endpoints. There are exactly two forms: a genus-three contracted component with a rational tail, or genus-one and genus-two components joined by a rational bridge. Let
\[
p:\overline{\mathcal M}_{3,1}\to\overline{\mathcal M}_3,\qquad
\xi:\overline{\mathcal M}_{1,1}\times\overline{\mathcal M}_{2,1}
 \to\overline{\mathcal M}_3.
\]
With each ordered pair \(i\ne j\) occurring once, the full auxiliary expression is
\[
\boxed{
Z_1(t;a)=
\sum_{i\ne j}E_{ij}\left[
 w_{ji}\,p_*V_{3,i|j}
 +\xi_*\bigl(V_{1,i|j}\boxtimes V_{2,j|i}\bigr)
 \right].}
\]
In the first term \(i\) is the genus-three endpoint. In the second, \(i\) is genus one and \(j\) is genus two. This is six tails and six bridges. The univalent unstable endpoint contributes \(w_{ji}\). Edge degree one has no nontrivial deck group, and distinct endpoint genera prohibit exchanging endpoints. There is **no extra factor \(1/2\)** in this indexing.

Stabilization must precede integration:
\[
p_*(1)=0,\qquad p_*(\psi^a)=\kappa_{a-1}\quad(a\ge1),\qquad
p_*(\psi^a p^*\alpha)=\kappa_{a-1}\alpha.
\]
In particular \(\kappa_0=2g-2=4\), and a tail monomial of degree \(a+\deg\alpha\) becomes codimension \(a-1+\deg\alpha\). The bridge is contracted to a node; the pushforward codimension is
\[
1+a+b+\deg\alpha+\deg\beta
\]
for \(\xi_*(\psi_1^a\alpha\boxtimes\psi_2^b\beta)\). Here \(\xi\) is the ordinary gluing pushforward. Stack stabilizers are already included in Chow pushforward; adding a manual elliptic stabilizer factor would double-count them.

The vertex expressions are polynomials after truncation in Chow degree: retain degrees through seven on \(\overline{\mathcal M}_{3,1}\), through one on \(\overline{\mathcal M}_{1,1}\), and through four on \(\overline{\mathcal M}_{2,1}\). Integrating a stable vertex first loses the lower-dimensional information needed here. These observations match the intended architecture of \`o3_unmarked_cycles.sage\`.

## 5. Auxiliary poles and exact Pixton certificates

Set \(a=\epsilon(0,1,3)\) and keep \(t\) invertible. The global twisted class exists without inverting auxiliary characters: expanding its Euler ratio in \(t^{-1}\) introduces no inverse powers of \(\epsilon\). The fixed-graph sum may contain such poles termwise, so every negative auxiliary coefficient vanishes **after summation in Chow**.

For each ordinary Chow codimension \(r\), homogeneity has the form
\[
[Z_D(t;\epsilon a)]_{\text{Chow degree }r}
=\sum_k\epsilon^k t^{2-r-k}C_{D,r,k}(a).
\]
Thus one may compute at \(t=1\) while retaining the Chow grading:
\[
C_{D,r,k}=0\quad(k<0),\qquad C_{D,r,0}=Z_{D,r}.
\]
The first equalities are pure tautological relations; the second gives the twisted class. Setting \(t=1\) is not the CJR nonequivariant specialization \(t=0\).

Repeat with \((0,2,5)\). An exact check of
\[
C_{D,r,0}(0,1,3)-C_{D,r,0}(0,2,5)=0
\]
is a useful independent implementation check. Two rays alone do not prove regularity in every direction; the geometric global-class argument supplies that proof. Conversely, checking only the finite parts is insufficient: every negative power must be tested, not silently discarded.

For a computational certificate, put both a pole and the Pixton relations in the same rational strata coordinates. If \(v\) is the pole vector and the rows of \(M\) are valid 3-spin relations (including the required tautological operations), save a rational vector \(u\) satisfying
\[
\boxed{uM=v.}
\]
This is a sufficient Chow-vanishing certificate by [J, Theorem 1.1]. PPZ alone establishes cohomological vanishing. General completeness of the Pixton relations is unnecessary for this implication. Failure to find \(u\) is inconclusive: it could signal an implementation error, an incomplete generated relation space, or a relation outside that space.

An additional useful certificate is the same equation for the difference between the two weight rays. Save the strata ordering and graph normalization as well as the matrix: identical vector lengths do not guarantee identical coordinates.

The target value
\[
\int_{\overline{\mathcal M}_3}Z_{1,6}=-\frac{1637}{80640}
\]
is a **task-supplied computational benchmark**, not a value independently established by this literature audit. No cited primary source was found asserting that particular number. Its agreement with the older scalar backend checks normalization, but does not certify \(Z_{1,3},Z_{1,4},Z_{1,5}\).

## 6. Punctured descendants and limits of the conclusions

The descendant at an infinity vertex is
\[
\psi_{\min}=c_1\bigl(\mathcal O(-e_{\min})\bigr),
\]
not a cotangent-line class at a contact marking [CJR-P, §8.1.1]. In the reduced CJR formula its factor is
\[
\frac{t}{-t-\psi_{\min}}
=-1+\frac{\psi_{\min}}t-\frac{\psi_{\min}^2}{t^2}
 +\frac{\psi_{\min}^3}{t^3}-\cdots .
\]
Ordinary stable-vertex flags instead carry \((w-\psi)^{-1}\). Interchanging these two descendants changes both the geometry and the signs.

Suppressing the uniform-degeneracy pushforwards only in notation, Theorem 8.1 gives the pattern
\[
\operatorname{Forget}_*(\psi_{\min,+}^{\,k}[\mathcal R_+]^{\rm red})
=\sum_j\frac{|c_j|}{r_j}
 \sum_{\ell=0}^{k-1}
 (\operatorname{ev}_j^*\psi_{\rm DF})^\ell
 \psi_{\min}^{\,k-1-\ell}[\mathcal R]^{\rm red}.
\]
The actual statement includes \(F_{1,*}\circ F_*^{\rm f}\); those morphisms must be restored in a class-level use. The divisor equation adds the curve-degree term and divisor evaluations. At \(k=0\) the displayed inner sum is empty. For stabilized effective classes, Theorem 8.4, equation (8.8), contains explicit boundary-divisor corrections. Ordinary GW string/dilaton manipulations cannot simply be substituted for these formulas.

There are three distinct outcomes:

1. **Established here:** the CJR aggregate identity in codimensions three through five (also six), the complete independent degree-zero formulas, and a class-valued degree-one localization expression with stabilization conventions.
2. **To be certified by the parent calculation:** the expanded degree-one finite classes, every auxiliary pole relation, and exact 3-spin membership and second-ray differences. This note does not claim those computations have finished.
3. **Not established:** reconstruction of individual positive-dimensional effective cycles, new relations independent of Pixton, or blanket tautologicality of effective cycles. [CJR-P, §9.8, Question 9.36] and [CJR-III, Proposition 10.9] require keeping these claims separate.

Finally, distinguish ambient curve degree \(D\) from codimension \(r\). This calculation gives lower-codimension relations at ambient degrees zero and one; it does not investigate every higher ambient degree. A scalar top-degree result cannot recover positive-dimensional classes. Pairing lower classes with a few complementary classes is only a numerical test unless a suitable nondegeneracy theorem is separately supplied. Exact relation-space certificates avoid that additional assumption.
