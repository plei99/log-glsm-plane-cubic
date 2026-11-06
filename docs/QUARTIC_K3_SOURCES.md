# Smooth quartic K3: ordinary GW theory and CJR localization sources

Research date: 2026-11-05. Work over the complex numbers with rational Chow groups. Fix a smooth quartic surface \(S\subset\mathbb P^3\), use the hypersurface bundle \(L=\mathcal O_{\mathbb P^3}(4)\), and consider connected, unmarked genus \(g\ge2\). The symbol \(D\) is ambient degree, \(\int_\beta H\), not the coefficient of a primitive K3 curve class. This note verifies mathematical inputs and derives constant-map checks; it does not verify any parent implementation or compute individual higher-genus effective cycles.

The ordinary K3 virtual class vanishes for every positive curve class, and also for degree zero in the stated genus range. Consequently the CJR infinity aggregate satisfies
\[
\boxed{B_{D,r}=-Z_{D,r}\quad\text{for }2g-2\le r\le3g-3,\ D\ge0.}
\]
For \(r>2g-2\), this is the negative-power localization relation. At \(r=2g-2\), it additionally uses the ordinary K3 vanishing. This does not assert that the entire equivariant CJR class is zero.

## 1. Primary sources and exact locators

**[CJR-III]** Q. Chen, F. Janda, Y. Ruan, *Structural Formulae in Logarithmic Gauged Linear Sigma Models III: The Localization Formula*, supplied private draft dated **October 17, 2025**. Source: [local PDF](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>), checked through the supplied text extraction `/tmp/log-glsm-audit/cjr-iii.txt`. Page numbers below are printed draft pages.

- Theorem 8.4, equation (8.28), and Remark 8.5, p. 84: localization in equivariant Chow; the left side has no negative powers of the GLSM parameter \(t\).
- §9.3, pp. 89–90: hypersurface specialization; equation (9.7), starting on p. 89, identifies the nonequivariant class with signed **ordinary** hypersurface GW theory. The stable zero-vertex Euler ratio and infinity factors are on p. 90.
- Example 9.2, pp. 90–92: constant-map consistency check.
- Remark 10.7, p. 94: negative-power relations. Proposition 10.9, p. 95: a conditional tautologicality statement.

The requested “Theorem 9.7, 9.3” locators mean **equation (9.7)** and **section 9.3** in this draft. They are not theorem numbers. The localization claims here retain the scope of this private version.

**[CJR-GLSM]** Q. Chen, F. Janda, Y. Ruan, *The logarithmic gauged linear sigma model*, Invent. Math. **225** (2021), 1077–1154; [arXiv:1906.04345v2](https://arxiv.org/abs/1906.04345v2), [PDF](https://arxiv.org/pdf/1906.04345v2). Proposition 4.1 and its proof, preprint pp. 34–35; also Corollary 1.10, p. 7. This is reference [19] used in draft equation (9.7). Remark 1.8, p. 6, says the canonical and reduced log-GLSM virtual classes have the same dimension. The comparison recovers the ordinary class of the zero locus. The preceding paragraph to Proposition 4.1 explains the union over curve classes on that zero locus with the same ambient pushforward.

**[MPT]** D. Maulik, R. Pandharipande, R. P. Thomas, with appendix by A. Pixton, *Curves on K3 surfaces and modular forms*, J. Topol. **3** (2010), 937–996; [arXiv:1001.2719](https://arxiv.org/abs/1001.2719), [author PDF](https://people.math.ethz.ch/~rahul/mpt.pdf). Introduction §0.1, preprint p. 2, explicitly gives ordinary vanishing for every nonzero effective class through the holomorphic symplectic trivial obstruction quotient, and distinguishes the reduced class of dimension \(g\) from ordinary dimension \(g-1\). No primitivity hypothesis is needed for this vanishing.

**[M]** D. Mumford, *Towards an Enumerative Geometry of the Moduli Space of Curves*, in *Arithmetic and Geometry II*, Progress in Math. **36** (1983), 271–328; [author-hosted scan](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf). §5, formula (5.2), Corollary (5.3), and (5.4), printed pp. 305–307. The algebraic GRR calculation yields vanishing of positive even Chern-character components, hence \(c(\mathbb E)c(\mathbb E^\vee)=1\) in rational Chow on the compactification. The subsequent discussion also gives a cohomological explanation; that explanation alone is not the justification for a Chow claim.

**[FP]** C. Faber, R. Pandharipande, *Hodge integrals and Gromov–Witten theory*, Invent. Math. **139** (2000), 173–199; [arXiv:math/9810173](https://arxiv.org/abs/math/9810173), [PDF](https://arxiv.org/pdf/math/9810173). Introduction equation (1) and the following paragraph, preprint p. 1: degree-zero moduli and obstruction bundle. Theorem 4, p. 4, and §4.3, pp. 20–21: the Bernoulli integral and the relation between \(\lambda_{g-1}^3\) and \(2\lambda_g\lambda_{g-1}\lambda_{g-2}\).

**[GP]** T. Graber, R. Pandharipande, *Localization of virtual classes*, Invent. Math. **135** (1999), 487–518; [arXiv:alg-geom/9708001](https://arxiv.org/abs/alg-geom/9708001), [PDF](https://arxiv.org/pdf/alg-geom/9708001). Equation (1) and §4, preprint pp. 13–18: projective-space graphs, automorphism/deck factors, moving tangent terms, and Hodge/node factors. The quartic edge and vertex formulas below are direct specializations and derivations, not formulas claimed to appear verbatim there.

**[CJR-P]** Q. Chen, F. Janda, Y. Ruan, *Punctured logarithmic R-maps*; [arXiv:2208.04519v2](https://arxiv.org/abs/2208.04519v2), [PDF](https://arxiv.org/pdf/2208.04519v2). The arXiv banner reads **21 April 2023**; the manuscript's first-page date reads **April 25, 2023**. §9, pp. 86–88, defines effective cycles and gives balancing (9.5), dimension (9.8), and Corollary 9.3. Definition 9.19 through Corollary 9.24, pp. 94–95, impose additional hypotheses for reduction to basic scalar invariants. §9.8, Question 9.36, p. 98, asks about tautologicality.

## 2. Ordinary K3 vanishing, including constant maps

### Positive degree

For \(\beta\ne0\), [MPT, §0.1] gives
\[
[\overline{\mathcal M}_{g}(S,\beta)]^{\rm vir}=0.
\]
The holomorphic symplectic form identifies tangent and cotangent directions; composition with the differential of a nonconstant map and the trace on the curve gives the trivial obstruction quotient. The quotient is precisely what reduced K3 GW theory removes. The assertion is a vanishing of the virtual cycle, not just its numerical degree. An empty moduli space is not required.

Since \(H|_S\) is ample, every effective class of positive ambient degree is nonzero. Thus each class in the ambient-degree sum vanishes separately. In particular, one must not identify \(H_2(S)\) with \(H_2(\mathbb P^3)\), or assume a quartic has only one curve class of a given degree.

### Degree zero: an independent obstruction calculation

A constant map has zero differential, so the preceding surjectivity argument cannot simply be reused. Instead [FP, introduction] gives
\[
\overline{\mathcal M}_{g}(S,0)=\overline{\mathcal M}_g\times S,
\qquad
[\overline{\mathcal M}_{g}(S,0)]^{\rm vir}
=c_{2g}(\mathbb E^\vee\boxtimes TS)
 \cap[\overline{\mathcal M}_g\times S].
\]
Here \(\mathbb E\) is the rank-\(g\) Hodge bundle and \(\lambda_i=c_i(\mathbb E)\). Let \(x_1,\ldots,x_g\) be its formal roots. The symplectic rank-two bundle \(TS\) has roots \(y,-y\), with \(c_2(TS)=-y^2\). Therefore
\[
\begin{aligned}
c_{2g}(\mathbb E^\vee\boxtimes TS)
 &=\prod_{i=1}^g(-x_i+y)(-x_i-y)\\
 &=\prod_{i=1}^g(x_i^2+c_2(TS))\\
 &=\lambda_g^2+
 c_2(TS)\bigl(\lambda_{g-1}^2-2\lambda_g\lambda_{g-2}\bigr).
\end{aligned}
\]
The last equality uses \(c_2(TS)^2=0\) on the surface and \(\lambda_0=1\). Mumford's relation implies both
\[
\lambda_g^2=0,
\qquad \lambda_{g-1}^2-2\lambda_g\lambda_{g-2}=0
\quad(g\ge2).
\]
Indeed these are the degree \(2g\) and \(2g-2\) coefficients of \(c(\mathbb E)c(\mathbb E^\vee)=1\). Thus
\[
\boxed{[\overline{\mathcal M}_{g}(S,0)]^{\rm vir}=0\quad(g\ge2).}
\]
Equivalently, every positive elementary symmetric polynomial in the \(x_i^2\) vanishes, leaving \(c_2(TS)^g=0\). This calculation proves cycle-level vanishing before stabilization or integration. It must not be extrapolated to genus one: on the stable marked genus-one constant-map space the same calculation leaves \(c_2(TS)\), which need not vanish. Unmarked genus-one constant maps are unstable.

The quartic Chern data can also be checked directly:
\[
c(TS)=\frac{(1+H)^4}{1+4H}\bigg|_S=1+6H^2|_S,
\qquad \int_Sc_2(TS)=6\cdot4=24.
\]
Only \(c_1(TS)=0\), the surface dimension, and Mumford's relation were needed for the vanishing.

## 3. CJR grading, constant term, and the infinity aggregate

Fix \(g\), and let \(s_D:\overline{\mathcal M}_g(\mathbb P^3,D)\to\overline{\mathcal M}_g\) stabilize the source curve. With universal curve and map \(\pi,f\), define the zero-locus contribution
\[
Z_D(t)=s_{D,*}\left(
 e_t\bigl((R\pi_*f^*L)^\vee\bigr)
 \cap[\overline{\mathcal M}_g(\mathbb P^3,D)]^{\rm vir}\right).
\]
The dual twist has fibre weight \(t-4H\). The Euler symbol means the multiplicative Euler ratio of the derived complex, including its \(R^1\) term. This is the Serre-duality form of the stable zero-vertex contribution in [CJR-III, §9.3, p. 90].

The grading follows independently from Riemann–Roch:
\[
\operatorname{vdim}\overline{\mathcal M}_g(\mathbb P^3,D)=4D,
\qquad \chi(f^*L)=4D+1-g.
\]
Twisting leaves total equivariant dimension \(g-1\). As \(\dim\overline{\mathcal M}_g=3g-3\), the total codimension is **\(2g-2\)**, and
\[
Z_D(t)=\sum_{r=0}^{3g-3} t^{2g-2-r}Z_{D,r},
\qquad Z_{D,r}\in A^r(\overline{\mathcal M}_g).
\]
The expansion is in inverse \(t\) with nilpotent ordinary Chow classes. It is not evaluation of individual rational graph terms at \(t=0\).

After curve stabilization, write the reduced CJR identity as
\[
L_D(t)=Z_D(t)+I_D(t),
\qquad
I_D(t)=\sum_{\Gamma:\,V_\infty(\Gamma)\ne\varnothing}
 \frac{\operatorname{st}_{\Gamma,*}(\text{CJR graph class})}
 {|\operatorname{Aut}\Gamma|},
\qquad
B_{D,r}=[t^{2g-2-r}]I_D(t).
\]
Theorem 8.4 and Remark 8.5 imply that \(L_D(t)\) is polynomial in \(t\). Hence
\[
B_{D,r}=-Z_{D,r}\quad\text{whenever }r>2g-2.
\]
At the threshold, equation (9.7) instead gives
\[
[t^0]L_D(t)=(-1)^{1-g+4D}
 \sum_{\beta:\,i_*\beta=D[\mathrm{line}]}
 \operatorname{st}_{S,*}[\overline{\mathcal M}_g(S,\beta)]^{\rm vir}.
\]
For \(D>0\), this is zero by positive-degree ordinary K3 vanishing. For \(D=0\), ampleness permits only the zero effective class, and the obstruction calculation above again gives zero. This proves the boxed aggregate identity at the beginning for all \(D\ge0\). The sign is \((-1)^{1-g}\); there is no extra sign to attach to \(B=-Z\).

In particular, genus two has threshold \(r=2\), with negative powers starting at \(r=3\). Genus three has threshold \(r=4\), with negative powers at \(r=5,6\). For \(r<2g-2\), positive equivariant powers may remain, and ordinary K3 vanishing alone does not make their coefficients zero.

### Two meanings of “reduced”

The left side here uses the **reduced CJR/log-GLSM virtual class**. Its compact specialization is the **ordinary K3 GW virtual class**, by [CJR-GLSM, Proposition 4.1] and [CJR-III, (9.7)]. The **reduced K3 GW class** removes the symplectic obstruction quotient and has dimension \(g\) for nonzero degree. It is a different theory. Replacing ordinary K3 theory by that class would change the grading to codimension \(2g-3\) after curve stabilization and would require a new comparison/localization argument. One cannot substitute it into (9.7).

## 4. Auxiliary torus localization in ambient degree one

This is a second localization, applied to \(Z_D\) on \(\mathbb P^3\). Its edges are coordinate lines and are distinct from the vertical CJR edges joining zero and infinity vertices. A general smooth quartic need not be invariant under this auxiliary torus; only the ambient Euler-twisted calculation needs that action.

Use precisely
\[
H|_{p_i}=a_i,\qquad
w_{ij}=a_i-a_j,\qquad q_i=t-4a_i,
\qquad i,j\in\{0,1,2,3\}.
\]
Define
\[
\Lambda_h^\vee(u)=\sum_{\ell=0}^h(-1)^\ell\lambda_\ell u^{h-\ell},
\qquad
\Lambda_h^+(u)=\sum_{\ell=0}^h\lambda_\ell u^{h-\ell}.
\]
For a stable genus-\(h\), valence-one contracted vertex over \(p_i\), attached to the line \(ij\), the contribution on \(\overline{\mathcal M}_{h,1}\) is
\[
\boxed{
V_{h,i|j}=
\frac{\displaystyle\prod_{k\ne i}\Lambda_h^\vee(a_i-a_k)}
 {\Lambda_h^+(q_i)(w_{ij}-\psi)}.}
\]
There are **three** tangent Hodge factors. To see all cancellations, the contracted component initially contributes
\[
\frac{\prod_{k\ne i}\Lambda_h^\vee(a_i-a_k)}
 {\prod_{k\ne i}(a_i-a_k)}
\quad\text{and}\quad
\frac{q_i}{\Lambda_h^+(q_i)}.
\]
Normalization at its unique node supplies the tangent Euler class and the inverse twist fibre \(q_i^{-1}\). The remaining node smoothing contributes \((w_{ij}-\psi)^{-1}\). Thus there is no additional tangent Euler denominator or numerator \(q_i\) in the boxed vertex.

On a degree-one line, the five weights of the dual-twisted section space of \(\mathcal O(4)\) are
\[
t-(4-m)a_i-ma_j,\qquad 0\le m\le4.
\]
The moving section space of the tangent bundle has weights \(w_{ij},w_{ji}\), and \(a_i-a_k,a_j-a_k\) for \(k\ne i,j\). The zero-weight tangent direction cancels the corresponding infinitesimal automorphism. Since \(w_{ij}w_{ji}=-w_{ij}^2\), the edge factor is
\[
\boxed{
E_{ij}=-\frac{\displaystyle\prod_{m=0}^4
 \bigl(t-(4-m)a_i-ma_j\bigr)}
 {(a_i-a_j)^2\displaystyle\prod_{k\notin\{i,j\}}
 (a_i-a_k)(a_j-a_k)}=E_{ji}.}
\]
This derives the requested minus sign from tangent weights. It is unrelated to the compact-target sign in (9.7). The degree-one deck group is trivial. An unmarked unstable endpoint at \(p_j\) contributes \(w_{ji}\), from the remaining moving automorphism of the rational tail.

For clarity, the complete degree-one graph assembly can be written without choosing different conventions for equal endpoint genera. Let \(p:\overline{\mathcal M}_{g,1}\to\overline{\mathcal M}_g\) forget the marking and \(\xi_h:\overline{\mathcal M}_{h,1}\times\overline{\mathcal M}_{g-h,1}\to\overline{\mathcal M}_g\) glue. Then
\[
Z_1(t;a)=\sum_{i\ne j}E_{ij}\left[
 w_{ji}\,p_*V_{g,i|j}
 +\frac12\sum_{h=1}^{g-1}
 \xi_{h,*}(V_{h,i|j}\boxtimes V_{g-h,j|i})\right].
\]
The half removes the duplicate indexing \((i,j,h)\leftrightarrow(j,i,g-h)\). It is not an extra generic elliptic-stack factor. In genus three one may instead keep only \(h=1\), all ordered \(i\ne j\), and omit the half. In genus two the \(h=1\) ordered-pair sum retains the half.

Pushforward must precede numerical integration: \(p_*1=0\), \(p_*\psi^a=\kappa_{a-1}\) for \(a\ge1\); a tail forgetful pushforward lowers codimension by one, while a bridge gluing pushforward raises it by one. Take the auxiliary nonequivariant limit only after summing the fixed loci. Distinct auxiliary weights or cancelled poles are computational devices, not independent values of the CJR infinity vertices.

## 5. Degree-zero ambient class and independent genus-two/three checks

Let \(H^4=0\), \(\int_{\mathbb P^3}H^3=1\), and set
\[
F_g(u)=\frac{u}{\Lambda_g^+(u)}.
\]
Mapping to a point and the dual twist give, respectively,
\[
[\overline{\mathcal M}_g(\mathbb P^3,0)]^{\rm vir}
=c_{3g}(\mathbb E^\vee\boxtimes T\mathbb P^3)
 \cap[\overline{\mathcal M}_g\times\mathbb P^3],
\]
\[
(R\pi_*f^*\mathcal O(4))^\vee
=(\mathcal O-\mathbb E)\boxtimes\mathcal O(-4),
\qquad e_t=\frac{t-4H}{\Lambda_g^+(t-4H)}.
\]
Thus the required constant-map integrand and its pushforward are
\[
\boxed{
Z_0^{(g)}(t)=[H^3]\left\{
 c_{3g}(\mathbb E^\vee\boxtimes T\mathbb P^3)
 \frac{t-4H}{\Lambda_g^+(t-4H)}\right\}.}
\]
The rank-three Euler sequence gives \(c_1=4H,c_2=6H^2,c_3=4H^3\). Therefore an independent root expression is
\[
c_{3g}(\mathbb E^\vee\boxtimes T\mathbb P^3)
=\prod_{i=1}^g(-x_i^3+4Hx_i^2-6H^2x_i+4H^3)
\pmod{H^4}.
\]

### Genus two

Direct expansion gives
\[
\begin{aligned}
c_6={}&\lambda_2^3-4H\lambda_1\lambda_2^2
 +H^2(6\lambda_1^2\lambda_2+4\lambda_2^2)\\
 &-H^3(4\lambda_1^3+12\lambda_1\lambda_2).
\end{aligned}
\]
Using \(\lambda_1^2=2\lambda_2\) and \(\lambda_2^2=0\), this reduces to
\[
c_6=-20H^3\lambda_1\lambda_2.
\]
The Hodge coefficient already has top degree three on \(\overline{\mathcal M}_2\), so only the leading term \(F_2(t)=t^{-1}+\cdots\) contributes. Hence
\[
\boxed{Z_0^{(2)}(t)=-20\lambda_1\lambda_2\,t^{-1}.}
\]
In the grading above, \(Z_{0,2}=0\), \(Z_{0,3}=-20\lambda_1\lambda_2\), and the corresponding aggregates are \(B_{0,2}=0\), \(B_{0,3}=20\lambda_1\lambda_2\).

### Genus three

Direct expansion gives
\[
\begin{aligned}
c_9={}&-\lambda_3^3+4H\lambda_2\lambda_3^2
 -H^2(4\lambda_1\lambda_3^2+6\lambda_2^2\lambda_3)\\
 &+H^3(12\lambda_1\lambda_2\lambda_3+4\lambda_2^3+4\lambda_3^2).
\end{aligned}
\]
Now use \(\lambda_1^2=2\lambda_2\), \(\lambda_2^2=2\lambda_1\lambda_3\), and \(\lambda_3^2=0\). This yields
\[
c_9=20H^3\lambda_1\lambda_2\lambda_3,
\qquad
\boxed{Z_0^{(3)}(t)=20\lambda_1\lambda_2\lambda_3\,t^{-2}.}
\]
Consequently
\[
Z_{0,4}=Z_{0,5}=0,\quad Z_{0,6}=20\lambda_1\lambda_2\lambda_3,
\qquad
B_{0,4}=B_{0,5}=0,\quad B_{0,6}=-20\lambda_1\lambda_2\lambda_3.
\]
These nonzero ambient twisted classes coexist with a zero ordinary K3 virtual class because infinity contributions cancel them.

For numerical checks, [FP, Theorem 4 and §4.3] gives
\[
\int_{\overline{\mathcal M}_g}\lambda_g\lambda_{g-1}\lambda_{g-2}
=\frac{|B_{2g}|\,|B_{2g-2}|}
 {2(2g)(2g-2)(2g-2)!}.
\]
Thus \(\int Z_{0,3}=-1/288\) in genus two and \(\int Z_{0,6}=1/72576\) in genus three. These are integrals of the indicated Laurent coefficients, not ordinary K3 GW invariants.

The expansions, symmetric-polynomial conversions, Mumford-ideal reductions, and exact rational numbers were checked independently with SymPy in the pinned Sage environment on OrbStack `agents`, using Bash and stdin scripts. No parent computational code was used or tested.

## 6. Quartic effective contacts and limits of reconstruction

For a connected effective punctured component with genus \(g\), ambient degree \(D\), and negative integral contacts \(c_i\le-1\), specialize [CJR-P, (9.5), p. 88] to \(N=\mathcal O(4)\):
\[
\boxed{\sum_i(c_i+1)=4D-(2g-2)\le0.}
\]
In particular \(4D\le2g-2\). Equality forces every contact to be \(-1\) (including the vacuous unmarked case). With every contact equal to \(-2\), their number is \(n=2g-2-4D\). In genus two only \(D=0\) is possible for such effective components; in genus three \(D=0,1\) are possible by balance, and \(D=1\) forces contacts \(-1\). These are necessary conditions, not existence or nonvanishing results. Inside a CJR graph, apply them with each component's own genus and degree, not automatically the total graph data.

The punctured convention uses negative \(c_i\). In [CJR-III, §9.3, p. 90], the vertical edge multiplicities are positive and satisfy
\[
\sum_{E\text{ at }V}c(E)=2g(V)-2+\operatorname{val}(V)-4D(V).
\]
For the corresponding punctures \(c_i=-c(E)\); rearranging gives exactly the boxed effective balance. Auxiliary projective-space edge degree one in §4 is a different decoration.

Specializing the reduced punctured dimension formula [CJR-P, (9.8)] gives
\[
\operatorname{vdim}_{\rm CJR,red}
=(g-1)+n+\sum_i(c_i+1)=1-g+4D+n.
\]
For all \(-2\) contacts it is \(g-1\), not zero. Also the expression required to be nonpositive in Assumption 9.22 is here
\[
(3-\dim\mathbb P^3+\operatorname{rk}L)(g-1)
 -\int_Dc_1(K_{\mathbb P^3}\otimes L)=g-1>0.
\]
Therefore the basic-scalar reconstruction hypotheses of Definition 9.19/Corollary 9.24 do **not** apply to the quartic in genus at least two. The ordinary compact K3 vanishing is not a negative-virtual-dimension hypothesis and cannot be inserted into Corollary 9.3 to conclude that each effective cycle vanishes.

Finally, \(B_{D,r}\) denotes the sum of **all** infinity-containing CJR contributions, including their gluing and descendant factors. An identity for this sum neither isolates each effective cycle nor proves its individual tautologicality. Infinity-containing also does not imply support on the boundary of \(\overline{\mathcal M}_g\), because rational components may disappear under stabilization. The explicit ambient degree-zero and degree-one expressions give tautological presentations of \(Z\), hence of the aggregates in the stated range; they do not resolve the general effective-cycle question in [CJR-P, §9.8].
