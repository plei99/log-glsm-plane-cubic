This calculation extends the projective and toric searches to complete intersections in the non-toric homogeneous spaces \(\operatorname{Gr}(2,4)\) and \(\operatorname{Gr}(2,5)\). The sample has fourteen targets. It computes connected, unmarked genera two and three, in Plücker degrees \(D=0,1,2\), at two auxiliary torus specializations. There are 42 target/degree pairs and 168 target/degree/genus/specialization evaluations. The ambient dimension is now allowed to be six: this is an explicit enlargement of the preceding toric dimension-four bound.

The scan retained **648 independent per-target relation vectors**. Their union has rank **206: 19 in genus two and 187 in genus three**. Every retained vector has an exact Pixton certificate. The 206-dimensional union is **equal to the previous bounded toric span**, not merely equal in dimension, and is contained in the previous 463-dimensional projective span. Explicit certificates verify both comparisons. Thus this bounded Grassmannian expansion adds **zero independent tautological directions** to either preceding span, while providing new ambient class formulas and CJR aggregate comparisons.

The implementation is [grassmannian_complete_intersections.sage](../grassmannian_complete_intersections.sage). It supplies the Grassmannian fixed-point, invariant-line, normal-bundle, and Plücker-character data to the existing GKM localization and curve-stabilization backend. It also computes ordinary Chern numbers by a separate Chern-root calculation. The [source audit](GRASSMANNIAN_COMPLETE_INTERSECTION_SOURCES.md) records the geometric references and the scope of the unpublished CJR draft. All saved relation vectors have exact rational coefficients and explicit certificates; no numerical rank tolerance is used.

**Ambient geometry.** A point of \(\operatorname{Gr}(2,n)\) is a two-dimensional vector subspace of \(\mathbb C^n\), equivalently a projective line in \(\mathbb P^{n-1}\). The tautological sequence is
\[
0\longrightarrow S\longrightarrow\mathcal O^{\oplus n}\longrightarrow Q\longrightarrow0,
\qquad T_X=S^\vee\otimes Q,
\qquad H=c_1(\det S^\vee).
\]
Thus \(\dim X=2n-4\), \(\operatorname{Pic}(X)=\mathbb ZH\), and \(-K_X=nH\). The Plücker embedding sends a basis \(u,v\) to \([u\wedge v]\). In coordinates its image is defined by
\[
p_{ij}p_{kl}-p_{ik}p_{jl}+p_{il}p_{jk}=0
\quad(i<j<k<l).
\]
For \(n=4\) this is the smooth quadric fourfold \(Q^4\subset\mathbb P^5\), with equation
\(p_{12}p_{34}-p_{13}p_{24}+p_{14}p_{23}=0\).
For \(n=5\), it is the smooth sixfold in \(\mathbb P^9\) cut out by the five quadratic Plücker relations, equivalently the five maximal Pfaffians of a generic skew-symmetric \(5\times5\) matrix of coordinates. Its degree is five.

These two Grassmannians are not toric varieties, and there is no toric fan describing them. The diagonal torus has effective rank \(n-1\), smaller than \(2n-4\). More intrinsically, a smooth projective toric variety of Picard rank one is projective space, whereas these Grassmannians are not projective spaces. The replacement for the fan in localization is the GKM graph described in [Liu–Sheshmani, Examples 2.3 and 2.5](https://sigma-journal.com/2017/048/sigma17-048.pdf).

**The complete intersections.** Let
\[
Y=Z(s)\subset X=\operatorname{Gr}(2,n),\qquad
E=\bigoplus_{a=1}^{c}\mathcal O_X(d_a),\quad d_a>0.
\]
The notation \(X[d_1,\ldots,d_c]\) below means a general smooth zero locus. Adjunction gives
\[
\dim Y=2n-4-c,\qquad -K_Y=(n-\sum_a d_a)H|_Y.
\]
All summands are ample and globally generated. The CY and K3 labels use the cohomological vanishings from the Koszul complex and Kodaira vanishing, in addition to trivial canonical bundle. For these two ambients, the sample exhausts positive Plücker degree lists, up to ordering, with target dimension two or three and \(\sum d_a\le n\). It does not exhaust Grassmannians or homogeneous vector-bundle zero loci.

| ID | Target | Geometric description |
|---|---|---|
| G1 | \(\operatorname{Gr}(2,4)[1]\) | Smooth quadric threefold \(Q^3\subset\mathbb P^4\), index 3 |
| G2 | \(\operatorname{Gr}(2,4)[2]\) | Intersection of two quadrics in \(\mathbb P^5\), index 2 and degree 4 |
| G3 | \(\operatorname{Gr}(2,4)[3]\) | Intersection \((2,3)\subset\mathbb P^5\), prime Fano threefold of genus 4 |
| G4 | \(\operatorname{Gr}(2,4)[4]\) | CY3 \((2,4)\subset\mathbb P^5\) |
| G5 | \(\operatorname{Gr}(2,4)[1,1]\) | Quadric surface \(\mathbb P^1\times\mathbb P^1\), degree-8 del Pezzo |
| G6 | \(\operatorname{Gr}(2,4)[1,2]\) | Degree-4 del Pezzo surface \((2,2)\subset\mathbb P^4\) |
| G7 | \(\operatorname{Gr}(2,4)[1,3]\) | Degree-6, genus-4 polarized K3 \((2,3)\subset\mathbb P^4\) |
| G8 | \(\operatorname{Gr}(2,4)[2,2]\) | Degree-8, genus-5 polarized K3 \((2,2,2)\subset\mathbb P^5\) |
| G9 | \(\operatorname{Gr}(2,5)[1,1,1]\) | Fano threefold \(V_5\subset\mathbb P^6\), index 2 and degree 5 |
| G10 | \(\operatorname{Gr}(2,5)[1,1,2]\) | Ordinary Gushel–Mukai threefold, degree 10 and genus 6 |
| G11 | \(\operatorname{Gr}(2,5)[1,1,3]\) | CY3 of degree 15 and Euler characteristic \(-150\) |
| G12 | \(\operatorname{Gr}(2,5)[1,2,2]\) | CY3 of degree 20 and Euler characteristic \(-120\) |
| G13 | \(\operatorname{Gr}(2,5)[1,1,1,1]\) | Degree-5 del Pezzo surface |
| G14 | \(\operatorname{Gr}(2,5)[1,1,1,2]\) | Degree-10, genus-6 polarized K3 surface |

The genus of a polarized K3 is \(1+H^2/2\); it is unrelated to the genus \(g=2,3\) of maps in the computation. G11 and G12 have Hodge numbers \((h^{1,1},h^{2,1})=(1,76),(1,61)\). Their Chern data agree with [Haghighat–Klemm, Appendix A](https://arxiv.org/abs/0802.2908).

If \(x,y\) are the Chern roots of \(S^\vee\), the independent Chern calculation uses
\[
c(T_X)=\frac{(1+x)^n(1+y)^n}{1-(x-y)^2},\qquad H=x+y,
\]
truncated to the required degree, and
\[
\int_{\operatorname{Gr}(2,n)}f(x,y)
=[x^{n-1}y^{n-2}]\,(x-y)f(x,y)
\]
for symmetric \(f\) of degree \(2n-4\). Hence
\[
\int_Y c_I(TY)=\int_X
\left(\prod_a d_aH\right)c_I\!\left(\frac{c(TX)}{\prod_a(1+d_aH)}\right).
\]
Every Gr(2,4) result is also checked by the corresponding projective complete-intersection formula. Ambient top Chern integrals are independently checked by fixed-point integration. These checks compare different constructions rather than reproducing the same polynomial routine.

**Fixed-point and edge data.** Write \(\lambda_1,\ldots,\lambda_n\) for coordinate torus characters. Fixed points are the coordinate two-planes \(p_I\), with \(|I|=2\), and
\[
T_{p_I}X=\bigoplus_{i\in I,j\notin I}\mathbb C_{\lambda_j-\lambda_i},
\qquad \mathcal O_X(d)|_{p_I}\text{ has weight }-d\sum_{i\in I}\lambda_i.
\]
Two fixed points are adjacent precisely when \(J=(I\setminus\{a\})\cup\{b\}\). Their invariant line consists of the planes \(V\) between the fixed one-space \(I\cap J\) and fixed three-space \(I\cup J\). At \(p_I\), its tangent weight is \(w=\lambda_b-\lambda_a\). This is the Johnson graph \(J(n,2)\): Gr(2,4) has six vertices and twelve edges; Gr(2,5) has ten vertices and thirty edges.

Along such a line,
\[
T_X|_\ell=\mathcal O(2)\oplus\mathcal O(1)^{\oplus(n-2)}
\oplus\mathcal O^{\oplus(n-3)}.
\]
This follows by splitting \(S|_\ell=\mathcal O(-1)\oplus\mathcal O\) and \(Q|_\ell=\mathcal O(1)\oplus\mathcal O^{n-3}\). It distinguishes the Grassmannian edge factor from a projective-space edge factor: the trivial normal directions must be retained. Their initial weights, and those of the degree-one directions, are stored explicitly. A normal line with degree \(a\) and initial character \(\mu\) contributes
\[
\prod_{r=0}^{ae}(\mu-rw/e)^{-1}
\]
on a degree-\(e\) cover; here \(a=0,1\). The moving tangent factor, twist, and separate deck factor are respectively
\[
\frac{(-1)^ee^{2e}}{(e!)^2w^{2e}},\qquad
\prod_{a=1}^{c}\prod_{r=0}^{ed_a}(t-L_{a,I}+rw/e),\qquad \frac1e.
\]
Vertices, flags, Hodge polynomials, graph automorphisms, and forgetful stabilization are handled by the already checked backend. The general GKM virtual-localization reference is [Liu–Sheshmani, Lemma 4.5 and Theorem 4.7](https://sigma-journal.com/2017/048/sigma17-048.pdf).

**What the CJR coefficients mean.** Keep the logarithmic parameter \(t\) separate from the auxiliary Grassmannian torus. For \(m=\dim Y\), \(\iota=n-\sum d_a\), define
\[
Z_D(t)=\operatorname{st}_*\left(e_t((R\pi_*f^*E)^\vee)
[\overline{\mathcal M}_g(X,D)]^{\rm vir}\right),\qquad
Z_D(t)=\sum_r t^{\alpha-r}Z_{D,r},
\]
\[
\alpha=m(g-1)-\iota D,\qquad
\sigma=(-1)^{c(1-g)+D\sum_a d_a}.
\]
Let \(G_{g,D}\) be the stabilized ordinary target class, summed over target curve classes with Plücker degree \(D\), and \(B_D\) the aggregate of all CJR terms involving infinity data. Then
\[
B_{D,r}=-Z_{D,r}\quad(r>\alpha),\qquad
B_{D,\alpha}=\sigma G_{g,D}-Z_{D,\alpha}.
\]
The second identity is used when \(\alpha\) lies in the Chow range. These are aggregate identities, not reconstructions of individual punctured vertices. For a higher-rank bundle the infinity target is \(\mathbb P_X(E^\vee)\), not a product of rank-one theories. The CJR construction allows a smooth projective ambient space; toricity is used only as one way to evaluate the ambient contribution. The supplied draft and its exact statements are documented in the source audit; the public compact comparison is [Chen–Janda–Ruan, Proposition 4.1](https://arxiv.org/abs/1906.04345v2).

For K3 targets G7, G8, G14, ordinary virtual classes vanish in the tested genera, including the separately checked constant classes. Thus in genus three \(\alpha=4\) and \(B_{D,r}=-Z_{D,r}\) for \(r=4,5,6\). These are not reduced K3 invariants; see [Maulik–Pandharipande–Thomas, §0.1](https://people.math.ethz.ch/~rahul/mpt.pdf).

For CY3 targets G4, G11, G12, \(\alpha=3g-3\), so there is no negative-\(t\) range in the unmarked pushforward. The top coefficient retains ordinary GW theory. The ambient finite coefficients are still useful geometric expressions, and auxiliary equivariance still supplies relations.

There is a stronger geometric simplification for degree one in Gr(2,5):
\[
\operatorname{vdim}\overline{\mathcal M}_3(\operatorname{Gr}(2,5),1)
=(1-3)(6-3)+5=-1.
\]
Its nonequivariant virtual class vanishes, hence the entire \(Z_1(t)\) vanishes, regardless of the defining bundle. This is checked in every codimension for all six Gr(2,5) targets. The same argument proves ambient constant-map vanishing for both ambients and genera. In total, the separate geometric verification checks 196 coefficient slots with negative ambient virtual dimension.

The same variety can occur with two ambient presentations, for example G4 as Gr(2,4)[4] and as (2,4) in P5. Its ordinary invariants agree. Its higher-genus ambient dual-Euler term \(Z\) need not agree: higher-genus quantum Lefschetz does not identify these two ambient expressions without correction. This is why the regression compares geometry and genus-zero ordinary values, rather than imposing equality of their higher-genus \(Z\) terms.

**How relations were certified.** Substitute \(\lambda_i=\epsilon b^{i-1}\) with \(b=4,5\). All relevant tangent, normal, and node denominators are checked for resonance at the degree-two bound. After summing fixed graphs, properness gives a regular auxiliary nonequivariant limit. Negative \(\epsilon\)-coefficients are therefore zero in Chow, and the finite coefficients agree in Chow between the two choices of \(b\). Their strata vectors need not be literally equal before imposing relations.

For every retained vector \(v\), the files contain rational coefficients \(u\) such that \(uM=v\), where \(M\) is the Pixton matrix in the same ordered strata basis. The export replays each identity. It also computes the union with the previous projective relation span; when contained, it saves explicit coefficients expressing each new union vector in that previous basis. The Chow validity of the Pixton system follows from [Janda, Theorem 1.1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf). This computation does not assume that every tautological relation is Pixton.

**Independent ordinary CY3 calculations.** Lines on \(\operatorname{Gr}(2,n)\) are parameterized by \(F(1,3;n)\). Over a flag \(A\subset B\), the line is \(\mathbb P(B/A)\), and the restriction of the defining degree-\(d\) equation is a section of
\[
\mathcal V_d=(A^\vee)^{\otimes d}\otimes\operatorname{Sym}^d(B/A)^\vee.
\]
The sum of these bundles has rank \(\sum_a(d_a+1)=3n-7=\dim F(1,3;n)\) for our CY3s. Restrictions of ambient sections generate it, so a general section has a transverse finite zero scheme. The integral
\[
L=\int_{F(1,3;n)}c_{\rm top}\left(\bigoplus_a\mathcal V_{d_a}\right)
\]
is evaluated by flag-variety localization, independently of stable-map graphs. It gives \(1280,540,400\) for G4, G11, G12. Transversality gives \(H^0(N_{\ell/Y})=0\); adjunction gives \(\deg N_{\ell/Y}=-2\); hence every counted line has normal bundle \(\mathcal O(-1)^{\oplus2}\). Every degree-one nonconstant map has one of these lines as image. The local-curve formula of [Faber–Pandharipande, Theorem 3](https://arxiv.org/abs/math/9810173) therefore determines the ordinary connected invariants
\[
N_{2,1}=L/240,\qquad N_{3,1}=L/6048.
\]
This includes contracted higher-genus components; it does not assert that there are embedded positive-genus degree-one curves.

In genus three this gives especially simple integrated CJR comparisons:
\[
\begin{array}{c|c|c|c|c}
Y&N_{3,1}&\int Z_{1,6}&\sigma&\int B_{1,6}\\\hline
G4&40/189&-163/6720&+1&14267/60480\\
G11&5/56&0&-1&-5/56\\
G12&25/378&0&-1&-25/378.
\end{array}
\]
For G11 and G12, all of the compact comparison is therefore supplied by the infinity aggregate. For comparison, the K3 G14 in degree two has
\[
\int Z_{2,6}=-8673817/193536,\qquad
\int B_{2,6}=8673817/193536,
\]
using ordinary K3 vanishing. These equalities concern integrated full aggregates; they are not values assigned to individual infinity vertices.

Genus-zero stable-map localization with convex quantum Lefschetz independently gives
\[
\begin{array}{c|r|r|r}
Y&N_{0,1}&N_{0,2}&N_{0,2}-N_{0,1}/8\\\hline
G4&1280&92448&92288\\
G11&540&25245/2&12555\\
G12&400&5590&5540.
\end{array}
\]
All six values agree at both auxiliary specializations. The dual-Euler localization value is converted to the ordinary convention by the comparison sign. Convexity applies because the ambient Grassmannian is homogeneous and every defining degree is positive; the relevant virtual comparison is [Chen–Janda–Webb, Theorem 1.2.3](https://arxiv.org/abs/1911.09825v2). The last column for G11 and G12 agrees with the genus-zero entries in [Haghighat–Klemm, Tables B.1 and B.2](https://arxiv.org/abs/0802.2908). Their higher-genus mirror predictions are not imported as proven ordinary values. In particular, degree-two higher-genus ordinary invariants remain unknown in this module; subtracting genus-zero multiple covers alone does not establish the needed higher-genus geometry.

At degree zero, for a threefold put \(C=\int_Y(c_3-c_1c_2)\). The compact classes are
\[
G_{2,0}=C\lambda_1\lambda_2,\qquad
G_{3,0}=-C\lambda_1\lambda_2\lambda_3.
\]
For the CY3s \(C=\chi(Y)\), so
\[
N_{2,0}=\chi(Y)/5760,\qquad N_{3,0}=-\chi(Y)/1451520.
\]
For a surface, the compact classes are \(c_1^2\lambda_2\) and \(c_1^2\lambda_1\lambda_3\). The ambient constant terms for our ambient dimensions four and six vanish after stabilization in the two tested genera; this is checked against every saved finite coefficient modulo Pixton, rather than conflated with the generally nonzero compact target constants.

| Target | Ambient | Bundle coefficients | Type | Euler characteristic |
|---|---|---|---|---:|
| G1 | gr2-4 | ((1,),) | Fano | 4 |
| G2 | gr2-4 | ((2,),) | Fano | 0 |
| G3 | gr2-4 | ((3,),) | Fano | -36 |
| G4 | gr2-4 | ((4,),) | CY3 | -176 |
| G5 | gr2-4 | ((1,), (1,)) | Fano | 4 |
| G6 | gr2-4 | ((1,), (2,)) | Fano | 8 |
| G7 | gr2-4 | ((1,), (3,)) | K3 | 24 |
| G8 | gr2-4 | ((2,), (2,)) | K3 | 24 |
| G9 | gr2-5 | ((1,), (1,), (1,)) | Fano | 4 |
| G10 | gr2-5 | ((1,), (1,), (2,)) | Fano | -16 |
| G11 | gr2-5 | ((1,), (1,), (3,)) | CY3 | -150 |
| G12 | gr2-5 | ((1,), (2,), (2,)) | CY3 | -120 |
| G13 | gr2-5 | ((1,), (1,), (1,), (1,)) | Fano | 7 |
| G14 | gr2-5 | ((1,), (1,), (1,), (2,)) | K3 | 24 |

Per-target independent localization relations (two auxiliary torus specializations):

| Target | Genus 2 ranks r=1,2,3 | Genus 3 ranks r=2,3,4,5,6 |
|---|---|---|
| G1 | 1, 5, 10 | 3, 3, 6, 6, 10 |
| G2 | 1, 6, 9 | 3, 3, 6, 6, 10 |
| G3 | 1, 6, 9 | 3, 3, 6, 6, 10 |
| G4 | 1, 6, 10 | 3, 3, 6, 6, 10 |
| G5 | 1, 6, 10 | 3, 3, 6, 6, 10 |
| G6 | 1, 6, 12 | 3, 5, 8, 10, 14 |
| G7 | 1, 6, 12 | 3, 5, 8, 10, 14 |
| G8 | 1, 6, 9 | 3, 3, 6, 6, 10 |
| G9 | 1, 6, 12 | 1, 3, 5, 7, 9 |
| G10 | 1, 6, 12 | 1, 3, 5, 7, 9 |
| G11 | 1, 6, 12 | 1, 3, 5, 7, 9 |
| G12 | 1, 6, 12 | 1, 3, 5, 7, 9 |
| G13 | 1, 6, 12 | 1, 3, 5, 7, 9 |
| G14 | 1, 6, 12 | 1, 3, 5, 7, 9 |

| Genus | Codimension | Grassmannian rank | Previous projective rank | Combined rank | Pixton rank |
|---|---|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 6 | 6 | 6 | 6 |
| 2 | 3 | 12 | 16 | 16 | 16 |
| 3 | 2 | 6 | 6 | 6 | 6 |
| 3 | 3 | 16 | 30 | 30 | 30 |
| 3 | 4 | 39 | 86 | 86 | 114 |
| 3 | 5 | 59 | 150 | 150 | 259 |
| 3 | 6 | 67 | 168 | 168 | 377 |

Exact twisted top integrals for G1 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 133/720 | -17/1792 |
| (2) | 7/144 | -48257/72576 |

Exact twisted top integrals for G2 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 91/288 | -661/30240 |
| (2) | 0 | -7451/7560 |

Exact twisted top integrals for G3 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 0 | -6791/241920 |
| (2) | -3/8 | -7/384 |

Exact twisted top integrals for G4 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | -22/45 | -163/6720 |
| (2) | -887/4 | -568469/362880 |

Exact twisted top integrals for G5 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 83/120 | -23389/725760 |
| (2) | 1/360 | -60611/8640 |

Exact twisted top integrals for G6 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 49/72 | -9157/181440 |
| (2) | -53/720 | -1173841/362880 |

Exact twisted top integrals for G7 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | -181/1440 | -41653/725760 |
| (2) | 2717/720 | -28099/40320 |

Exact twisted top integrals for G8 (gr2-4):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 11/30 | -17039/241920 |
| (2) | -73/30 | -326/405 |

Exact twisted top integrals for G9 (gr2-5):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 1285/384 | 0 |
| (2) | -11/128 | -25633969/193536 |

Exact twisted top integrals for G10 (gr2-5):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 295/144 | 0 |
| (2) | -217/96 | -616169/32256 |

Exact twisted top integrals for G11 (gr2-5):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | -25/32 | 0 |
| (2) | -1834/9 | -1062577/96768 |

Exact twisted top integrals for G12 (gr2-5):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 955/1152 | 0 |
| (2) | -32597/1152 | 10891/145152 |

Exact twisted top integrals for G13 (gr2-5):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 95/18 | 0 |
| (2) | -25/288 | -21215507/72576 |

Exact twisted top integrals for G14 (gr2-5):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0) | 0 | 0 |
| (1) | 3875/1152 | 0 |
| (2) | 14509/1152 | -8673817/193536 |

Ordinary genus-zero CY3 invariants, checked at both rays:

| Target | Mori tuple | N₀,β |
|---|---|---|
| G4 | (1,) | 1280 |
| G4 | (2,) | 92448 |
| G11 | (1,) | 540 |
| G11 | (2,) | 25245/2 |
| G12 | (1,) | 400 |
| G12 | (2,) | 5590 |

Known ordinary higher-genus CY3 numbers and integrated infinity aggregates:

| Target | Genus | Mori tuple | N_g,β | ∫B_β,3g−3 |
|---|---|---|---|---|
| G4 | 2 | (0,) | -11/360 | 11/360 |
| G4 | 3 | (0,) | 11/90720 | 11/90720 |
| G4 | 2 | (1,) | 16/3 | -218/45 |
| G4 | 3 | (1,) | 40/189 | 14267/60480 |
| G11 | 2 | (0,) | -5/192 | 5/192 |
| G11 | 3 | (0,) | 5/48384 | 5/48384 |
| G11 | 2 | (1,) | 9/4 | 97/32 |
| G11 | 3 | (1,) | 5/56 | -5/56 |
| G12 | 2 | (0,) | -1/48 | 1/48 |
| G12 | 3 | (0,) | 1/12096 | 1/12096 |
| G12 | 2 | (1,) | 5/3 | 965/1152 |
| G12 | 3 | (1,) | 25/378 | -25/378 |


**Reproduction and saved data.** Run from the project directory in the Linux Sage environment:

```bash
python -c 'from sage.all import *; load("tests/test_grassmannian_complete_intersections.sage")'
python calculations/grassmannian_ci_scan.py --workers 4
python calculations/grassmannian_ci_export.py
python calculations/grassmannian_ci_verify.py
```

For an individual calculation:

```python
load('grassmannian_complete_intersections.sage')
Y = GrassmannianCompleteIntersection(5, (1, 1, 3))
C = GrassmannianCompleteIntersectionCycles(Y, genus=3)
Z = C.cycle((1,), 6)
comparison = C.aggregate_relation((1,), 6)
ordinary_number = grassmannian_known_ordinary_gw(Y, 3, (1,))
```

The compact *number* in the last line is known, while the general compact *class* API does not silently reconstruct an unknown target class from an ambient expression. All-codimension formal terms and stabilized vectors are saved separately. The directory [results/grassmannian-complete-intersections/2026-11-06](../results/grassmannian-complete-intersections/2026-11-06) contains per-target coefficients, Pixton certificates, union and previous-span certificates, ordinary values, exact tables, logs, and source hashes. `independent-relations.txt` gives a readable basis and the associated strata dictionaries. `verification.json` records the checks, including fourteen genus-three degree-one direct Hodge integrations and six additional degree-two integrations for the CY3s.

All computations are bounded by genus, target list, and Plücker degree. The module does not yet accept arbitrary Gr(k,n), arbitrary homogeneous vector bundles, marked insertions, or degrees above two. It does not infer the individual punctured infinity invariants from their aggregate, and it preserves unknown ordinary GW data explicitly.
