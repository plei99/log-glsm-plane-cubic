The bounded toric expansion is complete. It treats twelve targets in five smooth toric Fano ambient spaces, connected unmarked genera two and three, and nonnegative Mori tuples whose coefficients sum to at most two. This gives 88 target/class pairs and 352 target/class/genus/torus-ray evaluations. The [new module](../toric_complete_intersections.sage) computes the stabilized dual-Euler classes, retains unknown ordinary GW classes, and exports exact rational certificates.

The scan retained 1306 independent per-target localization vectors, and checked 38 nonzero factor-exchange expressions. After removing repetitions, their union has rank **206: 19 in genus two and 187 in genus three**. Every row has an exact Pixton certificate. All 206 rows also lie in the previously computed 463-dimensional projective-family span, with explicit containment certificates. Thus this bounded expansion adds **zero new independent tautological directions** to that previous search. It does provide new toric class expressions, infinity-aggregate comparisons, and ordinary GW values.

The sample contains Fano surfaces and threefolds, three K3 surfaces, and three CY3s. Its ambient spaces are \(\mathbb P^1\times\mathbb P^2\), \(\mathbb P^1\times\mathbb P^3\), \(\mathbb P^2\times\mathbb P^2\), \(\mathbb F_1\times\mathbb P^1\), and \(\mathbb F_1\times\mathbb P^2\). All have dimension three or four. This is a selected sample; it does not classify all smooth toric Fanos or all their complete intersections. The code uses explicit toric fixed-point and invariant-curve data for these product constructions, rather than claiming to accept an arbitrary fan.

Write \(A,B\) for projective hyperplane classes. On \(\mathbb F_1\), let \(S\) be the negative section and \(F\) a fiber; let \(H\) be the hyperplane class of the other factor. The conventions are
\[
S^2=-SF,\quad F^2=0,\quad\int_{\mathbb F_1}SF=1,
\qquad -K_{\mathbb F_1}=2S+3F.
\]
The twelve examples are:

| ID | Ambient | Defining divisor(s) | Target |
|---|---|---|---|
| T1 | \(\mathbb P^1\times\mathbb P^2\) | \(A+B\) | Fano surface, degree 8 del Pezzo |
| T2 | \(\mathbb P^1\times\mathbb P^2\) | \(2A+3B\) | K3 |
| T3 | \(\mathbb P^1\times\mathbb P^3\) | \(A+2B\) | Fano threefold |
| T4 | \(\mathbb P^1\times\mathbb P^3\) | \(2A+4B\) | CY3 |
| T5 | \(\mathbb P^2\times\mathbb P^2\) | \(A+2B\) | Fano threefold |
| T6 | \(\mathbb P^2\times\mathbb P^2\) | \(3A+3B\) | CY3 |
| T7 | \(\mathbb P^2\times\mathbb P^2\) | \(A+B,\ A+B\) | Fano surface, degree 6 del Pezzo |
| T8 | \(\mathbb P^2\times\mathbb P^2\) | \(A+B,\ 2A+2B\) | K3 |
| T9 | \(\mathbb F_1\times\mathbb P^1\) | \(S+F+H\) | Fano surface, degree 7 del Pezzo |
| T10 | \(\mathbb F_1\times\mathbb P^1\) | \(2S+3F+2H\) | K3 |
| T11 | \(\mathbb F_1\times\mathbb P^2\) | \(S+F+H\) | Fano threefold |
| T12 | \(\mathbb F_1\times\mathbb P^2\) | \(2S+3F+3H\) | CY3 |

All are general smooth zero loci of globally generated bundles. Adjunction proves the Fano assertions: the residual anticanonical divisors are \(A+2B\) for T1 and T3, \(2A+B\) for T5, \(A+B\) for T7, \(S+2F+H\) for T9, and \(S+2F+2H\) for T11. They are ample. The remaining rows have trivial canonical bundle; the hypersurface exact sequence or the rank-two Koszul complex proves the needed \(H^i(\mathcal O_Y)\) vanishings. The [source audit](TORIC_COMPLETE_INTERSECTION_SOURCES.md) gives these arguments and the Chern calculations independently of localization.

The defining bundle \(S+F+H\) in T9 and T11 is nef and globally generated, but not ample: it has degree zero on \(S\times\{p\}\). The general CJR comparison and localization formulas still apply. The ampleness-based finiteness assertion in CJR III, Corollary 10.2, is not invoked for these two examples. Finiteness of this computation comes from the explicit Mori bound. The audited CJR input is the supplied October 17, 2025 draft, Theorem 8.4, Remark 8.5, and §9.4; the public compact comparison is [Chen–Janda–Ruan, Proposition 4.1](https://arxiv.org/abs/1906.04345v2).

The curve bound is vector-valued. For a product of projective spaces, \(\beta=(d_1,d_2)\) records degrees in the two line classes. For an \(\mathbb F_1\) product,
\[
\beta=d_SS+d_FF+d_H\ell,
\qquad (aS+bF+cH)\cdot\beta=(b-a)d_S+a d_F+c d_H.
\]
The permitted tuples are all nonnegative tuples with sum at most two: six tuples in rank two and ten in rank three, including zero. The invariant positive section has class \(S+F\). Thus a one-edge map to that section and a two-edge map to \(S\cup F\) both contribute to \((1,1,0)\). A double cover of the positive section has tuple \((2,2,0)\) and is outside the bound. The implementation groups graphs by
\[
\beta(\Gamma)=\sum_e d_e[C_e],
\qquad g=\sum_vg_v+b_1(\Gamma),
\]
so these cases are not confused with the sum of cover degrees. Fixed-point adjacency also comes from the toric invariant curves; arbitrary pairs of fixed points need not be adjacent.

For \(Y=Z(s)\subset X\), \(s\in H^0(X,E)\), put \(m=\dim Y\) and \(c=\operatorname{rk}E\). The ambient term and grading are
\[
Z_\beta(t)=\operatorname{st}_{X,*}
\left(e_t((R\pi_*f^*E)^\vee)
[\overline{\mathcal M}_g(X,\beta)]^{\mathrm{vir}}\right),
\]
\[
\alpha=m(g-1)-c_1(TY)\cdot\beta,
\qquad Z_\beta(t)=\sum_{r=0}^{3g-3}t^{\alpha-r}Z_{\beta,r},
\qquad \sigma=(-1)^{c(1-g)+c_1(E)\cdot\beta}.
\]
Let \(B_\beta\) be the pushforward of the full sum of CJR localization terms with infinity data, and let \(G_{g,\beta}\) be the sum of ordinary target virtual classes with ambient pushforward \(\beta\), stabilized to \(\overline{\mathcal M}_g\). Then
\[
\boxed{B_{\beta,r}=-Z_{\beta,r}\quad(r>\alpha)},
\qquad
\boxed{B_{\beta,\alpha}=\sigma G_{g,\beta}-Z_{\beta,\alpha}}.
\]
The second formula applies when \(\alpha\) is in the Chow range. These are aggregate identities; they do not isolate individual punctured effective vertices. For rank two, the infinity target remains \(\mathbb P_X(E^\vee)\). Only the zero-vertex Euler twist factors over the summands of \(E\), all with the same parameter \(t\).

In genus three, the Fano values of \(\alpha\) are
\[
\begin{array}{c|l}
T1&4-d_1-2d_2\\
T3&6-d_1-2d_2\\
T5&6-2d_1-d_2\\
T7&4-d_1-d_2\\
T9&4-d_S-d_F-d_H\\
T11&6-d_S-d_F-2d_H.
\end{array}
\]
For T2, T8, and T10, \(\alpha=4\); ordinary K3 virtual vanishing gives \(G_{3,\beta}=0\), including the degree-zero case checked below. Their CJR identity consequently applies to \(r=4,5,6\). This is ordinary K3 theory; reduced K3 invariants are different, as explained by [Maulik–Pandharipande–Thomas, §0.1](https://people.math.ethz.ch/~rahul/mpt.pdf). For T4, T6, and T12, \(\alpha=6\) at every degree. They have no negative-\(t\) range in this unmarked pushforward, and their top comparison retains ordinary GW theory.

The toric edge factors require more than substituting a different number of projective fixed points. Orient an invariant curve from \(i\) to \(j\), with tangent weight \(w\) at \(i\), and take a degree-\(e\) cover. A normal line of degree \(a\), with fiber weight \(\ell\) at \(i\), has weight \(\ell-aw\) at \(j\). Its inverse Euler ratio is
\[
\mathcal N(a,\ell,w,e)=
\begin{cases}
\displaystyle\prod_{r=0}^{ae}(\ell-rw/e)^{-1},&ae\ge0,\\[4pt]
\displaystyle\prod_{r=ae+1}^{-1}(\ell-rw/e),&ae<0.
\end{cases}
\]
The moving tangent factor is
\[
\mathcal T(w,e)=\frac{(-1)^e e^{2e}}{(e!)^2w^{2e}},
\]
and a nef twisting line \(L\) contributes
\[
\prod_{r=0}^{eL\cdot C}(t-L_i+rw/e).
\]
The deck factor \(1/e\) and graph automorphisms are included separately. These conventions agree with [Liu–Sheshmani, Lemma 4.5 and Theorem 4.7](https://sigma-journal.com/2017/048/sigma17-048.pdf).

For the exceptional section of \(\mathbb F_1\), \(a=-1\). A degree-one edge has no normal \(H^0\) or \(H^1\) contribution, while a double cover contributes \(\ell+w/2\) in the numerator. Replacing this by a positive-normal-degree denominator would change the answer. The positive section has normal degree \(+1\), and fibers have normal degree zero. Directions from the other product factor also have normal degree zero.

At a contracted vertex of genus \(h\), valence \(n\), tangent weights \(\tau_k\), and twisting weights \(q_j=t-L_{j,i}\), the factor before flag denominators is
\[
\left(\prod_k\tau_k\right)^{n-1}
\prod_k\Lambda_h^\vee(\tau_k)
\prod_jq_j^{1-n-2h}\Lambda_h^\vee(q_j).
\]
Each flag contributes \((w/e-\psi)^{-1}\). Mumford's relation converts the original inverse Hodge polynomial into this expression. A bivalent node contributes
\(\prod_k\tau_k/((w_1/e_1+w_2/e_2)\prod_jq_j)\), and an univalent endpoint contributes its flag weight. The implementation then applies gluing and the actual forgetful stabilization of temporary tail markings, using the previously checked curve-class backend.

The pure tautological relations come from setting the auxiliary torus parameters to \(\epsilon(1,b,b^2,\ldots)\), with \(b=4,5\). Every tangent, normal-cohomology, and node denominator is checked for resonance. The full proper pushforward is regular at \(\epsilon=0\), so negative Laurent coefficients vanish; its finite term is independent of the ray. The export also checks factor exchange for the symmetric bundles of T6, T7, and T8. These operations generate strata vectors before applying relations.

Every saved vector \(v\) has an exact rational certificate \(uM=v\), where \(M\) is the generated Pixton matrix in the same ordered strata basis. The validity in Chow is [Janda, Theorem 1.1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf), for the [Pandharipande–Pixton–Zvonkine system](https://arxiv.org/abs/1303.1043). The comparison with the preceding projective search is stronger than a comparison of ranks: each toric union row is explicitly written as a rational linear combination of that earlier union basis. No completeness conjecture about all tautological relations is used.

The following generated tables give every computed top integral, all genus-zero CY3 values, and every ordinary higher-genus value established below. The relation counts describe linear independence in the strata algebra before taking the Chow quotient. Counts from different targets must not be added to claim independent relations.

| Target | Ambient | Bundle coefficients | Type | Euler characteristic |
|---|---|---|---|---:|
| T1 | p1xp2 | ((1, 1),) | Fano | 4 |
| T2 | p1xp2 | ((2, 3),) | K3 | 24 |
| T3 | p1xp3 | ((1, 2),) | Fano | 4 |
| T4 | p1xp3 | ((2, 4),) | CY3 | -168 |
| T5 | p2xp2 | ((1, 2),) | Fano | 6 |
| T6 | p2xp2 | ((3, 3),) | CY3 | -162 |
| T7 | p2xp2 | ((1, 1), (1, 1)) | Fano | 6 |
| T8 | p2xp2 | ((1, 1), (2, 2)) | K3 | 24 |
| T9 | f1xp1 | ((1, 1, 1),) | Fano | 5 |
| T10 | f1xp1 | ((2, 3, 2),) | K3 | 24 |
| T11 | f1xp2 | ((1, 1, 1),) | Fano | 8 |
| T12 | f1xp2 | ((2, 3, 3),) | CY3 | -144 |

Per-target independent localization relations (factor-exchange checks are added to the union separately):

| Target | Genus 2 ranks r=1,2,3 | Genus 3 ranks r=2,3,4,5,6 |
|---|---|---|
| T1 | 1, 6, 12 | 5, 8, 22, 32, 44 |
| T2 | 1, 6, 12 | 5, 9, 21, 33, 43 |
| T3 | 1, 6, 12 | 6, 9, 16, 23, 33 |
| T4 | 1, 6, 12 | 6, 9, 16, 23, 33 |
| T5 | 1, 6, 12 | 2, 7, 15, 23, 33 |
| T6 | 1, 6, 12 | 2, 7, 15, 23, 33 |
| T7 | 1, 6, 12 | 2, 7, 15, 23, 33 |
| T8 | 1, 6, 12 | 2, 7, 15, 23, 33 |
| T9 | 1, 6, 12 | 2, 6, 19, 34, 47 |
| T10 | 1, 6, 12 | 2, 6, 19, 36, 53 |
| T11 | 1, 6, 12 | 1, 4, 10, 20, 34 |
| T12 | 1, 6, 12 | 1, 4, 10, 20, 34 |

| Genus | Codimension | Toric rank | Previous projective rank | Combined rank | Pixton rank |
|---|---|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 6 | 6 | 6 | 6 |
| 2 | 3 | 12 | 16 | 16 | 16 |
| 3 | 2 | 6 | 6 | 6 | 6 |
| 3 | 3 | 16 | 30 | 30 | 30 |
| 3 | 4 | 39 | 86 | 86 | 114 |
| 3 | 5 | 59 | 150 | 150 | 259 |
| 3 | 6 | 67 | 168 | 168 | 377 |

Exact twisted top integrals for T1 (p1xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | -1/320 | 1/80640 |
| (0,1) | 7/320 | 13/241920 |
| (1,0) | -7/1920 | -37/241920 |
| (0,2) | -1/192 | 901/483840 |
| (1,1) | 43/384 | -8237/96768 |
| (2,0) | -1/240 | -1/2688 |

Exact twisted top integrals for T2 (p1xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | -1/320 | 1/80640 |
| (0,1) | 11/160 | -121/32256 |
| (1,0) | 1/80 | 1/20160 |
| (0,2) | -33/160 | 121/3584 |
| (1,1) | 109/160 | -3167/161280 |
| (2,0) | 1/40 | 1/2016 |

Exact twisted top integrals for T3 (p1xp3):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | 0 | 0 |
| (0,1) | 23/180 | -53/20160 |
| (1,0) | -1/72 | 1/9072 |
| (0,2) | 0 | -79/3780 |
| (1,1) | 0 | -3889/9072 |
| (2,0) | -1/288 | 1/4032 |

Exact twisted top integrals for T4 (p1xp3):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | 0 | 0 |
| (0,1) | -37/180 | -19/3024 |
| (1,0) | -1/36 | 1/9072 |
| (0,2) | -229/18 | 2/63 |
| (1,1) | -338/45 | 41/1260 |
| (2,0) | -1/18 | 5/4536 |

Exact twisted top integrals for T5 (p2xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | 0 | 0 |
| (0,1) | 3/80 | 1/2560 |
| (1,0) | -1/640 | -1/10752 |
| (0,2) | 0 | 1/630 |
| (1,1) | 0 | -19941/35840 |
| (2,0) | -1/128 | -17/13440 |

Exact twisted top integrals for T6 (p2xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | 0 | 0 |
| (0,1) | 9/128 | 1/1280 |
| (1,0) | 9/128 | 1/1280 |
| (0,2) | -27/128 | -9/1280 |
| (1,1) | -6177/640 | 151/5120 |
| (2,0) | -27/128 | -9/1280 |

Exact twisted top integrals for T7 (p2xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | 0 | 0 |
| (0,1) | 11/640 | -1/26880 |
| (1,0) | 11/640 | -1/26880 |
| (0,2) | -1/320 | 1/5376 |
| (1,1) | 41/160 | -4145/1792 |
| (2,0) | -1/320 | 1/5376 |

Exact twisted top integrals for T8 (p2xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0) | 0 | 0 |
| (0,1) | 5/64 | 1/2688 |
| (1,0) | 5/64 | 1/2688 |
| (0,2) | -13/80 | 23/10752 |
| (1,1) | 39/40 | -196723/161280 |
| (2,0) | -13/80 | 23/10752 |

Exact twisted top integrals for T9 (f1xp1):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0,0) | -1/360 | 1/90720 |
| (0,0,1) | -1/288 | -13/90720 |
| (0,1,0) | -1/288 | -13/90720 |
| (1,0,0) | 7/1152 | 101/362880 |
| (0,0,2) | -11/2880 | -41/120960 |
| (0,1,1) | 47/1440 | -2689/725760 |
| (0,2,0) | -11/2880 | -41/120960 |
| (1,0,1) | 1/160 | -1021/2903040 |
| (1,1,0) | 67/2880 | 389/1451520 |
| (2,0,0) | 1/1440 | 1/2880 |

Exact twisted top integrals for T10 (f1xp1):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0,0) | -1/360 | 1/90720 |
| (0,0,1) | 1/90 | 1/22680 |
| (0,1,0) | 1/90 | 1/22680 |
| (1,0,0) | 0 | 25/290304 |
| (0,0,2) | 1/45 | 1/2268 |
| (0,1,1) | -2/9 | -191/7560 |
| (0,2,0) | 1/45 | 1/2268 |
| (1,0,1) | 0 | -37/483840 |
| (1,1,0) | 1/10 | -151/41472 |
| (2,0,0) | 0 | -367/725760 |

Exact twisted top integrals for T11 (f1xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0,0) | 0 | 0 |
| (0,0,1) | 1/960 | -1/40320 |
| (0,1,0) | -1/80 | 1/10080 |
| (1,0,0) | 7/640 | 0 |
| (0,0,2) | -11/1920 | -13/16128 |
| (0,1,1) | 63/640 | -53/3584 |
| (0,2,0) | -1/320 | 1/4480 |
| (1,0,1) | 1/96 | -53/15120 |
| (1,1,0) | 11/1920 | 1/16128 |
| (2,0,0) | 1/480 | 1/2240 |

Exact twisted top integrals for T12 (f1xp2):

| Mori tuple β | Genus 2 | Genus 3 |
|---|---|---|
| (0,0,0) | 0 | 0 |
| (0,0,1) | 1/16 | 1/1440 |
| (0,1,0) | -1/40 | 1/10080 |
| (1,0,0) | 7/640 | 0 |
| (0,0,2) | -3/16 | -1/160 |
| (0,1,1) | 51/80 | -17/1260 |
| (0,2,0) | -1/20 | 1/1008 |
| (1,0,1) | -1/160 | -19/16128 |
| (1,1,0) | 57/640 | 29/20160 |
| (2,0,0) | -7/320 | 1/2240 |

Ordinary genus-zero CY3 invariants, checked at both rays:

| Target | Mori tuple | N₀,β |
|---|---|---|
| T4 | (0, 1) | 640 |
| T4 | (1, 0) | 64 |
| T4 | (0, 2) | 10112 |
| T4 | (1, 1) | 6912 |
| T4 | (2, 0) | 8 |
| T6 | (0, 1) | 189 |
| T6 | (1, 0) | 189 |
| T6 | (0, 2) | 1701/8 |
| T6 | (1, 1) | 8262 |
| T6 | (2, 0) | 1701/8 |
| T12 | (0, 0, 1) | 168 |
| T12 | (0, 1, 0) | 54 |
| T12 | (1, 0, 0) | 9 |
| T12 | (0, 0, 2) | 189 |
| T12 | (0, 1, 1) | 1080 |
| T12 | (0, 2, 0) | 27/4 |
| T12 | (1, 0, 1) | 36 |
| T12 | (1, 1, 0) | 135 |
| T12 | (2, 0, 0) | 9/8 |

Known ordinary higher-genus CY3 numbers and integrated infinity aggregates:

| Target | Genus | Mori tuple | N_g,β | ∫B_β,3g−3 |
|---|---|---|---|---|
| T4 | 2 | (0, 0) | -7/240 | 7/240 |
| T4 | 3 | (0, 0) | 1/8640 | 1/8640 |
| T4 | 2 | (0, 1) | 8/3 | -443/180 |
| T4 | 3 | (0, 1) | 20/189 | 113/1008 |
| T4 | 2 | (1, 0) | 4/15 | -43/180 |
| T4 | 3 | (1, 0) | 2/189 | 95/9072 |
| T4 | 2 | (2, 0) | 8/15 | -43/90 |
| T4 | 3 | (2, 0) | 16/189 | 379/4536 |
| T6 | 2 | (0, 0) | -9/320 | 9/320 |
| T6 | 3 | (0, 0) | 1/8960 | 1/8960 |
| T6 | 2 | (0, 1) | 63/80 | 459/640 |
| T6 | 3 | (0, 1) | 1/32 | -41/1280 |
| T6 | 2 | (1, 0) | 63/80 | 459/640 |
| T6 | 3 | (1, 0) | 1/32 | -41/1280 |
| T12 | 2 | (0, 0, 0) | -1/40 | 1/40 |
| T12 | 3 | (0, 0, 0) | 1/10080 | 1/10080 |
| T12 | 2 | (0, 0, 1) | 7/10 | 51/80 |
| T12 | 3 | (0, 0, 1) | 1/36 | -41/1440 |
| T12 | 2 | (0, 1, 0) | 9/40 | -1/5 |
| T12 | 3 | (0, 1, 0) | 1/112 | 89/10080 |
| T12 | 2 | (1, 0, 0) | 3/80 | 17/640 |
| T12 | 3 | (1, 0, 0) | 1/672 | -1/672 |
| T12 | 2 | (0, 2, 0) | 9/20 | -2/5 |
| T12 | 3 | (0, 2, 0) | 1/14 | 71/1008 |
| T12 | 2 | (2, 0, 0) | 3/40 | -17/320 |
| T12 | 3 | (2, 0, 0) | 1/84 | 11/960 |


The Chern numbers are independent of the GW localization. They come from
\[
c(TY)=c(TX)\prod_j(1+L_j)^{-1}|_Y,
\qquad [Y]=\prod_jL_j,
\]
using the ordinary product Chow ring and \(S^2=-SF\). The Fano surfaces T1, T7, and T9 have \(\int c_1^2=8,6,7\) and Euler characteristics \(4,6,5\). Their constant classes are
\[
G_{2,0}=(8,6,7)\lambda_2,
\qquad G_{3,0}=(8,6,7)\lambda_1\lambda_3,
\]
with the entries assigned in that order. The K3 cases have \(c_1=0\), \(\int c_2=24\), and zero ordinary constant classes in these genera.

For threefolds put \(C_Y=\int(c_3-c_1c_2)\). For the Fano rows T3, T5, T11, the constants are \(-20,-18,-16\); for the CY3 rows T4, T6, T12, they are the Euler characteristics \(-168,-162,-144\). The constant-map obstruction bundle gives
\[
G_{2,0}=C_Y\lambda_1\lambda_2,
\qquad G_{3,0}=-C_Y\lambda_1\lambda_2\lambda_3.
\]
Thus the ordinary CY3 numbers are
\[
N_{2,0}=\chi(Y)/5760,
\qquad N_{3,0}=-\chi(Y)/1451520,
\]
by [Faber–Pandharipande, Theorem 4](https://arxiv.org/abs/math/9810173). The genus-three values for T4, T6, T12 are \(1/8640,1/8960,1/10080\).

For the ambient threefolds, \(\int(c_3(TX)-c_1(TX)c_2(TX))=-18\) on \(\mathbb P^1\times\mathbb P^2\), and \(-16\) on \(\mathbb F_1\times\mathbb P^1\). These independently predict the ambient constant terms
\[
Z_0(t)=(-18\text{ or }-16)\lambda_1\lambda_2t^{-1}\quad(g=2),
\qquad
Z_0(t)=(18\text{ or }16)\lambda_1\lambda_2\lambda_3t^{-2}\quad(g=3).
\]
Every other ambient in the sample is a fourfold, whose degree-zero stable-map virtual dimension is negative for \(g\ge2\); hence its ordinary nonequivariant ambient constant term is zero. This does not force the compact complete-intersection constant to vanish. All these statements were checked against the computed classes modulo the Pixton relations.

There are additional exact ordinary CY3 values. The source audit independently derives the counts 64 for T4 in \((1,0)\), 189 for each of T6's primitive classes, and 9 for T12 in \((1,0,0)\). The following further parameter-space calculations, carried out in this implementation, extend those checks to all seven primitive extremal classes.

For T4, a \(\mathbb P^1\) fiber lies in the hypersurface when three quartics on \(\mathbb P^3\) vanish, giving \(4^3=64\). Lines in the \(\mathbb P^3\) factor are parameterized by \(\mathbb P^1\times\operatorname{Gr}(2,4)\). Restricting the equation gives the globally generated rank-five bundle \(\mathcal O(2)\boxtimes\operatorname{Sym}^4U^\vee\). Its top Chern number is
\[
[hx^3y^2]\left((x-y)\prod_{r=0}^4(2h+(4-r)x+ry)\right)=640.
\]
Here \(h\) is the hyperplane class on \(\mathbb P^1\), and \(x,y\) are the roots of \(U^\vee\). The coefficient rule integrates the Grassmannian factor.

For T6, lines in one \(\mathbb P^2\) factor are parameterized by \(\operatorname{Gr}(2,3)\times\mathbb P^2\). One has \(\int_{\operatorname{Gr}(2,3)}c_2(\operatorname{Sym}^3U^\vee)=21\), and the other cubic factor contributes \(3^2\), giving \(21\cdot9=189\) in either class.

For T12, the negative-section fibers are cut out by two cubics on \(\mathbb P^2\), giving 9. For the fiber class of \(\mathbb F_1\to\mathbb P^1\), the parameter space is \(\mathbb P^1\times\mathbb P^2\). The restriction bundle splits as
\[
\mathcal O(1,3)\oplus\mathcal O(2,3)\oplus\mathcal O(3,3),
\]
since \(\pi_*\mathcal O_{\mathbb F_1}(2S+3F)=\mathcal O(1)\oplus\mathcal O(2)\oplus\mathcal O(3)\). Therefore its top Chern number is \(9(1+2+3)=54\). Finally, lines in the \(\mathbb P^2\) factor are parameterized by \(\mathbb F_1\times\operatorname{Gr}(2,3)\). Their count is \(21\int_{\mathbb F_1}(-K)^2=21\cdot8=168\).

In each of these seven classes, the parameter space exhausts the possible images. The bundle is generated by restrictions of defining equations, so a general equation has transverse isolated zeros. The derivative of this section is the map \(H^0(N_{C/X})\to H^0(E|_C)\). Transversality gives \(H^0(N_{C/Y})=0\); adjunction gives degree \(-2\), so \(N_{C/Y}=\mathcal O(-1)^{\oplus2}\). The ordinary genus-zero localization values independently reproduce all seven counts at both torus rays. The top-Chern interpretation uses the [regular-section formula](https://stacks.math.columbia.edu/tag/0FA8).

By [Faber–Pandharipande, Theorem 3](https://arxiv.org/abs/math/9810173), a fixed \((-1,-1)\) curve contributes
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3}.
\]
Consequently the primitive ordinary genus-three numbers are
\[
\begin{array}{c|c|c}
Y&\beta&N_{3,\beta}\\\hline
T4&(1,0)&2/189\\
T4&(0,1)&20/189\\
T6&(1,0),(0,1)&1/32\\
T12&(1,0,0)&1/672\\
T12&(0,1,0)&1/112\\
T12&(0,0,1)&1/36.
\end{array}
\]
For genus two divide the corresponding curve counts by 240; for genus three divide by 6048. These are ordinary connected GW numbers, including contracted higher-genus components, and they extend to smooth deformations of the general member.

The pure double classes \((2,0)\) for T4, and \((2,0,0),(0,2,0)\) for T12, also exhaust images on the same isolated curves. A zero degree in the other projections forces a map onto a single fiber or negative section. Distinct such curves are disjoint, so a connected map cannot combine them. Thus their ordinary values are the curve count times \(C(g,2)\). In genus three they are \(16/189,1/84,1/14\); in genus two they are \(8/15,3/40,9/20\). The genus-zero values are respectively \(8,9/8,27/4\), and agree with the independent convex-Lefschetz localization.

No such exhaustion has been asserted for doubles of a projective line class, where conics can contribute, or for mixed Mori classes. The report nevertheless gives ordinary genus-zero invariants in all nineteen positive CY3 tuples: nefness makes \(R^1\pi_*f^*E=0\) in genus zero, so [Chen–Janda–Webb, Theorem 1.2.3](https://arxiv.org/abs/1911.09825v2) applies. Higher-genus ordinary values in the remaining classes stay unknown. The genus-zero dual-Euler sum is converted with the sign \((-1)^{c+c_1(E)\cdot\beta}\); omitting it would, for example, change the T4 fiber count to \(-64\).

The known ordinary numbers give genuine numerical infinity-aggregate comparisons. For instance T12 in genus three has
\[
\int Z_{(1,0,0),6}=0,\qquad \sigma_{3,(1,0,0)}=-1,
\qquad \int B_{(1,0,0),6}=-1/672.
\]
For the exceptional double class,
\[
\int Z_{(2,0,0),6}=1/2240,\qquad \sigma_{3,(2,0,0)}=1,
\qquad \boxed{\int B_{(2,0,0),6}=1/84-1/2240=11/960.}
\]
These examples show why ordinary compact theory remains separate from the ambient twist. For the K3 T10, the mixed sector provides a class identity
\[
B_{(1,1,0),r}=-Z_{(1,1,0),r}\quad(r=4,5,6),
\qquad \int B_{(1,1,0),6}=151/41472.
\]
The full stabilized expressions in every codimension, beyond their top integrals, are preserved in the per-target `finite-classes-g*.json.gz` files.

The validation includes all twelve independent Chern rows, the Mori-class distinction for the positive section, seven primitive curve counts and three double-cover counts at two rays, recovery of the projective-space coefficient backend, constant-map class checks, and direct Hodge integration for sixteen selected higher-genus integrals. The last group includes both the exceptional double cover and the mixed \(S+F\) sector in both genera. The module also rejects a nonample CY surface construction whose vanishing first Chern class alone would not prove that it is K3: a product of two elliptic curves must not be mislabeled K3.

Reproduce the calculation with Sage's Python from the repository root:

```bash
python -c 'from sage.all import *; load("tests/test_toric_complete_intersections.sage")'
python calculations/toric_ci_scan.py --workers 3
python calculations/toric_ci_export.py
```

The environment used is Sage 10.9 on the native Linux OrbStack `agents` machine. From macOS, use `orb -m agents -w /home/patrick/agent-work/repos/research/test/log-glsm-plane-cubic /bin/bash -lc`, with `/home/patrick/agent-work/tools/mamba/envs/sage/bin` on `PATH`. The existing projective scan supplies the shared stabilization cache and Pixton matrices; when the cache is absent it is rebuilt. The export's comparison uses the previous projective union basis and checks that its strata ordering agrees.

For an individual computation:

```python
load('toric_complete_intersections.sage')
Y = toric_sample()[11]                       # T12, anticanonical F1 x P2
C = ToricCompleteIntersectionCycles(Y, 3)
C.cycle((1,1,0), 5)                         # full stabilized codimension-five class
C.top_integral((2,0,0))                     # 1/2240
C.aggregate_relation((2,0,0), 6)            # retains unknown compact class at class level
toric_known_ordinary_gw(Y, 3, (2,0,0))      # scalar refinement: 1/84
```

The scalar refinement is deliberately separate from the class-valued aggregate API. An unknown compact class is returned as `None`, not zero. The last example does not purport to compute an individual effective infinity vertex.

The [result directory](../results/toric-complete-intersections/2026-11-06) contains exact raw relation vectors and their origins, Pixton matrices and certificates, the [readable independent relation basis](../results/toric-complete-intersections/2026-11-06/independent-relations.txt), full finite classes, ordinary GW values, verification records, and source/artifact hashes. The `previous-span-g*-r*.json.gz` files explicitly express the toric basis in the earlier projective basis; `combined-relations-g*-r*.json.gz` certifies the combined span in Pixton's system. This search makes no claim about higher Mori degree, additional torus rays, markings, insertions, or other toric ambient families.
