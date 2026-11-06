# Genus four, up to two markings: exact calculations and bounded Chow relations

This pass retains all **43 ambient presentations**, including presentations of isomorphic compact targets. It uses genus four, zero through two labeled markings, and ambient degree zero through two (Mori coefficient sum at most two for the toric examples). Each mark carries `1`, `psi = tau_1(1)` on the **stable-map space**, or the ambient divisor `H`.

The numerical calculation covers every profile in this box. The Chow relation search is deliberately bounded to **codimensions 0–3** on each of Mbar(4,0), Mbar(4,1), and Mbar(4,2). It does not search every Chow degree up to dimensions 9, 10, and 11. Numerical top coefficients are nevertheless computed in those top dimensions. No individual punctured effective cycles are reconstructed.

## Definitions and the marked comparison

Let Y be the zero locus of a regular section of E on X, with m = dim Y and c = rank E. The ambient class beta includes the sum over compact curve classes mapping to it. The calculation uses the full derived dual Euler ratio with a **single** equivariant parameter t:
\[
\Theta(t)=e_t\bigl((R\pi_*f^*E)^\vee\bigr),\qquad
\mathcal Z_{\beta,\mathbf I}(t)=s_*\left(\Theta(t)\prod_i I_i\cap[\overline{\mathcal M}_{4,n}(X,\beta)]^{\rm vir}\right).
\]
Here s forgets the map and stabilizes. R^1 is retained. This is different from the ordinary compact GW class of Y.

Put a = number of psi or H insertions. Then
\[
\alpha_0=3m-(c_1(TX)-c_1(E))\cdot\beta,\quad
\mathcal Z_{\beta,\mathbf I}(t)=\sum_{r=0}^{9+n}t^{\alpha_0+a-r}\mathcal Z_{\beta,\mathbf I;r},\quad
\sigma_{4,\beta}=(-1)^{-3c+c_1(E)\cdot\beta}.
\]
In particular, the numerical top coefficient multiplies t^(alpha_0+a−9−n). Markings do not add another term to alpha_0 or the comparison sign.

The globally inserted CJR identity has the form
\[
\mathcal Z_{\beta,\mathbf I}(t)+\mathcal B_{\beta,\mathbf I}(t)
=\sigma_{4,\beta}\,\mathcal G^Y_{\beta,\mathbf I}.
\]
B denotes the **total infinity aggregate** after the correct stabilization, including its higher-rank projective-bundle geometry. At a negative t power, B = −Z. At t^0, B = sigma G − Z. A negative t power does **not** by itself imply Z = 0. Positive t powers are not labeled compact comparisons in the numerical files.

The private CJR III draft is used with compact-type zero-contact unit-sector markings. The evaluation class on this sector may be 1 or H. True stable-map descendants are inserted by pullback along the log-to-stable-map stabilization before localization. This avoids identifying CJR's log-domain cotangent line with a stable-map line. See the [source audit](GENUS_FOUR_MARKED_SOURCES.md), especially §§3–6, for the exact private and public references.

## Why all requested numerical markings reduce to one number

Let Z_beta be the unmarked numerical coefficient of the ambient twist and d_H = H·beta. The index complex pulls back when a marking is forgotten, so string, dilaton, and divisor equations apply to the full twist. They give
\[
\begin{array}{c|c}
\text{profile}&\text{numerical coefficient}\\\hline
\varnothing&Z_\beta\\
(1)&0\\
(\psi)&6Z_\beta\\
(H)&d_HZ_\beta\\
(1,1),(1,\psi),(\psi,1),(1,H),(H,1)&0\\
(\psi,\psi)&42Z_\beta\\
(\psi,H),(H,\psi)&7d_HZ_\beta\\
(H,H)&d_H^2Z_\beta.
\end{array}
\]
For example, forgetting the psi mark in (psi,H) gives 2g−2+1 = 7, rather than 6. The same formulas hold for ordinary compact numbers, replacing Z_beta by N_(4,beta). A zero multiplier gives a proved zero even when N_(4,beta) is unknown. There is no 1/n! because markings are labeled.

H is the standard hyperplane in projective space, the Plücker class on a Grassmannian, the sum of factor hyperplanes on a product, and S+2F+H_P on F_1 × P^q. The latter pairs with a Mori tuple as d_S+d_F+d_P. Thus d_H equals the coefficient sum for the chosen toric polarizations; the implementation still evaluates the actual divisor pairing.

These numerical formulas do not construct marked Chow classes. At class level,
\[
\psi_i^{\rm map}=s^*\psi_i^{\rm curve}+T_i,
\]
where T_i is the positive-degree rational-tail contraction divisor. The new marked localization retains that correction: a descendant at an unstable marked endpoint has weight −w; it is not replaced by a psi class on the stabilized curve. A unit-marked class is a forgetful pullback, and need not vanish even though its numerical integral does.

## Exact Hodge and localization algorithms

The new module is `genus_four_marked.sage`. Existing genus-two/three modules are unchanged. The rank-four Hodge algebra uses
\[
\lambda_1^2=2\lambda_2,\quad
\lambda_2^2=2\lambda_1\lambda_3-2\lambda_4,\quad
\lambda_3^2=2\lambda_2\lambda_4,\quad\lambda_4^2=0.
\]
In particular, the extra −2 lambda_4 term cannot be omitted. The general rewrite of a repeated factor is
\[
\lambda_i^2=2\sum_{k=1}^i(-1)^{k+1}\lambda_{i-k}\lambda_{i+k},
\qquad\lambda_0=1,\quad\lambda_j=0\ (j>g).
\]
It terminates because each replacement increases the sum of squared indices at fixed total degree. Hodge classes pull back from Mbar_4, so Hodge degree above nine vanishes even after adding marks.

Numerical Hodge integrals are computed by Mumford's GRR recursion rather than expanding large boundary classes:
\[
\operatorname{ch}_k(\mathbb E)=\frac{B_{k+1}}{(k+1)!}
\left[\kappa_k-\sum_i\psi_i^k+
\frac12\xi_*\sum_{a=0}^{k-1}(-1)^a\psi_+^a\psi_-^{k-1-a}\right]
\quad(k\text{ odd}).
\]
On a separating boundary, positive Chern characters add across the two components; on a nonseparating boundary they pull back from genus g−1. Expanding the remaining character factors reduces the integrals to lower character degree and pure psi intersections, evaluated by DVV. The factor 1/2 and both ordered branch partitions are retained.

Each fixed-locus integral is formed exactly over Q(epsilon), then expanded only to the finite jet needed for epsilon powers ≤0. If numerator and denominator have valuations u and v, the required precision is 1−u+v. This avoids forming a common denominator over all fixed graphs. The calculation is repeated on rays with bases four and five.

For Chow classes, real marking labels are preserved. Nonconstant edges are contracted; empty rational tails become temporary markings and are forgotten; marked tails become retained marking legs. Contracted rational components which become unstable are stabilized globally. A Hodge factor is pulled through forgetting before expansion, and psi monomials of degree below the forgetful relative dimension are discarded as zero. Graph monomials are then converted to explicit admcycles strata coordinates.

There is also an exact grading shortcut. Write alpha_X = dim(X)(g−1) − c_1(TX)·beta. The degree-j term of the Euler twist has ordinary codimension j≥0. At curve codimension r, each corresponding auxiliary term has epsilon degree alpha_X+a+j−r. If alpha_X+a exceeds the requested codimension bound, every relevant Laurent coefficient is absent, graph by graph. Such sectors are recorded as empty searches in that window, not as computationally certified new relations. This explains many zero rows for the higher-dimensional ambients.

## Compact constants, vanishing, and local curves

The independent genus-four constants are
\[
\int_{\overline{\mathcal M}_4}\lambda_2\lambda_3\lambda_4=\frac1{87091200},\qquad
N^Y_{4,0}=\frac{\int_Y(c_3-c_1c_2)}{87091200}\quad(\dim Y=3).
\]
For a CY3 this is chi(Y)/87091200. On a surface the unmarked constant-map class is (integral c_1^2) lambda_2 lambda_4; it is not a top number. On a fourfold, the current constant-map virtual dimension is negative. Ordinary positive-degree numbers for non-CY targets vanish by dimension unless that virtual dimension is zero; the degree-one cubic fourfold is an explicitly unknown exception in this pass.

All K3 presentations use **ordinary** theory. Positive degrees vanish by the ordinary K3 obstruction. In degree zero, the full obstruction polynomial is
\[
\lambda_4^2+c_2(TS)(\lambda_3^2-2\lambda_2\lambda_4)=0.
\]
Thus every requested ordinary K3 class vanishes, not merely its integral. The ambient twisted classes can still be nonzero.

For the elliptic cubic, ambient degrees one and two contain no compact maps: a nonconstant cover has ambient degree divisible by three. At degree zero, a single H insertion gives the nonzero compact class 3 lambda_4; with a psi at another marking it gives 3 lambda_4 psi. These classes integrate to zero in the present profiles. This is an explicit example of why scalar zero is not cycle zero.

For an isolated (−1,−1) rational curve, Faber–Pandharipande give
\[
C(4,k)=\frac{k^5}{172800}.
\]
Hence L isolated primitive lines contribute L/172800 and pure double covers contribute 32L/172800. The previously justified toric pure doubles are T4:(2,0), T12:(2,0,0), and T12:(0,2,0). Mixed toric classes and unsupported embedded-conic cases remain unknown. A genus-zero count alone is not used to infer their higher-genus values.

For Gr(2,7)[1^7], the previously audited incidence and obstruction argument exhausts degree-two images by 196 lines and 1225 smooth conics, each with normal bundle O(−1)⊕O(−1). The conic count retains the dependency on the verified genus-zero value N_(0,2)=2499/2. Therefore
\[
N_{4,0}=-\frac7{6220800},\qquad N_{4,1}=\frac{49}{43200},\qquad
N_{4,2}=\frac{1225+32\cdot196}{172800}=\frac{833}{19200}.
\]
Its ambient virtual dimension is −21+7D+n < 0 for D≤2 and n≤2, so all its nonequivariant ambient twisted classes vanish. The CJR infinity aggregate is the signed compact class. The signs in degrees 0,1,2 are −,+,−. In degree two the ordinary marked values are
\[
N(\psi)=\frac{833}{3200},\quad N(H)=\frac{833}{9600},\quad
N(\psi,\psi)=\frac{5831}{3200},\quad N(\psi,H)=\frac{5831}{9600},\quad
N(H,H)=\frac{833}{4800};
\]
unit-containing profiles are zero. Apply the negative comparison sign only when forming the degree-two infinity aggregate.

## What the relation search proves

Negative auxiliary epsilon coefficients vanish because the proper ambient pushforward has a regular nonequivariant limit. Differences between the finite coefficients on the two auxiliary rays also vanish. These are the candidate Chow relations. They are distinct from cancellation of negative powers of the CJR parameter t against infinity contributions.

The calculation removes linear dependence over Q within each target, and then across targets separately in every (n,r). It stores the strata basis, the actual generated Pixton matrix M, independent candidate rows V, and coefficients U, and replays
\[
\boxed{UM=V.}
\]
It also verifies that every stored per-target row lies in the certified combined span. Thus the reported bounded relations are inside the **implemented Pixton relation span**. This neither asserts completeness of that span for all tautological relations nor extends the computation to codimensions above three. Equality of scalar integrals was never used as a substitute for this certificate.

The numerical batch contains **181 sectors and 2353 ordered profile values**, checked on two rays. There are **18 unknown ordinary unmarked sectors**; their numerical marking reductions remain explicitly symbolic. The relation batch retains **234 per-target independent rows**, whose combined rank sum is **120**, counted separately by marking number and codimension.

| Markings n | Codimension r | Strata columns | Found rank | Pixton rank |
|---:|---:|---:|---:|---:|
| 0 | 2 | 17 | 4 | 4 |
| 0 | 3 | 73 | 15 | 41 |
| 1 | 2 | 39 | 7 | 9 |
| 1 | 3 | 210 | 30 | 117 |
| 2 | 2 | 95 | 11 | 19 |
| 2 | 3 | 673 | 53 | 354 |


## A concrete relation and its smooth-locus interpretation

The elliptic-cubic ambient twist in degree two gives a codimension-two relation from its epsilon^(−2) coefficient on the base-four ray. In the ordered basis S_0,...,S_16 in `strata-n0-r2.json`, its primitive integral coefficient vector is
\[
(-1920,180,-312,1260,-120,300,84,360,72,-660,-28,140,-28,20,52,20,1).
\]
In other words, the dot product of this vector with those **specified compactified strata classes** is zero. The first two basis classes are kappa_2 and kappa_1^2; all other classes are boundary-supported. The complete graph descriptions and exact certificate are exported in `sample-relation.md` and `sample-relation.json`.

Restricting to smooth curves therefore gives the particularly simple relation
\[
\boxed{3\kappa_1^2=32\kappa_2\quad\text{in }R^2(\mathcal M_4).}
\]
The full compactified equation includes the boundary terms above. Geometrically, kappa_a is the pushforward of the (a+1)-st power of the relative log canonical class of the universal curve. The equation constrains these characteristic classes in a family of smooth genus-four curves; when the family acquires nodes, the displayed boundary correction records the degeneration contribution. This is a recovered Pixton consequence, not a claim of a new relation outside the known ideal.

## All 43 examples

The tables below list the unmarked seed numbers. Every one- and two-mark value is obtained by the displayed universal factors and is also stored explicitly, with ordered marking labels, in the corresponding `*-numbers.json` file. Z is the ambient dual-Euler top coefficient; N is the ordinary compact number. Their meanings must not be interchanged. A question mark means unknown, not zero. The curve tuples and defining multidegrees use the previously documented ambient bases.

### 1. P^2[3] (p2-3; elliptic)

Target dimension 1; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 3, (n=0, r=3): 6, (n=1, r=2): 5, (n=1, r=3): 13, (n=2, r=2): 8, (n=2, r=3): 28.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 3 | -1 | -1/1612800 | 0 |
| (1,) | 1 | 3 | 1 | 5617/6451200 | 0 |
| (2,) | 2 | 3 | -1 | 92921/2457600 | 0 |

### 2. P^3[4] (p3-4; K3)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 2, (n=1, r=3): 5, (n=2, r=2): 3, (n=2, r=3): 11.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | -1 | -1/4354560 | 0 |
| (1,) | 1 | 6 | -1 | 10463/2073600 | 0 |
| (2,) | 2 | 6 | -1 | 21391/120960 | 0 |

### 3. P^3[2] (p3-2; Fano)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 2, (n=1, r=3): 5, (n=2, r=2): 3, (n=2, r=3): 11.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | -1 | -1/4354560 | 0 |
| (1,) | 1 | 4 | -1 | 15977/14515200 | 0 |
| (2,) | 2 | 2 | -1 | 1430689/3628800 | 0 |

### 4. P^3[3] (p3-3; Fano)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 2, (n=1, r=3): 5, (n=2, r=2): 3, (n=2, r=3): 11.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | -1 | -1/4354560 | 0 |
| (1,) | 1 | 5 | 1 | 9067/3096576 | 0 |
| (2,) | 2 | 4 | -1 | 16658659/58060800 | 0 |

### 5. P^4[2] (p4-2; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 1, (n=1, r=3): 5, (n=2, r=2): 1, (n=2, r=3): 7.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/4354560 |
| (1,) | 1 | 6 | -1 | 635/387072 | 0 |
| (2,) | 2 | 3 | -1 | 6785615/3096576 | 0 |

### 6. P^4[3] (p4-3; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 1, (n=1, r=3): 5, (n=2, r=2): 1, (n=2, r=3): 7.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/2903040 |
| (1,) | 1 | 7 | 1 | 45653/18579456 | 0 |
| (2,) | 2 | 5 | -1 | 5506027/6967296 | 0 |

### 7. P^4[4] (p4-4; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 1, (n=1, r=3): 5, (n=2, r=2): 1, (n=2, r=3): 7.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/1088640 |
| (1,) | 1 | 8 | -1 | 649/229376 | 0 |
| (2,) | 2 | 7 | -1 | -249335/6967296 | 0 |

### 8. P^5[2] (p5-2; Fano)

Target dimension 4; chosen H coefficients [1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 12 | -1 | 0 | 0 |
| (1,) | 1 | 8 | -1 | 521/737280 | 0 |
| (2,) | 2 | 4 | -1 | 219669283/29030400 | 0 |

### 9. P^5[3] (p5-3; Fano)

Target dimension 4; chosen H coefficients [1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 12 | -1 | 0 | 0 |
| (1,) | 1 | 9 | 1 | 521/737280 | ? |
| (2,) | 2 | 6 | -1 | 113652649/77414400 | 0 |

### 10. P^5[4] (p5-4; Fano)

Target dimension 4; chosen H coefficients [1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 12 | -1 | 0 | 0 |
| (1,) | 1 | 10 | -1 | 521/737280 | 0 |
| (2,) | 2 | 8 | -1 | 577939/1935360 | 0 |

### 11. P^5[5] (p5-5; Fano)

Target dimension 4; chosen H coefficients [1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 12 | -1 | 0 | 0 |
| (1,) | 1 | 11 | 1 | 521/737280 | 0 |
| (2,) | 2 | 10 | -1 | 63267893/580608 | 0 |

### 12. P^4[5] (p4-5; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks (n=0, r=2): 1, (n=0, r=3): 3, (n=1, r=2): 1, (n=1, r=3): 5, (n=2, r=2): 1, (n=2, r=3): 7.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/435456 |
| (1,) | 1 | 9 | 1 | 146963/55738368 | 115/6912 |
| (2,) | 2 | 9 | -1 | -25375061/27869184 | ? |

### 13. P^5[2,4] (p5-2-4; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | 1 | 0 | -11/5443200 |
| (1,) | 1 | 9 | 1 | 521/737280 | 1/135 |
| (2,) | 2 | 9 | 1 | -2530901/9676800 | ? |

### 14. P^5[3,3] (p5-3-3; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | 1 | 0 | -1/604800 |
| (1,) | 1 | 9 | 1 | 521/737280 | 39/6400 |
| (2,) | 2 | 9 | 1 | -6770959/38707200 | ? |

### 15. P^6[2,2,3] (p6-2-2-3; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/604800 |
| (1,) | 1 | 9 | 1 | 0 | 1/240 |
| (2,) | 2 | 9 | -1 | 34170049/16588800 | ? |

### 16. P^7[2,2,2,2] (p7-2-2-2-2; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | 1 | 0 | -1/680400 |
| (1,) | 1 | 9 | 1 | 0 | 2/675 |
| (2,) | 2 | 9 | 1 | 655885351/2721600 | ? |

### 17. P^1 × P^2 [(1,1)] (T1; Fano)

Target dimension 2; chosen H coefficients [1, 1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 6 | -1 | -1/4838400 | 0 |
| (0, 1) | 1 | 4 | 1 | -191/116121600 | 0 |
| (1, 0) | 1 | 5 | 1 | -53/9676800 | 0 |
| (0, 2) | 2 | 2 | -1 | -3601/19353600 | 0 |
| (1, 1) | 2 | 3 | -1 | 800791/92897280 | 0 |
| (2, 0) | 2 | 4 | -1 | -467/19353600 | 0 |

### 18. P^1 × P^2 [(2,3)] (T2; K3)

Target dimension 2; chosen H coefficients [1, 1]. Per-target relation ranks (n=0, r=3): 1, (n=1, r=3): 1, (n=2, r=3): 1.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 6 | -1 | -1/4838400 | 0 |
| (0, 1) | 1 | 6 | 1 | 989/38707200 | 0 |
| (1, 0) | 1 | 6 | -1 | -1/151200 | 0 |
| (0, 2) | 2 | 6 | -1 | 577/129024 | 0 |
| (1, 1) | 2 | 6 | 1 | 69757/1105920 | 0 |
| (2, 0) | 2 | 6 | -1 | -41/302400 | 0 |

### 19. P^1 × P^3 [(1,2)] (T3; Fano)

Target dimension 3; chosen H coefficients [1, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 9 | -1 | 0 | -1/4354560 |
| (0, 1) | 1 | 7 | -1 | -17/2073600 | 0 |
| (1, 0) | 1 | 8 | 1 | 0 | 0 |
| (0, 2) | 2 | 5 | -1 | 10021/725760 | 0 |
| (1, 1) | 2 | 6 | 1 | 2071/56700 | 0 |
| (2, 0) | 2 | 7 | -1 | -1/90720 | 0 |

### 20. P^1 × P^3 [(2,4)] (T4; CY3)

Target dimension 3; chosen H coefficients [1, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 9 | -1 | 0 | -1/518400 |
| (0, 1) | 1 | 9 | -1 | -5/290304 | 1/270 |
| (1, 0) | 1 | 9 | -1 | 0 | 1/2700 |
| (0, 2) | 2 | 9 | -1 | -2239/544320 | ? |
| (1, 1) | 2 | 9 | -1 | 49369/1814400 | ? |
| (2, 0) | 2 | 9 | -1 | -1/45360 | 8/675 |

### 21. P^2 × P^2 [(1,2)] (T5; Fano)

Target dimension 3; chosen H coefficients [1, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 9 | -1 | 0 | -1/4838400 |
| (0, 1) | 1 | 8 | -1 | 1/204800 | 0 |
| (1, 0) | 1 | 7 | 1 | 1/204800 | 0 |
| (0, 2) | 2 | 7 | -1 | 4093/6451200 | 0 |
| (1, 1) | 2 | 6 | 1 | 48371/921600 | 0 |
| (2, 0) | 2 | 5 | -1 | -113/921600 | 0 |

### 22. P^2 × P^2 [(3,3)] (T6; CY3)

Target dimension 3; chosen H coefficients [1, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 9 | -1 | 0 | -1/537600 |
| (0, 1) | 1 | 9 | 1 | 1/204800 | 7/6400 |
| (1, 0) | 1 | 9 | 1 | 1/204800 | 7/6400 |
| (0, 2) | 2 | 9 | -1 | 123/204800 | ? |
| (1, 1) | 2 | 9 | -1 | 578227/17203200 | ? |
| (2, 0) | 2 | 9 | -1 | 123/204800 | ? |

### 23. P^2 × P^2 [(1,1); (1,1)] (T7; Fano)

Target dimension 2; chosen H coefficients [1, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 6 | 1 | 0 | 0 |
| (0, 1) | 1 | 5 | 1 | 1/204800 | 0 |
| (1, 0) | 1 | 5 | 1 | 1/204800 | 0 |
| (0, 2) | 2 | 4 | 1 | 11/184320 | 0 |
| (1, 1) | 2 | 4 | 1 | 15491/92160 | 0 |
| (2, 0) | 2 | 4 | 1 | 11/184320 | 0 |

### 24. P^2 × P^2 [(1,1); (2,2)] (T8; K3)

Target dimension 2; chosen H coefficients [1, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0) | 0 | 6 | 1 | 0 | 0 |
| (0, 1) | 1 | 6 | -1 | 1/204800 | 0 |
| (1, 0) | 1 | 6 | -1 | 1/204800 | 0 |
| (0, 2) | 2 | 6 | 1 | 129/89600 | 0 |
| (1, 1) | 2 | 6 | 1 | 11902451/51609600 | 0 |
| (2, 0) | 2 | 6 | 1 | 129/89600 | 0 |

### 25. F_1 × P^1 [(1,1,1)] (T9; Fano)

Target dimension 2; chosen H coefficients [1, 2, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0, 0) | 0 | 6 | -1 | -1/5443200 | 0 |
| (0, 0, 1) | 1 | 5 | 1 | -109/21772800 | 0 |
| (0, 1, 0) | 1 | 5 | 1 | -109/21772800 | 0 |
| (1, 0, 0) | 1 | 5 | -1 | 15389/1393459200 | 0 |
| (0, 0, 2) | 2 | 4 | -1 | -11/518400 | 0 |
| (0, 1, 1) | 2 | 4 | -1 | -59/2419200 | 0 |
| (0, 2, 0) | 2 | 4 | -1 | -11/518400 | 0 |
| (1, 0, 1) | 2 | 4 | 1 | 19/23224320 | 0 |
| (1, 1, 0) | 2 | 4 | 1 | 281/69672960 | 0 |
| (2, 0, 0) | 2 | 4 | -1 | 1349/14515200 | 0 |

### 26. F_1 × P^1 [(2,3,2)] (T10; K3)

Target dimension 2; chosen H coefficients [1, 2, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0, 0) | 0 | 6 | -1 | -1/5443200 | 0 |
| (0, 0, 1) | 1 | 6 | -1 | -1/170100 | 0 |
| (0, 1, 0) | 1 | 6 | -1 | -1/170100 | 0 |
| (1, 0, 0) | 1 | 6 | 1 | 1847/348364800 | 0 |
| (0, 0, 2) | 2 | 6 | -1 | -41/340200 | 0 |
| (0, 1, 1) | 2 | 6 | -1 | 739/680400 | 0 |
| (0, 2, 0) | 2 | 6 | -1 | -41/340200 | 0 |
| (1, 0, 1) | 2 | 6 | 1 | 3263/348364800 | 0 |
| (1, 1, 0) | 2 | 6 | 1 | 379/12902400 | 0 |
| (2, 0, 0) | 2 | 6 | -1 | -4321/43545600 | 0 |

### 27. F_1 × P^2 [(1,1,1)] (T11; Fano)

Target dimension 3; chosen H coefficients [1, 2, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0, 0) | 0 | 9 | -1 | 0 | -1/5443200 |
| (0, 0, 1) | 1 | 7 | 1 | 1/230400 | 0 |
| (0, 1, 0) | 1 | 8 | 1 | 0 | 0 |
| (1, 0, 0) | 1 | 8 | -1 | 0 | 0 |
| (0, 0, 2) | 2 | 5 | -1 | -251/4147200 | 0 |
| (0, 1, 1) | 2 | 6 | -1 | -11/4838400 | 0 |
| (0, 2, 0) | 2 | 7 | -1 | -1/100800 | 0 |
| (1, 0, 1) | 2 | 6 | 1 | -173/1720320 | 0 |
| (1, 1, 0) | 2 | 7 | 1 | 1/161280 | 0 |
| (2, 0, 0) | 2 | 7 | -1 | 0 | 0 |

### 28. F_1 × P^2 [(2,3,3)] (T12; CY3)

Target dimension 3; chosen H coefficients [1, 2, 1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0, 0, 0) | 0 | 9 | -1 | 0 | -1/604800 |
| (0, 0, 1) | 1 | 9 | 1 | 1/230400 | 7/7200 |
| (0, 1, 0) | 1 | 9 | -1 | 0 | 1/3200 |
| (1, 0, 0) | 1 | 9 | 1 | 0 | 1/19200 |
| (0, 0, 2) | 2 | 9 | -1 | 41/76800 | ? |
| (0, 1, 1) | 2 | 9 | 1 | 23/179200 | ? |
| (0, 2, 0) | 2 | 9 | -1 | -1/50400 | 1/100 |
| (1, 0, 1) | 2 | 9 | -1 | -991/51609600 | ? |
| (1, 1, 0) | 2 | 9 | 1 | 1/161280 | ? |
| (2, 0, 0) | 2 | 9 | -1 | 0 | 1/600 |

### 29. Gr(2,4)[1] (G1; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/4354560 |
| (1,) | 1 | 6 | 1 | 607/9953280 | 0 |
| (2,) | 2 | 3 | -1 | 29744887/87091200 | 0 |

### 30. Gr(2,4)[2] (G2; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/3628800 |
| (1,) | 1 | 7 | -1 | 779/7741440 | 0 |
| (2,) | 2 | 5 | -1 | 39358051/34836480 | 0 |

### 31. Gr(2,4)[3] (G3; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/1451520 |
| (1,) | 1 | 8 | 1 | 887/6967296 | 0 |
| (2,) | 2 | 7 | -1 | 1432217/3628800 | 0 |

### 32. Gr(2,4)[4] (G4; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -11/5443200 |
| (1,) | 1 | 9 | -1 | 4913/34836480 | 1/135 |
| (2,) | 2 | 9 | -1 | 44533/7257600 | ? |

### 33. Gr(2,4)[1,1] (G5; Fano)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | 1 | 0 | 0 |
| (1,) | 1 | 4 | 1 | 607/4976640 | 0 |
| (2,) | 2 | 2 | 1 | 21156073/4838400 | 0 |

### 34. Gr(2,4)[1,2] (G6; Fano)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | 1 | 0 | 0 |
| (1,) | 1 | 5 | -1 | 563/3483648 | 0 |
| (2,) | 2 | 4 | 1 | 37826191/5806080 | 0 |

### 35. Gr(2,4)[1,3] (G7; K3)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | 1 | 0 | 0 |
| (1,) | 1 | 6 | 1 | 4373/23224320 | 0 |
| (2,) | 2 | 6 | 1 | 75006989/43545600 | 0 |

### 36. Gr(2,4)[2,2] (G8; K3)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | 1 | 0 | 0 |
| (1,) | 1 | 6 | 1 | 779/3870720 | 0 |
| (2,) | 2 | 6 | 1 | 12373745/1741824 | 0 |

### 37. Gr(2,5)[1,1,1] (G9; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/4354560 |
| (1,) | 1 | 7 | 1 | 0 | 0 |
| (2,) | 2 | 5 | -1 | 81343067/18579456 | 0 |

### 38. Gr(2,5)[1,1,2] (G10; Fano)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/2177280 |
| (1,) | 1 | 8 | -1 | 0 | 0 |
| (2,) | 2 | 7 | -1 | 107720609/27869184 | 0 |

### 39. Gr(2,5)[1,1,3] (G11; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/580608 |
| (1,) | 1 | 9 | 1 | 0 | 1/320 |
| (2,) | 2 | 9 | -1 | 10899433/6193152 | ? |

### 40. Gr(2,5)[1,2,2] (G12; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -1/725760 |
| (1,) | 1 | 9 | 1 | 0 | 1/432 |
| (2,) | 2 | 9 | -1 | 186853235/55738368 | ? |

### 41. Gr(2,5)[1,1,1,1] (G13; Fano)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | 1 | 0 | 0 |
| (1,) | 1 | 5 | 1 | 0 | 0 |
| (2,) | 2 | 4 | 1 | 81343067/13934592 | 0 |

### 42. Gr(2,5)[1,1,1,2] (G14; K3)

Target dimension 2; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 6 | 1 | 0 | 0 |
| (1,) | 1 | 6 | -1 | 0 | 0 |
| (2,) | 2 | 6 | 1 | 14132585/2654208 | 0 |

### 43. Gr(2,7)[1,1,1,1,1,1,1] (G15; CY3)

Target dimension 3; chosen H coefficients [1]. Per-target relation ranks zero in the declared window.

| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |
|---|---:|---:|---:|---|---|
| (0,) | 0 | 9 | -1 | 0 | -7/6220800 |
| (1,) | 1 | 9 | 1 | 0 | 49/43200 |
| (2,) | 2 | 9 | -1 | 0 | 833/19200 |

## Validation and reproducibility

`tests/test_genus_four_marked.sage` checks the GRR recursion against the lambda-g theorem and the top-Hodge-triple formula in genera two, three, and four; mixed Hodge values against the separate admcycles evaluator; three earlier Gr(2,7) twisted values; all twelve nonempty ordered marked profiles by direct genus-four degree-one localization; the earlier unmarked genus-three Chow implementation; and unit pullback, dilaton, and divisor identities as Chow classes. These checks include unstable marked endpoints and rational bubbles. The complete numerical batch checks all negative auxiliary poles and ray independence in 181 sectors. Matrix certificates are checked independently of numerical integration.

Run inside the pinned Sage Linux environment from the repository root:

```bash
sage tests/test_genus_four_marked.sage
python calculations/genus_four_marked_scan.py --workers 4
python calculations/genus_four_marked_relations.py --workers 4 --max-r 3
python calculations/genus_four_marked_certify.py --max-r 3
python calculations/genus_four_marked_export.py
python calculations/genus_four_marked_verify.py
python calculations/genus_four_marked_visual.py
```

The numerical batch's recorded wall time is in `numerical-summary.json`. Class scans are checkpointed per target and accept `--resume`. The first unoptimized class run was stopped; the completed run uses Hodge pullback and forgetful-dimension pruning. Its partial log is not a completed result. Every reported result comes from the completed per-target artifacts and the final replay.

## References

The [source audit](GENUS_FOUR_MARKED_SOURCES.md) gives page and theorem locators, separate proof dependencies, and the links to all target-specific audits. Main references are:

1. Chen–Janda–Ruan, private *Structural Formulae in Logarithmic Gauged Linear Sigma Models III*, October 17, 2025, Theorem 8.4, Remark 8.5, §9.4, and §10.1; the supplied draft is distinguished from published results.
2. [Chen–Janda–Ruan, The logarithmic gauged linear sigma model](https://arxiv.org/abs/1906.04345), Corollary 1.10 and Proposition 4.1; [Chen–Janda–Webb, Virtual cycles of stable (quasi-)maps with fields](https://arxiv.org/abs/1911.09825), Corollary 1.2.1 and Theorem 1.2.3.
3. [Coates–Givental, Quantum Riemann–Roch, Lefschetz and Serre](https://math.berkeley.edu/~giventh/papers/rls.pdf), Appendices 1–2, for forgetting compatibility, normalization, and descendant/ancestor corrections.
4. [Mumford, Towards an Enumerative Geometry of the Moduli Space of Curves](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf), §5, for Hodge GRR and c(E)c(E dual)=1.
5. [Faber–Pandharipande, Hodge integrals and Gromov–Witten theory](https://arxiv.org/abs/math/9810173), Theorems 3–4 and §4.3, for local curves and constant maps.
6. [Graber–Pandharipande, Localization of virtual classes](https://arxiv.org/abs/alg-geom/9708001), equation (1) and §4; the target-specific GKM factors are documented in the preceding projective, toric, and Grassmannian audits.
7. [Maulik–Pandharipande–Thomas, Curves on K3 surfaces and modular forms](https://arxiv.org/abs/1001.2719), §0.1, for ordinary positive-degree K3 vanishing.
8. [Janda, Relations on Mbar(g,n) via equivariant Gromov–Witten theory of P1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf), Theorem 1.1, for the Chow validity of the Pixton generators. The exact finite membership statements here additionally use the saved U M = V certificates.
