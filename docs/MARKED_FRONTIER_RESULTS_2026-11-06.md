# Higher markings, non-Fano toric ambients, and partial flags

The completed catalogue has **54 presentations**: the previous 43 plus five non-Fano toric presentations and six type-A flag presentations. Both genera 2 and 3 are covered. The ordinary curve-degree bound remains two; on Picard-rank-greater-than-one ambients this means Mori coefficient sum at most two.

Every ordered insertion profile through **five labeled marks** is represented numerically, with insertions 1, P=tau_1(1), and H=tau_0(H_X). Chow coefficients and relations are computed in **codimensions 0–2**. A separate probe studies the elliptic-cubic ambient calculation at (g,n,r)=(3,3,3). No exhaustive all-codimension search is claimed.

The [geometry and proof report](CONVEX_GKM_GEOMETRY_2026-11-06.md) gives every new ambient, fans for the Hirzebruch factors, defining bundles, Chern numbers, fixed-point/edge data, convexity checks, and the curve-count arguments. The calculations use ordinary K3 theory, not reduced K3 theory.

## The class being calculated

For a smooth ambient X and a regular section of E cutting out Y, the ambient class is
\[
Z_{g,\beta}(I;t)=\operatorname{st}_*\left(
e_t\!\left((R\pi_*f^*E)^\vee\right)
\prod_i\psi_i^{a_i}\operatorname{ev}_i^*(\gamma_i)
\cap[\overline{\mathcal M}_{g,n}(X,\beta)]^{\mathrm{vir}}\right),
\qquad I_i=\tau_{a_i}(\gamma_i).
\]
The Euler class is the ratio for the full derived complex, with a single invertible parameter t. Convexity is a genus-zero property; the higher-genus R¹ contribution remains in this expression. The psi classes in the insertions are stable-map descendants, with the required corrections when the source curve is stabilized.

We compute the fixed-graph sum with auxiliary ambient torus weights specialized to two generic rays, lambda=epsilon w. After summing graphs, negative powers of epsilon vanish, and the epsilon⁰ term is independent of the ray. Expanding both facts in the decorated strata algebra gives rational relation vectors. Equal-insertion label symmetries give further presentations of the same class. These are the direct localization rows used below. Scalar arithmetic normalizes the common invertible t to one; homogeneity restores the power \(t^{\alpha+a-r}\), where \(\alpha=\dim(Y)(g-1)-c_1(TY)\cdot\beta\) and a is the insertion degree. This normalization does not take the limit t→0.

Write L for the full equivariant GLSM pushforward. The comparison with compact theory has the form
\[
L=Z+B,\qquad [t^q]L=0\ (q<0),\qquad [t^0]L=\sigma G,
\qquad
\sigma=(-1)^{\operatorname{rk}(E)(1-g)+c_1(E)\cdot\beta}.
\]
Here G includes the sum over curve classes of Y mapping to beta, and B denotes the aggregate infinity contribution. Put A=alpha+a. In curve codimension r this determines B_r=−Z_r for r>A and B_A=sigma G−Z_A. Positive t-coefficients, corresponding to r<A, are not determined by this compact comparison. The supplied CJR III draft treats arbitrary-rank regular sections in §9.4; the coefficientwise statement and stabilization conventions are recorded in [GENUS_FOUR_MARKED_SOURCES.md, §4](GENUS_FOUR_MARKED_SOURCES.md). Its ample-bundle finiteness statement does not apply to all of the new nef bundles, so the Mori-degree cutoff here is imposed explicitly.

The present 54-presentation driver does not compute individual positive-dimensional punctured infinity cycles or assemble their class-valued CJR graph contributions. Its infinity output consists of numerical aggregates: at negative t-power it records −Z, at t⁰ it records sigma N−Z when the compact number N is known, and otherwise it leaves the entry unknown. Knowing N supplies an integral; the driver does not construct a Chow representative of G from it. Even when G vanishes as a class, the equality B=−Z determines the specified aggregate coefficient without identifying its individual infinity-cycle constituents.

## What changed at higher markings

Four-mark cases initially evaluated all 15 canonical multisets directly. Once the exact unit-forgetting identity was tested, remaining unit profiles were reconstructed from lower markings. All five-mark cases localize the six canonical profiles with no primary units and reconstruct the other fifteen. Exact label permutations recover all ordered profiles. Relations already obtained by full direct localization are retained. Thus the search spans include available raw localization relations and the explicitly declared functorial consequences; the number of raw profiles is not identical for every four-mark checkpoint.

Writing C for a homogeneous component of Z, the reconstruction formula is
\[
\boxed{C_{g,\beta}(I_1,\ldots,I_{n-1},1)
=\pi^*C_{g,\beta}(I_1,\ldots,I_{n-1})
+\sum_{I_i=P}(s_i)_*C_{g,\beta}(I_1,\ldots,1_i,\ldots,I_{n-1}).}
\]
Here pi forgets the new marking, and s_i attaches a degree-zero three-pointed rational tail carrying i and n. In codimension r, the first term uses codimension r and the summands use codimension r−1. The formula is restricted to the stated insertion alphabet: P is tau_1(1), not tau_1(H). It follows from the unit/string identities and the cotangent boundary formula. The positive-degree two-point tail correction with a primary unit vanishes by the fundamental-class axiom; this is why no extra descendant calibration term occurs for tau_1(1). The derived twisting index is unchanged under forgetting a degree-zero bubble.

The boundary correction is essential: direct four-mark tests find it nonzero. Direct five-mark fixed-locus calculations with one and two unit insertions agree with the reconstruction. The rational-tail linear map is also checked against multiplication by the boundary divisor on every divisor generator. Numerical string reduction alone is not used to establish these class identities.

The implementation uses canonical insertion multisets: at five marks there are 21 multisets, compared with 243 ordered profiles. Six contain no primary unit. Label permutations recover the ordered classes. For relation spans, adjacent transpositions generate the full symmetric-group closure; an independent row is expanded only once. A graph is also skipped when its boundary codimension, minus forgotten tails, plus the psi degree already on stable vertices exceeds the requested codimension. These are exact reductions. The GKM graph loop computes each graph's Mori degree once. Large jobs have separate atomic checkpoints for insertion-profile partitions.

For the final linear algebra, modular row selection only chooses a candidate independent subset. A rational echelon calculation and an exact matrix identity prove that it spans every input row. The Pixton certificates are then solved and replayed over the rationals. No claim of independence or membership rests on a finite-field calculation alone.

## A dimension ceiling in genus three

Put
\[
s_X=\dim(X)(g-1)-c_1(TX)\cdot\beta.
\]
A codimension-r term with a non-unit insertions and twisting characteristic degree j≥0 has auxiliary equivariant degree s_X+a+j−r. Therefore
\[
s_X+a>r\quad\Longrightarrow\quad C_r=0
\]
and there are no nonpositive auxiliary coefficients to extract from that profile. Checking all curve classes in all 54 presentations gives a maximum non-unit marking count of **four for genus three** in codimension at most two. At five marks every genus-three non-unit coefficient is zero. Unit profiles are supplied by the preceding formula. Further raw localization presentations may still yield redundant strata equations, but the marked classes in this window contain no additional independent non-unit input.

The analogous maximum in genus two is eleven, attained among the higher-dimensional projective presentations. This is a dimension bound, not a claim that every allowed coefficient is nonzero. The computation stops at five marks after finding the flag curve families below and verifying the stated relation spans. Increasing curve degree or codimension is a more direct way to enlarge the genus-three search.

This ceiling concerns codimensions at most two. Five-mark numerical descendant integrals lie in top curve codimension eight in genus two and eleven in genus three, outside this class window. They are computed separately from numerical localization and the descendant equations, and need not vanish.

The ceiling applies to the ambient classes Z. It is not a vanishing theorem for individual infinity cycles or a stopping criterion for a full infinity-cycle elimination search. For a CY3, A=3(g−1)+a, so codimensions r≤2 in genera two and three lie entirely in the positive-t range; the compact comparison alone does not determine B_r there.

## Relations and certificates

The certified rows in this pass come from ambient torus-pole cancellation, ray independence, insertion symmetries, and the stated tautological operations. They do not come from eliminating an independently computed family of infinity cycles. Thus the Pixton membership result below is restricted to these actual rows; it does not establish membership for all relations obtainable from the full CJR infinity theory.

For each (g,n,r), the files contain rational matrices U,M,V satisfying **UM=V**, where M is the generated Pixton matrix and V is an independent basis of the found span. The final verifier checks every stored per-target basis row against V as well. Four- and five-mark spans include all label choices of unit pullbacks and rational-tail pushforwards from the preceding marking level. These are known tautological operations; independence from unit pullbacks alone is not novelty relative to the whole tautological ideal.

Two [explicit four-mark equations](../results/marked-frontier/2026-11-06/sample-relations.md) give short representatives outside the comparison span of unit pullbacks and divisor products. The adjacent JSON file records their expressions in both the found basis and the Pixton generators. Every graph and coefficient is specified literally.

| Genus | Marks | Codimension | Found rank | Pixton rank | Strata columns |
|---:|---:|---:|---:|---:|---:|
| 2 | 0 | 1 | 1 | 1 | 3 |
| 2 | 0 | 2 | 6 | 6 | 8 |
| 2 | 1 | 1 | 1 | 1 | 4 |
| 2 | 1 | 2 | 12 | 12 | 17 |
| 2 | 2 | 1 | 1 | 1 | 7 |
| 2 | 2 | 2 | 23 | 24 | 38 |
| 2 | 3 | 1 | 1 | 1 | 13 |
| 2 | 3 | 2 | 51 | 56 | 100 |
| 3 | 0 | 1 | 0 | 0 | 3 |
| 3 | 0 | 2 | 6 | 6 | 13 |
| 3 | 1 | 1 | 0 | 0 | 5 |
| 3 | 1 | 2 | 11 | 11 | 27 |
| 3 | 2 | 1 | 0 | 0 | 9 |
| 3 | 2 | 2 | 21 | 22 | 64 |
| 3 | 3 | 1 | 0 | 0 | 17 |
| 3 | 3 | 2 | 47 | 52 | 170 |
| 2 | 4 | 1 | 1 | 1 | 25 |
| 2 | 4 | 2 | 123 | 140 | 284 |
| 3 | 4 | 1 | 0 | 0 | 33 |
| 3 | 4 | 2 | 114 | 132 | 478 |
| 2 | 5 | 1 | 1 | 1 | 49 |
| 2 | 5 | 2 | 301 | 350 | 824 |
| 3 | 5 | 1 | 0 | 0 | 65 |
| 3 | 5 | 2 | 270 | 335 | 1374 |

The higher-mark comparison is:

| g | n | r | Unit-pullback rank | Beyond unit pullbacks | Unit plus divisor-product rank | Beyond that span |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 4 | 2 | 113 | 10 | 121 | 2 |
| 3 | 4 | 2 | 105 | 9 | 105 | 9 |
| 2 | 5 | 2 | 290 | 11 | 298 | 3 |
| 3 | 5 | 2 | 270 | 0 | 270 | 0 |

The separate codimension-three probe found 117 independent rows in a generated Pixton space of rank 859, with 1269 strata columns. Its exact certificate is saved separately. The new ambient presentations did not enlarge the old n≤3 spans: this is checked in `baseline-extension.json`.

## New exact flag invariants

For Y=Fl(1,2;4)[(1,1),(1,2)], the two projections identify **33 isolated P¹ fibers**, **88 isolated line fibers**, and **5 isolated conic fibers**. Convex localization and the boundary-incidence argument give **634 isolated mixed-degree curves**. Their normal bundles are O(−1)². The proof excludes other images and intersecting components in the degree-two window; genus-zero numbers alone would not establish the higher-genus formulas.

Let \(c_g=|B_{2g}|/[2g(2g-2)!]\). Then
\[
N_{g,(k,0)}=33c_g k^{2g-3},\qquad
N_{g,(0,k)}=c_g\bigl(88k^{2g-3}+\mathbf1_{2\mid k}5(k/2)^{2g-3}\bigr),
\]
and \(N_{g,(1,1)}=634c_g\). The pure-direction statements hold for every k≥1 and g≥2 by the image classification; direct localization was checked only in the declared degree-two window.

| beta | N₂ | N₃ | Z₂, ambient twist | Z₃, ambient twist |
|---|---:|---:|---:|---:|
| (1,0) | 11/80 | 11/2016 | 0 | 0 |
| (0,1) | 11/30 | 11/756 | −1/24 | 0 |
| (2,0) | 11/40 | 11/252 | 0 | 0 |
| (0,2) | 181/240 | 709/6048 | 25/288 | 1/2016 |
| (1,1) | 317/120 | 317/3024 | −325/288 | −10811/72576 |

For k psi insertions and h H insertions, the numerical multiplier is
\[
d_H^h\prod_{j=0}^{k-1}(2g-2+h+j).
\]
Five psi insertions multiply by 720 in genus two and 6720 in genus three. Thus the (1,0) five-dilaton ordinary values are **99** and **110/3**. In mixed degree (1,1) they are **1902** and **6340/9**. Primary-unit profiles integrate to zero, and there is no 1/n! for labeled marks.

The CJR sign for PF4 is sigma=(-1)^(2d₁+3d₂), and B=sigma N−Z at the numerical t⁰ coefficient. The five-dilaton infinity aggregates in degree (1,1) are therefore **−2179/2** and **16015/54**. In pure degree (1,0), Z=0 and sigma=1, so the positive ordinary numbers are supplied entirely by the infinity aggregate. These are aggregate values, not assignments to individual infinity vertices.

For the previously studied CY3 linear section of Gr(2,7), degree two with five psi insertions gives ordinary values 4851 and 9310/3 in genera two and three. The corresponding ambient coefficients are −57687 and −389102483/432, and the infinity aggregates are 52836 and 390443123/432. The signs are respectively −1 and +1. These values use the previously audited ordinary seeds; they are also stored explicitly in the new numerical files.

In contrast, the anticanonical CY3 in F₂×P² contains P¹×E over a plane cubic E; the fibers have normal O⊕O(−2). They are not isolated (−1,−1) curves. Their unknown ordinary invariants are not replaced by the ambient twisted zeros.

## Counts, validation, and memory

The numerical files contain 194,376 ordered profile records in 534 genus/curve-class sectors. They reproduce 362 previous seeds. There are 65 unknown ordinary seed sectors; nonzero-multiplier marked values remain null.

The certificate replay checked 991 independent basis rows, counted separately over (g,n,r), and 13,309 per-target rows, in addition to the 117-row higher-codimension probe. At four marks, 71 target/genus cases localized every profile directly, and 37 used unit reconstruction. At five marks, all 108 cases used direct non-unit profiles followed by unit reconstruction.

The host has 24 GiB RAM; the guest exposes 19.57 GiB. The largest recorded worker used 5.17 GiB RSS. Minimum sampled guest available memory was 1.50 GiB. Guest swap-out counters changed by 0 pages and OOM-kill counters by 0. Host compressed-memory occupancy grew from 3.48 GiB to a sampled peak of 11.45 GiB; host swap-out traffic since the initial snapshot was 14.50 GiB. Host counters include other applications.

Memory constrained concurrency: the five-mark pool was stopped at atomic checkpoints when the large four-mark Grassmannian pool left about 4.5 GiB guest memory available. It resumed with two workers once available memory recovered to about 11 GiB. The largest five-mark Grassmannian jobs wait for the corresponding four-mark checkpoint. Later, several five-mark workers were temporarily suspended with their in-memory state intact when guest headroom reached about 2 GiB, then resumed after completed partitions released their caches and headroom recovered above 5 GiB. This scheduling cost wall time. No completed results were lost. RSS sums double-count shared pages and are not physical-memory totals. The guest kernel does not expose PSI; this absence is recorded. See `memory-summary.json` and the raw host/guest logs for limits of the measurements.

Direct four/five-mark localization, unit and dilaton Chow identities, section-map identities, old three-mark and genus-four regressions, new incidence Chern calculations, normal endpoint characters, and independent genus-zero counts all pass. A fresh full-profile NF4 run reproduces all 43 stored finite class vectors and the codimension-one/two relation rows literally, using the final graph loop and CLI. Separate NF4 and PF4 class replays also agree exactly. Source changes and reports are hashed in `provenance.json`; previous versions of the shared backends are preserved in `baseline/`. Interrupted logs are retained and are not treated as completed results.

## Every presentation and seed

The following tables give all unmarked inputs. The marking formula recovers every higher-mark numerical value; the explicit 364 ordered profiles per sector are stored in `numbers/*-n5.json.gz`. A question mark means an unknown ordinary invariant.

### G1: gr2-4 [(1,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/288 | 1/72576 |
| (1,) | 133/720 | -17/1792 | 0 | 0 |
| (2,) | 7/144 | -48257/72576 | 0 | 0 |

### G10: gr2-5 [(1,); (1,); (2,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/144 | 1/36288 |
| (1,) | 295/144 | 0 | 0 | 0 |
| (2,) | -217/96 | -616169/32256 | 0 | 0 |

### G11: gr2-5 [(1,); (1,); (3,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -5/192 | 5/48384 |
| (1,) | -25/32 | 0 | 9/4 | 5/56 |
| (2,) | -1834/9 | -1062577/96768 | ? | ? |

### G12: gr2-5 [(1,); (2,); (2,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/48 | 1/12096 |
| (1,) | 955/1152 | 0 | 5/3 | 25/378 |
| (2,) | -32597/1152 | 10891/145152 | ? | ? |

### G13: gr2-5 [(1,); (1,); (1,); (1,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 95/18 | 0 | 0 | 0 |
| (2,) | -25/288 | -21215507/72576 | 0 | 0 |

### G14: gr2-5 [(1,); (1,); (1,); (2,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 3875/1152 | 0 | 0 | 0 |
| (2,) | 14509/1152 | -8673817/193536 | 0 | 0 |

### G15: gr2-7 [(1,); (1,); (1,); (1,); (1,); (1,); (1,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -49/2880 | 7/103680 |
| (1,) | 1351/288 | 0 | 49/60 | 7/216 |
| (2,) | -19229/240 | -55586069/414720 | 539/80 | 133/288 |

### G2: gr2-4 [(2,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/240 | 1/60480 |
| (1,) | 91/288 | -661/30240 | 0 | 0 |
| (2,) | 0 | -7451/7560 | 0 | 0 |

### G3: gr2-4 [(3,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/96 | 1/24192 |
| (1,) | 0 | -6791/241920 | 0 | 0 |
| (2,) | -3/8 | -7/384 | 0 | 0 |

### G4: gr2-4 [(4,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -11/360 | 11/90720 |
| (1,) | -22/45 | -163/6720 | 16/3 | 40/189 |
| (2,) | -887/4 | -568469/362880 | ? | ? |

### G5: gr2-4 [(1,); (1,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 83/120 | -23389/725760 | 0 | 0 |
| (2,) | 1/360 | -60611/8640 | 0 | 0 |

### G6: gr2-4 [(1,); (2,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 49/72 | -9157/181440 | 0 | 0 |
| (2,) | -53/720 | -1173841/362880 | 0 | 0 |

### G7: gr2-4 [(1,); (3,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | -181/1440 | -41653/725760 | 0 | 0 |
| (2,) | 2717/720 | -28099/40320 | 0 | 0 |

### G8: gr2-4 [(2,); (2,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 11/30 | -17039/241920 | 0 | 0 |
| (2,) | -73/30 | -326/405 | 0 | 0 |

### G9: gr2-5 [(1,); (1,); (1,)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/288 | 1/72576 |
| (1,) | 1285/384 | 0 | 0 | 0 |
| (2,) | -11/128 | -25633969/193536 | 0 | 0 |

### NF1: f2xp1 [(1, 2, 0)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | -1/360 | 1/90720 | 0 | 0 |
| (0, 0, 1) | 0 | 1/9072 | 0 | 0 |
| (0, 1, 0) | -1/360 | -11/90720 | 0 | 0 |
| (1, 0, 0) | 0 | 0 | 0 | 0 |
| (0, 0, 2) | 0 | 0 | 0 | 0 |
| (0, 1, 1) | 0 | 17/22680 | 0 | 0 |
| (0, 2, 0) | -1/288 | -19/60480 | 0 | 0 |
| (1, 0, 1) | 0 | 0 | 0 | 0 |
| (1, 1, 0) | -1/360 | -11/90720 | 0 | 0 |
| (2, 0, 0) | 0 | 0 | 0 | 0 |

### NF2: f2xp1 [(2, 4, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | -1/360 | 1/90720 | 0 | 0 |
| (0, 0, 1) | 1/90 | 1/22680 | 0 | 0 |
| (0, 1, 0) | 1/90 | 1/22680 | 0 | 0 |
| (1, 0, 0) | 0 | 0 | 0 | 0 |
| (0, 0, 2) | 1/45 | 1/2268 | 0 | 0 |
| (0, 1, 1) | -2/9 | -191/7560 | 0 | 0 |
| (0, 2, 0) | 1/45 | 1/2268 | 0 | 0 |
| (1, 0, 1) | 0 | 0 | 0 | 0 |
| (1, 1, 0) | 1/90 | 1/22680 | 0 | 0 |
| (2, 0, 0) | 0 | 0 | 0 | 0 |

### NF3: f2xp2 [(2, 4, 3)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | 0 | 0 | -1/40 | 1/10080 |
| (0, 0, 1) | 1/16 | 1/1440 | ? | ? |
| (0, 1, 0) | -1/40 | 1/10080 | ? | ? |
| (1, 0, 0) | 0 | 0 | ? | ? |
| (0, 0, 2) | -3/16 | -1/160 | ? | ? |
| (0, 1, 1) | 51/80 | -17/1260 | ? | ? |
| (0, 2, 0) | -1/20 | 1/1008 | ? | ? |
| (1, 0, 1) | 0 | 0 | ? | ? |
| (1, 1, 0) | -1/40 | 1/10080 | ? | ? |
| (2, 0, 0) | 0 | 0 | ? | ? |

### NF4: f3xp1 [(1, 3, 0)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | -1/360 | 1/90720 | 0 | 0 |
| (0, 0, 1) | 0 | 23/181440 | 0 | 0 |
| (0, 1, 0) | -1/288 | -13/90720 | 0 | 0 |
| (1, 0, 0) | 0 | 0 | 0 | 0 |
| (0, 0, 2) | 0 | 0 | 0 | 0 |
| (0, 1, 1) | 0 | 7/8640 | 0 | 0 |
| (0, 2, 0) | -11/2880 | -41/120960 | 0 | 0 |
| (1, 0, 1) | 0 | 0 | 0 | 0 |
| (1, 1, 0) | 0 | 1/45360 | 0 | 0 |
| (2, 0, 0) | 0 | 0 | 0 | 0 |

### NF5: f3xp2 [(1, 3, 0)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | 0 | 0 | -1/320 | 1/80640 |
| (0, 0, 1) | 0 | -13/20160 | 0 | 0 |
| (0, 1, 0) | -1/80 | 1/10080 | 0 | 0 |
| (1, 0, 0) | 0 | 0 | 0 | 0 |
| (0, 0, 2) | 0 | 0 | 0 | 0 |
| (0, 1, 1) | 0 | 437/322560 | 0 | 0 |
| (0, 2, 0) | -1/320 | 1/4480 | 0 | 0 |
| (1, 0, 1) | 0 | 0 | 0 | 0 |
| (1, 1, 0) | 7/640 | 0 | 0 | 0 |
| (2, 0, 0) | 0 | 0 | 0 | 0 |

### PF1: flag3-1-2 [(1, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | -1/320 | 1/80640 | 0 | 0 |
| (0, 1) | -3/640 | -1/5376 | 0 | 0 |
| (1, 0) | -3/640 | -1/5376 | 0 | 0 |
| (0, 2) | -3/640 | -11/26880 | 0 | 0 |
| (1, 1) | 73/320 | -8117/161280 | 0 | 0 |
| (2, 0) | -3/640 | -11/26880 | 0 | 0 |

### PF2: flag3-1-2 [(2, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | -1/320 | 1/80640 | 0 | 0 |
| (0, 1) | 1/80 | 1/20160 | 0 | 0 |
| (1, 0) | 1/80 | 1/20160 | 0 | 0 |
| (0, 2) | 1/40 | 1/2016 | 0 | 0 |
| (1, 1) | -87/160 | -559/4480 | 0 | 0 |
| (2, 0) | 1/40 | 1/2016 | 0 | 0 |

### PF3: flag4-1-2 [(1, 0); (0, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -1/320 | 1/80640 |
| (0, 1) | -1/48 | 0 | 0 | 0 |
| (1, 0) | 0 | 0 | 0 | 0 |
| (0, 2) | -1/288 | 1/2016 | 0 | 0 |
| (1, 1) | 4111/1440 | -6743/60480 | 0 | 0 |
| (2, 0) | 0 | 0 | 0 | 0 |

### PF4: flag4-1-2 [(1, 1); (1, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -11/576 | 11/145152 |
| (0, 1) | -1/24 | 0 | 11/30 | 11/756 |
| (1, 0) | 0 | 0 | 11/80 | 11/2016 |
| (0, 2) | 25/288 | 1/2016 | 181/240 | 709/6048 |
| (1, 1) | -325/288 | -10811/72576 | 317/120 | 317/3024 |
| (2, 0) | 0 | 0 | 11/40 | 11/252 |

### PF5: flag4-1-3 [(1, 1); (1, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -1/240 | 1/60480 |
| (0, 1) | -1/24 | 0 | 0 | 0 |
| (1, 0) | -1/24 | 0 | 0 | 0 |
| (0, 2) | 0 | 5/4032 | 0 | 0 |
| (1, 1) | 161/120 | -42481/30240 | 0 | 0 |
| (2, 0) | 0 | 5/4032 | 0 | 0 |

### PF6: flag4-1-3 [(1, 1); (2, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -1/45 | 1/11340 |
| (0, 1) | -1/24 | 0 | ? | ? |
| (1, 0) | -1/24 | 0 | ? | ? |
| (0, 2) | 25/288 | 1/2016 | ? | ? |
| (1, 1) | 2057/720 | -22159/30240 | ? | ? |
| (2, 0) | 25/288 | 1/2016 | ? | ? |

### T1: p1xp2 [(1, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | -1/320 | 1/80640 | 0 | 0 |
| (0, 1) | 7/320 | 13/241920 | 0 | 0 |
| (1, 0) | -7/1920 | -37/241920 | 0 | 0 |
| (0, 2) | -1/192 | 901/483840 | 0 | 0 |
| (1, 1) | 43/384 | -8237/96768 | 0 | 0 |
| (2, 0) | -1/240 | -1/2688 | 0 | 0 |

### T10: f1xp1 [(2, 3, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | -1/360 | 1/90720 | 0 | 0 |
| (0, 0, 1) | 1/90 | 1/22680 | 0 | 0 |
| (0, 1, 0) | 1/90 | 1/22680 | 0 | 0 |
| (1, 0, 0) | 0 | 25/290304 | 0 | 0 |
| (0, 0, 2) | 1/45 | 1/2268 | 0 | 0 |
| (0, 1, 1) | -2/9 | -191/7560 | 0 | 0 |
| (0, 2, 0) | 1/45 | 1/2268 | 0 | 0 |
| (1, 0, 1) | 0 | -37/483840 | 0 | 0 |
| (1, 1, 0) | 1/10 | -151/41472 | 0 | 0 |
| (2, 0, 0) | 0 | -367/725760 | 0 | 0 |

### T11: f1xp2 [(1, 1, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | 0 | 0 | -1/360 | 1/90720 |
| (0, 0, 1) | 1/960 | -1/40320 | 0 | 0 |
| (0, 1, 0) | -1/80 | 1/10080 | 0 | 0 |
| (1, 0, 0) | 7/640 | 0 | 0 | 0 |
| (0, 0, 2) | -11/1920 | -13/16128 | 0 | 0 |
| (0, 1, 1) | 63/640 | -53/3584 | 0 | 0 |
| (0, 2, 0) | -1/320 | 1/4480 | 0 | 0 |
| (1, 0, 1) | 1/96 | -53/15120 | 0 | 0 |
| (1, 1, 0) | 11/1920 | 1/16128 | 0 | 0 |
| (2, 0, 0) | 1/480 | 1/2240 | 0 | 0 |

### T12: f1xp2 [(2, 3, 3)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | 0 | 0 | -1/40 | 1/10080 |
| (0, 0, 1) | 1/16 | 1/1440 | 7/10 | 1/36 |
| (0, 1, 0) | -1/40 | 1/10080 | 9/40 | 1/112 |
| (1, 0, 0) | 7/640 | 0 | 3/80 | 1/672 |
| (0, 0, 2) | -3/16 | -1/160 | ? | ? |
| (0, 1, 1) | 51/80 | -17/1260 | ? | ? |
| (0, 2, 0) | -1/20 | 1/1008 | 9/20 | 1/14 |
| (1, 0, 1) | -1/160 | -19/16128 | ? | ? |
| (1, 1, 0) | 57/640 | 29/20160 | ? | ? |
| (2, 0, 0) | -7/320 | 1/2240 | 3/40 | 1/84 |

### T2: p1xp2 [(2, 3)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | -1/320 | 1/80640 | 0 | 0 |
| (0, 1) | 11/160 | -121/32256 | 0 | 0 |
| (1, 0) | 1/80 | 1/20160 | 0 | 0 |
| (0, 2) | -33/160 | 121/3584 | 0 | 0 |
| (1, 1) | 109/160 | -3167/161280 | 0 | 0 |
| (2, 0) | 1/40 | 1/2016 | 0 | 0 |

### T3: p1xp3 [(1, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -1/288 | 1/72576 |
| (0, 1) | 23/180 | -53/20160 | 0 | 0 |
| (1, 0) | -1/72 | 1/9072 | 0 | 0 |
| (0, 2) | 0 | -79/3780 | 0 | 0 |
| (1, 1) | 0 | -3889/9072 | 0 | 0 |
| (2, 0) | -1/288 | 1/4032 | 0 | 0 |

### T4: p1xp3 [(2, 4)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -7/240 | 1/8640 |
| (0, 1) | -37/180 | -19/3024 | 8/3 | 20/189 |
| (1, 0) | -1/36 | 1/9072 | 4/15 | 2/189 |
| (0, 2) | -229/18 | 2/63 | ? | ? |
| (1, 1) | -338/45 | 41/1260 | ? | ? |
| (2, 0) | -1/18 | 5/4536 | 8/15 | 16/189 |

### T5: p2xp2 [(1, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -1/320 | 1/80640 |
| (0, 1) | 3/80 | 1/2560 | 0 | 0 |
| (1, 0) | -1/640 | -1/10752 | 0 | 0 |
| (0, 2) | 0 | 1/630 | 0 | 0 |
| (1, 1) | 0 | -19941/35840 | 0 | 0 |
| (2, 0) | -1/128 | -17/13440 | 0 | 0 |

### T6: p2xp2 [(3, 3)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | -9/320 | 1/8960 |
| (0, 1) | 9/128 | 1/1280 | 63/80 | 1/32 |
| (1, 0) | 9/128 | 1/1280 | 63/80 | 1/32 |
| (0, 2) | -27/128 | -9/1280 | ? | ? |
| (1, 1) | -6177/640 | 151/5120 | ? | ? |
| (2, 0) | -27/128 | -9/1280 | ? | ? |

### T7: p2xp2 [(1, 1); (1, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | 0 | 0 |
| (0, 1) | 11/640 | -1/26880 | 0 | 0 |
| (1, 0) | 11/640 | -1/26880 | 0 | 0 |
| (0, 2) | -1/320 | 1/5376 | 0 | 0 |
| (1, 1) | 41/160 | -4145/1792 | 0 | 0 |
| (2, 0) | -1/320 | 1/5376 | 0 | 0 |

### T8: p2xp2 [(1, 1); (2, 2)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0) | 0 | 0 | 0 | 0 |
| (0, 1) | 5/64 | 1/2688 | 0 | 0 |
| (1, 0) | 5/64 | 1/2688 | 0 | 0 |
| (0, 2) | -13/80 | 23/10752 | 0 | 0 |
| (1, 1) | 39/40 | -196723/161280 | 0 | 0 |
| (2, 0) | -13/80 | 23/10752 | 0 | 0 |

### T9: f1xp1 [(1, 1, 1)]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0, 0, 0) | -1/360 | 1/90720 | 0 | 0 |
| (0, 0, 1) | -1/288 | -13/90720 | 0 | 0 |
| (0, 1, 0) | -1/288 | -13/90720 | 0 | 0 |
| (1, 0, 0) | 7/1152 | 101/362880 | 0 | 0 |
| (0, 0, 2) | -11/2880 | -41/120960 | 0 | 0 |
| (0, 1, 1) | 47/1440 | -2689/725760 | 0 | 0 |
| (0, 2, 0) | -11/2880 | -41/120960 | 0 | 0 |
| (1, 0, 1) | 1/160 | -1021/2903040 | 0 | 0 |
| (1, 1, 0) | 67/2880 | 389/1451520 | 0 | 0 |
| (2, 0, 0) | 1/1440 | 1/2880 | 0 | 0 |

### p2-3: P^2[3]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | -1/320 | 1/40320 | 0 | 0 |
| (1,) | 29/640 | -1637/80640 | 0 | 0 |
| (2,) | -87/640 | 1637/8960 | 0 | 0 |

### p3-2: P^3[2]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | -1/288 | 1/72576 | 0 | 0 |
| (1,) | 17/120 | -179/6048 | 0 | 0 |
| (2,) | 0 | -491/3780 | 0 | 0 |

### p3-3: P^3[3]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | -1/288 | 1/72576 | 0 | 0 |
| (1,) | 1/72 | -21293/362880 | 0 | 0 |
| (2,) | 0 | 0 | 0 | 0 |

### p3-4: P^3[4]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | -1/288 | 1/72576 | 0 | 0 |
| (1,) | -17/72 | -491/7560 | 0 | 0 |
| (2,) | -67/36 | 5183/7560 | 0 | 0 |

### p4-2: P^4[2]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/288 | 1/72576 |
| (1,) | 53/192 | -22249/290304 | 0 | 0 |
| (2,) | 0 | -1541/3024 | 0 | 0 |

### p4-3: P^4[3]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/192 | 1/48384 |
| (1,) | 25/1152 | -1087/10752 | 0 | 0 |
| (2,) | 0 | 0 | 0 | 0 |

### p4-4: P^4[4]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/72 | 1/18144 |
| (1,) | 65/1152 | -21017/290304 | 0 | 0 |
| (2,) | -25/16 | -1525/13824 | 0 | 0 |

### p4-5: P^4[5]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -5/144 | 5/36288 |
| (1,) | 1645/1152 | -3265/290304 | 575/48 | 2875/6048 |
| (2,) | -351685/144 | -4368535/72576 | ? | ? |

### p5-2-4: P^5[2,4]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -11/360 | 11/90720 |
| (1,) | 467/480 | -12091/69120 | 16/3 | 40/189 |
| (2,) | 78323/144 | 948043/362880 | ? | ? |

### p5-2: P^5[2]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 343/720 | -949/6912 | 0 | 0 |
| (2,) | 0 | -10009/6480 | 0 | 0 |

### p5-3-3: P^5[3,3]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/40 | 1/10080 |
| (1,) | 49/64 | -1637/5760 | 351/80 | 39/224 |
| (2,) | 44011/160 | 323153/53760 | ? | ? |

### p5-3: P^5[3]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 1/32 | -737/5760 | 0 | 0 |
| (2,) | 0 | 0 | 0 | 0 |

### p5-4: P^5[4]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 29/96 | -121/1728 | 0 | ? |
| (2,) | 0 | -89/192 | 0 | 0 |

### p5-5: P^5[5]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | 0 | 0 |
| (1,) | 3509/576 | -397/6912 | ? | 0 |
| (2,) | 9361/4 | -136723687/82944 | 0 | ? |

### p6-2-2-3: P^6[2,2,3]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/40 | 1/10080 |
| (1,) | 469/480 | -609/1280 | 3 | 5/42 |
| (2,) | -42593/192 | -1152557/138240 | ? | ? |

### p7-2-2-2-2: P^7[2,2,2,2]

| beta | Z₂ | Z₃ | N₂ | N₃ |
|---|---|---|---|---|
| (0,) | 0 | 0 | -1/45 | 1/11340 |
| (1,) | 181/180 | -221/864 | 32/15 | 16/189 |
| (2,) | 16933/90 | -53243/22680 | ? | ? |

## Reproduction and references

Run inside the pinned Sage Linux environment, from the repository root. Completed checkpoints are reused. The initial four-mark full-profile checkpoints are retained; current code uses exact unit reconstruction for new runs.

```bash
sage tests/test_marked_frontier.sage
sage tests/test_convex_gkm_extensions.sage
sage tests/test_marked_forgetting.sage
python calculations/marked_frontier.py --markings 4 --workers 4
python calculations/marked_frontier.py --markings 5 --workers 2
python calculations/frontier_baseline.py
python calculations/frontier_certify.py --markings 4
python calculations/frontier_certify.py --markings 5
python calculations/frontier_numbers.py --max-markings 5
python calculations/frontier_verify_numbers.py
python calculations/frontier_codimension_probe.py
python calculations/frontier_validate_units.py
python calculations/frontier_gkm_replay.py
python calculations/frontier_verify.py
python calculations/frontier_examples.py
python calculations/frontier_resources.py
python calculations/frontier_report.py
```

For a clean checkout, the eleven added presentations also require the n=0,1,2,3 runs with `--indices 43 44 45 46 47 48 49 50 51 52 53` before the four-mark reconstruction. The previous forty-three n≤3 artifacts are inputs. A two-ray numerical seed check is carried forward for these previous targets; the new targets are recomputed on both rays.

To reproduce a full direct four-mark checkpoint, including its unit profiles, use `--include-unit-profiles --direct-only`. A fresh `MARKED_FRONTIER_OUT` environment variable selects a separate output directory. For example, `MARKED_FRONTIER_OUT=results/frontier-replay python calculations/marked_frontier.py --markings 4 --indices 46 --genera 2 --workers 1 --include-unit-profiles --direct-only` recomputes NF4 without reading its existing checkpoint. The saved `profiles` and `direct_profiles` fields specify which method was used for each original case. Full direct and reconstructed representatives need not give identical raw relation spans, so this method metadata is part of the result.

- [Geometry, convexity, fans, and curve-count proofs](CONVEX_GKM_GEOMETRY_2026-11-06.md), including the supplied CJR III draft locators and the nonample finiteness caveat.
- [Coates–Givental, Quantum Riemann–Roch, Lefschetz and Serre](https://math.berkeley.edu/~giventh/papers/rls.pdf), appendices 1–2, for twisting, forgetting, and ancestor/descendant conventions.
- [Faber–Pandharipande, Hodge integrals and Gromov–Witten theory](https://arxiv.org/abs/math/9810173), Theorems 3–4, for local rational curves and constant maps.
- [Graber–Pandharipande, Localization of virtual classes](https://arxiv.org/abs/alg-geom/9708001), and [Liu–Sheshmani, Equivariant Gromov–Witten invariants of algebraic GKM manifolds](https://sigma-journal.com/2017/048/sigma17-048.pdf), for the graph formulas, including negative normal degrees.
- [Janda, Relations via equivariant GW theory of P¹](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf), for Chow validity of the Pixton generators; [Janda, Frobenius manifolds near the discriminant](https://arxiv.org/abs/1505.03419), for the semisimple reconstruction comparison. The actual membership claims here are certified by the saved matrices.
