The new target is the general linear section
\[
Y=\operatorname{Gr}(2,7)\cap\mathbb P^{13}
\subset\mathbb P(\wedge^2\mathbb C^7)=\mathbb P^{20}.
\]
Equivalently, it is the zero locus of seven general sections of the Plücker line bundle, with \(E=\mathcal O(1)^{\oplus7}\). The computation treats connected unmarked stable maps, genera two and three, and Plücker degrees \(D=0,1,2\). The former fourteen-target outputs are preserved; the extension is stored in [results/grassmannian-27/2026-11-06](../results/grassmannian-27/2026-11-06).

The bounded scan gives **20 independent relation vectors: 19 in genus two (ranks 1, 6, 12 in codimensions 1, 2, 3) and one in genus three (codimension 6)**. Every vector has an exact Pixton certificate and an explicit expression in the previous 206-dimensional Grassmannian/toric span. Thus it adds zero independent directions to that span. The genus-three relation is the difference of finite degree-two top coefficients at the two auxiliary specializations; it has 377 nonzero coefficients in the saved strata basis. The final exact scan completed in approximately 185 seconds. The full ordinary numerical invariants in the stated genus/degree bound are now determined, including degree two, by the audited line-and-conic argument below.

**Geometry and independent Chern checks.** Write \(X=\operatorname{Gr}(2,7)\), \(H=c_1(S^\vee)\), and \(P=c_2(S^\vee)\). Since \(\dim X=10\) and \(-K_X=7H\), adjunction gives \(\dim Y=3\) and \(K_Y\cong\mathcal O_Y\). Bertini gives smoothness. The Koszul resolution, Kodaira vanishing, and Serre duality give \(h^0(\mathcal O_Y)=h^3(\mathcal O_Y)=1\), \(h^1(\mathcal O_Y)=h^2(\mathcal O_Y)=0\). Lefschetz gives \(h^{1,1}(Y)=1\). All seven positive Plücker degrees must equal one: the target-dimension condition requires seven equations, and the CY condition requires their degrees to sum to seven.

With roots \(x,y\) of \(S^\vee\),
\[
c(TX)=\frac{(1+x)^7(1+y)^7}{1-(x-y)^2},\qquad
c(TY)=\frac{c(TX)}{(1+H)^7}\bigg|_Y,
\]
which gives
\[
c_1(TY)=0,\qquad c_2(TY)=H^2+3P,\qquad c_3(TY)=-7HP.
\]
The coefficient formula
\[
\int_Xf(x,y)=[x^6y^5](x-y)f(x,y)
\]
gives \(\int_XH^{10}=42\) and \(\int_XH^8P=14\). Since \([Y]=H^7\),
\[
\boxed{H_Y^3=42,\qquad c_2(Y)\cdot H_Y=84,\qquad\chi(Y)=-98.}
\]
Thus \((h^{1,1},h^{2,1})=(1,50)\). The values agree with [Haghighat–Klemm, §5.4 and Appendix A](https://arxiv.org/abs/0802.2908). Here \(P\) is the second Chern class of the rank-two subbundle dual; it should not be confused with the second Chern class of the rank-five quotient in other conventions.

**Grassmannian localization.** The ambient is non-toric. Its diagonal torus fixed points are the 21 coordinate two-planes. Invariant lines connect pairs of two-subsets differing by one element, giving the Johnson graph \(J(7,2)\) with 105 edges. At a fixed plane \(I\), tangent weights are \(\lambda_j-\lambda_i\), for \(i\in I,j\notin I\), and the Plücker fiber has weight \(-\sum_{i\in I}\lambda_i\). Along every invariant line,
\[
TX|_\ell=\mathcal O(2)\oplus\mathcal O(1)^5\oplus\mathcal O^4.
\]
These follow from \(TX=S^\vee\otimes Q\), and specialize the general GKM formulas of [Liu–Sheshmani, Examples 2.3/2.5, Lemma 4.5, and Theorem 4.7](https://sigma-journal.com/2017/048/sigma17-048.pdf). The code retains the four degree-zero normal summands and the seven separate defining summands, using the same logarithmic parameter \(t\) for all of them.

The decorated stable-map graph counts are:

| Genus | Degree 0 | Degree 1 | Degree 2 |
|---|---:|---:|---:|
| 2 | 21 | 315 | 7035 |
| 3 | 21 | 420 | 11445 |

Degree-two graphs include one double-covered invariant line, a chain of two invariant lines, and two parallel edges forming a graph cycle. The last case has vertex genera summing to \(g-1\); omitting it would change the answer. Auxiliary weights are \(\lambda_i=\epsilon b^{i-1}\), with \(b=4,5\), and all denominators are checked for resonance. Graph summation and actual stabilization to \(\overline{\mathcal M}_g\) precede the auxiliary nonequivariant limit.

The larger ambient exposed a performance issue in adding rational functions with denominators from 21 fixed points. The final scan instead expands each graph summand through exponent zero and then sums its exact rational Laurent coefficients. For a summand \(f=\epsilon^v a(\epsilon)/b(\epsilon)\), with \(a(0)b(0)\ne0\), terms with \(v>0\) cannot contribute. Otherwise only the Taylor jet of \(a/b\) through degree \(-v\) is needed. The denominator has nonzero constant term, so this finite jet is computed exactly; there is no floating arithmetic or uncontrolled truncation. The code checks that the requested jet precision is available.

This method agrees coefficient by coefficient with the original rational-function method for Gr(2,4), genus three, degree two, and for Gr(2,7), genus two, degrees one and two, in every codimension and every nonpositive Laurent power. The independent regression is saved as `laurent-regression.log`; the original genus-two scan is retained as `original-genus-two.json`. The earlier fourteen-target default behavior remains unchanged. The original large rational-function scan was interrupted after its first genus-three specialization; `scan-graphwise.log` records the completed final scan.

**Which classes are being computed.** The ambient dual-Euler term is
\[
Z_D(t)=\operatorname{st}_{X,*}
\left(e_t((R\pi_*f^*E)^\vee)
[\overline{\mathcal M}_g(X,D)]^{\rm vir}\right)
=\sum_{r=0}^{3g-3}t^{3g-3-r}Z_{D,r}.
\]
For this CY3, \(\alpha=3g-3\) at every degree. There is no negative-\(t\) range after unmarked stabilization. If \(G_{g,D}\) denotes the ordinary target class and \(B_D\) the full CJR infinity aggregate, the top comparison is
\[
\boxed{B_{D,3g-3}=(-1)^{7(1-g+D)}G_{g,D}-Z_{D,3g-3}.}
\]
The signs for \(D=0,1,2\) are \((-,+,-)\) in genus two and \((+,-,+)\) in genus three. This uses the supplied CJR III draft, Theorem 8.4, Remark 8.5, and §9.4, under the same conventions as the [Grassmannian source audit](GRASSMANNIAN_COMPLETE_INTERSECTION_SOURCES.md). The public arbitrary-rank fields comparison is [Chen–Janda–Webb, Corollary 1.2.1](https://arxiv.org/abs/1911.09825v2). The infinity target is \(\mathbb P_X(E^\vee)\), with fiber \(\mathbb P^6\); no product of seven independent hypersurface infinity theories is substituted for it.

The ambient virtual dimension is
\[
\nu_X=(1-g)(10-3)+7D=7(D+1-g).
\]
It is negative at degree zero in both genera, and at genus three, degree one. Thus the entire nonequivariant \(Z_D(t)\) vanishes in those cases, even though ordinary target invariants can be nonzero. At \((g,D)=(2,1),(3,2)\), \(\nu_X=0\) and the twist has virtual rank zero; only the top stabilized coefficient can survive, and its integral equals the corresponding untwisted ambient GW number. These statements are checked against the computed finite coefficients.

**Ordinary genus-zero checks.** Lines in \(X\) are parameterized by \(F(1,3;7)\), of dimension fourteen. For a flag \(A\subset B\), the universal line is \(\mathbb P(B/A)\), and a Plücker hyperplane restricts to a section of the rank-two bundle \(A^\vee\otimes(B/A)^\vee\). Therefore
\[
L=\int_{F(1,3;7)}c_{14}\left((A^\vee\otimes(B/A)^\vee)^{\oplus7}\right)=196.
\]
This flag-variety intersection was evaluated at two independent choices of weights, separately from stable-map graph localization. Genus-zero convex quantum Lefschetz then gives
\[
\boxed{N_{0,1}=196,\qquad N_{0,2}=2499/2.}
\]
Both were computed at both auxiliary specializations. The degree-two invariant includes double covers of lines. Removing their contribution gives
\[
\boxed{N_{0,2}-N_{0,1}/8=1225.}
\]
These values agree with [Haghighat–Klemm, Table B.5](https://arxiv.org/abs/0802.2908). The localization computation, rather than that paper's higher-genus mirror predictions, supplies the ordinary genus-zero values here. Convex Lefschetz is [Chen–Janda–Webb, Theorem 1.2.3](https://arxiv.org/abs/1911.09825v2).

**Why degree two is also determined in higher genus.** This extension closes a geometric question that was left open in the preceding bounded report. The complete incidence and obstruction argument is in [GRASSMANNIAN_27_SOURCES.md](GRASSMANNIAN_27_SOURCES.md). It uses the genus-zero value computed above, not higher-genus mirror predictions.

First, the open Hilbert scheme of smooth embedded conics in the homogeneous ambient is smooth of dimension 21. Indeed, \(TX\) is globally generated, so the normal bundle of an embedded conic is globally generated and has no \(H^1\); its rank is nine and its degree is twelve. Seven Plücker equations restrict to a globally generated bundle of rank \(7\cdot3=21\) on this parameter space. General equations therefore cut out a finite transverse set of smooth conics. The parameter space need not be proper: the incidence is finite type, and a zero-dimensional finite-type fiber is finite.

Second, the family of two distinct meeting lines has dimension 20. A marked ambient line has dimension 15, and the evaluation fiber product over the ten-dimensional ambient has dimension \(15+15-10=20\). Their union is a nodal plane conic; each hyperplane imposes three independent restrictions. Seven equations impose rank 21, so a general target contains no such pair. The incidence of bad pairs and equations has dimension 146 inside the 147-dimensional equation parameter space. The diagonal of coincident lines is excluded from this dimension argument and is retained as multiple-cover geometry.

For each line or smooth conic in the general target, transversality gives \(H^0(N_{C/Y})=0\). Adjunction gives rank two and degree \(-2\), hence
\[
N_{C/Y}=\mathcal O_{\mathbb P^1}(-1)\oplus\mathcal O_{\mathbb P^1}(-1).
\]
For a conic these degrees are intrinsic to its normalization \(\mathbb P^1\); the Plücker bundle has degree two. A planar embedded double line cannot occur: it would require a surjection \(N_{L/Y}^{\vee}=\mathcal O(1)^2\to\mathcal O(-1)\). Stable-map double covers of reduced lines are different objects and are included.

Every connected stable map of total degree two consequently has image a line or a smooth conic. Distinct images could only be two meeting lines, already excluded. Contracted components cannot connect disjoint images. For a positive-degree stable map to one of these supports, \(H^0(f^*N_{C/Y})=0\): sections vanish on the positive-degree components and then on all contracted components by connectedness. The compatible normal obstruction theory identifies the ordinary virtual contribution with the local \((-1,-1)\)-curve integral, including contracted tails and reducible domains. No assertion that different conics are pairwise disjoint is needed.

Thus each smooth conic contributes one to \(N_{0,2}\), and each line contributes \(1/8\). There are exactly 1225 smooth conics. The all-genus local-curve formula of [Faber–Pandharipande, Theorem 3](https://arxiv.org/abs/math/9810173) is
\[
C(g,k)=\frac{|B_{2g}|}{2g(2g-2)!}k^{2g-3},\qquad
C(2,k)=k/240,\quad C(3,k)=k^3/6048.
\]
Therefore
\[
N_{g,1}=196C(g,1),\qquad
N_{g,2}=1225C(g,1)+196C(g,2),
\]
which proves
\[
\boxed{N_{2,1}=49/60,\quad N_{2,2}=539/80,\qquad
N_{3,1}=7/216,\quad N_{3,2}=133/288.}
\]
These are full connected ordinary scalar invariants. The proof applies to general smooth sections, and deformation invariance extends the resulting numbers within this smooth family. Only this newly audited case is added to the code's degree-two ordinary inputs; the earlier targets' recorded unknowns are not silently overwritten.


At degree zero, the constant-map formula gives
\[
N_{2,0}=\frac{\chi(Y)}{5760}=-\frac{49}{2880},\qquad
N_{3,0}=-\frac{\chi(Y)}{1451520}=\frac7{103680}.
\]
The corresponding compact classes are \(-98\lambda_1\lambda_2\) and \(98\lambda_1\lambda_2\lambda_3\). These must not be replaced by the vanishing ambient constant terms. The reference is [Faber–Pandharipande, Theorem 4 and §4.3](https://arxiv.org/abs/math/9810173).

All six integrated CJR comparisons are therefore explicit:

| Genus | Plücker degree | Ordinary N | Ambient twisted top ∫Z | Sign σ | Infinity aggregate ∫B |
|---|---:|---:|---:|---:|---:|
| 2 | 0 | -49/2880 | 0 | -1 | 49/2880 |
| 2 | 1 | 49/60 | 1351/288 | 1 | -5579/1440 |
| 2 | 2 | 539/80 | -19229/240 | -1 | 4403/60 |
| 3 | 0 | 7/103680 | 0 | 1 | 7/103680 |
| 3 | 1 | 7/216 | 0 | -1 | -7/216 |
| 3 | 2 | 133/288 | -55586069/414720 | 1 | 55777589/414720 |


**Relation certificates and scope.** The pure tautological relation vectors come from negative auxiliary \(\epsilon\)-coefficients and differences of finite parts between the two auxiliary specializations. This is distinct from deleting the ordinary target term in the CJR comparison. Each vector \(v\) has an exact rational certificate \(uM=v\), with \(M\) the Pixton matrix in the saved strata basis. The export replays every certificate and compares to the preceding 206-dimensional Grassmannian/toric span, saving explicit expressions in that basis whenever contained. [Janda, Theorem 1.1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf) supplies the Chow validity of the Pixton system; no completeness conjecture is used.

The mathematical identities involving \(B\) determine the full stabilized infinity aggregate, or its integral. They do not isolate individual effective punctured invariants. The ordinary values in the numerical tables are scalar invariants; the code does not silently reconstruct a full target Chow class from its integral.

**Reproduction.** In the Linux Sage environment, from the project directory:

```bash
python -c 'from sage.all import *; load("tests/test_grassmannian_27.sage")'
python -c 'from sage.all import *; load("tests/test_grassmannian_laurent.sage")'
python calculations/grassmannian_27_scan.py
python calculations/grassmannian_27_export.py
```

The target constructor is:

```python
load('grassmannian_complete_intersections.sage')
Y = grassmannian_seven_hyperplane_cy3()
C = GrassmannianCompleteIntersectionCycles(Y, genus=3)
Z = C.cycle((2,), 6)
comparison = C.aggregate_relation((2,), 6)
```

The output includes formal and stabilized coefficients in every codimension, top integrals, exact Pixton certificates, previous-span certificates, independent direct Hodge-integral checks, genus-zero cross-checks, test logs, and source/artifact hashes. The original fourteen-target sample and its result files remain unchanged. The new scope is still \(g=2,3\), unmarked connected maps, and \(D\le2\).
