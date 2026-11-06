# Genus-three relations in codimensions three, four, and five

Computed 5 November 2026, with rational coefficients throughout.

The calculation produces nine auxiliary-localization relations on
\(\overline{\mathcal M}_3\): two in codimension three, three in codimension
four, and four in codimension five. Every one has an exact certificate as
a linear combination of Pixton relations. It also computes the six
dual-Euler twisted classes \(Z_{D,r}\), for ambient curve degrees
\(D=0,1\) and codimensions \(r=3,4,5\). The Chen–Janda–Ruan localization
formula identifies their negatives with the corresponding sums of all
infinity-containing graph contributions.

Thus these are explicit additional equations obtained from the GLSM/O(3)
calculation. They are not relations outside Pixton's system. The calculation
does not reconstruct the individual positive-dimensional punctured effective
cycles occurring in the CJR sum.

The complete expanded formulas are in
[explicit-formulas.md](../results/genus3/lower-cycles-2026-11-05/explicit-formulas.md).
All nine pole equations, normalized to primitive integer coefficients, are
printed in [all-pole-relations.md](../results/genus3/lower-cycles-2026-11-05/all-pole-relations.md).
The source-by-source audit, including the unpublished draft's equation and
page numbers, is in
[GENUS_THREE_LOWER_CYCLES_SOURCES.md](GENUS_THREE_LOWER_CYCLES_SOURCES.md).

Write \(\lambda_i=c_i(\mathbb E_3)\) and
\(\kappa_i=p_*(\psi_1^{i+1})\) on \(\overline{\mathcal M}_3\), where
\(p:\overline{\mathcal M}_{3,1}\to\overline{\mathcal M}_3\) forgets the
marking. In particular, \(\kappa_0=4\). Let
\[
\xi:\overline{\mathcal M}_{1,1}\times\overline{\mathcal M}_{2,1}
\longrightarrow\overline{\mathcal M}_3
\]
glue the two marked points. On the first factor put \(a=\lambda_1=\psi_1\);
on the second put \(b_i=\lambda_i\) and \(\psi=\psi_1\).
The symbol \(\xi_*\) denotes the actual stack pushforward, without an
additional factor of \(1/2\). The two genera distinguish the components.

One short relation in \(A^3(\overline{\mathcal M}_3)\) is
\[
\boxed{
2\lambda_2\kappa_1-\lambda_1\kappa_2-4\lambda_1\lambda_2
+\xi_*\!\left[
 1\boxtimes(\psi b_1-2b_2)
 +a\boxtimes(\psi-b_1)
\right]=0.}
\tag{R3b}
\]
A second useful representative is
\[
\begin{aligned}
0={}&3\kappa_3-5\lambda_1\kappa_2+7\lambda_2\kappa_1
 -8\lambda_1\lambda_2-12\lambda_3\\
&+\xi_*\!\left[
 1\boxtimes(-5b_2+3\psi b_1-\psi^2)
 +a\boxtimes(2\psi-2b_1)
\right].
\end{aligned}\tag{R3a}
\]
These equations exhibit the boundary correction to a relation among global
Hodge and kappa classes. The correction is supported on the separating
divisor of type \(1+2\). This statement concerns the displayed pure
tautological equations; a CJR infinity-containing graph need not itself
stabilize into the boundary.

In codimension four the computation gives, among its three pole relations,
\[
\begin{aligned}
0={}&-11\kappa_4+17\lambda_1\kappa_3-24\lambda_2\kappa_2
 +7\lambda_1\lambda_2\kappa_1+11\lambda_3\kappa_1
 -32\lambda_1\lambda_3\\
&+\xi_*\!\left[
1\boxtimes(2\psi^3-5\psi^2b_1+9\psi b_2-4b_1b_2)
+a\boxtimes(\psi b_1-2b_2)
\right].
\end{aligned}\tag{R4}
\]
In codimension five one representative is
\[
\begin{aligned}
0={}&36\kappa_5-116\lambda_1\kappa_4+280\lambda_2\kappa_3
 -173\lambda_1\lambda_2\kappa_2-31\lambda_3\kappa_2\\
&+279\lambda_1\lambda_3\kappa_1-460\lambda_2\lambda_3\\
&+\xi_*\!\left[
1\boxtimes(6\psi^4-7b_1\psi^3-28b_2\psi^2+44b_1b_2\psi)
\right.\\[-2pt]
&\hspace{42mm}\left.
+a\boxtimes(8\psi^3-32b_1\psi^2+78b_2\psi-46b_1b_2)
\right].
\end{aligned}\tag{R5}
\]
Every term has the advertised codimension: \(\xi_*\) adds one. The four
displayed expressions were checked both against the original localization
poles and against their own exact Pixton certificates, so the simplification
was not accepted solely on the basis of a transcription.

Here is the relation to the CJR calculation. Let
\(s:\overline{\mathcal M}_{3,0}(\mathbb P^2,D)\to\overline{\mathcal M}_3\)
be stabilization, and let \(\pi,f\) be the universal curve and map. Use
the dual Euler convention from CJR III, §9.3:
\[
Z_D(t)=s_*\left[
e_t\!\left((R\pi_*f^*\mathcal O(3))^\vee\right)
\cap[\overline{\mathcal M}_{3,0}(\mathbb P^2,D)]^{\rm vir}
\right],\qquad Z_{D,r}=[t^{2-r}]Z_D(t).
\]
The Euler class is the ratio for the derived complex, not just its degree-zero
cohomology. A dual fibre has equivariant first Chern class \(t-3H\).
The virtual dimension is \(2+3D\), whereas the Euler ratio has degree
\(3D-2\). The pushforward therefore has total equivariant codimension two.
Consequently \(t^{-1},t^{-2},t^{-3}\) have coefficients in codimensions
three, four, five. The earlier scalar computation used \(t^{-4}\), whose
coefficient is top degree. This explains precisely which other “degrees”
are computed here: lower Chow codimensions, at ambient curve degrees zero
and one.

After stabilization, write CJR localization as
\[
L_D(t)=Z_D(t)+I_D(t),
\]
where \(I_D\) sums every graph with at least one infinity vertex, including
its automorphism factors. The left side is unlocalized and has no negative
powers of \(t\). Therefore
\[
\boxed{B_{D,r}:=[t^{2-r}]I_D(t)=-Z_{D,r},\qquad r=3,4,5.}
\]
This is Theorem 8.4, equation (8.28), and Remark 8.5 of the
[private CJR III draft](</Users/patrick/Library/Mobile Documents/com~apple~CloudDocs/texts/math/Chen-Janda-Ruan_Structural_Formulae_in_Logarithmic_Gauged_Linear_Sigma_Models_III_the_Localization_Formula.pdf>),
specialized using §9.3; Remark 10.7 describes its use for relations. The sign
in the draft's constant-term equation (9.7) is not an extra sign to insert
into this negative-coefficient identity.

At degree zero the complete answer has an independent short derivation.
Put \(h=H\) on \(\mathbb P^2\). The constant-map obstruction class is
\[
c_6(\mathbb E^\vee\boxtimes T_{\mathbb P^2})
=\lambda_3^2-3h\lambda_2\lambda_3
 +3h^2(\lambda_2^2+\lambda_1\lambda_3).
\]
Mumford's relation \(c(\mathbb E)c(\mathbb E^\vee)=1\) gives
\(\lambda_3^2=0\), \(\lambda_2^2=2\lambda_1\lambda_3\), and
\(\lambda_1^2=2\lambda_2\). Set
\[
F(u)=\frac{u}{u^3+\lambda_1u^2+\lambda_2u+\lambda_3}.
\]
The twist is \(F(t-3h)\), so integration over \(\mathbb P^2\) gives
\[
Z_0(t)=9\lambda_1\lambda_3F(t)+9\lambda_2\lambda_3F'(t)
=9\lambda_1\lambda_3t^{-2}-36\lambda_2\lambda_3t^{-3}
 +36\lambda_1\lambda_2\lambda_3t^{-4}.
\]
In particular the CJR aggregates are exactly
\[
\boxed{B_{0,3}=0,\qquad B_{0,4}=-9\lambda_1\lambda_3,
\qquad B_{0,5}=36\lambda_2\lambda_3.}
\]
The constant-map obstruction calculation and the Hodge identities are
justified by [Faber–Pandharipande](https://arxiv.org/abs/math/9810173),
introduction and §4.3, and [Mumford](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf),
§5. The latter yields the identities in rational Chow.

At degree one, for example, the finite class in codimension three is
\[
\begin{aligned}
Z_{1,3}={}&66\kappa_3-135\lambda_1\kappa_2+213\lambda_2\kappa_1
 -252\lambda_1\lambda_2\\
&+\xi_*\!\left[
1\boxtimes(-93b_2+27\psi b_1-6\psi^2)
+a\boxtimes(-12b_1+6\psi)
\right].
\end{aligned}
\]
Thus the degree-one CJR identity is \(B_{1,3}=-Z_{1,3}\). The equally
explicit \(Z_{1,4}\) and \(Z_{1,5}\) are printed in
[the complete formulas](../results/genus3/lower-cycles-2026-11-05/explicit-formulas.md).
These are representatives of classes: their coefficients can change by
adding the displayed Pixton relations, without changing the answers.

The new implementation retains stable-vertex Hodge and cotangent
polynomials as classes, instead of integrating each vertex immediately.
It applies [Graber–Pandharipande virtual localization](https://arxiv.org/abs/alg-geom/9708001),
§4, with the dual \(\mathcal O(3)\) Euler factors. There are three fixed
loci at degree zero. At degree one there are six genus-three tails and six
genus-one/genus-two bridges. For tails, stabilization uses
\[
p_*(1)=0,\qquad p_*(\psi^m p^*\alpha)=\kappa_{m-1}\alpha
\quad(m\ge1).
\]
For bridges it contracts the rational component and applies \(\xi_*\).
There is no degree-one deck factor and no exchange symmetry between
genus-one and genus-two vertices. The source audit gives the full vertex
and edge rational functions from which the coefficients were obtained.

The auxiliary \(\mathbb P^2\) torus has weights
\(\epsilon(0,1,3)\). This auxiliary parameter is different from CJR's
\(t\). In codimension \(r\), homogeneity reads
\[
[Z_D(t;\epsilon)]_r
=\sum_k\epsilon^k t^{2-r-k}C_{D,r,k}.
\]
The global twisted class exists without inverting auxiliary characters,
so its negative \(\epsilon\)-coefficients vanish in Chow. Degree-one
localization produces powers \((-2,-1)\) in codimension three,
\((-3,-2,-1)\) in codimension four, and \((-4,-3,-2,-1)\) in
codimension five. These are the nine relations. Their raw strata vectors
have ranks two, three, and four within their respective degrees. This is
linear independence of those expressions before quotienting by relations,
not independence from Pixton relations.

For a direct link between the output and the displayed formulas, write
\(C_{r,k}=C_{1,r,k}\), simplify \(\kappa_0=4\) and
\(\psi_1=\lambda_1\) on \(\overline{\mathcal M}_{1,1}\), and take
\[
R3a=\frac{36}{49}C_{3,-2},\qquad
R3b=\frac{216}{245}C_{3,-2}+\frac9{70}C_{3,-1},\qquad
R4=\frac9{35}C_{4,-1},\qquad R5=\frac97 C_{5,-1}.
\]
The finite parts \(C_{D,r,0}\) give the \(Z_{D,r}\). Repeating with
weights \(\epsilon(0,2,5)\) gives identical finite classes modulo Pixton
in all six cases. This is a useful second computation; geometric regularity,
rather than agreement on just two rays, justifies the general limit.

Pixton membership was checked over \(\mathbb Q\), not with floating point
arithmetic or intersection pairings. The three relation matrices have sizes
\(39\times40\), \(245\times121\), and \(970\times262\), and ranks
\(30,114,259\). These leave quotient dimensions \(10,7,3\) for the
generated relation spaces. For every pole vector \(v\) the program found
and saved a rational row vector \(u\) satisfying
\[
\boxed{uM=v.}
\]
It separately saved the analogous certificates for the six finite-part
differences and the four displayed formulas. The matrix columns were
converted to the same decorated-stratum coordinates used for the cycles;
the exact column ordering is saved alongside the matrices.

The 3-spin relations come from
[Pandharipande–Pixton–Zvonkine](https://arxiv.org/abs/1303.1043).
For the assertion in **Chow**, the needed result is
[Janda, Theorem 1.1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf).
The certificate \(uM=v\), together with that theorem, proves the claimed
membership and vanishing. No conjecture that Pixton's relations generate
all tautological relations is needed.

As top-degree normalization checks, the new class calculation reproduces
\[
\int_{\overline{\mathcal M}_3} Z_{0,6}=\frac1{40320},\qquad
\int_{\overline{\mathcal M}_3} Z_{1,6}=-\frac{1637}{80640}.
\]
The first also follows independently from
\(\int\lambda_1\lambda_2\lambda_3=1/1451520\), by the Bernoulli
formula in Faber–Pandharipande. The second agrees with the previous scalar
backend; it is a computed value, not a number claimed to be printed in
the cited literature. Because the two backends share fixed-locus factors,
this comparison checks stabilization and integration but is not an
independent derivation of every localization factor.

The implementation is
[o3_unmarked_cycles.sage](../o3_unmarked_cycles.sage). The calculation,
presentation checks, and saved-certificate replay are, respectively,
[genus3_lower_cycles.py](../calculations/genus3_lower_cycles.py),
[genus3_present_relations.py](../calculations/genus3_present_relations.py),
and [verify_genus3_lower_certificates.py](../calculations/verify_genus3_lower_certificates.py).
The original scalar backend's positive-dimensional guards remain in place.
The new cycle backend deliberately supports only unmarked genus three in
ambient degrees zero and one.

Run from the project's native Linux directory, with its Sage environment:

```bash
export PATH=/home/patrick/agent-work/tools/mamba/envs/sage/bin:$PATH
python calculations/genus3_lower_cycles.py
python calculations/genus3_present_relations.py
python calculations/verify_genus3_lower_certificates.py
```

From macOS, prefix the command with
`orb -m agents -w /home/patrick/agent-work/repos/research/test/log-glsm-plane-cubic /bin/bash -lc`.
The full calculation checks both top integrals, all nine auxiliary poles,
all six changes of torus weights, and the independent degree-zero formulas.
The presentation script verifies the four displayed formulas and their
origins as combinations of the poles. The replay script verifies nineteen
saved exact certificates and reconstructs all pole/displayed vectors from
their formal classes. Logs and machine-readable results are in
[lower-cycles-2026-11-05](../results/genus3/lower-cycles-2026-11-05/).

Finally, the identity for \(B_{D,r}\) is an identity for the entire
infinity-containing sum. It proves that sum is tautological by identifying
it with \(-Z_{D,r}\). It does not prove tautologicality of each constituent
or give its value separately. The distinction is required by CJR III,
Proposition 10.9, which assumes tautologicality of effective cycles, and by
the corresponding question in
[Punctured logarithmic R-maps, §9.8](https://arxiv.org/abs/2208.04519v2).
Likewise the punctured descendant \(\psi_{\min}\) is not an ordinary
marked-point cotangent class. The present calculation avoids identifying
those two classes.
