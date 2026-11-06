# Quartic K3: ordinary-GW vanishing and exact CJR relations

Run started 5 November 2026 and completed 6 November 2026 (Asia/Shanghai).

The completed search produced **344 independent tautological relation
vectors**, all certified inside the Pixton relation space, and **20 explicit
CJR aggregate identities**. The relation-vector ranks are
\[
\begin{aligned}
g=2:&\quad (r=1,2,3):\ (1,6,16),\\
g=3:&\quad (r=2,3,4,5,6):\ (6,28,72,98,117).
\end{aligned}
\]
In genus two these equal the full generated Pixton ranks in every
codimension. In genus three they fill that space in codimension two; the
generated Pixton ranks in codimensions three through six are
\(30,114,259,377\), so those degrees are only partially spanned by this
quartic search. The last ten degree-three rays added no new directions,
including all eight rays outside the initial quadratic family. This is an
observed stabilization of the finite search, not an exhaustion theorem.

This calculation replaces the plane cubic by a smooth quartic surface
\(S\subset\mathbb P^3\), and the ambient dual Euler twist by
\((R\pi_*f^*\mathcal O(4))^\vee\). It computes unmarked genera two and
three through ambient curve degree three. Ambient degree means
\(D=\int_\beta H\); it is not a choice of one primitive class in
\(H_2(S)\).
The hyperplane curve class on the quartic has ambient degree four and lies
beyond this finite scan.

The implementation is separate from the elliptic compact-target backend.
It supplies exact class-valued ambient localization, CJR aggregate
identities, and pure tautological relations from auxiliary torus
localization. It does not claim to reconstruct each punctured effective
cycle individually. The two main modules are
[quartic_k3.sage](../quartic_k3.sage) and
[quartic_k3_fixed_graphs.sage](../quartic_k3_fixed_graphs.sage).
The complete source audit is [QUARTIC_K3_SOURCES.md](QUARTIC_K3_SOURCES.md).

All positive-degree **ordinary** K3 GW virtual classes vanish. The
holomorphic symplectic form gives a trivial obstruction quotient for a
nonconstant map. Removing that quotient defines the distinct reduced K3
theory, whose invariants need not vanish. This is the ordinary/reduced
distinction in [Maulik–Pandharipande–Thomas, §0.1](https://people.math.ethz.ch/~rahul/mpt.pdf).

The reduced *CJR log-GLSM* class is not the reduced *K3 GW* class. The CJR
compact-target comparison recovers ordinary hypersurface theory: see
[Chen–Janda–Ruan, Proposition 4.1](https://arxiv.org/abs/1906.04345v2),
and equation (9.7) of the supplied private CJR III draft. Substituting
reduced K3 invariants into that comparison would change both the virtual
dimension and the mathematical theory. The new module explicitly refuses
to label reduced K3 invariants as zero.

There is also vanishing at ambient degree zero for every unmarked
\(g\ge2\). This requires a different argument because the differential
of a constant map is zero. The obstruction bundle is
\(\mathbb E^\vee\boxtimes TS\) on \(\overline{\mathcal M}_g\times S\).
Using the symplectic roots \(y,-y\) of \(TS\), its top Chern class is
\[
\begin{aligned}
c_{2g}(\mathbb E^\vee\boxtimes TS)
&=\prod_{i=1}^g(x_i^2+c_2(TS))\\
&=\lambda_g^2+c_2(TS)
  (\lambda_{g-1}^2-2\lambda_g\lambda_{g-2})=0.
\end{aligned}
\]
Here \(c_2(TS)^2=0\) by dimension, and the two Hodge coefficients vanish
by Mumford's relation \(c(\mathbb E)c(\mathbb E^\vee)=1\). This does not
extend to genus-one constant maps. The ingredients are
[Faber–Pandharipande, introduction](https://arxiv.org/abs/math/9810173)
and [Mumford, §5](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf).

Let \(s_D:\overline{\mathcal M}_{g,0}(\mathbb P^3,D)
\to\overline{\mathcal M}_g\) stabilize the domain. Define
\[
Z_D^{(g)}(t)=s_{D,*}\left(
e_t((R\pi_*f^*\mathcal O(4))^\vee)
\cap[\overline{\mathcal M}_{g,0}(\mathbb P^3,D)]^{\rm vir}\right),
\qquad Z_{D,r}^{(g)}=[t^{2g-2-r}]Z_D^{(g)}(t).
\]
The twist is an Euler ratio for a derived complex, with fibre weight
\(t-4H\). The ambient virtual dimension is \(4D\); the Euler ratio has
degree \(4D+1-g\). Thus the twisted class has total equivariant dimension
\(g-1\), or codimension \(2g-2\) on \(\overline{\mathcal M}_g\).

After stabilization the CJR localization identity has the form
\[
L_D^{(g)}(t)=Z_D^{(g)}(t)+I_D^{(g)}(t),
\]
where \(I_D\) sums **every** infinity-containing graph with its
automorphism and gluing factors. Its unlocalized left side has no negative
powers of \(t\). Moreover its constant term is
\[
[t^0]L_D^{(g)}(t)=(-1)^{1-g+4D}
\sum_{\beta:\,i_*\beta=D[\mathrm{line}]}
\operatorname{st}_{S,*}[\overline{\mathcal M}_g(S,\beta)]^{\rm vir}=0.
\]
Therefore, defining \(B_{D,r}^{(g)}=[t^{2g-2-r}]I_D^{(g)}(t)\),
\[
\boxed{B_{D,r}^{(g)}=-Z_{D,r}^{(g)},\qquad
2g-2\le r\le3g-3.}
\]
These statements use CJR III, Theorem 8.4, (8.28), Remark 8.5, §9.3,
and (9.7), with precise draft page references in the source audit. At
\(r=2g-2\) the equation uses K3 vanishing; at larger \(r\) it uses
negative-power cancellation. Positive powers of \(t\) are not forced to
vanish by this argument.

Consequently there are two aggregate identities per ambient degree in
genus two (codimensions two and three), and three in genus three
(codimensions four, five, six). The run through \(D=0,1,2,3\) therefore
gives twenty aggregate identities. All twenty right sides are explicit in
[cjr-aggregate-classes.md](../results/quartic-k3/2026-11-05/cjr-aggregate-classes.md).

The ambient degree-zero answers admit independent closed formulas:
\[
\boxed{
Z_0^{(2)}(t)=-20\lambda_1\lambda_2\,t^{-1},\qquad
Z_0^{(3)}(t)=20\lambda_1\lambda_2\lambda_3\,t^{-2}.}
\]
Indeed, for \(H^4=0\),
\[
Z_0^{(g)}(t)=[H^3]
\left\{\prod_{i=1}^g
(-x_i^3+4Hx_i^2-6H^2x_i+4H^3)
\frac{t-4H}{\Lambda_g^+(t-4H)}\right\},
\]
where \(\Lambda_g^+(u)=\sum_i\lambda_i u^{g-i}\).
The product is the constant-map obstruction Euler class for \(\mathbb P^3\).
Mumford reduction makes it \(-20H^3\lambda_1\lambda_2\) in genus two
and \(20H^3\lambda_1\lambda_2\lambda_3\) in genus three. Their Hodge
degrees are already top degree, so only the leading part of the Euler
ratio survives. In particular,
\[
B_{0,2}^{(2)}=0,\quad B_{0,3}^{(2)}=20\lambda_1\lambda_2,
\qquad
B_{0,4}^{(3)}=B_{0,5}^{(3)}=0,\quad
B_{0,6}^{(3)}=-20\lambda_1\lambda_2\lambda_3.
\]

The top-degree integrals of the ambient twisted coefficients, ordered by
ambient degree \(D=0,1,2,3\), are exactly
\[
\begin{aligned}
\left(\int Z_{D,3}^{(2)}\right)_{D=0}^3
&=\left(-\frac1{288},-\frac{17}{72},-\frac{67}{36},\frac{133}{2}\right),\\
\left(\int Z_{D,6}^{(3)}\right)_{D=0}^3
&=\left(\frac1{72576},-\frac{491}{7560},\frac{5183}{7560},-\frac{277}{1260}\right).
\end{aligned}
\]
The corresponding infinity-aggregate integrals are their negatives. These
are computed twisted/aggregate values, not nonzero ordinary K3 invariants.
The two degree-zero numbers also follow independently from the
Faber–Pandharipande Bernoulli formula for
\(\int\lambda_g\lambda_{g-1}\lambda_{g-2}\).

The genus-two constant-term classes have a much shorter presentation after
reduction by the certified relations:
\[
\begin{aligned}
Z_{0,2}^{(2)}&=0,\\
Z_{1,2}^{(2)}&=24\kappa_2+352\lambda_2,\\
Z_{2,2}^{(2)}&=2544\kappa_2-6848\lambda_2,\\
Z_{3,2}^{(2)}&=235296\kappa_2-880512\lambda_2.
\end{aligned}
\]
Thus each corresponding CJR infinity aggregate is the negative of the
displayed class. The two classes \(\kappa_2,\lambda_2\), together with
the Pixton rows, span all eight strata coordinates in codimension two.
The reduction coefficients and exact matrix certificates are saved in
[genus-two-constant-term.json](../results/quartic-k3/2026-11-05/genus-two-constant-term.json).

There is a second family of relations, involving only tautological classes.
These arise from an auxiliary torus on \(\mathbb P^3\); its scale
\(\epsilon\) is different from the CJR parameter \(t\). Write
\(H|_{p_i}=a_i=\epsilon A_i\), and
\(\Lambda_h^\vee(u)=\sum_\ell(-1)^\ell\lambda_\ell u^{h-\ell}\).
For a contracted stable vertex of genus \(h\), colour \(i\), and valence
\(v\), the contribution is
\[
\frac{(t-4a_i)^{1-v}}{\Lambda_h^+(t-4a_i)}
\left(\prod_{k\ne i}(a_i-a_k)\right)^{v-1}
\prod_{k\ne i}\Lambda_h^\vee(a_i-a_k)
\prod_{e\ni i}\frac1{w_{e,i}-\psi_{e,i}}.
\]
There are three tangent Hodge factors and four fixed-point colours. On a
degree-\(d\) edge joining colours \(i,j\), put \(\Delta=a_i-a_j\).
Its factor is
\[
\frac{(-1)^d d^{2d}}{(d!)^2\Delta^{2d}}
\frac{\displaystyle\prod_{m=0}^{4d}
\left(t-\frac{4d-m}{d}a_i-\frac m d a_j\right)}
{\displaystyle\prod_{k\ne i,j}\prod_{m=0}^{d}
\left(\frac m d a_i+\frac{d-m}{d}a_j-a_k\right)}.
\]
The graph weight also includes \(1/|\operatorname{Aut}\Gamma|\) and
\(\prod_e1/d_e\). An univalent unmarked rational endpoint contributes
its flag weight. A bivalent rational vertex contributes
\[
\frac{\prod_{k\ne i}(a_i-a_k)}
{(w_1+w_2)(t-4a_i)}.
\]
These are the projective-space virtual-localization factors of
[Graber–Pandharipande, §4](https://arxiv.org/abs/alg-geom/9708001),
with the explicit dual \(\mathcal O(4)\) Euler ratio.

The graph backend retains all stable-vertex classes. A rational tail is
replaced by a temporary marking; a bivalent rational chain is contracted;
the remaining components are glued; temporary markings are forgotten with
the actual tautological pushforward. For one tail this gives
\(p_*1=0\) and \(p_*\psi^m=\kappa_{m-1}\). For several tails the
implementation uses the full forgetful map, including its correction
terms. A trivalent genus-zero contracted component can become unstable
after the tails are forgotten, and that case is stabilized globally.
Parallel edges and self-loops are retained.

The degree-one calculation is implemented twice: directly as tails and
bridges, and by the general graph enumerator. Their stabilized classes
agree in every codimension in genera two and three. In genus two the
ordered-colour bridge sum has a half to remove duplicate indexing of equal
genus endpoints. This is distinct from a stack-stabilizer correction.

Mumford's relation also gives an exact computational simplification:
\[
\frac1{\Lambda_h^+(q)}=\frac{\Lambda_h^\vee(q)}{q^{2h}}.
\]
This avoids unnecessarily long formal inverse-Hodge expansions. Hodge
classes are pulled through forgetful maps where the projection formula
allows it, and permutations of temporary markings on one component are
collected together. Neither simplification discards a potentially nonzero
class.

The nonequivariant ambient twisted class exists without inverse powers of
the auxiliary parameter. At Chow codimension \(r\), homogeneity gives
\[
[Z_D^{(g)}(t;\epsilon A)]_r
=\sum_k\epsilon^k t^{2g-2-r-k}C_{g,D,r,k}(A).
\]
Thus every negative coefficient \(C_{g,D,r,k}\), \(k<0\), vanishes in
Chow. The finite part is \(Z_{D,r}^{(g)}\), so finite parts obtained from
two weight choices must agree. Both types of zero expressions are included
in the relation search. The geometric regularity supplies the justification
for arbitrary nonresonant weight choices; checking several rays alone is
not the proof of regularity.

For example, the genus-two codimension-one relation simplifies to
\[
\boxed{\kappa_1-2\lambda_1-
\frac12\xi_{1,1,*}(1\boxtimes1)=0.}
\]
Since \(\delta_1=\frac12\xi_{1,1,*}(1\boxtimes1)\), this says
\(\kappa_1=2\lambda_1+\delta_1\). Here the half converts the ordered
gluing pushforward into the boundary-divisor class; it is not an additional
factor on the already weighted localization graph sum.

A genus-three codimension-two relation is
\[
\boxed{\kappa_2-\lambda_1\kappa_1+4\lambda_2
=\xi_{1,2,*}\bigl(1\boxtimes(\psi-b_1)\bigr).}
\]
Here \(b_1=\lambda_1\) and \(\psi=\psi_1\) on
\(\overline{\mathcal M}_{2,1}\). The genus-one terms cancel using
\(\psi_1=\lambda_1\) on \(\overline{\mathcal M}_{1,1}\).
This relation comes from comparing finite auxiliary limits, not from
claiming that a positive CJR power vanishes. The complete independent
relation list, including all degree-two and degree-three graph terms, is
in [independent-relations.md](../results/quartic-k3/2026-11-05/independent-relations.md).
Its strata dictionary specifies the genera, flags, node pairings, and
decorations. The aggregate-class appendix additionally lists temporary
markings and their forgetful maps. The readable relation basis is row-reduced;
the original localization expressions and their origins are preserved in
losslessly compressed `final-relations-g*-r*.json.gz` files.

For each codimension the program expands the relations in the same
decorated-stratum coordinates as the Pixton matrix \(M\), retains an
independent set of vectors \(v\), and solves over \(\mathbb Q\) for
\[
uM=v.
\]
The saved rational coefficients are explicit Pixton membership certificates.
They are replayed in batches and the ranks of the retained vectors are
checked again. Independence here means independence as vectors of strata
expressions before imposing relations. Every such expression is zero after
passing to Chow; no independence from Pixton is asserted.

The 3-spin relations are those of
[Pandharipande–Pixton–Zvonkine](https://arxiv.org/abs/1303.1043).
Their Chow validity, needed here, is
[Janda, Theorem 1.1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf).
The certificates prove membership without any assumption that Pixton's
relations generate all tautological relations. Equality with the generated
Pixton space in a particular degree is established only when the computed
subspace has the same rank as that matrix.

The search uses 64 auxiliary rays at ambient degree one, eight nonresonant
rays at degree two, and sixteen at degree three. The degree-three list
includes eight rays outside the initial quadratic family. It covers every Chow
codimension on \(\overline{\mathcal M}_2\) and
\(\overline{\mathcal M}_3\), rather than only the CJR-negative range.
Redundant relations are discarded by exact elimination. This is a finite
search through ambient degree three, not an assertion that all possible
degrees, markings, insertions, or genera have been exhausted.

The reproducible commands, run from the project directory with Sage's
Python, are:

```bash
python calculations/quartic_k3_relations.py --rays 64
python calculations/quartic_k3_degree_two.py --degree 2
python calculations/quartic_k3_degree_two.py --degree 3
python calculations/quartic_k3_export.py
```

The pinned environment is Sage 10.9 on the native Linux OrbStack `agents`
machine. From macOS use
`orb -m agents -w /home/patrick/agent-work/repos/research/test/log-glsm-plane-cubic /bin/bash -lc`
and set `PATH=/home/patrick/agent-work/tools/mamba/envs/sage/bin:$PATH`.
The [test suite](../test_quartic_k3.sage) checks the independently derived
degree-zero classes and integrals, the ordinary/reduced distinction, the
quartic contact balance, the numerical normalization of the direct/general
graph formulas, and a second weight choice. The existing cubic fixed-locus
suite verifies that the shared enumerator extension preserves its behavior.

Finally, quartic effective vertices have the balance
\[
\sum_i(c_i+1)=4D-(2g-2),\qquad c_i\le-1,
\]
and reduced punctured virtual dimension \(1-g+4D+n\).
The all-\(-2\) profile has dimension \(g-1\), so it is not automatically
a scalar basic invariant. The hypotheses for basic-scalar reconstruction
in CJR's punctured theory fail in this K3 setting for \(g\ge2\).
Ordinary compact K3 vanishing therefore supplies equations among the
effective contributions; it does not set each effective contribution to
zero. The relevant sources are [CJR, Punctured logarithmic R-maps,
§9](https://arxiv.org/abs/2208.04519v2), especially (9.5), (9.8),
Definition 9.19, Assumption 9.22, and the tautologicality question in §9.8.
