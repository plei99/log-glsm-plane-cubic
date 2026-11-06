The corrected genus-three rerun produces exact elliptic Gromov–Witten
series, exact O(3)-twisted numbers, and a consistent system of relations
among logarithmic infinity vertices. It does **not** determine the five
one-contact infinity invariants. All calculations below use rational
arithmetic; no free parameter was assigned an arbitrary value.

Write

\[
 V_{g,D;(c_1,\ldots,c_n)}^k(\alpha_1,\ldots,\alpha_n)
 =\int_{[\mathcal M^\infty_{g,D;\mathbf c}]^{\rm red}}
       \psi_{\min}^k\prod_i\operatorname{ev}_i^*\alpha_i.
\]

Here D is ambient degree in P2, the contacts are negative, and the
insertions are aligned with the listed contacts. There are no ordinary
contact psi classes in this notation. H is the ambient hyperplane class.
The contact balance and dimension conditions are

\[
 \sum_i(-c_i-1)=2g-2-3D,\qquad
 k+\sum_i a_i=n+3D\quad(\alpha_i=H^{a_i}).
\]

For genus three this allows D=0 or D=1. With all contacts at most -2,
the profiles are (-5), (-4,-2), (-3,-3), (-3,-2,-2), (-2,-2,-2,-2)
at D=0, and (-2) at D=1. Assigning H powers in {0,1,2}, imposing the
dimension condition, and identifying equal-contact permutations gives
35 primary-contact candidates. This is a finite collection of candidates,
not a proof that arbitrary descendants and all (-1)-contact insertions
reduce to these 35.

The completed degree-at-most-one two-marking extension, retaining the first-stage
positive-degree equations, contains 1,411 equations. After substituting
408 determined values, 921 equations close and the remaining system is

\[
 A x=b,\quad A\in\mathbb Q^{490\times302},\qquad
 \operatorname{rank}A=248,\quad \dim\ker A=54.
\]

There are no remaining nonlinear rows in this finite system. An independent
full augmented-matrix reduction checked consistency, including dependent
rows. The saved certificate gives a particular solution p and a kernel
matrix K and verifies exactly

\[
 Ap=b,\qquad AK=0,\qquad x=p+K\lambda.
\]

Projection to the 35 candidates above has dimension 30: this finite
system imposes five independent affine relations on them. Their complete
coefficients are in
[primary-genus3-constraints.json](../results/genus3/corrected-2026-11-05/d1-extension/primary-genus3-constraints.json).
Projection to the five one-contact candidates has rank five. Thus none of
those five is fixed, and this particular collection of probes imposes no
relation on those five alone. This conclusion concerns the computed
truncation, not all relations available in the theory.

For a compact nontrivial consequence, put

\[
 A_3=V_{3,0;(-5)}^0(H),\quad
 B_3=V_{3,0;(-5)}^1(1),\quad
 C_3=V_{3,1;(-2)}^4(1).
\]

The exact system implies

\[
\begin{aligned}
 &V_{3,0;(-4,-2)}^2(1,1)
 +6V_{3,0;(-4,-2)}^1(H,1)
 +15V_{3,0;(-4,-2)}^0(H^2,1)\\
 &\quad+\frac1{16}C_3-\frac{125}{16}B_3-\frac{375}{8}A_3
 =\frac{113}{2304}.
\end{aligned}
\]

This is a consequence of the finite linear system, not a single graph's
contribution. It was checked independently against p and every column of K.
The degree-one extension also determines the lower-genus inputs
\(V_{2,0;(-3)}^0(H)=1/240\) and
\(V_{2,0;(-3)}^1(1)=-1/80\), which simplify this consequence.

The six-term equation has a more precise localization provenance. Let
R_D denote the residual (left-hand side minus right-hand side) of the
unmarked genus-three ambient-degree-D Chow probe at t^-4. After substituting
the determined values, the boxed equation in its integer-coefficient
normalization is -24 R_0 - R_1, modulo punctured axioms. The script
[genus3_relation_origin.py](../calculations/genus3_relation_origin.py)
checks this equality in every coordinate, including the constant, using
22 additional punctured axiom rows. Its explicit original-row certificate
is [boxed-equation-origin.json](../results/genus3/corrected-2026-11-05/d1-extension/boxed-equation-origin.json).
Thus the equation is a numerical consequence of cancellation of the t^-4
pole in equivariant localization, precisely the source of additional
relations discussed in CJR III, Remark 10.7. This does not by itself
identify a single boundary-divisor identity on a common moduli space.

If every inserted cycle in this six-term equation is stabilized and all
contact markings are forgotten, it becomes a relation of zero-cycles on
Mbar_3, hence has codimension six. The unmarked log-GLSM pushforward has
equivariant codimension two, so its t^-j coefficient has codimension 2+j;
in particular the t^-4 relation is in the top degree.

For this unmarked genus-three target the rational Chow ring is tautological
by [Faber's computation](https://annals.math.princeton.edu/1990/132-2/p05),
and its zero-cycles are determined by degree. This also agrees with
[Graber–Vakil's general tautological socle theorem](https://math.stanford.edu/~vakil/files/socle.pdf).
I directly computed the Pixton 3-spin relation matrix in codimension six:
2,202 rows, 378 strata coordinates, rank 377 over QQ. Thus the quotient by
these relations has dimension one. Its map to the actual top Chow group
is an isomorphism, so any decorated-strata representative of the pushed
six-term identity belongs to the Pixton relation space.

The normalization is explicit: integral(kappa_6)=1/82944. Consequently,
the scalar right-hand side 113/144 becomes 65088 kappa_6 as a zero-cycle
on Mbar_3, rather than a scalar multiple of its fundamental class.
This proves membership abstractly once the terms are represented in the
strata algebra; it does not provide geometric strata formulae for the
individual effective cycles, nor a coefficient vector for this particular
relation in the Pixton matrix. The original virtual-cycle information is
needed to compute those representatives. The conclusion here is specific
to the fully pushed top-degree identity, and is not a blanket assertion
about all log-GLSM pole-cancellation relations or marked effective cycles.
See [genus3_pixton_socle.py](../calculations/genus3_pixton_socle.py) and
[pixton-socle.json](../results/genus3/corrected-2026-11-05/d1-extension/pixton-socle.json).

The punctured string and divisor equations also give the simpler identities

\[
\begin{aligned}
 V_{3,0;(-5,-1)}^1(H,1)&=5A_3,\\
 V_{3,0;(-5,-1)}^1(1,H)&=5A_3,\\
 V_{3,0;(-5,-1)}^2(1,1)&=5B_3-15A_3,\\
 V_{3,0;(-5,-1)}^0(H,H)&=0,\\
 V_{3,0;(-5,-1)}^0(H^2,1)&=0.
\end{aligned}
\]

These specialize CJR's punctured equations (8.3)–(8.4), with
\(\psi_{DF}=-3H\). All five identities were checked against the same
exact affine certificate. The run determines 171 genus-three vertices,
all zero; these include (-1)-contact vertices, rather than any of the 35
candidates listed above. The complete list is in
[verified-values.json](../results/genus3/corrected-2026-11-05/d1-extension/verified-values.json).

For the elliptic curve E, let

\[
 F_{b_1,\ldots,b_n}(q)=\sum_{r\ge0}q^r
 \left\langle\prod_i\tau_{b_i}(\mathrm{pt})\right\rangle^E_{3,r},
 \qquad \sum_i b_i=4.
\]

These are connected invariants, with labeled markings and no extra n!
division. The degree r is intrinsic elliptic degree, so ambient degree is
3r. The point class lifts to H/3, not H2. In the paper normalization

\[
 E_{2m}(q)=-\frac{B_{2m}}{4m}
              +\sum_{r\ge1}\sigma_{2m-1}(r)q^r,
\]

the rerun of the Bloch–Okounkov input gives

\[
\begin{aligned}
 F_4&=\frac16E_2^3+\frac1{12}E_2E_4+\frac1{360}E_6,\\
 F_{3,1}&=-\frac83E_2^4-\frac23E_2^2E_4
              +\frac{25}{63}E_4^2+\frac{13}{180}E_2E_6,\\
 F_{2,2}&=-4E_2^4-\frac56E_2^2E_4
              +\frac{325}{504}E_4^2+\frac1{12}E_2E_6,\\
 F_{2,1,1}&=48E_2^5-\frac{40}{3}E_2^3E_4
              -\frac{25}{9}E_2E_4^2-\frac73E_2^2E_6
              +\frac{35}{36}E_4E_6.
\end{aligned}
\]

For example,

\[
 F_4=-\frac{31}{967680}+\frac{q}{1920}
        +\frac{81q^2}{640}+\frac{781q^3}{480}
        +\frac{16807q^4}{1920}+\frac{9801q^5}{320}
        +\frac{13421q^6}{160}+O(q^7).
\]

Also \([q^2]F_{3,1}=5/12\), \([q]F_{2,2}=1/576\), and
\([q^2]F_{2,1,1}=7/6\). These are known GW inputs to the reconstruction,
recomputed here; they are not newly solved individual infinity vertices.
The JSON artifact includes all four series through r=6. The four-marking
partition (1,1,1,1) was not completed in this rerun.

For O(3), write

\[
 T_{g,d}(H^a\psi^b)=
 \int_{[\overline M_{g,1}(\mathbb P^2,d)]^{\rm vir}}
 e_{\mathbb C^*}(R\pi_*f^*\mathcal O(3))
 \operatorname{ev}^*H^a\psi^b.
\]

The following numbers set the remaining twisting parameter t=1 after
taking the nonequivariant limit for the auxiliary base torus. Replacing
psi by the stabilized ancestor class gives

\[
\begin{aligned}
 T_{3,1}(H\bar\psi^4)&=-\frac{143}{40320},&
 T_{3,2}(H\bar\psi^4)&=-\frac{1823}{13440},\\
 T_{3,1}(H^2\bar\psi^3)&=-\frac1{840},&
 T_{3,2}(H^2\bar\psi^3)&=\frac{4553}{64512}.
\end{aligned}
\]

The ordinary descendants in degree one instead are

\[
 T_{3,1}(H\psi^4)=-\frac{115}{48384},\qquad
 T_{3,1}(H^2\psi^3)=-\frac{59}{60480}.
\]

CJR's dual Euler convention replaces the Euler class by
\(e((R\pi_*f^*\mathcal O(3))^\vee)\), and differs by
\((-1)^{3d+1-g}\). Thus the displayed degree-one numbers change sign
in the dual convention; degree-two numbers do not. The calibrated
Givental backend and the separate t-twisted modular implementation agree
on all four ancestors. For both degree-one descendants, calibrated and
direct fixed-locus localization agree as rational functions before the
base-torus limit, and yield the displayed exact values afterward. The
checks used base weights (0,epsilon,3epsilon).

The CJR equations use the corrected log cotangent factor

\[
 (t-3H)\left(\frac{3H-t}{c}\right)^b,
\]

from CJR III (8.21). Old solved values and old descendant equations were
not reused. The run imported 72 unaffected scalar primary equations from
the older checkpoint and regenerated the axioms. Five freshly compiled
complete localization relations, including the 53-term genus-three
one-point relation with right-hand side -31/967680, are recorded in
[selected-relations.json](../results/genus3/corrected-2026-11-05/selected-relations.json).

The completed calculation used at most two probe markings and powers
t^0 through t^-6. Its first stage included ambient probe degrees 0,1,2;
the completed second stage added degree-zero and degree-one probes. It
compiled in about 106 seconds, stopped at its first 120-second budget,
and finished its remaining solve pass on resumption in 4.1 seconds.
An attempted larger second stage also including degree-two probes was
stopped during compilation after more than 14 minutes, exceeding its
intended 600-second budget. Its incomplete equations
are not included in the rank claims above. Solver budgets are checked
between compilation/solve blocks, rather than interrupting an individual
symbolic compilation. Increasing the positive-degree probe coverage,
adding more markings, or supplying further independent lower-genus input
may strengthen the constraints; this run does not establish which of
those will determine the remaining values.

Reproduction uses SageMath 10.9 / Python 3.14.10 in the OrbStack agents
checkout. The scripts are
[genus3_corrected.py](../calculations/genus3_corrected.py),
[genus3_degree1_extension.py](../calculations/genus3_degree1_extension.py),
[genus3_analyze.py](../calculations/genus3_analyze.py),
[genus3_verify_conclusions.py](../calculations/genus3_verify_conclusions.py),
[genus3_known.py](../calculations/genus3_known.py), and
[genus3_o3.py](../calculations/genus3_o3.py).
The full exact linear certificate is
[linear-certificate.json](../results/genus3/corrected-2026-11-05/d1-extension/linear-certificate.json).
Run these from the log-glsm-plane-cubic repository with the Sage
environment's Python. The extension resumes the saved first-stage
checkpoint; the analysis takes its output directory via GENUS3_OUTPUT.

The draft references are Chen–Janda–Ruan, *Structural Formulae in
Logarithmic Gauged Linear Sigma Models I*, equation (1.10), and *III*,
equation (8.21) and Section 9.3, both October 17, 2025, plus their
*Punctured Logarithmic R-Maps*, Theorem 8.1, equations (8.3)–(8.4),
April 25, 2023. These are the files in the user's iCloud texts/math
directory identified in the preceding audit. The stationary elliptic
input uses the n-point formalism of
[Bloch–Okounkov](https://arxiv.org/abs/alg-geom/9712009) and its GW
interpretation in
[Okounkov–Pandharipande, Gromov–Witten theory, Hurwitz theory, and completed cycles](https://arxiv.org/abs/math/0204305).
