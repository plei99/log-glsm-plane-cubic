# Explicit four-mark equations

The equations below use the literal decorated-strata generators of the pinned admcycles basis. No extra automorphism factors are inserted. Each vector lies in the computed span and in the generated Pixton span, and lies outside the comparison span of unit pullbacks plus divisor products of the computed codimension-one relations. This last statement does not assert novelty relative to other tautological operations.

In the graph notation, the first list gives vertex genera, the second gives the half-edges at each vertex, and the final list pairs the half-edges forming nodes. Labels 1 through 4 are markings; all other labels are node half-edges. Polynomial subscripts index vertices starting at zero.

## Genus 2

The equation is $\sum_j c_j S_j=0$, with coefficients

```text
{233: 24, 267: -1}
```

S_233, coefficient 24:

```text
Graph :      [1, 1] [[6], [1, 2, 3, 4, 7]] [(6, 7)]
Polynomial : (kappa_1)_0
```

S_267, coefficient -1:

```text
Graph :      [0, 1] [[7, 8, 9], [1, 2, 3, 4, 10]] [(7, 8), (9, 10)]
Polynomial : 1
```

## Genus 3

The equation is $\sum_j c_j S_j=0$, with coefficients

```text
{28: 1, 88: -1}
```

S_28, coefficient 1:

```text
Graph :      [0, 3] [[2, 3, 4, 6], [1, 7]] [(6, 7)]
Polynomial : (kappa_1)_0
```

S_88, coefficient -1:

```text
Graph :      [0, 0, 3] [[2, 3, 7], [4, 8, 9], [1, 10]] [(7, 8), (9, 10)]
Polynomial : 1
```

## Geometric interpretation

The genus-two equation is
\[
24\,\xi_*(\kappa_1\boxtimes1)=\eta_*1.
\]
Here xi glues Mbar_(1,1) to Mbar_(1,5); all four external markings lie on the second factor. For eta, replace the first elliptic component by a rational component with a self-node. On Mbar_(1,1), kappa_1=lambda_1 and 12 lambda_1=delta_irr. The gluing map Mbar_(0,3)→Mbar_(1,1) has pushforward 2 delta_irr, because the two branches of the node are ordered. Thus its class is 24 kappa_1. The factor 24, rather than 12, records this stack and gluing convention. Equivalently, its integral is 1, while the integral of kappa_1 is 1/24.

The genus-three equation is
\[
\xi_*(\kappa_1\boxtimes1)=\eta_*1.
\]
Now xi glues Mbar_(0,4) to Mbar_(3,2): marks 2,3,4 lie on the rational component, and mark 1 lies on the genus-three component. The map eta splits that rational component at the boundary {2,3}|{4,*}. Since Mbar_(0,4) is P¹ and kappa_1 has degree one, it is the class of this boundary point. All graphs in this genus-three equation have trivial automorphism group.

These are familiar lower-genus identities pushed forward by gluing. Their position outside the particular comparison span explains why rank growth at higher markings need not produce a new tautological relation. The exact Pixton and found-span certificates establish membership without relying on a normalization match to a drawing.

For the Hodge-integral convention, see [Faber–Pandharipande, Hodge integrals and Gromov–Witten theory](https://arxiv.org/abs/math/9810173). The general strata and Pixton conventions are those of the pinned admcycles package and the matrices accompanying this file.
