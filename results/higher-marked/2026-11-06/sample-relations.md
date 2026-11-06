# Explicit three-mark codimension-two relations

Each equation uses the literal decorated-strata generators listed below, with no extra automorphism factors. These representatives are obtained by exact linear algebra from the computed relation span and unit pullbacks. Their independence here is relative to unit pullbacks alone; it is not a novelty claim relative to all tautological operations.

## Genus 2

The equation is `sum c_j S_j = 0`, with coefficients

```text
{20: 1, 29: -1}
```

**S_20, coefficient 1:**

```text
Graph :      [0, 2] [[1, 2, 3, 5], [6]] [(5, 6)]
Polynomial : (kappa_1)_0
```

**S_29, coefficient -1:**

```text
Graph :      [0, 0, 2] [[1, 2, 6], [3, 7, 8], [9]] [(6, 7), (8, 9)]
Polynomial : 1
```

## Genus 3

The equation is `sum c_j S_j = 0`, with coefficients

```text
{20: 1, 29: -1}
```

**S_20, coefficient 1:**

```text
Graph :      [0, 3] [[1, 2, 3, 5], [6]] [(5, 6)]
Polynomial : (kappa_1)_0
```

**S_29, coefficient -1:**

```text
Graph :      [0, 0, 3] [[1, 2, 6], [3, 7, 8], [9]] [(6, 7), (8, 9)]
Polynomial : 1
```

## Geometry of these representatives

In each genus the selected vector is S_20 − S_29. Let
\[
\xi:\overline{\mathcal M}_{0,4}\times\overline{\mathcal M}_{g,1}
\longrightarrow\overline{\mathcal M}_{g,3}
\]
attach the last rational marking to the genus-g marking. Then S_20 is xi_*(kappa_1 ⊠ 1). The other stratum is the chain of genera (0,0,g), with markings 1,2 on the first rational component and marking 3 on the middle one. There are no graph automorphism factors in these two examples.

The identity follows from kappa_1 = [a boundary point] on Mbar_(0,4) ≅ P^1 and the boundary splitting {1,2}|{3,*}. Thus these are known genus-zero gluing consequences, despite lying outside the span of unit pullbacks used in the comparison. No novelty relative to all Pixton consequences is asserted.
