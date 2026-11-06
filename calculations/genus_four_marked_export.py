"""Readable mathematics and a self-contained explorer for the completed batch."""
from log_glsm.provenance import source_hashes
import gzip,html,json,hashlib
from fractions import Fraction
from pathlib import Path

OUT=Path('results/genus-four-marked/2026-11-06')
DOC=Path('docs/GENUS_FOUR_MARKED_RESULTS_2026-11-06.md')


def geometry(t):
    if 'ambient_dimension' in t:return 'P^%s[%s]'%(t['ambient_dimension'],','.join(map(str,t['degrees'])))
    a=t['ambient']
    if 'n' in a:return 'Gr(2,%s)[%s]'%(a['n'],','.join(str(L[0]) for L in t['bundles']))
    ambient=' × '.join('F_1' if f=='F1' else 'P^'+str(f) for f in a['factors'])
    return ambient+' ['+'; '.join('('+','.join(map(str,L))+')' for L in t['bundles'])+']'


def main():
    numerical=json.loads((OUT/'numerical-summary.json').read_text())
    classes=json.loads((OUT/'class-summary.json').read_text())
    records=[json.loads((OUT/name).read_text()) for name in numerical['files']]
    for d in records:
        t=d['target'];slug=t.get('slug',t.get('id'));d['slug']=slug;d['geometry']=geometry(t)
        with gzip.open(OUT/(slug+'-classes-r3.json.gz'),'rt') as f:c=json.load(f)
        d['relation_ranks']=[{'n':r['markings'],'r':r['codimension'],'rank':len(r['rows'])} for r in c['relations']]
    prose=r'''# Genus four, up to two markings: exact calculations and bounded Chow relations

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

'''
    prose+=f"The numerical batch contains **{numerical['sectors']} sectors and {numerical['ordered_profile_values']} ordered profile values**, checked on two rays. There are **{numerical['unknown_ordinary_unmarked']} unknown ordinary unmarked sectors**; their numerical marking reductions remain explicitly symbolic. The relation batch retains **{classes['individual_target_rows']} per-target independent rows**, whose combined rank sum is **{classes['combined_rank_sum']}**, counted separately by marking number and codimension.\n\n"
    prose+='| Markings n | Codimension r | Strata columns | Found rank | Pixton rank |\n|---:|---:|---:|---:|---:|\n'
    for row in classes['degrees']:
        if row['codimension']>=2:prose+='| {markings} | {codimension} | {strata} | {combined_rank} | {pixton_rank} |\n'.format(**row)
    prose+=r'''

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

'''
    for i,d in enumerate(records,1):
        slug=d['slug'];t=d['target']
        prose+=f"### {i}. {d['geometry']} ({slug}; {d['kind']})\n\n"
        prose+=f"Target dimension {t['dimension']}; chosen H coefficients {d['polarization']}. "
        ranks=[x for x in d['relation_ranks'] if x['rank']]
        prose+='Per-target relation ranks '+(', '.join(f"(n={x['n']}, r={x['r']}): {x['rank']}" for x in ranks) if ranks else 'zero in the declared window')+'.\n\n'
        prose+='| Ambient class beta | H·beta | alpha_0 | Sign | Z | N |\n|---|---:|---:|---:|---|---|\n'
        for row in d['rows']:
            prose+=f"| {tuple(row['beta'])} | {row['H_degree']} | {row['alpha_0']} | {row['comparison_sign']} | {row['twisted_unmarked']} | {row['ordinary_unmarked'] if row['ordinary_unmarked'] is not None else '?'} |\n"
        prose+='\n'
    prose+=r'''## Validation and reproducibility

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
'''
    DOC.write_text(prose)
    # An independently readable, exact compactified sample.
    sample=[-1920,180,-312,1260,-120,300,84,360,72,-660,-28,140,-28,20,52,20,1]
    names=json.loads((OUT/'strata-n0-r2.json').read_text())
    sample_record={'genus':4,'markings':0,'codimension':2,'coefficients':sample,'strata':names,
        'origin':{'target':'p2-3','beta':[2],'base':4,'epsilon_power':-2},'smooth_locus':'3 kappa_1^2 = 32 kappa_2'}
    (OUT/'sample-relation.json').write_text(json.dumps(sample_record,indent=2)+'\n')
    sample_doc='# A compactified genus-four codimension-two relation\n\nThe equation is sum_j c_j S_j = 0. Each S_j is the literal admcycles generator below; no additional automorphism factor is inserted. The labeling and node pairs fully specify each graph.\n\n'
    for j,(co,name) in enumerate(zip(sample,names)):sample_doc+=f'**S_{j}, coefficient {co}:**\n\n```text\n{name}\n```\n\n'
    sample_doc+='Only S_0 = kappa_2 and S_1 = kappa_1^2 survive on the smooth locus. Thus 180 kappa_1^2 − 1920 kappa_2 = 0, equivalently 3 kappa_1^2 = 32 kappa_2.\n'
    (OUT/'sample-relation.md').write_text(sample_doc)
    (OUT/'explorer-data.json').write_text(json.dumps({'examples':records,'relations':classes},indent=2)+'\n')
    sources=['genus_four_marked.sage','tests/test_genus_four_marked.sage','calculations/genus_four_marked_scan.py',
        'calculations/genus_four_marked_relations.py','calculations/genus_four_marked_certify.py','calculations/genus_four_marked_export.py',
        'calculations/genus_four_marked_verify.py','calculations/genus_four_marked_visual.py','docs/GENUS_FOUR_MARKED_SOURCES.md',str(DOC)]
    (OUT/'provenance.json').write_text(json.dumps({'files':source_hashes(sources),
        'scope':{'g':4,'n_max':2,'degree_bound':2,'class_codimension_max':3},'ordinary_unknown_sectors':numerical['unknown_ordinary_unmarked']},indent=2)+'\n')
    print(DOC)


if __name__=='__main__':main()
