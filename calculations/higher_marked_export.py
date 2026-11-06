"""Detailed report and compact explorer data for the higher-marking pass."""
from log_glsm.provenance import source_hashes
import gzip,json,hashlib
from pathlib import Path
from fractions import Fraction
OUT=Path('results/higher-marked/2026-11-06')
DOC=Path('docs/GENUS_TWO_THREE_MARKED_RESULTS_2026-11-06.md')


def geometry(t):
    if 'ambient_dimension' in t:return 'P^%s[%s]'%(t['ambient_dimension'],','.join(map(str,t['degrees'])))
    a=t['ambient']
    if 'n' in a:return 'Gr(2,%s)[%s]'%(a['n'],','.join(str(L[0]) for L in t['bundles']))
    return ' × '.join('F_1' if x=='F1' else 'P^'+str(x) for x in a['factors'])+' ['+'; '.join(str(tuple(L)) for L in t['bundles'])+']'


def main():
    ns=json.loads((OUT/'numerical-summary.json').read_text());cs=json.loads((OUT/'class-summary.json').read_text())
    records=[]
    for path in sorted(OUT.glob('*-numbers.json')):
        d=json.loads(path.read_text());d['slug']=d['target'].get('slug',d['target'].get('id'));d['geometry']=geometry(d['target'])
        with gzip.open(OUT/(d['slug']+'-g'+str(d['genus'])+'-classes.json.gz'),'rt') as f:c=json.load(f)
        d['ranks']={(str(x['markings'])+','+str(x['codimension'])):len(x['rows']) for x in c['relations']}
        records.append(d)
    text=r'''# Genera two and three with up to three markings

This extends the existing 43 ambient presentations to **g = 2,3 and n = 0,1,2,3**. The degree bound remains two: ordinary projective/Plücker degree for rank-one ambients, and Mori coefficient sum for toric products. Each labeled marking carries exactly one of U = tau_0(1), P = tau_1(1), or H = tau_0(H_X). P is the true stable-map cotangent insertion.

The numerical pass covers all 40 ordered profiles, since 1+3+9+27=40. The class search covers **codimensions zero, one, and two** on every Mbar_(g,n) in this box. It is not an exhaustive search through the top dimension 3g−3+n. Top numerical coefficients are calculated separately.

No new geometric assumption has been used to turn an unknown compact invariant into a number. Existing constant-map, isolated-curve, K3, and Gr(2,7) proofs are carried forward with their stated dependencies.

## The new three-mark numbers

Let Z_(g,beta) be the unmarked integrated ambient dual-Euler coefficient, let d = H_X·beta, and let k and h count P and H insertions. Repeated dilaton equations, followed by primary divisor equations, give
\[
\boxed{\langle P^k,H^h\rangle^{\rm twist}_{g,\beta}
=d^h\prod_{j=0}^{k-1}(2g-2+h+j)\,Z_{g,\beta}.}
\]
The order of the labels does not change this numerical value. There is no division by n!. The formula also holds for the ordinary compact numbers N_(g,beta) with the twist omitted.

The new numerical profiles at n=3 are therefore:

| Profile | Genus two multiplier | Genus three multiplier |
|---|---:|---:|
| P,P,P | 24 | 120 |
| P,P,H, in any order | 12d | 30d |
| P,H,H, in any order | 4d² | 6d² |
| H,H,H | d³ | d³ |
| Any profile containing U | 0 | 0 |

For example, three dilatons give (2g−2)(2g−1)(2g), rather than (2g−2)^3. For P,P,H the first forgotten dilaton sees two other marks, and the second sees one, giving (2g)(2g−1)d. The additional factors are essential.

The zero statement follows from the string equation by induction on the number of P insertions. Removing a U replaces one P by U, preserving at least one unit until no descendants remain. The last unit then integrates to zero. This is a statement about numerical integrals, not about the vanishing of the marked Chow class.

H is the usual hyperplane or Plücker class. On projective products it is the sum of factor hyperplanes; on F_1 × P^q it is S+2F+H_P. In the documented Mori basis this pairs as d_S+d_F+d_P. Every output records the actual divisor coefficients and intersection number.

## What is being localized

For a regular complete intersection Y=Z(s) in X cut out by a bundle E of rank c, define
\[
\Theta(t)=e_t((R\pi_*f^*E)^\vee),\qquad
\mathcal Z_{g,\beta,\mathbf I}(t)
=s_*\left(\Theta(t)\prod_iI_i\cap[\overline{\mathcal M}_{g,n}(X,\beta)]^{\rm vir}\right).
\]
This is the full derived dual-Euler ratio, including R^1, with one common CJR parameter t. It is not the ordinary GW class of the compact target Y.

With m=dim Y and a the number of P or H insertions,
\[
\alpha_0=m(g-1)-c_1(TY)\cdot\beta,\qquad
\mathcal Z(t)=\sum_r t^{\alpha_0+a-r}\mathcal Z_r,
\qquad \sigma_{g,\beta}=(-1)^{c(1-g)+c_1(E)\cdot\beta}.
\]
The numerical coefficient has t power alpha_0+a−(3g−3+n). All potentially nonzero profiles here have a=n, so this is the same power alpha_0−(3g−3) as the unmarked seed. Unit profiles have zero coefficient even when their candidate power differs.

The CJR identity determines the total stabilized infinity aggregate B through Z+B=sigma G. At negative t powers, B=−Z; at t^0, B=sigma N−Z numerically when N is known. The output never sets Z to zero merely because its t power is negative, and it does not split B into individual effective cycles.

The marked conventions and the arbitrary-rank sign follow the supplied October 17, 2025 CJR III draft, Theorem 8.4, Remark 8.5, §9.4, and §10.1. The [earlier marked source audit](GENUS_FOUR_MARKED_SOURCES.md), §§3–6, states these conventions and the general forgetting formulas used here. The extension from two to three marks is a repeated application of those formulas, not a new comparison theorem.

## Marked classes and rational tails

The implementation preserves the distinction
\[
\psi_i^{\rm map}=s^*\psi_i^{\rm curve}+T_i,
\]
where T_i is the positive-degree rational-tail contraction divisor. At an unstable marked fixed endpoint, the map descendant contributes −w. At a stable vertex it is the local cotangent class. A marked tail becomes a retained marking after stabilization; an empty tail becomes a temporary marking which is subsequently forgotten. Rational components that destabilize are contracted globally.

Consequently, a class with unit markings is generally a nonzero forgetful pullback, and an H insertion is not a class pulled back from the moduli of curves. Multiplying an unmarked Chow class by 24 or 120 would not produce its three-descendant counterpart. Those multipliers apply only after integration.

There are 20 canonical insertion multisets through n=3. The class calculation evaluates these and uses exact permutations of the labeled strata basis to recover all 40 ordered profiles. For relations, it also closes under every marking permutation. Permutations fixing an insertion profile can provide additional relations between its raw strata representatives. Such symmetry consequences are recorded explicitly.

## Relation certificates and independence

The candidates are negative auxiliary epsilon coefficients, finite-part differences between the base-four and base-five weight rays, and the stated marking-symmetry consequences. They are converted to documented decorated-strata coordinates before taking ranks.

For each (g,n,r), M is the generated Pixton matrix, V is the independent span found across targets, and U is an exact rational solution of
\[
\boxed{UM=V.}
\]
The verifier also checks that every stored per-target candidate basis row lies in the span of V. Scalar integration and agreement on two rays are not used as substitutes for a Chow certificate.

The comparison with unit pullbacks uses all label choices of the map Mbar_(g,n) → Mbar_(g,n−1), applied to the found relation space at the preceding marking level. If B is that inherited span, the number of found directions beyond it is rank(V+B)−rank(B). These directions need not be new to tautological geometry: psi products, clutching, or other known operations may also generate them, and the explicit Pixton certificate places them in the known relation system.

Three marks reach the setting of the classical Belorousski–Pandharipande codimension-two relation on Mbar_(2,3). That is a useful geometric benchmark. We do not identify a computed row with the authors' particular diagram normalization without a separate comparison; the exported strata vectors and certificates specify exactly which equations were obtained.

'''
    text+=f"The completed numerical pass contains **{ns['sectors']} genus-and-curve-class sectors** and **{ns['ordered_profiles']} ordered ambient numerical coefficients**. The class search retains **{cs['per_target_independent_rows']} per-target independent rows**, with a combined rank sum of **{cs['combined_rank_sum']}**, counted separately by (g,n,r).\n\n"
    text+='| g | n | Codimension | Found rank | Pixton rank | Strata columns |\n|---:|---:|---:|---:|---:|---:|\n'
    for r in cs['degrees']:
        if r['codimension']:
            text+='| {genus} | {markings} | {codimension} | {found_rank} | {pixton_rank} | {strata} |\n'.format(**r)
    text+='\nThe three-mark comparison with unit pullbacks is:\n\n| g | Codimension | Found rank | Unit-pullback rank | Directions beyond unit pullbacks |\n|---:|---:|---:|---:|---:|\n'
    for r in cs['unit_pullback_comparison']:
        if r['markings']==3:text+='| {genus} | {codimension} | {found_rank} | {unit_pullback_rank} | {found_directions_beyond_unit_pullbacks} |\n'.format(**r)
    text+=r'''

An explicit equation from each three-mark span is
\[
\boxed{\xi_*(\kappa_1\boxtimes 1)=\eta_*1\quad\text{in }R^2(\overline{\mathcal M}_{g,3}),\qquad g=2,3.}
\]
Here xi glues the last marking of Mbar_(0,4), whose other markings are 1,2,3, to the marking of Mbar_(g,1). The map eta is the boundary chain with genera (0,0,g): the first rational component carries marks 1 and 2, the middle carries mark 3, and the genus-g component is unmarked apart from its node. Both graphs have trivial automorphism group. Geometrically, this is the pushforward of kappa_1 = [a boundary point] on Mbar_(0,4) ≅ P^1. It illustrates why a relation can lie beyond unit pullbacks of our two-mark span while still being a familiar genus-zero gluing consequence. The exported vectors verify that this equation lies in the computed span in each genus.

The [explicit sample equations](../results/higher-marked/2026-11-06/sample-relations.md) give literal decorated graphs, and the adjacent JSON file includes their rational Pixton certificates and expressions in the found-plus-unit span. All 181 basis rows, not merely these two samples, are available in the `certificate-g*-n*-r*.json` files.

## Exact Gr(2,7) examples

For Y=Gr(2,7)[1^7], the previously audited ordinary degree-two values are
\[
N_{2,2}=\frac{539}{80},\qquad N_{3,2}=\frac{133}{288}.
\]
Here d=2. The three-mark ordinary numbers are

| Profile | Genus two | Genus three |
|---|---:|---:|
| P,P,P | 1617/10 | 665/12 |
| P,P,H | 1617/10 | 665/24 |
| P,H,H | 539/5 | 133/12 |
| H,H,H | 539/10 | 133/36 |

All permutations have the same numerical value; any profile with U is zero. These are ordinary GW numbers, before the CJR comparison sign. The signs at degree two are negative in genus two and positive in genus three.

The independently rerun ambient twisted seeds are instead
\[
Z_{2,2}=-\frac{19229}{240},\qquad
Z_{3,2}=-\frac{55586069}{414720}.
\]
For example, the corresponding three-dilaton ambient coefficients are −19229/10 and −55586069/3456. The discrepancy from the compact numbers is supplied by the CJR infinity aggregate.

For this same three-dilaton profile, applying the comparison signs gives the exact numerical infinity aggregates
\[
B_{2,2}(P,P,P)=\frac{8806}{5},\qquad
B_{3,2}(P,P,P)=\frac{55777589}{3456}.
\]
These are total stabilized aggregates, not values assigned to individual infinity vertices.

The compact inputs retain the prior proof dependency: 196 isolated lines and 1225 isolated conics, with the conic count extracted from the verified genus-zero value 2499/2. The incidence and normal-obstruction arguments in [the Gr(2,7) source audit](GRASSMANNIAN_27_SOURCES.md) show that degree-two images are exhausted by these curves. The local factor is |B_(2g)| k^(2g−3)/(2g(2g−2)!), as in Faber–Pandharipande.

## Unknown compact terms and all examples

There are **37 unknown ordinary unmarked sectors** across the two genera. Their nonzero-multiplier marked values remain unknown. Unit-containing profiles are nevertheless proved zero by the marking equations. The complete files distinguish `null` from the exact number zero. All ambient twisted coefficients were computed.

Ordinary K3 theory vanishes in the stated range, including constant-map classes by the separately audited obstruction calculation. For threefold constants, N_(g,0) = (−1)^g integral(c_3−c_1c_2) integral(lambda_g lambda_(g−1) lambda_(g−2)); for CY3s this becomes chi/5760 in genus two and −chi/1451520 in genus three. Positive-degree compact Fano integrals vanish by virtual dimension unless that dimension is zero; exceptional dimension-zero cases without an independent calculation remain unknown. These numerical statements do not erase nonzero classes in lower codimension.

The following tables give the seed numbers for every presentation and curve class. Apply the displayed marking formula to recover any profile, or read the explicit ordered values in the `*-numbers.json` files. Question marks mean unknown. Compact targets with different ambient presentations are kept separate because their twists and infinity aggregates differ.

'''
    by_slug={}
    for d in records:by_slug.setdefault(d['slug'],{})[d['genus']]=d
    for slug,pair in by_slug.items():
        a,b=pair[2],pair[3]
        text+=f"### {slug}: {a['geometry']}\n\n"
        text+=f"Polarization coefficients: {a['polarization']}. Three-mark codimension-two relation ranks: {a['ranks']['3,2']} in genus two and {b['ranks']['3,2']} in genus three.\n\n"
        text+='| beta | Z₂ | Z₃ | N₂ | N₃ |\n|---|---|---|---|---|\n'
        for x,y in zip(a['rows'],b['rows']):
            assert x['beta']==y['beta']
            text+='| %s | %s | %s | %s | %s |\n'%(tuple(x['beta']),x['twisted_unmarked'],y['twisted_unmarked'],
                x['ordinary_unmarked'] if x['ordinary_unmarked'] is not None else '?',y['ordinary_unmarked'] if y['ordinary_unmarked'] is not None else '?')
        text+='\n'
    text+=r'''## Code and validation

`higher_marked_complete_intersections.sage` exposes the genus-two/three catalogue, ordinary inputs, canonical profiles, and permutation maps. The shared backend in `genus_four_marked.sage` now accepts three marks. Its old n≤2 behavior is retained, and its exact previous source is archived as `baseline-genus-four-backend.sage.txt` with the prior provenance hash.

The tests independently check three-mark fixed-locus integrals against string/dilaton/divisor reductions in both genera, unit pullback and dilaton identities as Chow classes, label permutations, and the known Gr(2,7) values. The complete prior genus-four regression suite also passes. Every numerical sector is checked on two auxiliary rays, with all negative auxiliary poles cancelling.

The source formulas and Hodge recursion are unchanged: full derived dual Euler twists, Mumford reduction, exact numerical GRR recursion, and finite Laurent jets. The larger marking count is handled by actual graph stabilization. The numerical marking formulas do not stand in for that class calculation.

Commands, from the repository root inside the pinned Sage Linux environment:

```bash
sage tests/test_higher_marked.sage
sage tests/test_genus_four_marked.sage
python calculations/higher_marked_numbers.py --workers 4 --resume
python calculations/higher_marked_classes.py --workers 8 --resume
python calculations/higher_marked_certify.py
python calculations/higher_marked_export.py
python calculations/higher_marked_verify.py
python calculations/higher_marked_samples.py
python calculations/higher_marked_visual.py
```

Class checkpoints are atomic and locked per target/genus. A completed checkpoint is reused; remove the corresponding file deliberately to recompute that case. The large Gr(2,7), genus-two calculation can instead be split with `python calculations/higher_marked_split.py --index 42 --genus 2 --workers 8`, then assembled with the same command and `--merge-only`. Its profile spans are merged by exact rational elimination. Run the merge after other writers for that case have stopped. The logs include interrupted runs and completed continuations. Reported results come from completed target artifacts and final certificate replay, not partial logs.

## References

1. [Coates–Givental, Quantum Riemann–Roch, Lefschetz and Serre](https://math.berkeley.edu/~giventh/papers/rls.pdf), Appendices 1–2: forgetting compatibility of the twisting index, node normalization, and descendant/ancestor comparison. The numerical three-mark formulas above are deductions by repeated string, dilaton, and divisor equations.
2. [Pandharipande, Rational curves on hypersurfaces](https://numdam.org/item/SB_1997-1998__40__307_0.pdf), §1.2, for these marking equations and their higher-genus form.
3. Chen–Janda–Ruan, supplied private CJR III draft, October 17, 2025, Theorem 8.4, Remark 8.5, §9.4, §10.1; precise local references are in [the marked source audit](GENUS_FOUR_MARKED_SOURCES.md). For the public comparison, see [CJR, The logarithmic gauged linear sigma model](https://arxiv.org/abs/1906.04345), Proposition 4.1, and [Chen–Janda–Webb](https://arxiv.org/abs/1911.09825), Corollary 1.2.1.
4. [Graber–Pandharipande, Localization of virtual classes](https://arxiv.org/abs/alg-geom/9708001), for the graph calculation; [Mumford, Towards an Enumerative Geometry](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1983b--EnumGeomModuli-NC.pdf), §5, for Hodge GRR.
5. [Faber–Pandharipande, Hodge integrals and Gromov–Witten theory](https://arxiv.org/abs/math/9810173), Theorems 3–4, for the local rational-curve and constant-map formulas.
6. [Janda, Relations on Mbar(g,n) via equivariant GW theory of P1](https://content.algebraicgeometry.nl/2017-3/2017-3-018.pdf), Theorem 1.1: Chow validity of the Pixton generators. The saved rational matrices separately certify membership of the particular relations found here.
7. [Belorousski–Pandharipande, A descendent relation in genus 2](https://arxiv.org/abs/math/9803072), Ann. Scuola Norm. Sup. Pisa 29 (2000), 171–191: the classical codimension-two three-point relation, cited as context rather than an unverified normalization match.
'''
    DOC.write_text(text)
    compact=[{'id':d['slug'],'genus':d['genus'],'name':d['geometry'],
        'rows':[[r['beta'],r['H_degree'],r['alpha_0'],r['sign'],r['twisted_unmarked'],r['ordinary_unmarked']] for r in d['rows']]}
        for d in records]
    (OUT/'explorer-data.json').write_text(json.dumps({'examples':compact,'relations':cs,'numbers':ns},indent=2)+'\n')
    sources=['genus_four_marked.sage','higher_marked_complete_intersections.sage','tests/test_higher_marked.sage']
    sources += [str(p) for p in Path('calculations').glob('higher_marked_*.py')]+[str(DOC)]
    (OUT/'provenance.json').write_text(json.dumps({'files':source_hashes(sources),
        'baseline_genus_four_sha256':hashlib.sha256((OUT/'baseline-genus-four-backend.sage.txt').read_bytes()).hexdigest(),
        'scope':{'genera':[2,3],'max_markings':3,'curve_bound':2,'max_chow_codimension':2}},indent=2)+'\n')
    print(DOC)


if __name__=='__main__':main()
