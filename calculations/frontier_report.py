"""Assemble the completed marking/ambient report and compact visualization data."""
from log_glsm.provenance import source_hashes
import gzip,json,hashlib
from pathlib import Path
OUT=Path('results/marked-frontier/2026-11-06')
DOC=Path('docs/MARKED_FRONTIER_RESULTS_2026-11-06.md')
read=lambda p:json.loads(p.read_text())
numbers=read(OUT/'numerical-summary.json');base=read(OUT/'baseline-extension.json')
four=read(OUT/'summary-n4-r2.json');five=read(OUT/'summary-n5-r2.json')
assert four['complete'] and five['complete']
verification=read(OUT/'verification.json');memory=read(OUT/'memory-summary.json')
probe=read(OUT/'codimension-three-probe.json')
spaces=base['spaces']+four['spaces']+five['spaces'];records=[];methods={4:{},5:{}}
for p in sorted((OUT/'numbers').glob('*-n5.json.gz')):
    with gzip.open(p,'rt') as f:d=json.load(f)
    t=d['target'];slug=t.get('slug',t.get('id'));g=d['genus']
    if 'ambient_dimension' in t:name='P^%s[%s]'%(t['ambient_dimension'],','.join(map(str,t['degrees'])))
    else:name=t['ambient']['slug']+' ['+'; '.join(str(tuple(L)) for L in t['bundles'])+']'
    for n in (4,5):
        with gzip.open(OUT/('n%s-r2'%n)/(slug+'-g%s.json.gz'%g),'rt') as f:c=json.load(f)
        key='unit reconstruction' if c.get('unit_profiles_reconstructed') else 'all profiles by localization'
        methods[n][key]=methods[n].get(key,0)+1
    records.append({'id':slug,'name':name,'genus':g,'target':t,'rows':[{k:v for k,v in x.items() if k!='profiles'} for x in d['rows']]})
text=r'''# Higher markings, non-Fano toric ambients, and partial flags

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
'''
for s in spaces:
    if s['r']:text+='| {g} | {n} | {r} | {found} | {pixton} | {strata} |\n'.format(**s)
text+='\nThe higher-mark comparison is:\n\n| g | n | r | Unit-pullback rank | Beyond unit pullbacks | Unit plus divisor-product rank | Beyond that span |\n|---:|---:|---:|---:|---:|---:|---:|\n'
for s in four['spaces']+five['spaces']:
    if s['r']==2:text+='| {g} | {n} | {r} | {unit_pullback} | {beyond_unit_pullback} | {unit_and_divisor_products} | {beyond_unit_and_divisor_products} |\n'.format(**s)
text+=f"\nThe separate codimension-three probe found {probe['found_rank']} independent rows in a generated Pixton space of rank {probe['pixton_rank']}, with {probe['strata']} strata columns. Its exact certificate is saved separately. The new ambient presentations did not enlarge the old n≤3 spans: this is checked in `baseline-extension.json`.\n"
text+=r'''
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

'''
text+=f"The numerical files contain {numbers['ordered_profiles']:,} ordered profile records in {numbers['sectors']} genus/curve-class sectors. They reproduce {numbers['previous_seeds_matched']} previous seeds. There are {numbers['unknown_ordinary_sectors']} unknown ordinary seed sectors; nonzero-multiplier marked values remain null.\n\n"
text+=f"The certificate replay checked {verification['pixton_basis_rows_replayed']} independent basis rows, counted separately over (g,n,r), and {verification['per_target_rows_replayed']:,} per-target rows, in addition to the {probe['found_rank']}-row higher-codimension probe. At four marks, {methods[4].get('all profiles by localization',0)} target/genus cases localized every profile directly, and {methods[4].get('unit reconstruction',0)} used unit reconstruction. At five marks, all {methods[5].get('unit reconstruction',0)} cases used direct non-unit profiles followed by unit reconstruction.\n\n"
text+=f"The host has {memory['physical_gib']:.0f} GiB RAM; the guest exposes {memory['guest_total_gib']:.2f} GiB. The largest recorded worker used {memory['peak_worker_gib']:.2f} GiB RSS. Minimum sampled guest available memory was {memory['minimum_guest_available_gib']:.2f} GiB. Guest swap-out counters changed by {memory['guest_swapout_page_change']} pages and OOM-kill counters by {memory['guest_oom_kill_change']}. Host compressed-memory occupancy grew from {memory['initial_compressor_gib']:.2f} GiB to a sampled peak of {memory['peak_host_compressor_gib']:.2f} GiB; host swap-out traffic since the initial snapshot was {memory['host_swapout_gib_since_initial_snapshot']:.2f} GiB. Host counters include other applications.\n\n"
text+=r'''Memory constrained concurrency: the five-mark pool was stopped at atomic checkpoints when the large four-mark Grassmannian pool left about 4.5 GiB guest memory available. It resumed with two workers once available memory recovered to about 11 GiB. The largest five-mark Grassmannian jobs wait for the corresponding four-mark checkpoint. Later, several five-mark workers were temporarily suspended with their in-memory state intact when guest headroom reached about 2 GiB, then resumed after completed partitions released their caches and headroom recovered above 5 GiB. This scheduling cost wall time. No completed results were lost. RSS sums double-count shared pages and are not physical-memory totals. The guest kernel does not expose PSI; this absence is recorded. See `memory-summary.json` and the raw host/guest logs for limits of the measurements.

Direct four/five-mark localization, unit and dilaton Chow identities, section-map identities, old three-mark and genus-four regressions, new incidence Chern calculations, normal endpoint characters, and independent genus-zero counts all pass. A fresh full-profile NF4 run reproduces all 43 stored finite class vectors and the codimension-one/two relation rows literally, using the final graph loop and CLI. Separate NF4 and PF4 class replays also agree exactly. Source changes and reports are hashed in `provenance.json`; previous versions of the shared backends are preserved in `baseline/`. Interrupted logs are retained and are not treated as completed results.

## Every presentation and seed

The following tables give all unmarked inputs. The marking formula recovers every higher-mark numerical value; the explicit 364 ordered profiles per sector are stored in `numbers/*-n5.json.gz`. A question mark means an unknown ordinary invariant.

'''
groups={}
for d in records:groups.setdefault(d['id'],{})[d['genus']]=d
for slug,pair in groups.items():
    a,b=pair[2],pair[3];text+='### '+slug+': '+a['name']+'\n\n'
    text+='| beta | Z₂ | Z₃ | N₂ | N₃ |\n|---|---|---|---|---|\n'
    for x,y in zip(a['rows'],b['rows']):
        assert x['beta']==y['beta']
        text+='| %s | %s | %s | %s | %s |\n'%(tuple(x['beta']),x['twisted_unmarked'],y['twisted_unmarked'],x['ordinary_unmarked'] or '?',y['ordinary_unmarked'] or '?')
    text+='\n'
text+=r'''## Reproduction and references

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
'''
DOC.write_text(text)
(OUT/'explorer-data.json').write_text(json.dumps({'records':records,'spaces':spaces,'numbers':numbers,'memory':memory,'probe':{k:v for k,v in probe.items() if k not in ('U','V')}},indent=2)+'\n')
sources=[str(p) for p in Path('.').glob('*.sage')]
sources += [str(p) for p in Path('calculations').glob('frontier_*.py')]+['calculations/marked_frontier.py',str(DOC),'docs/CONVEX_GKM_GEOMETRY_2026-11-06.md']
(OUT/'provenance.json').write_text(json.dumps({'files':source_hashes(sources),
    'scope':{'presentations':54,'genera':[2,3],'max_marks':5,'curve_bound':2,'chow_codimensions':[0,1,2],
             'additional_probe':{'target':'elliptic cubic','g':3,'n':3,'r':3}}},indent=2)+'\n')
print(DOC)
