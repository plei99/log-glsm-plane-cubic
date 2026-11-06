"""Export short, exactly certified four-mark relations beyond the comparison span."""
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_certify import suffix,previous,pullback_matrix,exact_span,pixton

records=[]
text='# Explicit four-mark equations\n\n'
text+='The equations below use the literal decorated-strata generators of the pinned admcycles basis. No extra automorphism factors are inserted. Each vector lies in the computed span and in the generated Pixton span, and lies outside the comparison span of unit pullbacks plus divisor products of the computed codimension-one relations. This last statement does not assert novelty relative to other tautological operations.\n\n'
text+='In the graph notation, the first list gives vertex genera, the second gives the half-edges at each vertex, and the final list pairs the half-edges forming nodes. Labels 1 through 4 are markings; all other labels are node half-edges. Polynomial subscripts index vertices starting at zero.\n\n'
for g in (2,3):
    tag=suffix(g,4,2);V=previous(g,4,2);M=pixton(g,4,2);rows=[]
    for v in (previous(g,3,2)*pullback_matrix(g,4,2)).rows():
        for i in range(1,5):
            mapping=list(range(1,5));mapping[i-1],mapping[3]=4,i
            rows.append(run.common.permute_vector(v,g,4,2,tuple(mapping)))
    R=run.common.converter(g).adm.TautologicalRing(g,4);divisors=R.generators(1)
    for v in previous(g,4,1).rows():
        relation=sum((c*divisors[j] for j,c in enumerate(v) if c),R.zero())
        rows.extend(vector(QQ,(relation*d).vector(2)) for d in divisors)
    B=exact_span(matrix(QQ,rows));pivots=B.pivots()
    outside=lambda v:v!=vector(QQ,[v[j] for j in pivots])*B
    v=min((v for v in V.rows() if outside(v)),key=lambda v:(len(v.nonzero_positions()),sum(len(str(c)) for c in v if c)))
    v*=lcm(c.denominator() for c in v);v/=gcd(ZZ(c) for c in v)
    if next(c for c in v if c)<0:v=-v
    u=M.transpose().solve_right(v);q=V.transpose().solve_right(v)
    assert u*M==v and q*V==v and outside(v)
    names=run.read(run.OUT/('strata-'+tag+'.json'))
    record={'g':g,'n':4,'r':2,'vector':run.common.sparse(v),
            'pixton_certificate':run.common.sparse(u),'found_certificate':run.common.sparse(q),
            'comparison_rank':B.nrows(),'outside_comparison_span':True,
            'strata':{str(i):names[i] for i in v.nonzero_positions()}}
    records.append(record)
    text+='## Genus %s\n\nThe equation is $\\sum_j c_j S_j=0$, with coefficients\n\n```text\n%s\n```\n\n'%(g,{i:ZZ(c) for i,c in enumerate(v) if c})
    for i,c in enumerate(v):
        if c:text+='S_%s, coefficient %s:\n\n```text\n%s\n```\n\n'%(i,c,names[i])
assert records[0]['vector']=={'233':'24','267':'-1'}
assert records[1]['vector']=={'28':'1','88':'-1'}
text+=r'''## Geometric interpretation

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
'''
run.dump(run.OUT/'sample-relations.json',records)
(run.OUT/'sample-relations.md').write_text(text)
print([(d['g'],len(d['vector'])) for d in records])
