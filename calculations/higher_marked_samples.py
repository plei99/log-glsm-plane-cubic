"""Export explicit three-mark relations beyond unit-pullback spans."""
import json
from pathlib import Path
from sage.all import *
OUT=Path('results/higher-marked/2026-11-06')
def sparse(v):return {str(i):str(x) for i,x in enumerate(v) if x}
records=[];text='# Explicit three-mark codimension-two relations\n\n'
text+='Each equation uses the literal decorated-strata generators listed below, with no extra automorphism factors. These representatives are obtained by exact linear algebra from the computed relation span and unit pullbacks. Their independence here is relative to unit pullbacks alone; it is not a novelty claim relative to all tautological operations.\n\n'
for g in (2,3):
    suffix='g%s-n3-r2'%g;V=load(str(OUT/('relations-'+suffix+'.sobj')))
    B=load(str(OUT/('unit-pullbacks-'+suffix+'.sobj')));M=load(str(OUT/('pixton-'+suffix+'.sobj')))
    span=V.stack(B).row_space().basis_matrix();unit=B.row_space()
    choices=[v for v in span.rows() if v not in unit]
    if not choices:continue
    v=min(choices,key=lambda x:(sum(bool(c) for c in x),sum(len(str(c)) for c in x)))
    v*=lcm(c.denominator() for c in v);v/=gcd(ZZ(c) for c in v)
    first=next(c for c in v if c)
    if first<0:v=-v
    u=M.transpose().solve_right(v);q=V.stack(B).transpose().solve_right(v)
    assert u*M==v and q*V.stack(B)==v and v not in unit
    names=json.loads((OUT/('strata-'+suffix+'.json')).read_text())
    record={'genus':g,'markings':3,'codimension':2,'vector':sparse(v),'pixton_certificate':sparse(u),
        'found_plus_unit_certificate':sparse(q),'found_rows':V.nrows(),'unit_rows':B.nrows(),
        'strata':{str(i):names[i] for i,c in enumerate(v) if c},'outside_unit_pullback_span':True}
    records.append(record);text+='## Genus %s\n\n'%g
    text+='The equation is `sum c_j S_j = 0`, with coefficients\n\n```text\n'+str({i:ZZ(c) for i,c in enumerate(v) if c})+'\n```\n\n'
    for i,c in enumerate(v):
        if c:text+='**S_%s, coefficient %s:**\n\n```text\n%s\n```\n\n'%(i,c,names[i])
text+=r'''## Geometry of these representatives

In each genus the selected vector is S_20 − S_29. Let
\[
\xi:\overline{\mathcal M}_{0,4}\times\overline{\mathcal M}_{g,1}
\longrightarrow\overline{\mathcal M}_{g,3}
\]
attach the last rational marking to the genus-g marking. Then S_20 is xi_*(kappa_1 ⊠ 1). The other stratum is the chain of genera (0,0,g), with markings 1,2 on the first rational component and marking 3 on the middle one. There are no graph automorphism factors in these two examples.

The identity follows from kappa_1 = [a boundary point] on Mbar_(0,4) ≅ P^1 and the boundary splitting {1,2}|{3,*}. Thus these are known genus-zero gluing consequences, despite lying outside the span of unit pullbacks used in the comparison. No novelty relative to all Pixton consequences is asserted.
'''
(OUT/'sample-relations.json').write_text(json.dumps(records,indent=2)+'\n')
(OUT/'sample-relations.md').write_text(text)
print([(r['genus'],len(r['vector'])) for r in records])
