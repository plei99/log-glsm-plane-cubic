"""Replay Grassmannian certificates and compare the prior projective span."""
import gzip,hashlib,json
from pathlib import Path
from sage.all import *

from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersectionCycles,
    grassmannian_known_ordinary_gw,
    grassmannian_line_count,
    grassmannian_sample,
)
from log_glsm.targets.toric import (
    toric_degree_tuples,
    toric_genus_zero_number,
)
OUT=Path('results/grassmannian-complete-intersections/2026-11-06')
OLD=Path('results/complete-intersections/2026-11-06')

def read(p):
    p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def dump(p,data,compressed=False):
    raw=(json.dumps(data,indent=None if compressed else 2,default=int)+'\n').encode()
    Path(str(p)+('.gz' if compressed else '')).write_bytes(gzip.compress(raw,mtime=0) if compressed else raw)

targets=grassmannian_sample();catalog=[read(OUT/t.slug/'summary.json') for t in targets]
assert all(s['mori_bound']==2 and s['weight_bases']==[4,5] and len(s['genera'])==2 for s in catalog)
MATS={(g,r):load(str(OUT/('pixton-g%s-r%s.sobj'%(g,r)))) for g in (2,3) for r in range(3*g-2)}
rows={key:[] for key in MATS};certificates=0;symmetries=[];finite={}
for t,s in zip(targets,catalog):
    for gs in s['genera']:
        g=gs['genus'];data=read(OUT/t.slug/('finite-classes-g%s.json.gz'%g))
        assert len(data)==len(toric_degree_tuples(t.ambient))*(3*g-2)
        finite[t.slug,g]={(tuple(x['beta']),x['codimension']):x for x in data}
        for d in gs['degrees']:
            r=d['codimension'];M=MATS[g,r];records=read(OUT/t.slug/('relations-g%s-r%s.json.gz'%(g,r)))
            assert len(records)==d['relation_rank']
            if not records:continue
            V=matrix(QQ,[rec['strata_vector'] for rec in records]);U=matrix(QQ,len(records),M.nrows(),sparse=True)
            for i,rec in enumerate(records):
                for j,co in rec['pixton_weights']:U[i,j]=QQ(co)
            assert U*M==V and V.rank()==len(records)
            rows[g,r].extend(V.rows());certificates+=len(records)
    print('Replayed',t.slug,flush=True)


union=[];readable=['All expressions equal zero in rational Chow. S_i refers to the dictionary for its genus and codimension.\n']
for (g,r),M in MATS.items():
    V=matrix(QQ,rows[g,r]) if rows[g,r] else matrix(QQ,0,M.ncols());B=V.row_space().basis_matrix()
    records=[]
    if B.nrows():
        U=M.transpose().solve_right(B.transpose()).transpose();assert U*M==B
        for v,u in zip(B.rows(),U.rows()):records.append({'strata_vector':[str(c) for c in v],
            'pixton_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
    dump(OUT/('union-relations-g%s-r%s.json'%(g,r)),records,True)
    assert read(OLD/('strata-g%s-r%s.json'%(g,r)))==read(OUT/('strata-g%s-r%s.json'%(g,r)))
    old=read(OLD/('union-relations-g%s-r%s.json.gz'%(g,r)))
    A=matrix(QQ,[x['strata_vector'] for x in old]) if old else matrix(QQ,0,M.ncols())
    combined=A.stack(B).row_space().basis_matrix();added=combined.nrows()-A.nrows()
    if not added and B.nrows():
        W=A.transpose().solve_right(B.transpose()).transpose();assert W*A==B
        dump(OUT/('previous-span-g%s-r%s.json'%(g,r)),
             {'genus':g,'codimension':r,'weights':[[[int(i),str(c)] for i,c in enumerate(w) if c] for w in W.rows()]},True)
    # Save a certificate for the combined span, rather than only a dimension.
    combined_records=[]
    if combined.nrows():
        C=M.transpose().solve_right(combined.transpose()).transpose();assert C*M==combined
        for v,u in zip(combined.rows(),C.rows()):combined_records.append({'strata_vector':[str(c) for c in v],
            'pixton_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
    dump(OUT/('combined-relations-g%s-r%s.json'%(g,r)),combined_records,True)
    union.append({'genus':g,'codimension':r,'rank':B.nrows(),'pixton_rank':M.rank(),
                  'previous_rank':A.nrows(),'combined_rank':combined.nrows(),'additional_rank':added})
    readable.append('\nGenus %s, codimension %s\n'%(g,r))
    for i,s in enumerate(read(OUT/('strata-g%s-r%s.json'%(g,r)))):readable.append('S_%s = %s\n'%(i,s))
    for i,v in enumerate(B.rows()):readable.append('R_%s: '%i+' + '.join('(%s)*S_%s'%(c,j) for j,c in enumerate(v) if c)+' = 0\n')
    print('Union',g,r,B.nrows(),'previous',A.nrows(),'combined',combined.nrows(),flush=True)
(OUT/'independent-relations.txt').write_text(''.join(readable))
dump(OUT/'union-summary.json',{'individual_certificates':certificates,
    'degrees':union,'total_union_rank':sum(d['rank'] for d in union),
    'total_combined_rank':sum(d['combined_rank'] for d in union),'additional_rank':sum(d['additional_rank'] for d in union)})

# Direct Hodge integration is independent of the Pixton quotient integration.
checks=[]
for t,s in zip(targets,catalog):
    b=GrassmannianCompleteIntersectionCycles(t,3)
    actual=b.top_integral((1,));expected=QQ(s['genera'][1]['top_twisted_integrals']['1'])
    assert actual==expected,(t.slug,1)
    checks.append({'target':t.slug,'genus':3,'beta':[1],'integral':str(actual)})
for i in (3,10,11):
    t=targets[i]
    for g in (2,3):
        b=GrassmannianCompleteIntersectionCycles(t,g)
        actual=b.top_integral((2,));expected=QQ(catalog[i]['genera'][g-2]['top_twisted_integrals']['2'])
        assert actual==expected,(t.slug,g,2)
        checks.append({'target':t.slug,'genus':g,'beta':[2],'integral':str(actual)})
print('Direct Hodge checks passed',len(checks),flush=True)

ordinary=[];genus_zero=[]
for t,s in zip(targets,catalog):
    if t.kind!='CY3':continue
    for beta in toric_degree_tuples(t.ambient):
        name=','.join(map(str,beta))
        if any(beta):
            a=toric_genus_zero_number(t,beta,4);b=toric_genus_zero_number(t,beta,5);assert a==b
            genus_zero.append({'target':t.slug,'beta':list(beta),'ordinary_gw':str(a),
                'multiple_cover_subtracted':str(a if beta==(1,) else a-grassmannian_line_count(t)/8)})
        for gs in s['genera']:
            g=gs['genus'];value=grassmannian_known_ordinary_gw(t,g,beta)
            if value is not None:
                Z=QQ(gs['top_twisted_integrals'][name]);sign=t.comparison_sign(g,beta)
                ordinary.append({'target':t.slug,'genus':g,'beta':list(beta),'ordinary_gw':str(value),
                    'comparison_sign':str(sign),'twisted_top':str(Z),'infinity_aggregate_integral':str(sign*value-Z)})
dump(OUT/'ordinary-cy3-gw.json',ordinary);dump(OUT/'genus-zero-cy3-gw.json',genus_zero)
dump(OUT/'verification.json',{'individual_certificates_replayed':certificates,
    'union_rank':sum(d['rank'] for d in union),
    'previous_span_certificates_replayed':sum(d['rank'] for d in union if d['additional_rank']==0),
    'combined_rank':sum(d['combined_rank'] for d in union),'direct_hodge_checks':checks,
    'genus_zero_two_ray_checks':len(genus_zero)})

lines=['| Target | Ambient | Bundle coefficients | Type | Euler characteristic |\n|---|---|---|---|---:|\n']
for t in targets:lines.append('| %s | %s | %s | %s | %s |\n'%(t.slug,t.ambient.slug,str(t.bundles),t.kind,t.chern_numbers()['euler']))
lines+=['\nPer-target independent localization relations (two auxiliary torus specializations):\n\n',
        '| Target | Genus 2 ranks r=1,2,3 | Genus 3 ranks r=2,3,4,5,6 |\n|---|---|---|\n']
for t,s in zip(targets,catalog):
    a,b=s['genera'];lines.append('| %s | %s | %s |\n'%(t.slug,', '.join(str(x['relation_rank']) for x in a['degrees'][1:]),', '.join(str(x['relation_rank']) for x in b['degrees'][2:])))
lines+=['\n| Genus | Codimension | Grassmannian rank | Previous projective rank | Combined rank | Pixton rank |\n|---|---|---:|---:|---:|---:|\n']
for d in union:
    if d['pixton_rank']:lines.append('| %(genus)s | %(codimension)s | %(rank)s | %(previous_rank)s | %(combined_rank)s | %(pixton_rank)s |\n'%d)
for t,s in zip(targets,catalog):
    lines+=['\nExact twisted top integrals for %s (%s):\n\n'%(t.slug,t.ambient.slug),'| Mori tuple β | Genus 2 | Genus 3 |\n|---|---|---|\n']
    for beta in toric_degree_tuples(t.ambient):
        name=','.join(map(str,beta));lines.append('| (%s) | %s | %s |\n'%(name,s['genera'][0]['top_twisted_integrals'][name],s['genera'][1]['top_twisted_integrals'][name]))
lines+=['\nOrdinary genus-zero CY3 invariants, checked at both rays:\n\n| Target | Mori tuple | N₀,β |\n|---|---|---|\n']
for d in genus_zero:lines.append('| %s | %s | %s |\n'%(d['target'],str(tuple(d['beta'])),d['ordinary_gw']))
lines+=['\nKnown ordinary higher-genus CY3 numbers and integrated infinity aggregates:\n\n| Target | Genus | Mori tuple | N_g,β | ∫B_β,3g−3 |\n|---|---|---|---|---|\n']
for d in ordinary:lines.append('| %s | %s | %s | %s | %s |\n'%(d['target'],d['genus'],str(tuple(d['beta'])),d['ordinary_gw'],d['infinity_aggregate_integral']))
(OUT/'tables.md').write_text(''.join(lines))
print('VERIFIED AND EXPORTED',sum(d['rank'] for d in union),'Grassmannian relations;',sum(d['additional_rank'] for d in union),'additional directions',flush=True)
