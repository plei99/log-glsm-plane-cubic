"""Replay Gr(2,7) certificates, compare the earlier span, and integrate directly."""
import gzip,json
from pathlib import Path
from sage.all import *

from log_glsm.core.hodge import (
    AdmcyclesHodgeIntegralBackend,
)
from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersectionCycles,
    grassmannian_known_ordinary_gw,
    grassmannian_line_count,
    grassmannian_seven_hyperplane_cy3,
)
from log_glsm.targets.toric import (
    toric_genus_zero_number,
)
OUT=Path('results/grassmannian-27/2026-11-06')
OLD=Path('results/grassmannian-complete-intersections/2026-11-06')

def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_text())

def dump(path,data,compressed=False):
    raw=(json.dumps(data,indent=None if compressed else 2,default=int)+'\n').encode()
    Path(str(path)+('.gz' if compressed else '')).write_bytes(gzip.compress(raw,mtime=0) if compressed else raw)

t=grassmannian_seven_hyperplane_cy3();summary=read(OUT/'G15/summary.json')
assert len(summary['genera'])==2
degrees=[];certificates=0;vanishings=0;zero_dimensional_slots=0;finite={};readable=[]
for gs in summary['genera']:
    g=gs['genus'];finite[g]=read(OUT/'G15'/('finite-classes-g%s.json.gz'%g))
    for info in gs['degrees']:
        r=info['codimension'];M=load(str(OUT/('pixton-g%s-r%s.sobj'%(g,r))))
        strata=read(OUT/('strata-g%s-r%s.json'%(g,r)))
        assert strata==read(OLD/('strata-g%s-r%s.json'%(g,r)))
        records=read(OUT/'G15'/('relations-g%s-r%s.json.gz'%(g,r)))
        V=matrix(QQ,[rec['strata_vector'] for rec in records]) if records else matrix(QQ,0,M.ncols())
        U=matrix(QQ,len(records),M.nrows())
        for i,rec in enumerate(records):
            for j,c in rec['pixton_weights']:U[i,j]=QQ(c)
        assert U*M==V and V.rank()==info['relation_rank']==len(records)
        old=read(OLD/('union-relations-g%s-r%s.json.gz'%(g,r)))
        A=matrix(QQ,[rec['strata_vector'] for rec in old]) if old else matrix(QQ,0,M.ncols())
        added=A.stack(V).rank()-A.rank()
        if added==0 and V.nrows():
            W=A.transpose().solve_right(V.transpose()).transpose();assert W*A==V
            dump(OUT/('previous-span-g%s-r%s.json'%(g,r)),
                 {'weights':[[[int(i),str(c)] for i,c in enumerate(w) if c] for w in W.rows()]},True)
        degrees.append({'genus':g,'codimension':r,'rank':V.nrows(),'previous_rank':A.rank(),
                        'additional_rank':added,'pixton_rank':M.rank()})
        certificates+=len(records)
        readable.append('\nGenus %s, codimension %s\n'%(g,r))
        readable.extend('S_%s = %s\n'%(i,s) for i,s in enumerate(strata))
        readable.extend('R_%s: '%i+' + '.join('(%s)*S_%s'%(c,j) for j,c in enumerate(v) if c)+' = 0\n' for i,v in enumerate(V.rows()))
    for rec in finite[g]:
        D=rec['beta'][0]
        if 7*(D+1-g)<0:
            v=vector(QQ,rec['strata_vector']);r=rec['codimension']
            M=load(str(OUT/('pixton-g%s-r%s.sobj'%(g,r))))
            assert v in M.row_space();vanishings+=1
        if D+1-g==0 and rec['codimension']<3*g-3:
            v=vector(QQ,rec['strata_vector']);r=rec['codimension']
            M=load(str(OUT/('pixton-g%s-r%s.sobj'%(g,r))))
            assert v in M.row_space();zero_dimensional_slots+=1
print('Certificates replayed',certificates,'additional rank',sum(d['additional_rank'] for d in degrees),flush=True)
dump(OUT/'union-summary.json',{'degrees':degrees,'total_rank':certificates,'additional_rank':sum(d['additional_rank'] for d in degrees)})
(OUT/'independent-relations.txt').write_text(''.join(readable))

# Integrate the saved formal fixed-locus terms directly. This does not use
# their stabilized strata vectors or the Pixton integration functional.
hodge=AdmcyclesHodgeIntegralBackend();checks=[]
for g in (2,3):
    for D in (1,2):
        rec=next(rec for rec in finite[g] if rec['beta']==[D] and rec['codimension']==3*g-3)
        actual=QQ.zero()
        for term in rec['terms']:
            kind,data,powers,lambdas=term['key'];assert kind=='graph'
            actual+=QQ(term['coefficient'])*prod(hodge.integral(gg,tuple(p),tuple(l))
                        for gg,p,l in zip(data[0],powers,lambdas))
        expected=QQ(summary['genera'][g-2]['top_twisted_integrals'][str(D)])
        assert actual==expected
        checks.append({'genus':g,'degree':D,'top_twisted_integral':str(actual)})
        print('Direct Hodge check',g,D,actual,flush=True)
b=GrassmannianCompleteIntersectionCycles(t,2);b.hodge=hodge
assert b.top_integral((1,))==QQ(summary['genera'][0]['top_twisted_integrals']['1'])

genus_zero=[]
for D in (1,2):
    a=toric_genus_zero_number(t,(D,),4);b=toric_genus_zero_number(t,(D,),5);assert a==b
    genus_zero.append({'degree':D,'ordinary_gw':str(a),'multiple_cover_subtracted':str(a if D==1 else a-grassmannian_line_count(t)/8)})
dump(OUT/'genus-zero.json',genus_zero)
ordinary=[]
for g in (2,3):
    for D in (0,1,2):
        value=grassmannian_known_ordinary_gw(t,g,(D,))
        Z=QQ(summary['genera'][g-2]['top_twisted_integrals'][str(D)]);sign=t.comparison_sign(g,(D,))
        ordinary.append({'genus':g,'degree':D,'ordinary_gw':None if value is None else str(value),
                         'twisted_top':str(Z),'comparison_sign':str(sign),
                         'infinity_aggregate':None if value is None else str(sign*value-Z)})
dump(OUT/'ordinary-gw.json',ordinary)
dump(OUT/'verification.json',{'certificates_replayed':certificates,'previous_span_certificates':sum(d['rank'] for d in degrees if d['additional_rank']==0),
                             'negative_virtual_dimension_slots':vanishings,'zero_dimensional_lower_slots':zero_dimensional_slots,
                             'direct_hodge_checks':checks,'genus_zero_two_ray_checks':2})
print('VERIFIED AND EXPORTED Gr(2,7)[1^7]',flush=True)
