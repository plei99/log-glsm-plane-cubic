"""Check the previous toric span and geometric virtual-dimension vanishings."""
import gzip,json
from pathlib import Path
from sage.all import *

OUT=Path('results/grassmannian-complete-intersections/2026-11-06')
TORIC=Path('results/toric-complete-intersections/2026-11-06')


def read(p):
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())


degrees=[]
for g in (2,3):
    for r in range(3*g-2):
        name='strata-g%s-r%s.json'%(g,r)
        strata=read(OUT/name);assert strata==read(TORIC/name)
        name='union-relations-g%s-r%s.json.gz'%(g,r)
        def mat(folder):
            records=read(folder/name)
            return matrix(QQ,[v['strata_vector'] for v in records]) if records else matrix(QQ,0,len(strata))
        A=mat(TORIC);B=mat(OUT);combined=A.stack(B).rank()
        equal=(A.row_space()==B.row_space())
        if combined==A.rank() and B.nrows():
            W=A.transpose().solve_right(B.transpose()).transpose();assert W*A==B
            data={'genus':g,'codimension':r,'weights':[[[int(i),str(c)] for i,c in enumerate(w) if c] for w in W.rows()]}
            (OUT/('toric-span-g%s-r%s.json.gz'%(g,r))).write_bytes(gzip.compress(json.dumps(data).encode(),mtime=0))
        degrees.append({'genus':g,'codimension':r,'toric_rank':A.rank(),'grassmannian_rank':B.rank(),
                        'combined_rank':combined,'equal_spans':equal})

vanishing=[]
for s in read(OUT/'catalog.json'):
    target=s['target'];n=target['ambient']['n'];N=target['ambient']['dimension']
    for g in (2,3):
        for rec in read(OUT/target['id']/('finite-classes-g%s.json.gz'%g)):
            D=rec['beta'][0];vdim=(1-g)*(N-3)+n*D
            if vdim<0:
                v=vector(QQ,rec['strata_vector'])
                if v:
                    M=load(str(OUT/('pixton-g%s-r%s.sobj'%(g,rec['codimension']))))
                    assert v in M.row_space()
                vanishing.append({'target':target['id'],'genus':g,'degree':D,
                                  'codimension':rec['codimension'],'ambient_virtual_dimension':vdim})

result={'degrees':degrees,'equal_to_previous_toric_span':all(d['equal_spans'] for d in degrees),
        'virtual_dimension_vanishing_checks':len(vanishing),'vanishing_slots':vanishing}
(OUT/'geometric-verification.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
print('Equal to previous toric span:',result['equal_to_previous_toric_span'])
print('Verified negative ambient virtual dimension:',len(vanishing),'slots')
