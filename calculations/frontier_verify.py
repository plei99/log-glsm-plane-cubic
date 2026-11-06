"""Replay every stored basis row, exact Pixton certificates, and scope metadata."""
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_certify import suffix,OLD
total=0;basis=0;spaces=[]
for g in (2,3):
    for n in range(6):
        for r in range(3):
            tag=suffix(g,n,r);M=load(str(run.OUT/('pixton-'+tag+'.sobj')))
            V=load(str(run.OUT/('relations-'+tag+'.sobj')));U=load(str(run.OUT/('certificate-'+tag+'.sobj')))
            assert U*M==V and V.rank()==V.nrows();basis+=V.nrows();pivots=V.pivots();count=0
            for i,t in enumerate(run.targets()):
                path=OLD/(t.slug+'-g%s-classes.json.gz'%g) if n<=3 and i<43 else run.result_path(i,g,n,2)
                d=run.read(path)
                if n>=4:assert d['profiles']==[list(p) for p in run.profiles(n)]
                item=next(s for s in d['relations'] if s['codimension']==r and s.get('markings',n)==n)
                rows=[vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in v.items()}) for v in item['rows']]
                if rows:
                    A=matrix(QQ,rows);assert A.matrix_from_columns(pivots)*V==A,(t.slug,g,n,r)
                count+=len(rows)
                if g==3 and n==5:
                    assert not item['rows']
                    assert not any('1' not in x['profile'] for x in d['finite_canonical_classes'])
            total+=count;spaces.append({'g':g,'n':n,'r':r,'rank':V.nrows(),'per_target_rows':count})
            print('VERIFIED',tag,V.nrows(),count,flush=True)
probe=run.read(run.OUT/'codimension-three-probe.json');M=load(str(run.OUT/'pixton-g3-n3-r3.sobj'))
V=matrix(QQ,[vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in row.items()}) for row in probe['V']])
U=matrix(QQ,[vector(QQ,M.nrows(),{int(j):QQ(c) for j,c in row.items()}) for row in probe['U']])
assert U*M==V and V.rank()==probe['found_rank']
run.dump(run.OUT/'verification.json',{'pixton_basis_rows_replayed':basis,'per_target_rows_replayed':total,
    'spaces':spaces,'higher_codimension_probe_rows':V.nrows(),'all_certificates_exact':True})
