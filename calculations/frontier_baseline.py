"""Add the eleven new presentations to the already certified n<=3 spans."""
import shutil
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_certify import exact_span,suffix,OLD
summary=[]
for g in (2,3):
    for n in range(4):
        data=[run.read(run.result_path(i,g,n,2)) for i in range(43,54)]
        for r in range(3):
            tag=suffix(g,n,r);M=load(str(OLD/('pixton-'+tag+'.sobj')))
            prior=load(str(OLD/('relations-'+tag+'.sobj')));rows=list(prior.rows())
            for d in data:
                item=next(s for s in d['relations'] if s['codimension']==r)
                rows.extend(vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in v.items()}) for v in item['rows'])
            V=exact_span(matrix(QQ,rows)) if rows else matrix(QQ,0,M.ncols())
            U=M.transpose().solve_right(V.transpose()).transpose() if V.nrows() else matrix(QQ,0,M.nrows())
            assert U*M==V
            save(V,str(run.OUT/('relations-'+tag+'.sobj')));save(U,str(run.OUT/('certificate-'+tag+'.sobj')))
            shutil.copy2(OLD/('pixton-'+tag+'.sobj'),run.OUT/('pixton-'+tag+'.sobj'))
            shutil.copy2(OLD/('strata-'+tag+'.json'),run.OUT/('strata-'+tag+'.json'))
            run.dump(run.OUT/('certificate-'+tag+'.json'),{'equation':'U M = V','complete':True,
                'V':[run.common.sparse(v) for v in V.rows()],'U':[run.common.sparse(v) for v in U.rows()]})
            item={'g':g,'n':n,'r':r,'old_found':prior.nrows(),'found':V.nrows(),'pixton':M.rank(),
                  'strata':M.ncols(),'new_ambient_rows':len(rows)-prior.nrows(),'certificate_verified':True}
            summary.append(item);print(item,flush=True)
run.dump(run.OUT/'baseline-extension.json',{'presentations':54,'spaces':summary})
