"""Reuse audited seeds, compute new ambient seeds, expand all labeled profiles."""
from log_glsm.cycles.marked import ambient_polarization, marked_multiplier
from log_glsm.targets.projective import ProjectiveCompleteIntersection
from log_glsm.targets.convex_gkm import convex_gkm_ordinary_seed
from log_glsm.cycles.higher_marked import higher_marked_profiles
from log_glsm.core.hodge_grr import GenusFourHodgeIntegrals
import argparse,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from sage.all import *
from calculations import marked_frontier as run
OLD=run.Path('results/higher-marked/2026-11-06')


def scan(index,g,max_n):
    t=run.targets()[index];path=run.OUT/'numbers'/(t.slug+'-g%s-n%s.json.gz'%(g,max_n))
    if path.exists():return str(path),'reused'
    start=time.time();old=run.read(OLD/(t.slug+'-g%s-numbers.json'%g)) if index<43 else None
    rows=[];h=GenusFourHodgeIntegrals()
    backends=[run.common.higher_marked_backend(t,g,base,h) for base in (4,5)] if old is None else None
    for beta in run.common.genus_four_degrees(t):
        if old is not None:
            seed=next(row for row in old['rows'] if tuple(row['beta'])==beta)
            z=QQ(seed['twisted_unmarked']);ordinary=None if seed['ordinary_unmarked'] is None else QQ(seed['ordinary_unmarked'])
        else:
            values=[b.numerical_localization(beta) for b in backends]
            assert all(not v['poles'] for v in values) and values[0]['finite']==values[1]['finite'],(t.slug,g,beta,values)
            z=values[0]['finite'];ordinary=convex_gkm_ordinary_seed(t,g,beta)
        d=beta[0] if isinstance(t,ProjectiveCompleteIntersection) else t.ambient.pairing(ambient_polarization(t),beta)
        degree=beta[0] if isinstance(t,ProjectiveCompleteIntersection) else beta
        alpha=t.alpha(g,degree);sign=t.comparison_sign(g,degree);profiles=[]
        for p in higher_marked_profiles(max_markings=max_n):
            m=marked_multiplier(g,d,p);Z=m*z;N=QQ.zero() if not m else None if ordinary is None else m*ordinary
            power=alpha+sum(x!='1' for x in p)-(3*g-3+len(p))
            B=-Z if power<0 else sign*N-Z if power==0 and N is not None else None
            profiles.append({'insertions':p,'multiplier':str(m),'twisted':str(Z),'ordinary':None if N is None else str(N),
                             'top_t_power':int(power),'infinity_aggregate':None if B is None else str(B)})
        rows.append({'beta':list(map(int,beta)),'H_degree':int(d),'alpha':int(alpha),'sign':int(sign),
                     'twisted_unmarked':str(z),'ordinary_unmarked':None if ordinary is None else str(ordinary),'profiles':profiles})
        print('NUMBER',t.slug,g,beta,z,ordinary,flush=True)
    run.dump(path,{'target':t.record(),'genus':g,'max_markings':max_n,'rows':rows,'seconds':time.time()-start,
                   'seed_source':'previous two-ray audited artifact' if old else 'new two-ray direct localization'})
    return str(path),time.time()-start


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-markings',type=int,default=5);p.add_argument('--workers',type=int,default=2);a=p.parse_args()
    with ProcessPoolExecutor(max_workers=a.workers,max_tasks_per_child=4) as pool:
        for f in as_completed([pool.submit(scan,i,g,a.max_markings) for i in range(54) for g in (2,3)]):print('DONE',f.result(),flush=True)
