"""All ordered profiles through three markings in genera two and three."""
import argparse,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from sage.all import *

from log_glsm.core.hodge_grr import (
    GenusFourHodgeIntegrals,
)
from log_glsm.cycles.higher_marked import (
    higher_marked_backend,
    higher_marked_ordinary_seed,
    higher_marked_profiles,
)
from log_glsm.cycles.marked import (
    ambient_polarization,
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_degrees,
    genus_four_targets,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
)
OUT=Path('results/higher-marked/2026-11-06')


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,default=lambda x:int(x) if x in ZZ else str(x))+'\n')


def scan(index,g):
    target=genus_four_targets()[index];start=time.time();h=GenusFourHodgeIntegrals()
    backends=[higher_marked_backend(target,g,base,h) for base in (4,5)];rows=[]
    for beta in genus_four_degrees(target):
        values=[b.numerical_localization(beta) for b in backends]
        assert all(not v['poles'] for v in values),(target.slug,g,beta,values)
        assert values[0]['finite']==values[1]['finite'],(target.slug,g,beta)
        z=values[0]['finite'];ordinary=higher_marked_ordinary_seed(target,g,beta)
        dH=beta[0] if isinstance(target,ProjectiveCompleteIntersection) else target.ambient.pairing(ambient_polarization(target),beta)
        degree=beta[0] if isinstance(target,ProjectiveCompleteIntersection) else beta
        alpha=target.alpha(g,degree);sign=target.comparison_sign(g,degree);profiles=[]
        for profile in higher_marked_profiles():
            multiplier=marked_multiplier(g,dH,profile);a=sum(x!='1' for x in profile)
            power=alpha+a-(3*g-3+len(profile));value=multiplier*z
            compact=QQ.zero() if not multiplier else None if ordinary is None else multiplier*ordinary
            infinity=-value if power<0 else sign*compact-value if power==0 and compact is not None else None
            profiles.append({'insertions':profile,'multiplier':str(multiplier),'twisted':str(value),
                'ordinary':None if compact is None else str(compact),'top_t_power':int(power),
                'infinity_aggregate':None if infinity is None else str(infinity)})
        rows.append({'beta':list(map(int,beta)),'H_degree':int(dH),'alpha_0':int(alpha),'sign':str(sign),
            'twisted_unmarked':str(z),'ordinary_unmarked':None if ordinary is None else str(ordinary),
            'graphs_unmarked':values[0]['graphs'],'profiles':profiles})
        print(target.slug,'g',g,'beta',beta,'Z',z,flush=True)
    result={'target':target.record(),'genus':g,'markings':[0,1,2,3],'curve_bound':2,
        'polarization':list(map(int,ambient_polarization(target))),'rays':[4,5],'all_poles_cancel':True,
        'ray_finite_parts_agree':True,'rows':rows,'seconds':time.time()-start}
    dump(OUT/(target.slug+'-g'+str(g)+'-numbers.json'),result)
    return target.slug,g,result['seconds']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--resume',action='store_true');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    jobs=[(i,g) for i,t in enumerate(genus_four_targets()) for g in (2,3)
          if not args.resume or not (OUT/(t.slug+'-g'+str(g)+'-numbers.json')).exists()]
    with ProcessPoolExecutor(max_workers=args.workers,max_tasks_per_child=8) as pool:
        for future in as_completed([pool.submit(scan,*j) for j in jobs]):print('DONE',future.result(),flush=True)
    records=[json.loads((OUT/(t.slug+'-g'+str(g)+'-numbers.json')).read_text()) for t in genus_four_targets() for g in (2,3)]
    summary={'presentations':43,'genera':[2,3],'markings':[0,1,2,3],'curve_bound':2,
        'sectors':sum(len(d['rows']) for d in records),'ordered_profiles':sum(len(r['profiles']) for d in records for r in d['rows']),
        'unknown_ordinary_sectors':sum(r['ordinary_unmarked'] is None for d in records for r in d['rows'])}
    dump(OUT/'numerical-summary.json',summary);print(summary,flush=True)
