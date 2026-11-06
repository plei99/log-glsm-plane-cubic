"""Exact genus-four numerical batch, retaining every ordered marking profile."""
import argparse, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
from sage.all import *

from log_glsm.core.hodge_grr import (
    GenusFourHodgeIntegrals,
)
from log_glsm.cycles.marked import (
    ambient_polarization,
    genus_four_backend,
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_degrees,
    genus_four_kind,
    genus_four_targets,
    known_genus_four_ordinary,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
)

OUT=Path('results/genus-four-marked/2026-11-06')
PROFILES=[()] + [(x,) for x in ('1','psi','H')] + list(product(('1','psi','H'),repeat=2))


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,default=lambda x:int(x) if x in ZZ else str(x))+'\n')


def scan(index):
    target=genus_four_targets()[index];start=time.time()
    hodge=GenusFourHodgeIntegrals();backends=[genus_four_backend(target,base,hodge) for base in (4,5)]
    p=ambient_polarization(target);rows=[]
    for beta in genus_four_degrees(target):
        reference=None
        for base,b in zip((4,5),backends):
            value=b.numerical_localization(beta)
            assert not value['poles'],(target.slug,beta,base,value)
            if reference is None:reference=value
            else:assert value['finite']==reference['finite'],(target.slug,beta,'ray discrepancy')
        z=reference['finite'];ordinary=known_genus_four_ordinary(target,beta)
        dH=beta[0] if isinstance(target,ProjectiveCompleteIntersection) else target.ambient.pairing(p,beta)
        degree=beta[0] if isinstance(target,ProjectiveCompleteIntersection) else beta
        alpha=target.alpha(4,degree);sign=target.comparison_sign(4,degree)
        records=[]
        for profile in PROFILES:
            multiplier=marked_multiplier(4,dH,profile);inserted=sum(x!='1' for x in profile)
            power=alpha+inserted-9-len(profile);value=multiplier*z
            compact=QQ.zero() if not multiplier else None if ordinary is None else multiplier*ordinary
            infinity=(-value if power<0 else sign*compact-value if power==0 and compact is not None else None)
            records.append({'insertions':profile,'multiplier':str(multiplier),'twisted':str(value),
                'ordinary':None if compact is None else str(compact),'top_t_power':int(power),
                'infinity_aggregate':None if infinity is None else str(infinity)})
        rows.append({'beta':list(map(int,beta)),'H_degree':int(dH),'alpha_0':int(alpha),'comparison_sign':str(sign),
            'graphs_unmarked':reference['graphs'],'twisted_unmarked':str(z),
            'ordinary_unmarked':None if ordinary is None else str(ordinary),'profiles':records})
        print(target.slug,tuple(beta),'Z',z,'ordinary',ordinary,flush=True)
    result={'target':target.record(),'kind':genus_four_kind(target),'genus':4,'markings':[0,1,2],
        'polarization':list(map(int,p)),'curve_bound':2,'rays':[4,5],
        'method':'exact fixed-locus integration; string/dilaton/divisor for ordered profiles',
        'scalar_checks':{'all_negative_auxiliary_poles_zero':True,'finite_parts_equal_on_two_rays':True},
        'pixton_membership':'not inferred from scalar checks','rows':rows,'seconds':time.time()-start}
    dump(OUT/(target.slug+'-numbers.json'),result)
    return {'slug':target.slug,'sectors':len(rows),'profiles':len(rows)*len(PROFILES),'seconds':result['seconds']}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True);start=time.time()
    indices=[i for i,t in enumerate(genus_four_targets()) if not args.resume or not (OUT/(t.slug+'-numbers.json')).exists()]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(scan,i) for i in indices]):print('FINISHED',future.result(),flush=True)
    results=[json.loads((OUT/(t.slug+'-numbers.json')).read_text()) for t in genus_four_targets()]
    summary={'presentations':len(results),'genus':4,'markings':[0,1,2],'curve_bound':2,
        'sectors':sum(len(r['rows']) for r in results),'ordered_profile_values':sum(len(row['profiles']) for r in results for row in r['rows']),
        'unknown_ordinary_unmarked':sum(row['ordinary_unmarked'] is None for r in results for row in r['rows']),
        'rays':[4,5],'seconds_this_run':time.time()-start,'files':[r['target'].get('slug',r['target'].get('id'))+'-numbers.json' for r in results]}
    dump(OUT/'numerical-summary.json',summary);print(summary,flush=True)
