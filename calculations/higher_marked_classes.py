"""Genus 2/3, n<=3, degree<=2, codimension<=2 marked Chow relations."""
import argparse,gzip,json,time,gc,fcntl
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from sage.all import *
from itertools import permutations
from log_glsm.cycles.higher_marked import (
    higher_marked_backend,
    higher_marked_profiles,
    marked_profile_permutations,
)
from log_glsm.cycles.marked import (
    GenusFourMixin,
)
from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.catalogue import (
    genus_four_degrees,
    genus_four_targets,
)
from log_glsm.targets.projective import (
    CompleteIntersectionCycles,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersectionCycles,
)
OUT=Path('results/higher-marked/2026-11-06')
CONVERTERS={};VECTORS={};ACTIONS={}


def sparse(v):return {str(i):str(x) for i,x in enumerate(v) if x}


def converter(g):
    if g not in CONVERTERS:CONVERTERS[g]=higher_marked_backend(genus_four_targets()[0],g)
    return CONVERTERS[g]


def vector_of(terms,g,n,r):
    b=converter(g);length=len(b.adm.TautologicalRing(g,n).generators(r));result=vector(QQ,length)
    for key,co in terms.items():
        cache=(g,r,key)
        if cache not in VECTORS:VECTORS[cache]=vector(QQ,b.marked_monomial_class(key).vector(r))
        result+=co*VECTORS[cache]
    return result


def permute_vector(v,g,n,r,mapping):
    if mapping==tuple(range(1,n+1)):return v
    key=(g,n,r,mapping)
    if key not in ACTIONS:
        R=converter(g).adm.TautologicalRing(g,n);action=[]
        for term in R.generators(r):
            renamed=vector(QQ,term.rename_legs({i+1:j for i,j in enumerate(mapping)},inplace=False).vector(r))
            support=[j for j,c in enumerate(renamed) if c]
            assert len(support)==1 and renamed[support[0]]==1
            action.append(support[0])
        ACTIONS[key]=action
    result=vector(QQ,len(v))
    for i,j in enumerate(ACTIONS[key]):result[j]+=v[i]
    return result


def scan(index,g,max_r=2,profile_subset=None,part=None):
    target=genus_four_targets()[index];start=time.time();states={};reference={};finite=[]
    path=(OUT/(target.slug+'-g'+str(g)+'-classes.json.gz') if part is None else
          OUT/'parts'/(target.slug+'-g'+str(g))/('part%s.json.gz'%part))
    path.parent.mkdir(parents=True,exist_ok=True)
    lock_handle=open(path.with_suffix('.lock'),'w')
    fcntl.flock(lock_handle,fcntl.LOCK_EX)
    # A heavy target can be scheduled independently while a resumed pool is
    # still working through earlier targets. Completed atomic checkpoints win.
    if path.exists():
        with gzip.open(path,'rt') as f:existing=json.load(f)
        if existing['max_codimension']==max_r:return target.slug,g,'checkpoint already complete'
    for n in range(4):
        for r in range(max_r+1):states[n,r]={'echelon':{},'rows':[],'origins':[],'candidates':0}
    def consider(v,n,r,origin):
        state=states[n,r];state['candidates']+=1;remainder=copy(v)
        for pivot,row in sorted(state['echelon'].items()):
            if remainder[pivot]:remainder-=remainder[pivot]*row
        if not remainder:return
        pivot=next(i for i,x in enumerate(remainder) if x)
        state['echelon'][pivot]=remainder/remainder[pivot];state['rows'].append(sparse(v));state['origins'].append(origin)
    def consider_orbit(v,n,r,origin):
        if not v:return
        for mapping in permutations(range(1,n+1)):
            consider(permute_vector(v,g,n,r,mapping),n,r,dict(origin,marking_permutation=mapping))
    for base in (4,5):
        b=higher_marked_backend(target,g,base)
        for beta in genus_four_degrees(target):
            for profile in (higher_marked_profiles(canonical=True) if profile_subset is None else profile_subset):
                n=len(profile);coefficients=b.marked_coefficients(beta,profile,max_r)
                for r,powers in coefficients.items():
                    for power,terms in powers.items():
                        if power<0:
                            consider_orbit(vector_of(terms,g,n,r),n,r,{'kind':'auxiliary_pole','beta':list(map(int,beta)),
                                'canonical_profile':profile,'base':base,'epsilon_power':power})
                    v=vector_of(powers.get(0,{}),g,n,r)
                    if base==4:
                        reference[beta,profile,r]=v
                        for ordered,mapping in marked_profile_permutations(profile):
                            w=permute_vector(v,g,n,r,mapping)
                            if w:finite.append({'beta':list(map(int,beta)),'insertions':ordered,'codimension':r,'vector':sparse(w)})
                        # A permutation fixing the profile represents the same
                        # class, though the raw strata representatives may differ.
                        for mapping in permutations(range(1,n+1)):
                            if all(profile[old-1]==profile[new-1] for old,new in enumerate(mapping,1)):
                                consider_orbit(permute_vector(v,g,n,r,mapping)-v,n,r,
                                    {'kind':'equal_profile_symmetry','beta':list(map(int,beta)),'canonical_profile':profile,
                                     'base':base,'stabilizer_permutation':mapping})
                    else:consider_orbit(v-reference[beta,profile,r],n,r,
                        {'kind':'finite_part_difference','beta':list(map(int,beta)),'canonical_profile':profile,'base':5,'reference_base':4})
            print(target.slug,'g',g,'base',base,'beta',beta,'ranks',[(n,r,len(s['rows'])) for (n,r),s in states.items() if s['rows']],flush=True)
    result={'target':target.record(),'genus':g,'markings':[0,1,2,3],'curve_bound':2,'max_codimension':max_r,
        'rays':[4,5],'profile_method':'20 canonical profiles; exact labeled permutation action; all 40 ordered profiles recovered',
        'profile_subset':profile_subset,'split_part':part,
        'finite_classes':finite,'relations':[{'markings':n,'codimension':r,'candidates':s['candidates'],
            'rows':s['rows'],'origins':s['origins']} for (n,r),s in states.items()],
        'seconds':time.time()-start,'pixton_membership':'requires separate certificate replay'}
    temporary=path.with_suffix(path.suffix+'.tmp')
    with gzip.open(temporary,'wt') as f:json.dump(result,f,default=lambda x:int(x) if x in ZZ else str(x))
    temporary.replace(path)
    # These caches retain target-specific rational functions, not the reusable
    # strata vectors or the marking actions. Release them between targets.
    for cls in (GenusFourMixin,CompleteIntersectionCycles,ToricCompleteIntersectionCycles,QuarticK3GraphCycles):
        for name in ('marked_vertex_terms','_marked_hodge_vertex','_hodge_vertex','_general_edge','integrated_vertex'):
            method=cls.__dict__.get(name)
            if hasattr(method,'cache_clear'):method.cache_clear()
    gc.collect()
    return target.slug,g,result['seconds']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--indices',type=int,nargs='*')
    p.add_argument('--genera',type=int,nargs='*',default=[2,3]);p.add_argument('--resume',action='store_true');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);indices=range(43) if args.indices is None else args.indices
    jobs=[(i,g) for i in indices for g in args.genera if not args.resume or not (OUT/(genus_four_targets()[i].slug+'-g'+str(g)+'-classes.json.gz')).exists()]
    with ProcessPoolExecutor(max_workers=args.workers,max_tasks_per_child=6) as pool:
        for future in as_completed([pool.submit(scan,*j) for j in jobs]):print('DONE',future.result(),flush=True)
