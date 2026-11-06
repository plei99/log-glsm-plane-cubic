"""Marked Chow search in codimensions 0..3, separate from top numbers."""
import argparse, gzip, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
from sage.all import *

from log_glsm.cycles.marked import (
    ambient_polarization,
    genus_four_backend,
)
from log_glsm.targets.catalogue import (
    genus_four_degrees,
    genus_four_targets,
)
OUT=Path('results/genus-four-marked/2026-11-06')
PROFILES=[()] + [(x,) for x in ('1','psi','H')] + list(product(('1','psi','H'),repeat=2))
CONVERTER=None;VECTORS={}


def sparse(v):return {str(i):str(x) for i,x in enumerate(v) if x}


def vector_of(terms,n,r):
    global CONVERTER
    if CONVERTER is None:CONVERTER=genus_four_backend(genus_four_targets()[0])
    length=len(CONVERTER.adm.TautologicalRing(4,n).generators(r));result=vector(QQ,length)
    for key,co in terms.items():
        cache=(r,key)
        if cache not in VECTORS:VECTORS[cache]=vector(QQ,CONVERTER.marked_monomial_class(key).vector(r))
        result+=co*VECTORS[cache]
    return result


def scan(index,max_r):
    target=genus_four_targets()[index];start=time.time();states={};finite=[];reference={}
    for n in range(3):
        for r in range(max_r+1):states[n,r]={'echelon':{},'rows':[],'origins':[],'candidates':0}
    def consider(v,n,r,origin,terms):
        state=states[n,r];state['candidates']+=1;remainder=copy(v)
        for pivot,row in sorted(state['echelon'].items()):
            if remainder[pivot]:remainder-=remainder[pivot]*row
        if not remainder:return
        pivot=next(i for i,x in enumerate(remainder) if x)
        state['echelon'][pivot]=remainder/remainder[pivot];state['rows'].append(sparse(v))
        state['origins'].append(dict(origin,graph_terms=None if terms is None else [{'key':key,'coefficient':str(co)} for key,co in terms.items()]))
    for base in (4,5):
        b=genus_four_backend(target,base)
        for beta in genus_four_degrees(target):
            for profile in PROFILES:
                n=len(profile);coefficients=b.marked_coefficients(beta,profile,max_r)
                for r,powers in coefficients.items():
                    for power,terms in powers.items():
                        if power<0:
                            consider(vector_of(terms,n,r),n,r,{'kind':'auxiliary_pole','beta':list(map(int,beta)),
                                'insertions':profile,'base':base,'epsilon_power':power},terms)
                    terms=powers.get(0,{});v=vector_of(terms,n,r)
                    if base==4:
                        reference[beta,profile,r]=v
                        if v:finite.append({'beta':list(map(int,beta)),'insertions':profile,'codimension':r,'vector':sparse(v)})
                    else:consider(v-reference[beta,profile,r],n,r,{'kind':'finite_part_difference','beta':list(map(int,beta)),
                        'insertions':profile,'base':5,'reference_base':4},None)
            print(target.slug,'classes base',base,'beta',beta,'ranks',[(n,r,len(s['rows'])) for (n,r),s in states.items() if s['rows']],flush=True)
    result={'target':target.record(),'genus':4,'markings':[0,1,2],'max_codimension':max_r,'curve_bound':2,
        'profiles':PROFILES,'rays':[4,5],'polarization':list(map(int,ambient_polarization(target))),
        'finite_classes':finite,'relations':[{'markings':n,'codimension':r,'candidates':s['candidates'],
            'rows':s['rows'],'origins':s['origins']} for (n,r),s in states.items()],
        'seconds':time.time()-start,'pixton_membership':'awaiting separate exact matrix certificate'}
    path=OUT/(target.slug+'-classes-r'+str(max_r)+'.json.gz')
    with gzip.open(path,'wt') as f:json.dump(result,f,default=lambda x:int(x) if x in ZZ else str(x))
    return {'slug':target.slug,'seconds':result['seconds'],'ranks':[(n,r,len(s['rows'])) for (n,r),s in states.items() if s['rows']]}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--max-r',type=int,default=3)
    p.add_argument('--indices',type=int,nargs='*');p.add_argument('--resume',action='store_true');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    indices=range(43) if args.indices is None else args.indices
    indices=[i for i in indices if not args.resume or not (OUT/(genus_four_targets()[i].slug+'-classes-r'+str(args.max_r)+'.json.gz')).exists()]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(scan,i,args.max_r) for i in indices]):print('FINISHED',future.result(),flush=True)
