"""Exact bounded Grassmannian search; shared stabilization and explicit Pixton certificates."""
import argparse,gc,gzip,json,multiprocessing,shutil,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from types import SimpleNamespace
from sage.all import *
from calculations import projective_ci_scan as common
from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersection,
    GrassmannianCompleteIntersectionCycles,
    grassmannian_sample,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersectionCycles,
    toric_degree_tuples,
)
OUT=Path('results/grassmannian-complete-intersections/2026-11-06')
dump=common.dump


def reset_caches():
    for cls in (GrassmannianCompleteIntersectionCycles,ToricCompleteIntersectionCycles,QuarticK3GraphCycles,QuarticK3Cycles):
        for method in cls.__dict__.values():
            if hasattr(method,'cache_clear'):method.cache_clear()
    common.reset_caches();gc.collect()


def scan_target(spec):
    target=GrassmannianCompleteIntersection(spec['ambient']['n'],spec['degrees'],spec['id'])
    folder=OUT/target.slug;folder.mkdir(exist_ok=True)
    betas=toric_degree_tuples(target.ambient);start=time.time()
    summary={'target':target.record(),'mori_bound':2,'weight_bases':[4,5],'genera':[]}
    for g in (2,3):
        states={r:{'rows':[],'origins':[],'echelon':{},'candidates':0} for r in range(3*g-2)}
        reference={};finite=[];counts={}
        def as_vector(terms,r):
            v=common.vector_of(g,terms)
            return v if v is not None else vector(QQ,common.MODELS[g]['degrees'][r]['matrix'].ncols())
        def consider(v,r,origin):
            state=states[r];state['candidates']+=1;remainder=copy(v)
            for pivot,row in sorted(state['echelon'].items()):
                if remainder[pivot]:remainder-=remainder[pivot]*row
            if not remainder:return
            pivot=next(i for i,x in enumerate(remainder) if x)
            state['echelon'][pivot]=remainder/remainder[pivot]
            state['rows'].append(v);state['origins'].append(origin)
        for base in (4,5):
            b=GrassmannianCompleteIntersectionCycles(target,g,base)
            for beta in betas:
                name=','.join(map(str,beta));counts[name]=len(b.graphs(beta))
                for r in states:
                    coeffs=b.coefficients(beta,r)
                    for pole,terms in coeffs.items():
                        if pole<0:consider(as_vector(terms,r),r,{'kind':'auxiliary_pole','beta':list(beta),'epsilon_power':pole,'base':base})
                    terms=coeffs.get(0,{});v=as_vector(terms,r)
                    if base==4:
                        reference[beta,r]=v;power=target.t_power(g,beta,r)
                        compact=target.compact_degree_zero_terms(g) if not any(beta) and power==0 else None
                        known=power<0 or (power==0 and (target.kind=='K3' or not any(beta)))
                        finite.append({'beta':list(beta),'codimension':r,'t_power':int(power),
                            'kind':'negative_pole' if power<0 else 'compact_comparison' if power==0 else 'positive_power',
                            'compact_known':known,'comparison_sign':str(target.comparison_sign(g,beta)),
                            'compact_terms':None if compact is None else [{'key':k,'coefficient':str(c)} for k,c in compact.items()],
                            'terms':[{'key':k,'coefficient':str(c)} for k,c in terms.items()],
                            'strata_vector':[str(c) for c in v]})
                    else:consider(v-reference[beta,r],r,{'kind':'finite_part_difference','beta':list(beta),'base':base,'reference_base':4})
                print(target.slug,'g',g,'beta',beta,'base',base,'graphs',counts[name],
                      'ranks',[len(states[r]['rows']) for r in states],flush=True)
        model=common.MODELS[g];zero=tuple([0]*target.ambient.picard_rank)
        for r in states:
            expected={}
            if target.ambient.dimension==3 and r==3*g-3:
                X=target.ambient;c=[X.homogeneous(X.total_chern,k) for k in range(4)]
                co=X.integral(c[3]-c[1]*c[2])*(-1)**g
                expected={('interior',None,(1,2) if g==2 else (1,2,3)):co}
            assert reference[zero,r]-as_vector(expected,r) in model['degrees'][r]['matrix'].row_space(),('ambient constants',target.slug,g,r)
        top={','.join(map(str,beta)):str(reference[beta,3*g-3]*model['integral']) for beta in betas}
        gs={'genus':g,'top_twisted_integrals':top,'fixed_graph_counts':counts,'degrees':[],
            'negative_t_slots':sum(x['t_power']<0 for x in finite),
            'k3_zero_t_slots':sum(x['t_power']==0 and target.kind=='K3' for x in finite)}
        dump(folder/('finite-classes-g%s.json'%g),finite,compressed=True)
        for r,state in states.items():
            M=model['degrees'][r]['matrix'];records=[]
            if state['rows']:
                V=matrix(QQ,state['rows']);U=M.transpose().solve_right(V.transpose()).transpose()
                assert U*M==V and V.rank()==len(state['rows']),('Pixton membership',target.slug,g,r)
                for v,u,origin in zip(V.rows(),U.rows(),state['origins']):
                    records.append({'origin':origin,'strata_vector':[str(c) for c in v],
                                    'pixton_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
            dump(folder/('relations-g%s-r%s.json'%(g,r)),records,compressed=True)
            gs['degrees'].append({'codimension':r,'relation_rank':len(records),'candidates':state['candidates'],
                                 'pixton_rank':int(M.rank()),'strata_columns':M.ncols()})
        summary['genera'].append(gs);summary['seconds']=time.time()-start;dump(folder/'summary.json',summary)
        print('CERTIFIED',target.slug,'g',g,'ranks',[d['relation_rank'] for d in gs['degrees']],'top',top,flush=True)
        reset_caches()
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--targets',nargs='*');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    common.OUT=OUT;common.ARGS=SimpleNamespace(max_degree=2)
    old=Path('results/complete-intersections/2026-11-06')
    for name in ('stabilization-vectors.sobj','stabilization-cache.json'):
        if (old/name).exists() and not (OUT/name).exists():shutil.copy2(old/name,OUT/name)
    common.prepare()
    targets=grassmannian_sample()
    if args.targets:targets=tuple(t for t in targets if t.slug in args.targets)
    if not targets:raise ValueError('no matching targets')
    catalog=[]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('fork')) as pool:
        futures=[pool.submit(scan_target,t.record()) for t in targets]
        for f in as_completed(futures):
            catalog.append(f.result());dump(OUT/'catalog.json',sorted(catalog,key=lambda s:int(s['target']['id'][1:])))
    print('COMPLETE',len(catalog),'targets',flush=True)
