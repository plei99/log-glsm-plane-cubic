"""Exact Fano/CY complete-intersection scan, with shared stabilization vectors."""
import argparse,gc,gzip,hashlib,json,multiprocessing,shutil,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from log_glsm.provenance import source_fingerprint
from sage.all import *
from math import factorial
from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.projective import (
    CompleteIntersectionCycles,
    ProjectiveCompleteIntersection,
    complete_intersection_weights,
    requested_complete_intersections,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
)

OUT=Path('results/complete-intersections/2026-11-06')
VECTORS={}
MODELS={}
ARGS=None


def dump(path,data,compressed=False):
    raw=(json.dumps(data,separators=(',',':'),default=int)+'\n').encode()
    if compressed:
        path=Path(str(path)+'.gz');raw=gzip.compress(raw,mtime=0)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=Path(str(path)+'.tmp');temporary.write_bytes(raw);temporary.replace(path)


def reset_caches():
    for cls in (CompleteIntersectionCycles,QuarticK3GraphCycles,QuarticK3Cycles):
        for method in cls.__dict__.values():
            if hasattr(method,'cache_clear'):method.cache_clear()
    gc.collect()


def prepare():
    _load_vendored_admcycles()
    from admcycles import DR
    from admcycles.moduli import get_moduli
    for g in (2,3):
        b=QuarticK3GraphCycles(g,complete_intersection_weights(4))
        model={'backend':b,'degrees':{}}
        for r in range(3*g-2):
            source=Path('results/quartic-k3/2026-11-05')/('pixton-g%s-r%s.sobj'%(g,r))
            if source.exists():
                M=load(str(source))
                old=json.loads(source.with_name('strata-g%s-r%s.json'%(g,r)).read_text())
                assert old==[str(s) for s in b.R.generators(r)]
            else:
                raw=DR.rels_matrix(g,r,0,symm=0,moduli_type=get_moduli('st',DRpy=True),usespin=True)
                M=matrix(QQ,[DR.convert_vector_to_monomial_basis(row,g,r,(),
                    moduli_type=get_moduli('st',DRpy=True)) for row in raw.rows()])
            save(M,str(OUT/('pixton-g%s-r%s.sobj'%(g,r))))
            dump(OUT/('strata-g%s-r%s.json'%(g,r)),[str(s) for s in b.R.generators(r)])
            model['degrees'][r]={'matrix':M,'echelon':M.row_space().basis_matrix()}
        r=3*g-3;M=model['degrees'][r]['matrix']
        assert M.rank()==M.ncols()-1, 'top tautological quotient must be one-dimensional'
        kappa=vector(QQ,b.R.kappa(r).vector(r))
        integral=M.stack(matrix(QQ,[kappa])).solve_right(
            vector(QQ,[0]*M.nrows()+[QQ(1)/(24**g*factorial(g))]))
        assert M*integral==0 and kappa*integral==QQ(1)/(24**g*factorial(g))
        model['integral']=integral;MODELS[g]=model

    fingerprint=source_fingerprint()
    cache=OUT/'stabilization-vectors.sobj';meta=OUT/'stabilization-cache.json'
    if cache.exists() and meta.exists() and json.loads(meta.read_text())['fingerprint']==fingerprint:
        VECTORS.update(load(str(cache)));print('Loaded',len(VECTORS),'stabilization vectors',flush=True)
        return
    for g in (2,3):
        b=MODELS[g]['backend'];keys=set()
        for d in range(ARGS.max_degree+1):keys.update(b.rational_terms(d))
        for index,key in enumerate(sorted(keys,key=lambda k:(b.codimension(k),str(k)))):
            r=b.codimension(key)
            VECTORS[g,key]=vector(QQ,b.monomial_class(key).vector(r))
            if index%500==0:print('Stabilization cache genus',g,index,'/',len(keys),flush=True)
    save(VECTORS,str(cache));dump(meta,{'fingerprint':fingerprint,'vectors':len(VECTORS)})
    reset_caches()


def vector_of(g,terms):
    if not terms:return None
    r=QuarticK3GraphCycles.codimension(next(iter(terms)))
    model=MODELS[g];v=vector(QQ,model['degrees'][r]['matrix'].ncols())
    for key,co in terms.items():
        if (g,key) not in VECTORS:
            VECTORS[g,key]=vector(QQ,model['backend'].monomial_class(key).vector(r))
        v+=co*VECTORS[g,key]
    return v


def scan_target(spec):
    target=ProjectiveCompleteIntersection(spec['ambient_dimension'],spec['degrees'],spec['name'])
    folder=OUT/target.slug;folder.mkdir(exist_ok=True)
    started=time.time();summary={'target':target.record(),'max_degree':ARGS.max_degree,
                                'weight_bases':ARGS.bases,'genera':[]}
    for g in (2,3):
        states={r:{'rows':[],'origins':[],'echelon':{},'candidates':0,'growth':[]}
                for r in range(3*g-2)}
        reference={};finite=[];top={};counts={}
        def as_vector(terms,r):
            value=vector_of(g,terms)
            return value if value is not None else vector(QQ,MODELS[g]['degrees'][r]['matrix'].ncols())
        def consider(v,r,origin):
            state=states[r];state['candidates']+=1
            remainder=copy(v)
            for pivot,row in sorted(state['echelon'].items()):
                if remainder[pivot]:remainder-=remainder[pivot]*row
            if not remainder:return
            pivot=next(i for i,x in enumerate(remainder) if x)
            state['echelon'][pivot]=remainder/remainder[pivot]
            state['rows'].append(v);state['origins'].append(origin)

        for base_index,base in enumerate(ARGS.bases):
            weights=complete_intersection_weights(target.N+1,base)
            backend=CompleteIntersectionCycles(target,g,weights)
            for d in range(ARGS.max_degree+1):
                counts[d]=len(backend.graphs(d))
                for r,state in states.items():
                    coeffs=backend.coefficients(d,r)
                    for power,terms in coeffs.items():
                        if power<0:consider(as_vector(terms,r),r,
                            {'kind':'auxiliary_pole','ambient_degree':d,'epsilon_power':power,
                             'weight_base':base})
                    terms=coeffs.get(0,{})
                    v=as_vector(terms,r)
                    if base_index==0:
                        reference[d,r]=v
                        power=target.t_power(g,d,r)
                        compact=target.compact_degree_zero_terms(g) if d==0 and power==0 else None
                        finite.append({'ambient_degree':d,'codimension':r,'t_power':int(power),
                            'kind':('negative_pole' if power<0 else 'compact_comparison' if power==0 else 'positive_power'),
                            'compact_known':bool(power<0 or (power==0 and d==0)),'comparison_sign':str(target.comparison_sign(g,d)),
                            'compact_terms':None if compact is None else
                                [{'key':k,'coefficient':str(c)} for k,c in compact.items()],
                            'terms':[{'key':k,'coefficient':str(c)} for k,c in terms.items()],
                            'strata_vector':[str(c) for c in v]})
                    else:
                        consider(v-reference[d,r],r,{'kind':'finite_part_difference',
                            'ambient_degree':d,'weight_base':base,'reference_base':ARGS.bases[0]})
                print(target.slug,'g',g,'D',d,'base',base,'graphs',counts[d],
                      'ranks',[len(states[r]['rows']) for r in states],flush=True)
            for r,state in states.items():state['growth'].append(len(state['rows']))
        model=MODELS[g]
        for d in range(ARGS.max_degree+1):top[d]=str(reference[d,3*g-3]*model['integral'])
        # Independent ambient constant-map obstruction calculation.
        for r in states:
            expected={}
            if target.N==3 and r==3*g-3:
                expected={('interior',None,((1,2) if g==2 else (1,2,3))):(-1)**(g+1)*20}
            difference=reference[0,r]-as_vector(expected,r)
            assert difference in model['degrees'][r]['matrix'].row_space(),('constant maps',target.slug,g,r)

        gs={'genus':g,'top_twisted_integrals':top,'fixed_graph_counts':counts,
            'negative_cjr_identities':sum(rec['kind']=='negative_pole' for rec in finite),
            'constant_comparisons':sum(rec['kind']=='compact_comparison' for rec in finite),
            'degrees':[]}
        dump(folder/('finite-classes-g%s.json'%g),finite,compressed=True)
        for r,state in states.items():
            M=model['degrees'][r]['matrix'];records=[]
            if state['rows']:
                V=matrix(QQ,state['rows'])
                U=M.transpose().solve_right(V.transpose()).transpose()
                assert U*M==V and V.rank()==len(state['rows']),('Pixton membership',target.slug,g,r)
                for v,u,origin in zip(V.rows(),U.rows(),state['origins']):
                    records.append({'origin':origin,'strata_vector':[str(c) for c in v],
                        'pixton_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
            dump(folder/('relations-g%s-r%s.json'%(g,r)),records,compressed=True)
            gs['degrees'].append({'codimension':r,'candidates':state['candidates'],
                'relation_rank':len(records),'rank_by_base':state['growth'],
                'pixton_rank':int(M.rank()),'strata_columns':int(M.ncols())})
        summary['genera'].append(gs)
        summary['seconds']=time.time()-started;dump(folder/'summary.json',summary)
        print('CERTIFIED',target.slug,'g',g,'ranks',
              [x['relation_rank'] for x in gs['degrees']],'top',top,flush=True)
        reset_caches()
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--max-degree',type=int,choices=(1,2,3),default=3)
    parser.add_argument('--bases',type=int,nargs='+',default=[4,5])
    parser.add_argument('--targets',nargs='*')
    ARGS=parser.parse_args()
    if any(b<4 for b in ARGS.bases):raise ValueError('bases must be at least four to avoid cover resonances')
    OUT.mkdir(parents=True,exist_ok=True)
    targets=requested_complete_intersections()
    if ARGS.targets:targets=tuple(t for t in targets if t.slug in ARGS.targets)
    if not targets:raise ValueError('no matching targets')
    prepare()
    catalog=[]
    with ProcessPoolExecutor(max_workers=ARGS.workers,mp_context=multiprocessing.get_context('fork')) as pool:
        futures={pool.submit(scan_target,t.record()):t.slug for t in targets}
        for future in as_completed(futures):
            catalog.append(future.result())
            dump(OUT/'catalog.json',sorted(catalog,key=lambda x:x['target']['slug']))
    print('COMPLETE',len(catalog),'targets',flush=True)
