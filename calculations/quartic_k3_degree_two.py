"""Extend the quartic search successively to degrees two and three."""
import argparse,gzip,json,time
from pathlib import Path
from sage.all import *

from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
    QuarticK3Theory,
)
out=Path('results/quartic-k3/2026-11-05')
parser=argparse.ArgumentParser()
parser.add_argument('--degree',type=int,choices=(2,3),default=2)
parser.add_argument('--start-ray',type=int,default=0,
                    help='continue a completed degree-three scan from its saved basis')
args=parser.parse_args()
degree=args.degree
suffix='degree-two' if degree==2 else 'degree-three'
rays=([(0,1,3,7)]+[(0,1,q,q*q+q+1) for q in range(3,10)] if degree==2
      else [(0,1,q,q*q) for q in range(5,13)])
if degree==3:
    # Add rays off the one-parameter quadratic family. Reject all midpoint
    # and 1:2 resonances that occur on covers of degree at most three.
    for p,q in ((4,13),(4,17),(5,19),(5,23),(6,29),(7,31),(8,37),(9,41)):
        weights=(0,1,p,q)
        assert all((d-m)*x+m*y!=d*z for x in weights for y in weights for z in weights
                   if len({x,y,z})==3 for d in (2,3) for m in range(1,d))
        rays.append(weights)
if args.start_ray and degree!=3:raise ValueError('continuation is only supported for degree three')

def dump(name,data):
    (out/name).write_text(json.dumps(data,indent=2,default=int)+'\n')

def terms_record(terms):
    return [{'key':k,'coefficient':str(c)} for k,c in terms.items() if c]

summary={'ambient_degree':degree,'rays':rays,'genera':[]}
if args.start_ray:
    summary['previous_scan']=json.loads((out/(suffix+'-summary.json')).read_text())
    summary['start_ray']=args.start_ray
for g in (2,3):
    started=time.time()
    b=QuarticK3GraphCycles(g,rays[0]); states={}
    for r in range(3*g-2):
        prefix='relations' if degree==2 else ('final-relations' if args.start_ray else 'combined-relations')
        path=out/('%s-g%s-r%s.json'%(prefix,g,r))
        records=json.loads(path.read_text() if path.exists() else gzip.decompress(
            Path(str(path)+'.gz').read_bytes()).decode())
        M=load(str(out/('pixton-g%s-r%s.sobj'%(g,r))))
        state={'matrix':M,'records':records,'echelon':{},'monomials':{},
               'base_rank':len(records),'candidates':0,'growth':[]}
        def insert(v):
            remainder=copy(v)
            for pivot,row in sorted(state['echelon'].items()):
                if remainder[pivot]:remainder-=remainder[pivot]*row
            if not remainder:return False
            pivot=next(i for i,x in enumerate(remainder) if x)
            state['echelon'][pivot]=remainder/remainder[pivot]
            return True
        for record in records:assert insert(vector(QQ,record['strata_vector']))
        states[r]=state

    def vector_of(terms,r):
        state=states[r]; value=vector(QQ,state['matrix'].ncols())
        for key,co in terms.items():
            if key not in state['monomials']:
                state['monomials'][key]=vector(QQ,b.monomial_class(key).vector(r))
            value+=co*state['monomials'][key]
        return value

    def consider(terms,r,origin):
        state=states[r];state['candidates']+=1
        v=vector_of(terms,r);remainder=copy(v)
        for pivot,row in sorted(state['echelon'].items()):
            if remainder[pivot]:remainder-=remainder[pivot]*row
        if not remainder:return
        pivot=next(i for i,x in enumerate(remainder) if x)
        state['echelon'][pivot]=remainder/remainder[pivot]
        state['records'].append({'origin':origin,'terms':terms_record(terms),
                                 'strata_vector':[str(c) for c in v]})

    finite=[];reference={}
    for r in states:
        reference[r]=b.coefficients(degree,r).get(0,{})
        finite.append({'genus':g,'ambient_degree':degree,'codimension':r,
                       't_power':int(QuarticK3Theory.t_power(g,r)),
                       'is_cjr_relation':r>=2*g-2,'terms':terms_record(reference[r])})
    dump('finite-classes-degree%s-g%s.json'%(degree,g),finite)
    print('genus',g,'degree',degree,'graphs',len(b.graphs(degree)),flush=True)
    # The general graph implementation must reproduce the separate direct
    # degree-one tail/bridge formula, including the symmetric genus-two bridge.
    if degree==2:
        direct=QuarticK3Cycles(g,rays[0])
        for r,state in states.items():
            diff=vector_of(b.coefficients(1,r).get(0,{}),r)-vector_of(direct.coefficients(1,r).get(0,{}),r)
            assert diff in state['matrix'].row_space()
            print('compared genus',g,'codimension',r,flush=True)
        print('genus',g,'direct degree-one comparison passed',flush=True)
    for index,weights in enumerate(rays):
        if index<args.start_ray:continue
        current=b if index==0 else QuarticK3GraphCycles(g,weights)
        for r in states:
            for power,terms in current.coefficients(degree,r).items():
                if power<0:consider(terms,r,{'kind':'pole','degree':degree,
                                           'ray':weights,'epsilon_power':power})
            difference=dict(current.coefficients(degree,r).get(0,{}))
            for key,co in reference[r].items():difference[key]=difference.get(key,0)-co
            difference={k:c for k,c in difference.items() if c}
            if difference:consider(difference,r,{'kind':'finite_part_difference','degree':degree,
                                                 'ray':weights,'reference_ray':rays[0]})
            states[r]['growth'].append(len(states[r]['records']))
            if index==0:print('processed genus',g,'degree',degree,'codimension',r,flush=True)
        print('genus',g,'ray',index+1,'ranks',
              {r:len(s['records']) for r,s in states.items()},flush=True)
    gs={'genus':g,'fixed_graphs':len(b.graphs(degree)),
        'top_integral':str(b.top_integral(degree)),'degrees':[]}
    for r,state in states.items():
        records=state['records']; M=state['matrix']
        new=records[state['base_rank']:]
        if new:
            V=matrix(QQ,[record['strata_vector'] for record in new])
            U=M.transpose().solve_right(V.transpose()).transpose()
            assert U*M==V
            for record,u in zip(new,U.rows()):
                record['pixton_weights']=[[int(i),str(c)] for i,c in enumerate(u) if c]
        prefix='combined-relations' if degree==2 else 'final-relations'
        dump('%s-g%s-r%s.json'%(prefix,g,r),records)
        gs['degrees'].append({'codimension':r,'candidates':state['candidates'],
            'previous_rank':state['base_rank'],'combined_rank':len(records),
            'rank_by_ray':state['growth']})
        print('certified genus',g,'codimension',r,'new relations',len(new),flush=True)
    gs['seconds']=time.time()-started
    summary['genera'].append(gs);dump(suffix+'-summary.json',summary)
    print('completed genus',g,'top value',gs['top_integral'],'seconds',gs['seconds'],flush=True)
