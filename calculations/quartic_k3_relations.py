"""Search auxiliary weights, retain independent relations, certify over QQ."""
import argparse
import json
import time
from pathlib import Path
from sage.all import *


from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
    QuarticK3Theory,
)
_load_vendored_admcycles()
from admcycles import DR
from admcycles.moduli import get_moduli

parser=argparse.ArgumentParser()
parser.add_argument('--rays',type=int,default=64)
args=parser.parse_args()
out=Path('results/quartic-k3/2026-11-05')
out.mkdir(parents=True,exist_ok=True)


def terms_record(terms):
    return [{'key':key,'coefficient':str(co)} for key,co in terms.items() if co]


def dump(name,data):
    (out/name).write_text(json.dumps(data,indent=2,default=int)+'\n')


def ray(index):
    # Distinct integral weights; degree-one denominators only need distinctness.
    pairs=[(p,q) for q in range(4,40) for p in range(2,q)]
    if index==0:return (0,1,3,7)
    return (0,1)+[pair for pair in pairs if pair!=(3,7)][index-1]


summary={'scope':{'genera':[2,3],'ambient_degrees':[0,1],
                  'auxiliary_rays':args.rays,'theory':'ordinary K3; dual Euler O(4)'},
         'genera':[]}
for g in (2,3):
    started=time.time()
    backend=QuarticK3Cycles(g,ray(0))
    states={}
    for r in range(3*g-2):
        columns=len(backend.R.generators(r))
        raw=DR.rels_matrix(g,r,0,symm=0,moduli_type=get_moduli('st',DRpy=True),usespin=True)
        rows=[DR.convert_vector_to_monomial_basis(row,g,r,(),
              moduli_type=get_moduli('st',DRpy=True)) for row in raw.rows()]
        M=matrix(QQ,rows) if rows else matrix(QQ,0,columns)
        assert M.ncols()==columns
        states[r]={'matrix':M,'basis':[],'echelon':{},'records':[],
                   'monomials':{},'candidates':0,'zero_vectors':0,'growth':[]}
        save(M,str(out/('pixton-g%s-r%s.sobj'%(g,r))))
        dump('strata-g%s-r%s.json'%(g,r),[str(s) for s in backend.R.generators(r)])

    def vector_of(terms,r):
        state=states[r]
        value=vector(QQ,state['matrix'].ncols())
        for key,co in terms.items():
            if key not in state['monomials']:
                state['monomials'][key]=vector(QQ,backend.monomial_class(key).vector(r))
            value+=co*state['monomials'][key]
        return value

    def consider(terms,r,origin):
        state=states[r]
        state['candidates']+=1
        v=vector_of(terms,r)
        if not v:
            state['zero_vectors']+=1
            return
        remainder=copy(v)
        for pivot,row in sorted(state['echelon'].items()):
            if remainder[pivot]:remainder-=remainder[pivot]*row
        if not remainder:return
        pivot=next(i for i,x in enumerate(remainder) if x)
        state['echelon'][pivot]=remainder/remainder[pivot]
        state['basis'].append(v)
        state['records'].append({'origin':origin,'terms':terms_record(terms),
                                 'strata_vector':[str(x) for x in v]})

    reference={}
    finite=[]
    for degree in (0,1):
        for r in range(3*g-2):
            reference[degree,r]=backend.coefficients(degree,r).get(0,{})
            finite.append({'genus':g,'ambient_degree':degree,'codimension':r,
                           't_power':int(QuarticK3Theory.t_power(g,r)),
                           'is_cjr_relation':r>=2*g-2,
                           'terms':terms_record(reference[degree,r])})
    dump('finite-classes-g%s.json'%g,finite)
    for index in range(args.rays):
        b=backend if index==0 else QuarticK3Cycles(g,ray(index))
        for degree in (0,1):
            for r in states:
                for power,terms in b.coefficients(degree,r).items():
                    if power<0:
                        consider(terms,r,{'kind':'pole','degree':degree,
                                         'ray':ray(index),'epsilon_power':power})
                difference=dict(b.coefficients(degree,r).get(0,{}))
                for key,co in reference[degree,r].items():
                    difference[key]=difference.get(key,0)-co
                difference={k:c for k,c in difference.items() if c}
                if difference:
                    consider(difference,r,{'kind':'finite_part_difference','degree':degree,
                                           'ray':ray(index),'reference_ray':ray(0)})
        for r,state in states.items():state['growth'].append(len(state['basis']))
        if index%8==0 or index==args.rays-1:
            print('genus',g,'ray',index+1,'ranks',
                  {r:len(s['basis']) for r,s in states.items()},flush=True)

    top={d:str(backend.top_integral(d)) for d in (0,1)}
    assert QQ(top[0])==({2:-QQ(1)/288,3:QQ(1)/72576}[g])
    genus_summary={'genus':g,'top_integrals':top,'degrees':[]}
    for r,state in states.items():
        M=state['matrix']
        for record,v in zip(state['records'],state['basis']):
            u=M.transpose().solve_right(v)
            assert u*M==v,('Pixton certificate failed',g,r)
            record['pixton_weights']=[[int(i),str(c)] for i,c in enumerate(u) if c]
        # Independent source derivation of the degree-zero class.
        expected=backend.R.zero()
        if r==3*g-3:
            expected=((-1)**(g+1)*20)*backend.hodge_class(g,0,(g-2,g-1,g))
        diff=vector_of(reference[0,r],r)-vector(QQ,expected.vector(r))
        assert diff in M.row_space(),('constant-map formula',g,r)
        dump('relations-g%s-r%s.json'%(g,r),state['records'])
        genus_summary['degrees'].append({'codimension':r,'candidates':state['candidates'],
            'zero_strata_vectors':state['zero_vectors'],'independent_relations':len(state['basis']),
            'rank_by_ray':state['growth'],'pixton_rows':int(M.nrows()),
            'strata_columns':int(M.ncols()),'pixton_rank':int(M.rank())})
    genus_summary['seconds']=time.time()-started
    summary['genera'].append(genus_summary)
    dump('summary.json',summary)
    print('Certified genus',g,'top values',top,'seconds',genus_summary['seconds'],flush=True)
