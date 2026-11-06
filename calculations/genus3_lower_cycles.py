"""Compute class-valued genus-three localization and exact Pixton certificates."""
import json,time
from pathlib import Path
from sage.all import *

from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.cycles.plane import (
    O3UnmarkedCycleBackend,
)
_load_vendored_admcycles()
from admcycles import DR
from admcycles.admcycles import generating_indices
from admcycles.moduli import get_moduli

out=Path('results/genus3/lower-cycles-2026-11-05')
out.mkdir(parents=True,exist_ok=True)
b=O3UnmarkedCycleBackend((0,1,3))
alternate=O3UnmarkedCycleBackend((0,2,5))
records=[]
for degree in (0,1):
    for r in (3,4,5,6):
        records.append(b.record(degree,r))
(out/'cycle-formulas.json').write_text(json.dumps(records,indent=2,default=int)+'\n')
print('Saved exact cycle formulas',flush=True)
top={d:str(b.top_integral(d)) for d in (0,1)}
assert QQ(top[0])==QQ(1)/40320
assert QQ(top[1])==-QQ(1637)/80640
print('Top-degree cross-checks',top,flush=True)

summary={'top_integrals':top,'degrees':[]}
for r in (3,4,5):
    started=time.time()
    print('Computing codimension',r,flush=True)
    raw=DR.rels_matrix(3,r,0,symm=0,moduli_type=get_moduli('st',DRpy=True),usespin=True)
    # Convert from Pixton's coordinates to admcycles' decorated-monomial coordinates.
    M=matrix(QQ,[DR.convert_vector_to_monomial_basis(row,3,r,(),
        moduli_type=get_moduli('st',DRpy=True)) for row in raw.rows()])
    rowspace=M.row_space()
    pivots=rowspace.basis_matrix().pivots()
    print('Pixton matrix',M.nrows(),M.ncols(),'rank',len(pivots),flush=True)
    generators=generating_indices(3,0,r,FZ=True,FZmethod='3spin')
    # Store exact generators and the actual relation matrix for reproducing membership checks.
    save(M,str(out/('pixton-codimension-%s.sobj'%r)))
    (out/('strata-order-codimension-%s.json'%r)).write_text(json.dumps(
        [{'index':i,'stratum':str(s)} for i,s in enumerate(b.R.generators(r))],
        indent=2,default=int)+'\n')
    classes={}
    for degree in (0,1):
        classes[degree]=b.cycle(degree,r)
        print('Converted finite class',degree,r,flush=True)
    certificates=[]
    path_certificates=[]
    for degree in (0,1):
        for power,terms in sorted(b.coefficients(degree,r).items()):
            if power>=0:continue
            cycle=b.cycle(degree,r,power)
            v=vector(QQ,cycle.vector(r))
            assert v in rowspace,('uncancelled auxiliary pole',degree,r,power)
            weights=M.transpose().solve_right(v)
            assert weights*M==v
            certificates.append({'kind':'auxiliary_pole','ambient_degree':degree,
                'epsilon_power':int(power),'formal_terms':[{'key':k,'coefficient':str(c)} for k,c in terms.items()],
                'strata_vector':[str(x) for x in v],
                'pixton_weights':[[int(i),str(w)] for i,w in enumerate(weights) if w]})
            print('Certified pole',degree,r,power,flush=True)
        diff=classes[degree]-alternate.cycle(degree,r)
        v=vector(QQ,diff.vector(r))
        assert v in rowspace,('path dependent limit',degree,r)
        weights=M.transpose().solve_right(v)
        assert weights*M==v
        path_certificates.append({'ambient_degree':degree,
            'strata_vector':[str(x) for x in v],
            'pixton_weights':[[int(i),str(w)] for i,w in enumerate(weights) if w]})
        print('Path independence',degree,r,flush=True)
    # The degree-zero answer has an independent mapping-to-a-point derivation.
    expected={3:b.R.zero(),4:9*b.hodge_class(3,0,(1,3)),
              5:-36*b.hodge_class(3,0,(2,3))}[r]
    assert vector(QQ,(classes[0]-expected).vector(r)) in rowspace
    candidates=matrix(QQ,[c['strata_vector'] for c in certificates])
    entry={'codimension':r,'t_power':2-r,'pixton_rows':int(M.nrows()),
        'strata_dimension':int(M.ncols()),'pixton_rank':int(M.rank()),
        'quotient_dimension':int(M.ncols()-M.rank()),
        'pole_relations':len(certificates),'pole_relation_rank':int(candidates.rank()),
        'finite_class_vectors':{str(d):[str(x) for x in classes[d].vector(r)] for d in (0,1)},
        'finite_class_basis_vectors':{str(d):[str(x) for x in classes[d].basis_vector(r)] for d in (0,1)},
        'basis_indices':[int(x) for x in generators],
        'path_independence_verified':True,'seconds':time.time()-started}
    (out/('pixton-certificates-codimension-%s.json'%r)).write_text(json.dumps(certificates,indent=2,default=int)+'\n')
    (out/('path-certificates-codimension-%s.json'%r)).write_text(json.dumps(path_certificates,indent=2,default=int)+'\n')
    summary['degrees'].append(entry)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=int)+'\n')
    print('Completed codimension',r,'quotient dimension',entry['quotient_dimension'],
          'pole rank',entry['pole_relation_rank'],'seconds',entry['seconds'],flush=True)
