"""Audit every available direct four-mark unit coefficient against forgetting."""
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_forgetting import unit_coefficients
from calculations.frontier_certify import pixton,exact_span
verified=0;records=[]
for g in (2,3):
    spaces={r:exact_span(pixton(g,4,r)) for r in range(3)}
    for i,t in enumerate(run.targets()):
        path=run.result_path(i,g,4,2)
        if not path.exists():continue
        data=run.read(path)
        if data.get('unit_profiles_reconstructed') or len(data['profiles'])!=15:continue
        predicted=unit_coefficients(i,g,4,2)
        key=lambda x:(tuple(x['beta']),tuple(x['profile']),x['codimension'])
        left={key(x):x['vector'] for x in data['finite_canonical_classes'] if '1' in x['profile']}
        right={key(x):x['vector'] for x in predicted}
        for r in range(3):
            rows=[];V=spaces[r]
            for k in left.keys()|right.keys():
                if k[2]!=r:continue
                v=lambda d:vector(QQ,V.ncols(),{int(j):QQ(c) for j,c in d.items()})
                rows.append(v(left.get(k,{}))-v(right.get(k,{})))
            if rows:
                A=matrix(QQ,rows);assert A.matrix_from_columns(V.pivots())*V==A,(t.slug,g,r)
            verified+=len(rows)
        records.append((t.slug,g));print('UNIT AUDIT',t.slug,g,verified,flush=True)
run.dump(run.OUT/'unit-validation.json',{'direct_target_genus_cases':len(records),
         'nonzero_or_nontrivial_representatives_compared':verified,'cases':records,'all_passed':True})
