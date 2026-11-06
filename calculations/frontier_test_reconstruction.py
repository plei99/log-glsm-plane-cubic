"""Compare reconstructed five-mark classes with direct fixed-locus calculations."""
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_forgetting import complete
from calculations.frontier_certify import pixton
complete(0,2,5,2)
data=run.read(run.result_path(0,2,5,2))
stored={(tuple(x['beta']),tuple(x['profile']),x['codimension']):x['vector'] for x in data['finite_canonical_classes']}
b=run.common.higher_marked_backend(run.targets()[0],2)
for p in (('1','psi','H','H','H'),('1','1','psi','H','H')):
    co=b.marked_coefficients((2,),p,2)
    for r in range(3):
        v=run.common.vector_of(co[r].get(0,{}),2,5,r)
        expected=vector(QQ,len(v),{int(j):QQ(c) for j,c in stored.get(((2,),p,r),{}).items()})
        assert v-expected in pixton(2,5,r).row_space()
    print('Direct versus reconstructed',p,'passed',flush=True)
