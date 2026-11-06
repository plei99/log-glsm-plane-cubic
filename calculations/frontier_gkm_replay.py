"""Replay marked GKM classes against the earlier direct checkpoints."""
from sage.all import *
from calculations import marked_frontier as run
cases=((46,(1,1,0),('1',)*4),(51,(0,2),('1','psi','psi','H')))
for index,beta,profile in cases:
    target=run.targets()[index];data=run.read(run.result_path(index,2,4,2))
    assert not data.get('unit_profiles_reconstructed')
    expected={x['codimension']:x['vector'] for x in data['finite_canonical_classes']
              if tuple(x['beta'])==beta and tuple(x['profile'])==profile}
    backend=run.common.higher_marked_backend(target,2,4)
    coefficients=backend.marked_coefficients(beta,profile,2)
    nonzero=0
    for r in range(3):
        v=run.common.vector_of(coefficients[r].get(0,{}),2,4,r)
        assert run.common.sparse(v)==expected.get(r,{}),(target.slug,beta,profile,r)
        nonzero+=bool(v)
    assert nonzero
    print('REPLAY',target.slug,beta,profile,nonzero,'nonzero coefficients; exact vectors match',flush=True)
