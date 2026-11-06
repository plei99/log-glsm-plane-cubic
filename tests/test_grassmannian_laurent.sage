"""The two exact evaluation strategies must give every requested coefficient."""

from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersection,
    GrassmannianCompleteIntersectionCycles,
)
from time import perf_counter

for n,g,D,degrees in ((4,3,2,(4,)),(7,2,1,(1,)*7),(7,2,2,(1,)*7)):
    target=GrassmannianCompleteIntersection(n,degrees)
    fast=GrassmannianCompleteIntersectionCycles(target,g,graphwise_series=True)
    reference=GrassmannianCompleteIntersectionCycles(target,g,graphwise_series=False)
    clean=lambda cs:{p:{k:c for k,c in v.items() if c} for p,v in cs.items() if any(v.values())}
    start=perf_counter()
    actual=[clean(fast.coefficients((D,),r)) for r in range(3*g-2)]
    fast_seconds=perf_counter()-start;start=perf_counter()
    expected=[clean(reference.coefficients((D,),r)) for r in range(3*g-2)]
    reference_seconds=perf_counter()-start
    assert actual==expected,(n,g,D)
    print('ALL COEFFICIENTS MATCH',n,g,D,'graphwise',fast_seconds,'rational',reference_seconds,flush=True)
