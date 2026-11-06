"""A declared higher-codimension probe, separate from the all-example window."""
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_certify import exact_span,pixton
print(run.scan(0,3,3,3),flush=True)
data=run.read(run.result_path(0,3,3,3));item=next(s for s in data['relations'] if s['codimension']==3)
M=pixton(3,3,3)
A=matrix(QQ,[vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in row.items()}) for row in item['rows']])
V=exact_span(A);U=M.transpose().solve_right(V.transpose()).transpose();assert U*M==V
run.dump(run.OUT/'codimension-three-probe.json',{'target':run.targets()[0].slug,'genus':3,'markings':3,
    'codimension':3,'curve_bound':2,'found_rank':V.nrows(),'pixton_rank':M.rank(),'strata':M.ncols(),
    'certificate_verified':True,'V':[run.common.sparse(v) for v in V.rows()],'U':[run.common.sparse(v) for v in U.rows()]})
print('PROBE',V.nrows(),M.rank(),M.ncols(),flush=True)
