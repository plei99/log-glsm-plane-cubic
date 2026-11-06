"""More marks: independent dilaton values and actual Chow forgetting."""

from log_glsm.cycles.higher_marked import (
    higher_marked_backend,
    higher_marked_profiles,
)
from log_glsm.cycles.marked import (
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
assert len(higher_marked_profiles(max_markings=4))==121
assert len(higher_marked_profiles(max_markings=5))==364
for g,expected4,expected5 in ((2,120,720),(3,840,6720)):
    assert marked_multiplier(g,1,('psi',)*4)==expected4
    assert marked_multiplier(g,1,('psi',)*5)==expected5
    b=higher_marked_backend(genus_four_targets()[0],g)
    z=b.numerical_localization((1,))['finite']
    for profile in (('psi',)*4,('H',)*4,('1','psi','H','H'),('psi',)*5):
        v=b.numerical_localization((1,),profile)
        assert not v['poles'] and v['finite']==marked_multiplier(g,1,profile)*z,(g,profile)
    print('Direct four/five-mark numbers',g,'passed',flush=True)
b=higher_marked_backend(genus_four_targets()[0],2)
def C(profile,r):
    co=b.marked_coefficients((1,),profile,r)[r].get(0,{})
    R=b.adm.TautologicalRing(2,len(profile))
    return sum((a*b.marked_monomial_class(k) for k,a in co.items()),R.zero())
for r in (1,2):
    assert (C(('H','H','H','1'),r)-C(('H','H','H'),r).forgetful_pullback([4])).is_zero()
assert (C(('H','H','H','psi'),2).forgetful_pushforward([4])-5*C(('H','H','H'),1)).is_zero()
print('Four-mark unit/dilaton Chow identities passed',flush=True)
import sys,gzip,json
from calculations import marked_frontier
from tempfile import TemporaryDirectory
from pathlib import Path
_output = TemporaryDirectory()
marked_frontier.OUT = Path(_output.name)
from calculations.marked_frontier import scan,read,result_path
scan(0,2,3)
fresh=read(result_path(0,2,3,2))
with gzip.open('tests/fixtures/p2-3-g2-classes.json.gz','rt') as f:old=json.load(f)
for r in range(3):
    size=len(b.adm.TautologicalRing(2,3).generators(r))
    def span(rows):return matrix(QQ,[vector(QQ,size,{int(i):QQ(c) for i,c in row.items()}) for row in rows]).row_space() if rows else VectorSpace(QQ,size).zero_subspace()
    previous=next(s for s in old['relations'] if s['markings']==3 and s['codimension']==r)
    current=next(s for s in fresh['relations'] if s['codimension']==r)
    assert span(previous['rows'])==span(current['rows'])
print('Generator-based orbit closure reproduces previous full-permutation spans',flush=True)
