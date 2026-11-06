"""The extra unit carries a rational-tail correction for each psi insertion."""

from log_glsm.cycles.marked import (
    genus_four_backend,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
b=genus_four_backend(genus_four_targets()[0],genus=2)
def C(profile,r):
    data=b.marked_coefficients((1,),profile,r)[r].get(0,{})
    R=b.adm.TautologicalRing(2,len(profile))
    return sum((a*b.marked_monomial_class(k) for k,a in data.items()),R.zero())
M=load('tests/fixtures/pixton-g2-n4-r2.sobj')
for p in (('psi','H','H'),('psi','psi','H'),('1','psi','H')):
    actual=C(p+('1',),2)
    predicted=C(p,2).forgetful_pullback([4])
    correction=b.adm.TautologicalRing(2,4).zero()
    for i,x in enumerate(p,1):
        if x=='psi':
            lower=list(p);lower[i-1]='1'
            correction+=b.adm.sepbdiv(0,[i,4],2,4)*C(tuple(lower),1).forgetful_pullback([4])
    assert vector(QQ,(actual-predicted-correction).vector(2)) in M.row_space()
    print(p,'boundary correction nonzero:',vector(QQ,correction.vector(2)) not in M.row_space(),flush=True)
print('Descendant unit-forgetting identities passed',flush=True)
import sys
from calculations.frontier_forgetting import section_matrix
D=section_matrix(2,4,2,1);W=M.row_space()
for generator,row in zip(b.adm.TautologicalRing(2,3).generators(1),D.rows()):
    expected=b.adm.sepbdiv(0,[1,4],2,4)*generator.forgetful_pullback([4])
    assert vector(QQ,expected.vector(2))-row in W
print('Rational-tail section matrix agrees with divisor multiplication on every generator',flush=True)
