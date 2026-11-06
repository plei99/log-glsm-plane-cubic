"""Independent Hodge formulas, legacy values, and direct marked localization."""
from itertools import product
from log_glsm.core.hodge import (
    AdmcyclesHodgeIntegralBackend,
)
from log_glsm.core.hodge_grr import (
    GenusFourHodgeIntegrals,
    mumford_monomial,
)
from log_glsm.cycles.marked import (
    genus_four_backend,
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
    known_genus_four_ordinary,
)
from log_glsm.targets.grassmannian import (
    grassmannian_seven_hyperplane_cy3,
)
from log_glsm.targets.projective import (
    CompleteIntersectionCycles,
    ProjectiveCompleteIntersection,
)
import time
start=time.time();h=GenusFourHodgeIntegrals()
assert len(genus_four_targets())==43
assert dict(mumford_monomial(4,(2,2)))=={(1,3):2,(4,):-2}
assert dict(mumford_monomial(4,(3,3)))=={(2,4):2}
assert not mumford_monomial(4,(4,4))
assert not mumford_monomial(4,(1,2,3,4))
# Test the GRR recursion against closed formulas, bypassing its shortcuts.
for g in (2,3,4):
    for indices,powers,expected in (
        ((g,),(2*g-2,),h.lambda_g_constant(g)),
        (h.top_hodge_triple_indices(g),(),h.top_hodge_triple_constant(g))):
        actual=sum(co*h.character_integral(g,powers,chars) for chars,co in h._lambda_characters(indices))
        assert actual==expected,(g,indices,actual,expected)
print('GRR versus lambda-g and Faber-Pandharipande: passed',flush=True)
old=AdmcyclesHodgeIntegralBackend()
for case in ((4,(9,),(1,)),(4,(8,),(2,)),(3,(3,2),(1,2)),(2,(2,2),(1,))):
    assert h.integral(*case)==old.integral(*case),case
print('Mixed Hodge integrals versus independent admcycles evaluation: passed',flush=True)
# Previously verified independent numbers, using a new numerical algorithm.
t=grassmannian_seven_hyperplane_cy3()
for g,d,expected in ((2,1,QQ(1351)/288),(2,2,-QQ(19229)/240),(3,2,-QQ(55586069)/414720)):
    result=genus_four_backend(t,hodge=h,genus=g).numerical_localization((d,))
    assert result['finite']==expected and not result['poles'],(g,d,result)
assert known_genus_four_ordinary(t,(0,))==-QQ(7)/6220800
assert known_genus_four_ordinary(t,(1,))==QQ(49)/43200
assert known_genus_four_ordinary(t,(2,))==QQ(833)/19200
print('Independent legacy localization and ordinary genus-four values: passed',flush=True)
# Include unstable marked endpoints and contracted rational bubbles.
t=genus_four_targets()[0];b=genus_four_backend(t,hodge=h)
z=b.numerical_localization((1,))['finite']
for n in (1,2):
    for profile in product(('1','psi','H'),repeat=n):
        result=b.numerical_localization((1,),profile)
        assert not result['poles'] and result['finite']==marked_multiplier(4,1,profile)*z,(profile,result,z)
print('All twelve ordered marked profiles by direct localization: passed',flush=True)
# Recover the independent unmarked surface constant as an actual Chow class.
t=ProjectiveCompleteIntersection(3,(2,));b=genus_four_backend(t,hodge=h)
coeffs=b.marked_coefficients((0,),('H',),max_codimension=10)
cycle=sum((co*b.marked_monomial_class(key) for key,co in coeffs[10].get(0,{}).items()),b.adm.TautologicalRing(4,1).zero())
assert cycle.evaluate()==0  # divisor equation at beta=0
print('Marked degree-zero Chow stabilization: passed',flush=True)
# Compare the new stabilization with the independently existing g=3 module.
t=ProjectiveCompleteIntersection(3,(4,));new=genus_four_backend(t,genus=3);legacy=CompleteIntersectionCycles(t,3)
c=new.marked_coefficients((2,),(),3)
for r in range(4):
    actual=sum((co*new.marked_monomial_class(key) for key,co in c[r].get(0,{}).items()),new.R.zero())
    assert (actual-legacy.cycle(2,r)).is_zero(),r
print('Unmarked stabilization versus prior genus-three classes: passed',flush=True)
# Functorial marking identities are class statements, not just top integrals.
b=genus_four_backend(ProjectiveCompleteIntersection(2,(3,)),genus=2)
def C(profile,r):
    terms=b.marked_coefficients((1,),profile,r)[r].get(0,{})
    return sum((co*b.marked_monomial_class(k) for k,co in terms.items()),b.adm.TautologicalRing(2,len(profile)).zero())
for r in (1,2):
    assert (C(('1',),r)-C((),r).forgetful_pullback([1])).is_zero()
    assert (C(('psi',),r+1).forgetful_pushforward([1])-2*C((),r)).is_zero()
    assert (C(('H',),r+1).forgetful_pushforward([1])-C((),r)).is_zero()
print('Unit pullback, dilaton, and divisor equations as Chow identities: passed',flush=True)
print('TOTAL seconds',time.time()-start,flush=True)
