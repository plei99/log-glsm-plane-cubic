"""Three-mark descendants and labeled Chow functoriality."""

from log_glsm.core.hodge_grr import (
    GenusFourHodgeIntegrals,
)
from log_glsm.cycles.higher_marked import (
    higher_marked_backend,
    higher_marked_ordinary_seed,
    higher_marked_profiles,
    marked_profile_permutations,
)
from log_glsm.cycles.marked import (
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
from log_glsm.targets.grassmannian import (
    grassmannian_seven_hyperplane_cy3,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
)
assert len(higher_marked_profiles())==40
assert len(higher_marked_profiles(canonical=True))==20
h=GenusFourHodgeIntegrals()
for g,expected in ((2,(24,12,4)),(3,(120,30,6))):
    assert marked_multiplier(g,1,('psi','psi','psi'))==expected[0]
    assert marked_multiplier(g,1,('psi','psi','H'))==expected[1]
    assert marked_multiplier(g,1,('psi','H','H'))==expected[2]
    b=higher_marked_backend(genus_four_targets()[0],g,hodge=h)
    z=b.numerical_localization((1,))['finite']
    for profile in higher_marked_profiles(canonical=True):
        if len(profile)!=3:continue
        value=b.numerical_localization((1,),profile)
        assert not value['poles'] and value['finite']==marked_multiplier(g,1,profile)*z,(g,profile,value)
    print('Direct genus',g,'three-mark localization passed',flush=True)
# All labeled placements of a repeated profile are represented exactly once.
p=marked_profile_permutations(('1','psi','psi'))
assert {x[0] for x in p}=={('1','psi','psi'),('psi','1','psi'),('psi','psi','1')}
for renamed,mapping in p:
    assert all(renamed[new-1]==('1','psi','psi')[old] for old,new in enumerate(mapping))
# Real marked stabilization: unit pullback from two to three marks.
b=higher_marked_backend(ProjectiveCompleteIntersection(2,(3,)),2)
def C(profile,r):
    coeffs=b.marked_coefficients((1,),profile,r)[r].get(0,{})
    return sum((co*b.marked_monomial_class(key) for key,co in coeffs.items()),b.adm.TautologicalRing(2,len(profile)).zero())
for r in (1,2):
    assert (C(('H','H','1'),r)-C(('H','H'),r).forgetful_pullback([3])).is_zero()
    assert (C(('H','H','psi'),r+1).forgetful_pushforward([3])-4*C(('H','H'),r)).is_zero()
print('Three-mark unit and dilaton Chow identities passed',flush=True)
left=b.marked_coefficients((1,),('1','psi','H'),2)
right=b.marked_coefficients((1,),('H','1','psi'),2)
R=b.adm.TautologicalRing(2,3)
for r in range(3):
    for power in set(left[r])|set(right[r]):
        lc=sum((co*b.marked_monomial_class(key) for key,co in left[r].get(power,{}).items()),R.zero())
        rc=sum((co*b.marked_monomial_class(key) for key,co in right[r].get(power,{}).items()),R.zero())
        assert (lc.rename_legs({1:2,2:3,3:1},inplace=False)-rc).is_zero(),(r,power)
print('Direct labeled Chow permutation check passed',flush=True)
t=grassmannian_seven_hyperplane_cy3()
assert higher_marked_ordinary_seed(t,2,(2,))==QQ(539)/80
assert higher_marked_ordinary_seed(t,3,(2,))==QQ(133)/288
assert marked_multiplier(2,2,('psi','psi','psi'))*QQ(539)/80==QQ(1617)/10
assert marked_multiplier(3,2,('psi','psi','psi'))*QQ(133)/288==QQ(665)/12
print('Higher marked tests passed',flush=True)
