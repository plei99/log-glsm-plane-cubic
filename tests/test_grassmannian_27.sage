"""Geometry, signs, and independent low-degree counts of Gr(2,7)[1^7]."""

from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersection,
    grassmannian_known_ordinary_gw,
    grassmannian_line_count,
    grassmannian_seven_hyperplane_cy3,
)
from log_glsm.targets.toric import (
    toric_genus_zero_number,
)

t=grassmannian_seven_hyperplane_cy3();X=t.ambient;H=X.divisor((1,))
assert t.dimension==3 and t.kind=='CY3'
assert X.record()['fixed_points']==21 and X.record()['invariant_curves']==105
assert t.record()['pluecker_degree']==42
assert t.chern_numbers()=={'c1_cubed':0,'c1_c2':0,'euler':-98,'C':-98}
assert X.integral(t.fundamental*H*t.chern[2])==84
x,y=X.P.gens();P=X.chow(x*y)
assert t.chern[2]==H**2+3*P and t.chern[3]==-7*H*P
for weights in (range(7),(0,1,4,10,22,46,94)):
    assert grassmannian_line_count(t,weights)==196
for base in (4,5):
    assert toric_genus_zero_number(t,(1,),base)==196
    assert toric_genus_zero_number(t,(2,),base)==QQ(2499)/2
for g in (2,3):
    assert t.alpha(g,(2,))==3*g-3
    assert [t.comparison_sign(g,(d,)) for d in (0,1,2)]==([-1,1,-1] if g==2 else [1,-1,1])
assert grassmannian_known_ordinary_gw(t,2,(0,))==-QQ(49)/2880
assert grassmannian_known_ordinary_gw(t,3,(0,))==QQ(7)/103680
assert grassmannian_known_ordinary_gw(t,2,(1,))==QQ(49)/60
assert grassmannian_known_ordinary_gw(t,3,(1,))==QQ(7)/216
assert grassmannian_known_ordinary_gw(t,2,(2,))==QQ(539)/80
assert grassmannian_known_ordinary_gw(t,3,(2,))==QQ(133)/288
renamed=GrassmannianCompleteIntersection(7,(1,)*7,'renamed')
assert grassmannian_known_ordinary_gw(renamed,3,(2,))==QQ(133)/288
assert grassmannian_known_ordinary_gw(t,3,(3,)) is None
print('Gr(2,7) Chern, independent flag-count, genus-zero, and comparison-sign checks passed')
