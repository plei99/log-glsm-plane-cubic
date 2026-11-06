"""Independent enumerative checks and boundaries between ordinary and twisted GW."""

from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.projective import (
    CompleteIntersectionCycles,
    ProjectiveCompleteIntersection,
    complete_intersection_weights,
    cy3_line_count_grassmannian,
    cy3_line_count_localization,
    cy3_ordinary_gw_number,
    requested_complete_intersections,
)

targets=requested_complete_intersections()
assert len(targets)==14 and len({t.slug for t in targets})==14
assert [t.euler_characteristic() for t in targets]==[4,9,4,-6,-56,6,27,188,825,-200,-176,-144,-144,-128]
assert [t.index for t in targets]==[2,1,3,2,1,4,3,2,1,0,0,0,0,0]

# These enumerative values come from Grassmannian geometry / the cited literature,
# independently of the higher-genus fixed-graph backend.
for t,lines in zip(targets[-5:],(2875,1280,1053,720,512)):
    assert cy3_line_count_grassmannian(t)==lines
    assert cy3_line_count_localization(t)==lines
    assert cy3_line_count_localization(t,complete_intersection_weights(t.N+1,5))==lines
    for g in (2,3):
        assert all(t.t_power(g,D,r)>=0 for D in range(4) for r in range(3*g-2))

for t,co in zip(targets[:2],(8,3)):
    assert t.compact_degree_zero_terms(2)=={('interior',None,(2,)):co}
    assert t.compact_degree_zero_terms(3)=={('interior',None,(1,3)):co}
for t,co in zip(targets[2:5],(-20,-30,-80)):
    assert t.compact_degree_zero_terms(2)=={('interior',None,(1,2)):co}
    assert t.compact_degree_zero_terms(3)=={('interior',None,(1,2,3)):-co}
for t in targets[5:9]:
    assert t.compact_degree_zero_terms(2)=={} and t.compact_degree_zero_terms(3)=={}

# The new parameterized backend must recover the independently checked quartic
# formulas, including nonconstant maps, in every codimension.
quartic=ProjectiveCompleteIntersection(3,(4,))
for g in (2,3):
    old=QuarticK3GraphCycles(g,(0,1,4,5))
    new=CompleteIntersectionCycles(quartic,g,(0,1,4,5))
    for D in (0,1):
        for r in range(3*g-2):
            assert old.coefficients(D,r)==new.coefficients(D,r),(g,D,r)

# Ordinary CY theory is not K3 theory: a nonconstant compact class stays unknown.
cy=CompleteIntersectionCycles(targets[-1],3)
record=cy.aggregate_relation(1,6)
assert record['compact_known'] is False and record['compact_class'] is None
assert record['infinity_aggregate'] is None
assert record['comparison_sign']==1
try:
    cy.aggregate_relation(1,5)
except ValueError:
    pass
else:
    raise AssertionError('positive t powers must not become vanishing relations')
assert targets[9].comparison_sign(3,1)==-1
assert targets[10].comparison_sign(3,1)==1
assert targets[12].comparison_sign(3,1)==-1
print('complete-intersection enumerative, grading, constant-map, and quartic-regression checks passed')

quintic=ProjectiveCompleteIntersection(4,(5,))
assert cy3_ordinary_gw_number(quintic,2,0)==-QQ(5)/144
assert cy3_ordinary_gw_number(quintic,3,0)==QQ(5)/36288
assert cy3_ordinary_gw_number(quintic,2,1)==QQ(575)/48
assert cy3_ordinary_gw_number(quintic,3,1)==QQ(2875)/6048
try:
    cy3_ordinary_gw_number(quintic,3,2)
except NotImplementedError:
    pass
else:
    raise AssertionError('unknown degree-two ordinary GW theory must not be replaced by zero')
print('ordinary CY3 degree-zero and degree-one GW checks passed')

negative=CompleteIntersectionCycles(ProjectiveCompleteIntersection(3,(2,)),3).aggregate_relation(0,6)
assert negative['compact_known'] and not negative['compact_class']
assert negative['infinity_aggregate']+negative['zero_vertex']==0
print('known and unknown infinity-aggregate API checks passed')

# Published-in-repo quartic regression values test higher-cover edges and the
# degree-three stabilization cases through the public numerical interface.
quartic_covers=CompleteIntersectionCycles(ProjectiveCompleteIntersection(3,(4,)),2)
assert quartic_covers.top_integral(2)==-QQ(67)/36
assert quartic_covers.top_integral(3)==QQ(133)/2
print('degree-two and degree-three quartic specialization checks passed')
