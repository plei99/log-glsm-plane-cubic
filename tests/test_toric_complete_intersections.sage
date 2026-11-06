"""Observable geometry checks for the bounded toric extension."""

from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersection,
    ToricCompleteIntersectionCycles,
    ToricProductAmbient,
    toric_degree_tuples,
    toric_fixed_graphs,
    toric_genus_zero_number,
    toric_known_ordinary_gw,
    toric_primitive_curve_counts,
    toric_sample,
)
targets=toric_sample()
assert len(targets)==12
assert [t.chern_numbers()['euler'] for t in targets]==[4,24,4,-168,6,-162,6,24,5,24,8,-144]
assert [targets[i].chern_numbers()['c1_squared'] for i in (0,1,6,7,8,9)]==[8,0,6,0,7,0]
assert [targets[i].chern_numbers()['C'] for i in (2,3,4,5,10,11)]==[-20,-168,-18,-162,-16,-144]
assert sum(len(toric_degree_tuples(t.ambient)) for t in targets)==88

# Independent curve geometry / local (-1,-1) multiple-cover checks.
for index,beta,expected in ((3,(1,0),64),(5,(1,0),189),(5,(0,1),189),
                            (11,(1,0,0),9),(3,(2,0),8),(11,(2,0,0),QQ(9)/8)):
    for base in (4,5):assert toric_genus_zero_number(targets[index],beta,base)==expected
assert toric_known_ordinary_gw(targets[11],3,(2,0,0))==QQ(1)/84
assert toric_known_ordinary_gw(targets[5],3,(2,0)) is None
renamed=ToricCompleteIntersection(targets[3].ambient,targets[3].bundles,'renamed')
assert toric_known_ordinary_gw(renamed,3,(1,0))==QQ(2)/189

# S+F must include both the positive section and a broken S--F map.
graphs=toric_fixed_graphs(targets[8].ambient,3)
assert {sum(e for _,_,e in graph.edges) for graph in graphs[(1,1,0)]}=={1,2}
assert all(graph.mori_degree()==beta for beta,items in graphs.items() for graph in items)
assert any(any(e==2 and targets[8].ambient.edges[graph.vertex_colours[i],graph.vertex_colours[j]]['class']==(1,0,0)
                   for i,j,e in graph.edges) for graph in graphs[(2,0,0)])

# Recover the independent projective backend through the public coefficient API.
Y=ToricCompleteIntersection(ToricProductAmbient((3,)),((4,),))
for g,Ds in ((2,(0,1,2)),(3,(0,1))):
    new=ToricCompleteIntersectionCycles(Y,g)
    old=QuarticK3GraphCycles(g,(0,1,4,16))
    for D in Ds:
        for r in range(3*g-2):assert new.coefficients((D,),r)==old.coefficients(D,r),(g,D,r)

k3=ToricCompleteIntersectionCycles(targets[1],3)
assert k3.aggregate_relation((1,0),4)['compact_known']
assert k3.aggregate_relation((1,0),4)['compact_class']==0
cy=ToricCompleteIntersectionCycles(targets[3],3)
assert cy.aggregate_relation((0,1),6)['compact_known'] is False
try:cy.aggregate_relation((0,1),5)
except ValueError:pass
else:raise AssertionError('positive t powers are not automatic CJR relations')

# A product of elliptic curves has c1=0 but is not a K3; reject this extension
# until its CY-surface classification is explicitly implemented.
try:ToricCompleteIntersection(ToricProductAmbient(('F1',2)),((2,3,0),(0,0,3)))
except NotImplementedError:pass
else:raise AssertionError('c1=0 alone was incorrectly used to identify a K3 surface')
print('toric Chern, multidegree, curve-count, multiple-cover, projective-regression, and theory-boundary checks passed')

# Further independent parameter-space calculations in the results report.
for index,beta,expected in ((3,(0,1),640),(11,(0,1,0),54),(11,(0,0,1),168),(11,(0,2,0),QQ(27)/4)):
    for base in (4,5):assert toric_genus_zero_number(targets[index],beta,base)==expected
assert toric_primitive_curve_counts(targets[3])=={(1,0):64,(0,1):640}
assert toric_known_ordinary_gw(targets[11],3,(0,2,0))==QQ(1)/14
print('additional primitive-curve and pure-fiber multiple-cover checks passed')

try:ToricCompleteIntersection(ToricProductAmbient((1,2)),((0,0),))
except ValueError:pass
else:raise AssertionError('a general nonzero section of a trivial line bundle has empty zero locus')
print('nonempty-target input guard passed')
