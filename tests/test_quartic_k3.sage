"""Independent constant-map checks and ordinary/reduced theory boundaries."""

from log_glsm.cycles.quartic import (
    P3FixedLocusGraphEnumerator,
    QuarticK3GraphCycles,
)
from log_glsm.localization.plane_graphs import (
    P2FixedLocusGraphEnumerator,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
    QuarticK3Theory,
)

for genus,indices,co,integral in (
        (2,(1,2),-20,-QQ(1)/288),
        (3,(1,2,3),20,QQ(1)/72576)):
    backend=QuarticK3Cycles(genus)
    for r in range(backend.dimension):
        assert not backend.coefficients(0,r).get(0,{})
    assert backend.coefficients(0,backend.dimension).get(0,{}) == {
        ('interior',None,indices):QQ(co)}
    assert backend.top_integral(0)==integral
    assert QuarticK3Theory.compact_virtual_class(genus,0)==0
    assert QuarticK3Theory.compact_virtual_class(genus,7)==0
    assert QuarticK3Theory.t_power(genus,backend.dimension)==1-genus
    try:
        backend.aggregate_relation(0,2*genus-3)
    except ValueError:
        pass
    else:
        raise AssertionError('positive equivariant powers must not be labelled vanishing relations')

try:
    QuarticK3Theory.compact_virtual_class(3,1,reduced=True)
except NotImplementedError:
    pass
else:
    raise AssertionError('reduced K3 invariants were incorrectly declared zero')

assert QuarticK3Theory.balanced(3,1,(-1,-1))
assert not QuarticK3Theory.balanced(3,1,(-2,))
assert QuarticK3Theory.balanced(3,0,(-5,))

# Equal genus bridge factors and multiple-cover deck factors must agree with
# the direct, independently assembled degree-one tail/bridge expression.
direct=QuarticK3Cycles(2)
general=QuarticK3GraphCycles(2)
assert direct.top_integral(1)==general.top_integral(1)
assert general.top_integral(2)==QuarticK3GraphCycles(2,(0,1,4,21)).top_integral(2)
assert P3FixedLocusGraphEnumerator(0,0,1).graph_count()==6
assert P2FixedLocusGraphEnumerator(0,0,1).graph_count()==3

print('quartic K3 constant-map, theory, balance, and graph-normalization tests passed')
