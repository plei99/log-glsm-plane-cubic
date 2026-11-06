"""Tests for universal plane-cubic CJR contribution factors."""


from log_glsm.core.conventions import (
    PlaneCubicCoefficientRing,
)
from log_glsm.infinity.factors import (
    EffectivePolynomial,
    P2Class,
    PlaneCubicGraphFactors,
    p2_dual_power,
)
from log_glsm.infinity.solver import (
    EffectiveVertex,
)



def run_tests():
    rings = PlaneCubicCoefficientRing(10)
    factors = PlaneCubicGraphFactors(rings)
    t = rings.t

    assert factors.w().coefficients == (t, -3, 0)
    assert factors.edge(1).coefficients == (1/t, 3/t^2, 9/t^3)
    assert factors.edge(2).coefficients == (2/t^2, 12/t^3, 54/t^4)
    assert factors.zero_nonspecial(2).coefficients == (
        t^2/2, -3*t, QQ(9)/2
    )
    assert factors.nonspecial_edge_pair(2) == P2Class.one(rings)
    assert factors.zero_marked() == factors.w()
    assert factors.zero_marked(2, 0) == factors.w()
    # CJR III (8.21): the entire cotangent weight is divided by c.
    # For c=2: -w^2/2 and w^3/4, with w=t-3H and H^3=0.
    assert factors.zero_marked(2, 1).coefficients == (
        -t^2/2, 3*t, -QQ(9)/2
    )
    assert factors.zero_marked(2, 2).coefficients == (
        t^3/4, -QQ(9)*t^2/4, QQ(27)*t/4
    )
    assert factors.zero_nodal(2, 2) == factors.w()
    assert factors.stable_flag_descendant(2, 0).coefficients == (
        2/t, 6/t^2, 18/t^3
    )
    assert factors.stable_flag_descendant(2, 1).coefficients == (
        4/t^2, 24/t^3, 108/t^4
    )
    assert [factors.infinity_descendant_coefficient(k) for k in range(4)] == [
        -1, 1/t, -1/t^2, 1/t^3
    ]
    assert [p2_dual_power(k) for k in range(3)] == [2, 1, 0]

    v1 = EffectiveVertex(1, 0, (-1,))
    v2 = EffectiveVertex(2, 0, (-3,))
    left = EffectivePolynomial(rings.full_field, {(v1,): 2, tuple(): 3})
    right = EffectivePolynomial(rings.full_field, {(v2,): t})
    product_value = left * right
    assert product_value.terms[(v1, v2)] == 2*t
    assert product_value.terms[(v2,)] == 3*t
    assert product_value.constant_term() == 0


run_tests()
print("all universal CJR graph-factor tests passed")
