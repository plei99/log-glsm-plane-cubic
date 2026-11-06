"""Independent closed Hodge values for the dual Euler O(3) evaluator."""

from math import factorial
from log_glsm.core.conventions import (
    PlaneCubicCoefficientRing,
)
from log_glsm.localization.plane import (
    P2FixedLocusEvaluator,
)
from log_glsm.localization.twisted_plane import (
    TwistedInsertion,
    TwistedZeroVertexRequest,
)

rings = PlaneCubicCoefficientRing()
backend = P2FixedLocusEvaluator(rings, lift_strategy="standard")
lambdas = rings.base_weights()
t = rings.t

# With psi^(3g-2) on Mbar_(g,1), only the degree-zero Hodge term survives.
# The one-point psi theorem gives 1/(24^g*g!), independently of the
# characteristic-polynomial construction and the inverse-Hodge shortcut.
for g in range(1, 5):
    request = TwistedZeroVertexRequest(g, 0, (TwistedInsertion(2, 3*g-2),))
    expected = sum(
        lambdas[i]^2
        * prod(lambdas[i]-lambdas[j] for j in range(3) if j != i)^(g-1)
        * (t-3*lambdas[i])^(1-g)
        for i in range(3)
    ) / (24^g * factorial(g))
    assert backend.evaluate(request) == expected

# Mapping-to-a-point obstruction e(E^vee tensor TP2), multiplied by
# (t-3H)/(t-3H+lambda_1), and integral lambda_1=1/24.
for h_power, expected in ((0, -1/(8*t)), (1, -QQ(1)/8), (2, 0)):
    request = TwistedZeroVertexRequest(1, 0, (TwistedInsertion(h_power),))
    assert rings.nonequivariant_base_limit(backend.evaluate(request)) == expected

# Genus-zero twist has virtual rank one and Euler class t-3H.
request = TwistedZeroVertexRequest(0, 0, (
    TwistedInsertion(1), TwistedInsertion(0), TwistedInsertion(0)))
assert rings.nonequivariant_base_limit(backend.evaluate(request)) == -3

print("all independent O(3) Hodge specialization tests passed")
