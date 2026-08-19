r"""Exact universal CJR factors for a plane cubic hypersurface.

All classes on ``P^2`` are represented in the basis ``1,H,H^2`` and reduced
modulo ``H^3``.  Coefficients lie in the shared exact equivariant field from
``log_glsm_conventions.sage``.
"""

load("log_glsm_conventions.sage")


class P2Class(SageObject):
    """An equivariant coefficient vector ``a0+a1*H+a2*H^2``."""

    def __init__(self, rings, coefficients=(0, 0, 0)):
        if not isinstance(rings, PlaneCubicCoefficientRing):
            raise TypeError("rings must be a PlaneCubicCoefficientRing")
        if len(coefficients) != 3:
            raise ValueError("a P2 class needs three coefficients")
        self.rings = rings
        self.coefficients = tuple(rings.full(value) for value in coefficients)

    @classmethod
    def one(cls, rings):
        return cls(rings, (1, 0, 0))

    @classmethod
    def H(cls, rings):
        return cls(rings, (0, 1, 0))

    @classmethod
    def H2(cls, rings):
        return cls(rings, (0, 0, 1))

    def __getitem__(self, power):
        return self.coefficients[power]

    def __add__(self, other):
        self._require_same_ring(other)
        return P2Class(self.rings, tuple(a + b for a, b in zip(
            self.coefficients, other.coefficients
        )))

    def __sub__(self, other):
        self._require_same_ring(other)
        return P2Class(self.rings, tuple(a - b for a, b in zip(
            self.coefficients, other.coefficients
        )))

    def __neg__(self):
        return self.scale(-1)

    def __mul__(self, other):
        if not isinstance(other, P2Class):
            return self.scale(other)
        self._require_same_ring(other)
        answer = [self.rings.full(0)] * 3
        for left_power, left in enumerate(self.coefficients):
            for right_power, right in enumerate(other.coefficients):
                if left_power + right_power <= 2:
                    answer[left_power + right_power] += left * right
        return P2Class(self.rings, answer)

    def __rmul__(self, other):
        return self * other

    def scale(self, scalar):
        scalar = self.rings.full(scalar)
        return P2Class(self.rings, tuple(scalar * value for value in self.coefficients))

    def integrate(self):
        """Use ``int_P2 H^2=1``."""
        return self.coefficients[2]

    def pair(self, other):
        return (self * other).integrate()

    def specialize_base_weights(self, weights=(0, 1, 3)):
        return tuple(
            self.rings.specialize_base_weights(value, weights)
            for value in self.coefficients
        )

    def _require_same_ring(self, other):
        if not isinstance(other, P2Class) or other.rings is not self.rings:
            raise TypeError("P2 classes must use the same coefficient-ring instance")

    def __eq__(self, other):
        return (isinstance(other, P2Class)
                and self.rings is other.rings
                and self.coefficients == other.coefficients)

    def _repr_(self):
        return "P2Class(%s + (%s)H + (%s)H^2)" % self.coefficients


class PlaneCubicGraphFactors(SageObject):
    """CJR III, Section 9.3, specialized to ``(P2,O(3))``."""

    def __init__(self, rings=None):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.t = self.rings.t
        self._infinity_descendant_coefficients = {}

    def w(self):
        """Return ``t-H_infinity=t-3H`` in the hyperplane basis of ``P^2``."""
        return P2Class(self.rings, (self.t, -3, 0))

    def inverse_w_power(self, power):
        r"""Return ``(t-3H)^(-power)`` modulo ``H^3``."""
        power = ZZ(power)
        if power <= 0:
            raise ValueError("the inverse power must be positive")
        return P2Class(self.rings, (
            self.t ** (-power),
            3 * power * self.t ** (-power - 1),
            9 * binomial(power + 1, 2) * self.t ** (-power - 2),
        ))

    def edge(self, contact):
        r"""The edge factor ``c^c/(c! (t-3H)^c)``."""
        contact = ZZ(contact)
        if contact <= 0:
            raise ValueError("an edge contact must be positive")
        return self.inverse_w_power(contact).scale(
            QQ(contact ** contact) / factorial(contact)
        )

    def zero_nonspecial(self, contact):
        r"""The unstable nonspecial factor ``(t-3H)^2/c``."""
        contact = ZZ(contact)
        if contact <= 1:
            raise ValueError("a nonspecial tail requires contact greater than one")
        return P2Class(self.rings, (
            self.t ** 2 / contact,
            -6 * self.t / contact,
            QQ(9) / contact,
        ))

    def zero_marked(self, contact=None, psi_power=0):
        r"""Return the unstable marked-zero contribution.

        CJR III, equation (8.21), contributes the ordinary insertion times

        .. math::

            (t-H_\infty)\left(H_\infty/c-t\right)^b

        for a marking carrying ``psi^b`` on an edge of contact ``c``.  Here
        ``H_infinity=3H`` for the plane cubic.  The default arguments retain
        the primary factor ``t-3H`` used by older callers.
        """
        psi_power = ZZ(psi_power)
        if psi_power < 0:
            raise ValueError("a marked descendant power must be nonnegative")
        if contact is None:
            if psi_power:
                raise ValueError("a positive marked descendant needs a contact")
            return self.w()

        contact = ZZ(contact)
        if contact <= 0:
            raise ValueError("a marked edge contact must be positive")
        descendant_weight = P2Class(
            self.rings, (-self.t, QQ(3) / contact, 0)
        )
        answer = self.w()
        for _ in range(psi_power):
            answer *= descendant_weight
        return answer

    def zero_nodal(self, contact_left, contact_right):
        r"""The unstable nodal factor after identifying its evaluations.

        The CJR expression is

            w_left*w_right / (w_left/c_left + w_right/c_right).

        Both evaluations agree on the contracted unstable component, hence
        ``w_left=w_right=t-H``.
        """
        contact_left = ZZ(contact_left)
        contact_right = ZZ(contact_right)
        if contact_left <= 0 or contact_right <= 0:
            raise ValueError("nodal contacts must be positive")
        scale = QQ(contact_left * contact_right) / (contact_left + contact_right)
        return self.w().scale(scale)

    def stable_flag_descendant(self, contact, psi_power):
        r"""Coefficient of ``psi^k`` in ``1/((t-H)/c-psi)``."""
        contact = ZZ(contact)
        psi_power = ZZ(psi_power)
        if contact <= 0 or psi_power < 0:
            raise ValueError("contact must be positive and psi power nonnegative")
        return self.inverse_w_power(psi_power + 1).scale(
            contact ** (psi_power + 1)
        )

    def infinity_descendant_coefficient(self, psi_min_power):
        r"""Coefficient multiplying an effective invariant.

        CJR contributes ``t/(-t-psi_min)``.  Therefore the coefficient of
        ``psi_min^k`` is ``(-1)^(k+1)t^(-k)``.
        """
        psi_min_power = ZZ(psi_min_power)
        if psi_min_power < 0:
            raise ValueError("psi_min power must be nonnegative")
        if psi_min_power not in self._infinity_descendant_coefficients:
            self._infinity_descendant_coefficients[psi_min_power] = \
                self.rings.full(
                    (-1) ** (psi_min_power + 1)
                    * self.t ** (-psi_min_power)
                )
        return self._infinity_descendant_coefficients[psi_min_power]

    def nonspecial_edge_pair(self, contact):
        """Multiply an edge by its nonspecial unstable zero leaf."""
        return self.edge(contact) * self.zero_nonspecial(contact)


def p2_dual_power(power):
    """The basis dual to ``H^power`` under ``int_P2``."""
    power = ZZ(power)
    if power < 0 or power > 2:
        raise ValueError("a P2 basis power lies between zero and two")
    return 2 - power


class EffectivePolynomial(SageObject):
    r"""A finite exact sum of coefficients times effective-vertex monomials."""

    def __init__(self, coefficient_field, terms=None):
        self.coefficient_field = coefficient_field
        self.terms = {}
        for factors, coefficient in (terms or {}).items():
            self.add_term(coefficient, factors)

    def _monomial_key(self, factors):
        return tuple(sorted(tuple(factors), key=lambda item: item.order_key()))

    def add_term(self, coefficient, factors=()):
        coefficient = self.coefficient_field(coefficient)
        if not coefficient:
            return
        key = self._monomial_key(factors)
        self.terms[key] = self.terms.get(key, self.coefficient_field(0)) + coefficient
        if not self.terms[key]:
            del self.terms[key]

    def __add__(self, other):
        if not isinstance(other, EffectivePolynomial) \
                or other.coefficient_field is not self.coefficient_field:
            raise TypeError("effective polynomials must use the same field")
        answer = EffectivePolynomial(self.coefficient_field, self.terms)
        for factors, coefficient in other.terms.items():
            answer.add_term(coefficient, factors)
        return answer

    def __mul__(self, other):
        if not isinstance(other, EffectivePolynomial) \
                or other.coefficient_field is not self.coefficient_field:
            raise TypeError("effective polynomials must use the same field")
        answer = EffectivePolynomial(self.coefficient_field)
        for left_factors, left_coefficient in self.terms.items():
            for right_factors, right_coefficient in other.terms.items():
                answer.add_term(
                    left_coefficient * right_coefficient,
                    left_factors + right_factors,
                )
        return answer

    def scale(self, scalar):
        scalar = self.coefficient_field(scalar)
        return EffectivePolynomial(
            self.coefficient_field,
            {factors: scalar * coefficient
             for factors, coefficient in self.terms.items()},
        )

    def constant_term(self):
        return self.terms.get(tuple(), self.coefficient_field(0))

    def nonconstant_terms(self):
        return tuple(
            (coefficient, factors)
            for factors, coefficient in sorted(
                self.terms.items(), key=lambda item: tuple(
                    factor.order_key() for factor in item[0]
                )
            ) if factors
        )

    def _repr_(self):
        if not self.terms:
            return "0"
        return " + ".join(
            "(%s)*%s" % (
                coefficient,
                "*".join(repr(factor) for factor in factors) if factors else "1",
            )
            for factors, coefficient in self.terms.items()
        )
