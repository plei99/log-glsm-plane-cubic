r"""core / tautological_polynomial."""

from sage.all import *



__all__ = ['TautologicalPolynomial']


class TautologicalPolynomial(SageObject):
    r"""A truncated polynomial in vertex psi and lambda classes.

    A term key is ``(psi_powers, lambda_indices)``.  Coefficients lie in the
    full equivariant fraction field, while tautological codimension is kept
    separately so rational functions in torus weights never reach admcycles.
    """

    def __init__(self, field, marking_count, dimension, terms=None):
        self.field = field
        self.marking_count = ZZ(marking_count)
        self.dimension = ZZ(dimension)
        self.terms = {}
        for key, coefficient in (terms or {}).items():
            psi_powers, lambda_indices = key
            psi_powers = tuple(ZZ(value) for value in psi_powers)
            lambda_indices = tuple(sorted(
                ZZ(value) for value in lambda_indices if value
            ))
            if len(psi_powers) != self.marking_count:
                raise ValueError("a tautological monomial has the wrong arity")
            if sum(psi_powers) + sum(lambda_indices) > self.dimension:
                continue
            coefficient = self.field(coefficient)
            if coefficient:
                self.terms[(psi_powers, lambda_indices)] = coefficient

    @classmethod
    def one(cls, field, marking_count, dimension):
        return cls(field, marking_count, dimension, {
            ((ZZ.zero(),) * marking_count, ()): field.one(),
        })

    @classmethod
    def monomial(cls, field, marking_count, dimension,
                 psi_powers=None, lambda_indices=(), coefficient=1):
        psi_powers = ((ZZ.zero(),) * marking_count if psi_powers is None
                      else tuple(psi_powers))
        return cls(field, marking_count, dimension, {
            (psi_powers, tuple(lambda_indices)): coefficient,
        })

    def _compatible(self, other):
        return (isinstance(other, TautologicalPolynomial)
                and self.field == other.field
                and self.marking_count == other.marking_count
                and self.dimension == other.dimension)

    def __add__(self, other):
        if not self._compatible(other):
            raise TypeError("incompatible tautological polynomials")
        terms = dict(self.terms)
        for key, coefficient in other.terms.items():
            terms[key] = terms.get(key, self.field.zero()) + coefficient
            if not terms[key]:
                del terms[key]
        return TautologicalPolynomial(
            self.field, self.marking_count, self.dimension, terms
        )

    def __neg__(self):
        return self.scaled(-1)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not self._compatible(other):
            raise TypeError("incompatible tautological polynomials")
        terms = {}
        for (left_psi, left_lambda), left_coefficient in self.terms.items():
            for (right_psi, right_lambda), right_coefficient in other.terms.items():
                psi_powers = tuple(
                    left + right for left, right in zip(left_psi, right_psi)
                )
                lambda_indices = tuple(sorted(left_lambda + right_lambda))
                if sum(psi_powers) + sum(lambda_indices) > self.dimension:
                    continue
                key = psi_powers, lambda_indices
                terms[key] = (terms.get(key, self.field.zero())
                              + left_coefficient * right_coefficient)
        return TautologicalPolynomial(
            self.field, self.marking_count, self.dimension, terms
        )

    def __pow__(self, exponent):
        exponent = ZZ(exponent)
        if exponent < 0:
            raise ValueError("tautological powers must be nonnegative")
        answer = TautologicalPolynomial.one(
            self.field, self.marking_count, self.dimension
        )
        factor = self
        while exponent:
            if exponent % 2:
                answer = answer * factor
            exponent //= 2
            if exponent:
                factor = factor * factor
        return answer

    def scaled(self, scalar):
        scalar = self.field(scalar)
        return TautologicalPolynomial(
            self.field, self.marking_count, self.dimension,
            {key: scalar * coefficient
             for key, coefficient in self.terms.items()},
        )

    def integrate(self, genus, hodge_backend):
        answer = self.field.zero()
        for (psi_powers, lambda_indices), coefficient in self.terms.items():
            answer += coefficient * hodge_backend.integral(
                genus, psi_powers, lambda_indices
            )
        return self.field(answer)
