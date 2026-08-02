r"""Exact intersection-number backends needed by twisted fixed loci.

Pure psi intersections are implemented in every genus by the DVV recursion,
and the lambda-g theorem supplies all integrals with a single top Hodge
class.  ``AdmcyclesHodgeIntegralBackend`` extends this to arbitrary products
of lambda and psi classes by evaluating them in the tautological ring.
"""

load("log_glsm_conventions.sage")

import os
import sys


def _load_vendored_admcycles():
    """Import admcycles, preferring the workspace-local vendored package."""
    candidates = (
        os.path.join(os.getcwd(), "vendor"),
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "vendor"),
    )
    for candidate in candidates:
        if os.path.isdir(os.path.join(candidate, "admcycles")) \
                and candidate not in sys.path:
            sys.path.insert(0, candidate)
    try:
        import admcycles
        return admcycles
    except ImportError as error:
        raise UnsupportedGeometryError(
            "admcycles is required for general Hodge monomials; install it "
            "or keep the vendored package in ./vendor"
        ) from error


def _odd_double_factorial(number):
    number = ZZ(number)
    if number <= -1:
        return ZZ.one()
    return prod(ZZ(k) for k in range(1, number + 1, 2))


class HodgeIntegralRequest(SageObject):
    """A hashable unsupported/supported Hodge integral specification."""

    def __init__(self, genus, psi_powers, lambda_indices=()):
        self.genus = ZZ(genus)
        self.psi_powers = tuple(ZZ(value) for value in psi_powers)
        self.lambda_indices = tuple(sorted(ZZ(value) for value in lambda_indices))
        if self.genus < 0 or any(value < 0 for value in self.psi_powers):
            raise ValueError("genus and psi powers must be nonnegative")
        if any(value < 0 or value > self.genus for value in self.lambda_indices):
            raise ValueError("lambda indices must lie between zero and genus")

    def signature(self):
        return self.genus, self.psi_powers, self.lambda_indices

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, HodgeIntegralRequest) and self.signature() == other.signature()

    def _repr_(self):
        return "HodgeIntegral(g=%s,psi=%s,lambda=%s)" % self.signature()


class PsiIntersectionBackend(SageObject):
    """Witten--Kontsevich intersections by exact DVV recursion."""

    def __init__(self):
        self._cache = {}

    def integral(self, genus, degrees):
        genus = ZZ(genus)
        degrees = tuple(sorted((ZZ(value) for value in degrees), reverse=True))
        key = genus, degrees
        if key in self._cache:
            return self._cache[key]
        value = self._integral_uncached(genus, degrees)
        self._cache[key] = QQ(value)
        return self._cache[key]

    def _integral_uncached(self, genus, degrees):
        marking_count = len(degrees)
        if genus < 0 or any(value < 0 for value in degrees):
            return QQ.zero()
        if 2 * genus - 2 + marking_count <= 0:
            return QQ.zero()
        if sum(degrees) != 3 * genus - 3 + marking_count:
            return QQ.zero()
        if genus == 0 and degrees == (0, 0, 0):
            return QQ.one()
        if genus == 1 and degrees == (1,):
            return QQ(1) / 24

        first = degrees[0]
        rest = degrees[1:]
        if first == 0:
            # String equation; this case is mostly the genus-zero base.
            return sum(
                self.integral(genus, rest[:j] + (rest[j] - 1,) + rest[j + 1:])
                for j in range(len(rest)) if rest[j] > 0
            )

        denominator = _odd_double_factorial(2 * first + 1)
        merge = QQ.zero()
        for j, degree in enumerate(rest):
            coefficient = QQ(
                _odd_double_factorial(2 * first + 2 * degree - 1)
            ) / _odd_double_factorial(2 * degree - 1)
            merged_degrees = rest[:j] + (first + degree - 1,) + rest[j + 1:]
            merge += coefficient * self.integral(genus, merged_degrees)

        boundary = QQ.zero()
        if first >= 2:
            for left_degree in range(first - 1):
                right_degree = first - 2 - left_degree
                coefficient = (
                    _odd_double_factorial(2 * left_degree + 1)
                    * _odd_double_factorial(2 * right_degree + 1)
                )
                nonseparating = self.integral(
                    genus - 1, (left_degree, right_degree) + rest
                )
                separating = QQ.zero()
                rest_count = len(rest)
                for mask in range(1 << rest_count):
                    left_rest = tuple(
                        rest[j] for j in range(rest_count) if mask & (1 << j)
                    )
                    right_rest = tuple(
                        rest[j] for j in range(rest_count) if not mask & (1 << j)
                    )
                    for left_genus in range(genus + 1):
                        separating += self.integral(
                            left_genus, (left_degree,) + left_rest
                        ) * self.integral(
                            genus - left_genus, (right_degree,) + right_rest
                        )
                boundary += coefficient * (nonseparating + separating)

        return (merge + boundary / 2) / denominator


class HodgeIntegralBackend(SageObject):
    """Pure-psi plus lambda-g exact intersection backend."""

    def __init__(self):
        self.psi = PsiIntersectionBackend()
        self._cache = {}

    @staticmethod
    def lambda_g_constant(genus):
        genus = ZZ(genus)
        if genus < 0:
            raise ValueError("genus must be nonnegative")
        if genus == 0:
            return QQ.one()
        return QQ(
            (2 ** (2 * genus - 1) - 1) * abs(bernoulli(2 * genus))
        ) / (2 ** (2 * genus - 1) * factorial(2 * genus))

    def integral(self, request_or_genus, psi_powers=None, lambda_indices=()):
        if isinstance(request_or_genus, HodgeIntegralRequest):
            request = request_or_genus
        else:
            request = HodgeIntegralRequest(
                request_or_genus, psi_powers or (), lambda_indices
            )
        if request in self._cache:
            return self._cache[request]

        nonzero_lambda = tuple(index for index in request.lambda_indices if index)
        if not nonzero_lambda:
            value = self.psi.integral(request.genus, request.psi_powers)
        elif nonzero_lambda == (request.genus,):
            expected = 2 * request.genus - 3 + len(request.psi_powers)
            if sum(request.psi_powers) != expected:
                value = QQ.zero()
            else:
                multinomial = factorial(expected) / prod(
                    factorial(value) for value in request.psi_powers
                )
                value = QQ(multinomial) * self.lambda_g_constant(request.genus)
        else:
            raise UnsupportedGeometryError(
                "general Hodge monomial %s is not implemented; "
                "pure psi and one lambda_g factor are exact" % (request,)
            )
        self._cache[request] = QQ(value)
        return self._cache[request]


class AdmcyclesHodgeIntegralBackend(HodgeIntegralBackend):
    r"""All tautological Hodge/psi monomials via admcycles ``evaluate()``."""

    def __init__(self):
        super().__init__()
        self.admcycles = _load_vendored_admcycles()

    def integral(self, request_or_genus, psi_powers=None, lambda_indices=()):
        if isinstance(request_or_genus, HodgeIntegralRequest):
            request = request_or_genus
        else:
            request = HodgeIntegralRequest(
                request_or_genus, psi_powers or (), lambda_indices
            )
        if request in self._cache:
            return self._cache[request]

        genus = request.genus
        markings = len(request.psi_powers)
        if 2 * genus - 2 + markings <= 0:
            value = QQ.zero()
        else:
            codimension = sum(request.psi_powers) + sum(request.lambda_indices)
            dimension = 3 * genus - 3 + markings
            if codimension != dimension:
                value = QQ.zero()
            else:
                tautological_class = self.admcycles.fundclass(genus, markings)
                for marking, power in enumerate(request.psi_powers, start=1):
                    if power:
                        tautological_class *= self.admcycles.psiclass(
                            marking, genus, markings
                        ) ** int(power)
                for index in request.lambda_indices:
                    if index:
                        tautological_class *= self.admcycles.lambdaclass(
                            int(index), genus, markings
                        )
                value = QQ(tautological_class.evaluate())
        self._cache[request] = value
        return value
