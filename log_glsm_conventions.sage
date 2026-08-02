r"""Shared conventions and finite truncation for plane-cubic log GLSM.

The module is intentionally independent of the graph enumerator and of the
dynamic-programming solver.  It supplies immutable probe data, one exact
equivariant coefficient-ring factory, the plane-cubic degree conversion, and
dimension-driven bounds for descendant/Laurent expansions.
"""


class UnsupportedGeometryError(NotImplementedError):
    """A requested invariant lies outside an exact geometric backend."""


class TruncationError(ValueError):
    """A requested series truncation is not provably sufficient."""


class EllipticInsertion(SageObject):
    r"""An ordinary marking insertion on the cubic with a psi power.

    The point class on the cubic is represented after push-forward by
    ``H/3`` on ``P^2``, since ``int_E H=3``.  Keeping this conversion here
    avoids the common but incorrect use of the ambient point class ``H^2``.
    """

    _KINDS = ("unit", "point")

    def __init__(self, kind, psi_power=0):
        self.kind = str(kind)
        self.psi_power = ZZ(psi_power)
        if self.kind not in self._KINDS:
            raise ValueError("elliptic insertions are 'unit' or 'point'")
        if self.psi_power < 0:
            raise ValueError("a psi power must be nonnegative")

    @property
    def elliptic_codimension(self):
        return ZZ(1 if self.kind == "point" else 0)

    @property
    def total_codimension(self):
        return self.elliptic_codimension + self.psi_power

    def ambient_class(self):
        """Return ``(coefficient, H_power)`` for integration after push-forward."""
        if self.kind == "unit":
            return QQ.one(), ZZ.zero()
        return QQ(1) / 3, ZZ.one()

    def signature(self):
        return self.kind, self.psi_power

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, EllipticInsertion) and self.signature() == other.signature()

    def _repr_(self):
        if self.psi_power:
            return "%s*psi^%s" % (self.kind, self.psi_power)
        return self.kind


class ProbeSpec(SageObject):
    r"""Discrete data for one connected compact-type elliptic GW probe."""

    def __init__(self, genus, ambient_degree, insertions=(), connected=True,
                 label=""):
        self.genus = ZZ(genus)
        self.ambient_degree = ZZ(ambient_degree)
        self.insertions = tuple(
            insertion if isinstance(insertion, EllipticInsertion)
            else EllipticInsertion(*insertion)
            for insertion in insertions
        )
        self.connected = bool(connected)
        self.label = str(label)
        if self.genus < 0 or self.ambient_degree < 0:
            raise ValueError("genus and ambient degree must be nonnegative")

    @classmethod
    def stationary(cls, genus, ambient_degree, descendants, label=""):
        return cls(
            genus, ambient_degree,
            tuple(EllipticInsertion("point", power) for power in descendants),
            connected=True, label=label,
        )

    @property
    def marking_count(self):
        return ZZ(len(self.insertions))

    @property
    def intrinsic_degree(self):
        if self.ambient_degree % 3:
            raise ValueError(
                "a curve on the plane cubic has ambient degree divisible by three"
            )
        return self.ambient_degree // 3

    @property
    def virtual_dimension(self):
        # Complex virtual dimension of stable maps to a Calabi--Yau curve.
        return 2 * self.genus - 2 + self.marking_count

    @property
    def insertion_codimension(self):
        return ZZ(sum(item.total_codimension for item in self.insertions))

    def dimension_defect(self):
        return self.virtual_dimension - self.insertion_codimension

    def is_dimension_zero(self):
        return self.dimension_defect() == 0

    def stationary_descendants(self):
        if any(item.kind != "point" for item in self.insertions):
            raise UnsupportedGeometryError(
                "Bloch--Okounkov input currently supports stationary point insertions"
            )
        return tuple(item.psi_power for item in self.insertions)

    def signature(self):
        return (
            self.genus, self.ambient_degree,
            tuple(item.signature() for item in self.insertions),
            self.connected, self.label,
        )

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, ProbeSpec) and self.signature() == other.signature()

    def _repr_(self):
        return "ProbeSpec(g=%s,d=%s,insertions=%s,connected=%s%s)" % (
            self.genus, self.ambient_degree, self.insertions, self.connected,
            ",label=%r" % self.label if self.label else "",
        )


class PlaneCubicCoefficientRing(SageObject):
    r"""The unique exact equivariant ring used by localization components."""

    def __init__(self, laurent_precision=12):
        self.laurent_precision = ZZ(laurent_precision)
        if self.laurent_precision < 1:
            raise ValueError("Laurent precision must be positive")
        self.base_polynomial = PolynomialRing(QQ, names=("lambda0", "lambda1", "lambda2"))
        self.lambda0, self.lambda1, self.lambda2 = self.base_polynomial.gens()
        self.base_field = self.base_polynomial.fraction_field()
        self.full_polynomial = PolynomialRing(
            QQ, names=("lambda0", "lambda1", "lambda2", "t")
        )
        self.full_field = self.full_polynomial.fraction_field()
        self.full_lambda0, self.full_lambda1, self.full_lambda2, self.t = \
            self.full_polynomial.gens()
        self.laurent_ring = LaurentSeriesRing(
            self.base_field, "t", default_prec=self.laurent_precision
        )
        self.laurent_t = self.laurent_ring.gen()

    def full(self, value):
        return self.full_field(value)

    def base_weights(self):
        return self.full_lambda0, self.full_lambda1, self.full_lambda2

    def to_laurent(self, value, precision=None):
        """Expand an exact rational function at ``t=0``."""
        precision = self.laurent_precision if precision is None else ZZ(precision)
        if precision < 1:
            raise ValueError("Laurent precision must be positive")
        value = self.full_field(value)
        numerator = self.full_polynomial(value.numerator())
        denominator = self.full_polynomial(value.denominator())
        base_lambdas = self.base_field.gens()
        arguments = base_lambdas + (self.laurent_t,)
        series = self.laurent_ring(numerator(*arguments)) / self.laurent_ring(
            denominator(*arguments)
        )
        return series.add_bigoh(precision)

    def specialize_base_weights(self, value, weights=(0, 1, 3)):
        """Specialize only the auxiliary ``P^2`` torus weights."""
        if len(weights) != 3 or len(set(weights)) != 3:
            raise ValueError("base weights must be three distinct values")
        value = self.full_field(value)
        return value(
            QQ(weights[0]), QQ(weights[1]), QQ(weights[2]), self.t
        )


class PlaneCubicDimension(SageObject):
    """Exact dimension and denominator-expansion bounds."""

    @staticmethod
    def elliptic_virtual_dimension(genus, markings):
        return ZZ(2 * ZZ(genus) - 2 + ZZ(markings))

    @staticmethod
    def plane_virtual_dimension(genus, degree, markings):
        # vdim Mbar_(g,n)(P2,d) = g-1+3d+n.
        return ZZ(genus) - 1 + 3 * ZZ(degree) + ZZ(markings)

    @staticmethod
    def o3_euler_virtual_rank(genus, degree):
        # chi(f^* O(3)) = 3d+1-g.
        return 3 * ZZ(degree) + 1 - ZZ(genus)

    @classmethod
    def twisted_zero_virtual_dimension(cls, genus, degree, markings):
        return (
            cls.plane_virtual_dimension(genus, degree, markings)
            - cls.o3_euler_virtual_rank(genus, degree)
        )

    @staticmethod
    def infinity_balance_defect(genus, degree, contacts):
        contacts = tuple(ZZ(c) for c in contacts)
        if any(c >= 0 for c in contacts):
            raise ValueError("infinity contacts must be negative")
        return 3 * ZZ(degree) - (2 * ZZ(genus) - 2) \
            - sum(c + 1 for c in contacts)

    @staticmethod
    def denominator_power_bound(complex_dimension, inserted_codimension=0):
        r"""Maximum psi power needed from ``1/(w-psi)``.

        A term with psi power above the remaining complex dimension vanishes.
        """
        remaining = ZZ(complex_dimension) - ZZ(inserted_codimension)
        return max(ZZ(-1), remaining)

    @staticmethod
    def psi_min_power_bound(genus, valence, insertion_codimension=0):
        # The effective moduli class may lower this bound, but never raises it
        # above the dimension of Mbar_(g,valence).
        dimension = 3 * ZZ(genus) - 3 + ZZ(valence)
        return PlaneCubicDimension.denominator_power_bound(
            dimension, insertion_codimension
        )

    @staticmethod
    def require_finite_bound(bound, context="expansion"):
        bound = ZZ(bound)
        if bound < 0:
            raise TruncationError("%s vanishes by dimension" % context)
        return bound


def plane_cubic_reduced_sign(genus, ambient_degree):
    """The sign ``(-1)^(1-g+3d)`` in the reduced-class comparison."""
    return ZZ(-1) ** (1 - ZZ(genus) + 3 * ZZ(ambient_degree))

