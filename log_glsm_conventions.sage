r"""Shared conventions and finite truncation for plane-cubic log GLSM.

The module is intentionally independent of the graph enumerator and of the
dynamic-programming solver.  It supplies immutable probe data, one exact
equivariant coefficient-ring factory, the plane-cubic degree conversion, and
dimension-driven bounds for descendant/Laurent expansions.
"""

import os


NONEQUIVARIANT_BASE_PATH = "lambda=(0,epsilon,epsilon^2)"
NONEQUIVARIANT_CACHE_TAG = "nonequivariant-quadratic"


class UnsupportedGeometryError(NotImplementedError):
    """A requested invariant lies outside an exact geometric backend."""


class TruncationError(ValueError):
    """A requested series truncation is not provably sufficient."""


def specialized_zero_vertex_cache_path(path, specialization):
    """Return a convention-tagged cache path without double-tagging it."""
    if path is None or specialization is None:
        return path
    root, extension = os.path.splitext(str(path))
    if specialization == "nonequivariant":
        # The cache tag is part of the mathematical convention: individual
        # fixed loci are evaluated along the nonresonant path recorded above
        # before the auxiliary parameter tends to zero.
        tag = NONEQUIVARIANT_CACHE_TAG
    else:
        tag = "-".join(
            str(QQ(value)).replace("/", "_over_")
            for value in specialization
        )
    suffix = ".weights-%s" % tag
    if root.endswith(suffix):
        return str(path)
    return "%s%s%s" % (root, suffix, extension)


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

    def to_record(self):
        """Return a JSON-friendly exact checkpoint record."""
        return {
            "kind": self.kind,
            "psi_power": int(self.psi_power),
        }

    @classmethod
    def from_record(cls, record):
        return cls(record["kind"], record.get("psi_power", 0))

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, EllipticInsertion) and self.signature() == other.signature()

    def _repr_(self):
        if self.psi_power:
            return "%s*psi^%s" % (self.kind, self.psi_power)
        return self.kind


class ProbeSpec(SageObject):
    r"""Discrete data for one connected compact-type probe.

    ``psi_convention`` records where ordinary marking cotangent lines live.
    ``"log"`` means the cotangent line on the coarse universal log-GLSM
    curve, which is the class occurring in CJR III, Theorem 8.4 and equation
    (8.21).  ``"stabilized"`` means the stable-map cotangent line after the
    morphism ``st``; Bloch--Okounkov computes this latter GW descendant.

    The two classes differ on a fixed locus when stabilization contracts a
    marked zero tail.  Primary insertions do not see the distinction.
    """

    _PSI_CONVENTIONS = ("log", "stabilized")

    def __init__(self, genus, ambient_degree, insertions=(), connected=True,
                 label="", psi_convention="stabilized"):
        self.genus = ZZ(genus)
        self.ambient_degree = ZZ(ambient_degree)
        self.insertions = tuple(
            insertion if isinstance(insertion, EllipticInsertion)
            else EllipticInsertion(*insertion)
            for insertion in insertions
        )
        self.connected = bool(connected)
        self.label = str(label)
        self.psi_convention = str(psi_convention)
        if self.genus < 0 or self.ambient_degree < 0:
            raise ValueError("genus and ambient degree must be nonnegative")
        if self.psi_convention not in self._PSI_CONVENTIONS:
            raise ValueError(
                "psi_convention must be 'log' or 'stabilized'"
            )

    @classmethod
    def stationary(cls, genus, ambient_degree, descendants, label=""):
        return cls(
            genus, ambient_degree,
            tuple(EllipticInsertion("point", power) for power in descendants),
            connected=True, label=label, psi_convention="stabilized",
        )

    @property
    def has_descendants(self):
        return any(insertion.psi_power for insertion in self.insertions)

    def with_psi_convention(self, psi_convention):
        """Return the same discrete probe with an explicit psi convention."""
        return ProbeSpec(
            self.genus,
            self.ambient_degree,
            self.insertions,
            connected=self.connected,
            label=self.label,
            psi_convention=psi_convention,
        )

    @property
    def marking_count(self):
        return ZZ(len(self.insertions))

    def is_stable_map_type(self):
        r"""Whether the underlying connected stable-map data are stable.

        Positive ambient degree stabilizes the map.  In degree zero this is
        the usual pointed-curve condition ``2g-2+n>0``.  Keeping this check
        on the shared probe object prevents the exceptional unmarked
        constant genus-one type from reaching either localization or the
        Bloch--Okounkov wrapper.
        """
        return (
            self.ambient_degree > 0
            or 2 * self.genus - 2 + self.marking_count > 0
        )

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
            self.connected, self.label, self.psi_convention,
        )

    def to_record(self):
        """Return a JSON-friendly exact checkpoint record."""
        return {
            "genus": int(self.genus),
            "ambient_degree": int(self.ambient_degree),
            "insertions": [item.to_record() for item in self.insertions],
            "connected": bool(self.connected),
            "label": self.label,
            "psi_convention": self.psi_convention,
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["genus"],
            record["ambient_degree"],
            tuple(EllipticInsertion.from_record(item)
                  for item in record.get("insertions", ())),
            connected=record.get("connected", True),
            label=record.get("label", ""),
            psi_convention=record.get("psi_convention", "stabilized"),
        )

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, ProbeSpec) and self.signature() == other.signature()

    def _repr_(self):
        return "ProbeSpec(g=%s,d=%s,insertions=%s,connected=%s,psi=%s%s)" % (
            self.genus, self.ambient_degree, self.insertions, self.connected,
            self.psi_convention,
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
        # Generic Laurent-series division over the three-variable rational
        # function field repeatedly normalises very large expressions.  For
        # coefficient extraction, regard t as the sole polynomial variable
        # and extend the quotient recursively instead.
        self.t_polynomial = PolynomialRing(self.base_field, "t_coefficient")
        self.t_coefficient = self.t_polynomial.gen()
        self._coefficient_expansions = {}
        self._infinity_coefficient_expansions = {}
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

    def _coefficient_expansion(self, value):
        """Return cached data for the exact Laurent recurrence at ``t=0``."""
        value = self.full_field(value)
        cached = self._coefficient_expansions.get(value)
        if cached is not None:
            return cached

        if not value:
            cached = {
                "zero": True,
                "shift": ZZ(0),
                "numerator": self.t_polynomial.zero(),
                "denominator": self.t_polynomial.one(),
                "coefficients": [],
            }
            self._coefficient_expansions[value] = cached
            return cached

        numerator = self.full_polynomial(value.numerator())
        denominator = self.full_polynomial(value.denominator())
        arguments = self.base_field.gens() + (self.t_coefficient,)
        numerator_t = self.t_polynomial(numerator(*arguments))
        denominator_t = self.t_polynomial(denominator(*arguments))
        numerator_valuation = ZZ(numerator_t.valuation())
        denominator_valuation = ZZ(denominator_t.valuation())
        cached = {
            "zero": False,
            "shift": numerator_valuation - denominator_valuation,
            "numerator": numerator_t,
            "denominator": denominator_t,
            "numerator_valuation": numerator_valuation,
            "denominator_valuation": denominator_valuation,
            "coefficients": [],
        }
        self._coefficient_expansions[value] = cached
        return cached

    def laurent_coefficient(self, value, power=0):
        r"""Return ``[t^power] value`` by an exact univariate recurrence.

        This has no precision cutoff.  It removes the t-adic valuations of
        numerator and denominator, then solves the quotient identity one
        coefficient at a time.  Previously requested coefficients are kept.
        """
        power = ZZ(power)
        expansion = self._coefficient_expansion(value)
        if expansion["zero"]:
            return self.base_field.zero()

        target = power - expansion["shift"]
        if target < 0:
            return self.base_field.zero()

        numerator = expansion["numerator"]
        denominator = expansion["denominator"]
        numerator_valuation = expansion["numerator_valuation"]
        denominator_valuation = expansion["denominator_valuation"]
        denominator_constant = denominator[denominator_valuation]
        coefficients = expansion["coefficients"]

        while len(coefficients) <= target:
            degree = ZZ(len(coefficients))
            right_hand_side = numerator[numerator_valuation + degree]
            upper = min(degree, denominator.degree() - denominator_valuation)
            for offset in range(1, upper + 1):
                right_hand_side -= denominator[denominator_valuation + offset] * \
                    coefficients[degree - offset]
            coefficients.append(self.base_field(right_hand_side / denominator_constant))
        return coefficients[target]

    def _infinity_coefficient_expansion(self, value):
        r"""Return cached recurrence data after substituting ``u=1/t``."""
        value = self.full_field(value)
        cached = self._infinity_coefficient_expansions.get(value)
        if cached is not None:
            return cached

        if not value:
            cached = {
                "zero": True,
                "shift": ZZ(0),
                "numerator": self.t_polynomial.zero(),
                "denominator": self.t_polynomial.one(),
                "numerator_degree": ZZ(0),
                "denominator_degree": ZZ(0),
                "coefficients": [],
            }
            self._infinity_coefficient_expansions[value] = cached
            return cached

        numerator = self.full_polynomial(value.numerator())
        denominator = self.full_polynomial(value.denominator())
        arguments = self.base_field.gens() + (self.t_coefficient,)
        numerator_t = self.t_polynomial(numerator(*arguments))
        denominator_t = self.t_polynomial(denominator(*arguments))
        numerator_degree = ZZ(numerator_t.degree())
        denominator_degree = ZZ(denominator_t.degree())
        cached = {
            "zero": False,
            # P(1/u)/Q(1/u) starts with u^(deg(Q)-deg(P)).
            "shift": denominator_degree - numerator_degree,
            "numerator": numerator_t,
            "denominator": denominator_t,
            "numerator_degree": numerator_degree,
            "denominator_degree": denominator_degree,
            "coefficients": [],
        }
        self._infinity_coefficient_expansions[value] = cached
        return cached

    def laurent_coefficient_at_infinity(self, value, power=0):
        r"""Return ``[t^power] value`` in the expansion at ``t=infinity``.

        Equivalently, substitute ``u=1/t`` and extract ``[u^(-power)]`` at
        ``u=0``.  CJR R-torus localization uses this expansion: edge and
        infinity-node denominators are formal series in ``t^(-1)``.
        """
        power = ZZ(power)
        expansion = self._infinity_coefficient_expansion(value)
        if expansion["zero"]:
            return self.base_field.zero()

        target = -power - expansion["shift"]
        if target < 0:
            return self.base_field.zero()

        numerator = expansion["numerator"]
        denominator = expansion["denominator"]
        numerator_degree = expansion["numerator_degree"]
        denominator_degree = expansion["denominator_degree"]
        denominator_constant = denominator[denominator_degree]
        coefficients = expansion["coefficients"]

        while len(coefficients) <= target:
            degree = ZZ(len(coefficients))
            numerator_index = numerator_degree - degree
            right_hand_side = (
                numerator[numerator_index]
                if numerator_index >= 0 else self.base_field.zero()
            )
            upper = min(degree, denominator_degree)
            for offset in range(1, upper + 1):
                right_hand_side -= denominator[denominator_degree - offset] * \
                    coefficients[degree - offset]
            coefficients.append(self.base_field(
                right_hand_side / denominator_constant
            ))
        return coefficients[target]

    def specialize_base_weights(self, value, weights=(0, 1, 3)):
        """Specialize only the auxiliary ``P^2`` torus weights."""
        if len(weights) != 3 or len(set(weights)) != 3:
            raise ValueError("base weights must be three distinct values")
        value = self.full_field(value)
        return value(
            QQ(weights[0]), QQ(weights[1]), QQ(weights[2]), self.t
        )

    def nonequivariant_base_limit(self, value):
        r"""Set all auxiliary ``P^2`` weights to zero after fixed-locus summing.

        The input must already be the complete auxiliary-torus localization
        sum for one twisted invariant.  Individual fixed-locus summands have
        poles at this specialization, while their normalized sum has a
        removable limit.
        """
        value = self.full_field(value)
        return self.full_field(value(0, 0, 0, self.t))


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

    @classmethod
    def infinity_reduced_virtual_dimension(cls, genus, degree, contacts):
        r"""Complex reduced virtual dimension of a plane-cubic infinity vertex.

        CJR II, (9.9), specializes for ``X=P^2`` and ``E=O(3)`` to

            vdim_red = vdim Mbar_(g,n)(Z,beta_X) + sum_i(c_i+1),

        where ``Z`` is the elliptic cubic.  For balanced data this is the
        particularly simple value ``n + 3*degree``.
        """
        genus = ZZ(genus)
        degree = ZZ(degree)
        contacts = tuple(ZZ(c) for c in contacts)
        if genus < 0 or degree < 0:
            raise ValueError("genus and infinity degree must be nonnegative")
        if any(c >= 0 for c in contacts):
            raise ValueError("infinity contacts must be negative")
        return cls.elliptic_virtual_dimension(genus, len(contacts)) \
            + sum(c + 1 for c in contacts)

    @classmethod
    def infinity_required_psi_min_power(
            cls, genus, degree, contacts, insertion_codimension=0):
        r"""Return the unique ``psi_min`` power that produces a number.

        A negative result means that the evaluation insertions already
        exceed the reduced virtual dimension, so the invariant vanishes.
        """
        return cls.infinity_reduced_virtual_dimension(
            genus, degree, contacts
        ) - ZZ(insertion_codimension)

    @classmethod
    def infinity_dimension_defect(
            cls, genus, degree, contacts, psi_min=0,
            insertion_codimension=0):
        return cls.infinity_required_psi_min_power(
            genus, degree, contacts, insertion_codimension
        ) - ZZ(psi_min)

    @staticmethod
    def denominator_power_bound(complex_dimension, inserted_codimension=0):
        r"""Maximum psi power needed from ``1/(w-psi)``.

        A term with psi power above the remaining complex dimension vanishes.
        """
        remaining = ZZ(complex_dimension) - ZZ(inserted_codimension)
        return max(ZZ(-1), remaining)

    @staticmethod
    def psi_min_power_bound(genus, valence, insertion_codimension=0,
                            ambient_degree=0):
        r"""Exact ``psi_min`` bound for balanced plane-cubic infinity data.

        This compatibility helper retains the historical valence-based API.
        On a balanced infinity vertex the reduced virtual dimension is
        ``valence + 3*ambient_degree``.  New code with the contact profile in
        hand should use :meth:`infinity_required_psi_min_power` instead.
        """
        genus = ZZ(genus)
        valence = ZZ(valence)
        ambient_degree = ZZ(ambient_degree)
        if genus < 0 or valence < 0 or ambient_degree < 0:
            raise ValueError("genus, valence, and degree must be nonnegative")
        dimension = valence + 3 * ambient_degree
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
