r"""Exact supported backends for ``O(3)``-twisted ``P^2`` theory.

There are three deliberately distinct interfaces:

* ``TwistedZeroVertexBackend`` evaluates stable zero vertices.  It is exact
  for constant genus-zero vertices and accepts a registry of externally
  proven values.  Unsupported higher-genus Hodge data raises a diagnostic.
* ``O3TwistedIBackend`` evaluates every fixed-point coefficient of the
  conventional Euler-twisted genus-zero I-function.  Its fibre convention is
  ``s+3H``.  The CJR stable-zero class used elsewhere is the *dual* Euler
  twist with fibre weight ``t-3H``; these are related by dualization but must
  not be interchanged term by term.
* ``ResummedO3TwistedTheory`` provides the proven flat-coordinate primitive
  stationary blocks used by the existing plane-cubic reconstruction.

Keeping these interfaces separate prevents a resummed generating-series
identity from being misidentified with an individual localization vertex.
The complete individual-vertex implementation lives in the companion module
``o3_fixed_locus_graphs.sage`` as ``FullTwistedZeroVertexBackend``; this file
retains the lightweight backend so unsupported-geometry behavior can still be
tested or externally populated through the registry.
"""

load("log_glsm_conventions.sage")
load("hodge_integrals.sage")


class TwistedInsertion(SageObject):
    r"""A scaled ``H^power psi^psi_power`` insertion on a zero vertex.

    ``lift_kind`` records an optional sparse lift for the auxiliary
    ``(C*)^3`` localization on ``P^2``:

    * ``standard`` uses ``H`` or ``H^2``;
    * ``vanish`` uses ``H-lambda_i`` and is available for ``H``;
    * ``support`` uses ``prod_(j != i)(H-lambda_j)`` and is available for
      ``H^2``.

    All three specialize to the same ordinary cohomology class when the
    auxiliary base weights vanish.
    """

    _LIFT_KINDS = ("standard", "vanish", "support")

    def __init__(self, h_power=0, psi_power=0, scale=1, label="",
                 lift_kind="standard", lift_index=None):
        self.h_power = ZZ(h_power)
        self.psi_power = ZZ(psi_power)
        self.scale = QQ(scale)
        self.label = str(label)
        self.lift_kind = str(lift_kind)
        self.lift_index = (None if lift_index is None else ZZ(lift_index))
        if self.h_power < 0 or self.h_power > 2:
            raise ValueError("a P2 cohomology power lies between zero and two")
        if self.psi_power < 0:
            raise ValueError("a psi power must be nonnegative")
        if self.lift_kind not in self._LIFT_KINDS:
            raise ValueError("lift kind must be standard, vanish, or support")
        if self.lift_kind == "standard":
            if self.lift_index is not None:
                raise ValueError("a standard lift has no fixed-point index")
        else:
            if self.lift_index not in (0, 1, 2):
                raise ValueError("a sparse lift index must be 0, 1, or 2")
            if self.lift_kind == "vanish" and self.h_power != 1:
                raise ValueError("only H admits the H-lambda_i lift")
            if self.lift_kind == "support" and self.h_power != 2:
                raise ValueError("only H^2 admits a fixed-point-supported lift")

    @classmethod
    def from_elliptic(cls, insertion):
        if not isinstance(insertion, EllipticInsertion):
            raise TypeError("expected an EllipticInsertion")
        scale, h_power = insertion.ambient_class()
        return cls(h_power, insertion.psi_power, scale, label=insertion.kind + "_E")

    @classmethod
    def H_vanishing_at(cls, fixed_point, psi_power=0, scale=1, label=""):
        return cls(
            1, psi_power, scale, label,
            lift_kind="vanish", lift_index=fixed_point,
        )

    @classmethod
    def point_supported_at(cls, fixed_point, psi_power=0, scale=1,
                           label=""):
        return cls(
            2, psi_power, scale, label,
            lift_kind="support", lift_index=fixed_point,
        )

    def with_sparse_lift(self, fixed_point):
        """Return the same insertion with its natural sparse lift."""
        if self.lift_kind != "standard" or self.h_power == 0:
            return self
        lift_kind = "vanish" if self.h_power == 1 else "support"
        return TwistedInsertion(
            self.h_power, self.psi_power, self.scale, self.label,
            lift_kind=lift_kind, lift_index=fixed_point,
        )

    def with_scale(self, scale):
        """Return the same cohomology/descendant insertion with a new scale."""
        return TwistedInsertion(
            self.h_power, self.psi_power, scale, self.label,
            lift_kind=self.lift_kind, lift_index=self.lift_index,
        )

    def allowed_fixed_points(self):
        """Fixed points where the chosen lift can have nonzero restriction."""
        if self.lift_kind == "support":
            return (self.lift_index,)
        if self.lift_kind == "vanish":
            return tuple(index for index in range(3)
                         if index != self.lift_index)
        return (ZZ(0), ZZ(1), ZZ(2))

    def restriction(self, fixed_point, lambdas):
        """Return the scaled fixed-point restriction of the chosen lift."""
        fixed_point = ZZ(fixed_point)
        if fixed_point not in (0, 1, 2) or len(lambdas) != 3:
            raise ValueError("expected one of the three P2 fixed points")
        if self.h_power == 0:
            value = 1
        elif self.lift_kind == "standard":
            value = lambdas[fixed_point] ** self.h_power
        elif self.lift_kind == "vanish":
            value = lambdas[fixed_point] - lambdas[self.lift_index]
        else:
            if fixed_point != self.lift_index:
                value = 0
            else:
                value = prod(
                    lambdas[fixed_point] - lambdas[other]
                    for other in range(3) if other != fixed_point
                )
        return self.scale * value

    @property
    def codimension(self):
        return self.h_power + self.psi_power

    def signature(self):
        return (
            self.h_power, self.psi_power, self.scale, self.label,
            self.lift_kind, self.lift_index,
        )

    def localization_signature(self):
        r"""Return the mathematical insertion, omitting provenance labels.

        Edge labels only record which CJR flag created an insertion.  They do
        not change the corresponding class on ``Mbar(P2)`` and therefore must
        not split the expensive fixed-locus cache.
        """
        return (
            self.h_power, self.psi_power, self.scale,
            self.lift_kind, self.lift_index,
        )

    def localization_record(self):
        return {
            "h_power": int(self.h_power),
            "psi_power": int(self.psi_power),
            "scale": str(self.scale),
            "lift_kind": self.lift_kind,
            "lift_index": (
                None if self.lift_index is None else int(self.lift_index)
            ),
        }

    def to_record(self):
        return {
            "h_power": int(self.h_power),
            "psi_power": int(self.psi_power),
            "scale": str(self.scale),
            "label": self.label,
            "lift_kind": self.lift_kind,
            "lift_index": (
                None if self.lift_index is None else int(self.lift_index)
            ),
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record.get("h_power", 0),
            record.get("psi_power", 0),
            QQ(record.get("scale", "1")),
            record.get("label", ""),
            lift_kind=record.get("lift_kind", "standard"),
            lift_index=record.get("lift_index"),
        )

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, TwistedInsertion) and self.signature() == other.signature()

    def _repr_(self):
        lift = ""
        if self.lift_kind != "standard":
            lift = ",lift=%s@%s" % (self.lift_kind, self.lift_index)
        return "TwistedInsertion(H^%s psi^%s scale=%s%s%s)" % (
            self.h_power, self.psi_power, self.scale,
            ",%s" % self.label if self.label else "",
            lift,
        )


class TwistedZeroVertexRequest(SageObject):
    """Canonical input for one stable zero-level twisted invariant."""

    def __init__(self, genus, degree, insertions=(), label=""):
        self.genus = ZZ(genus)
        self.degree = ZZ(degree)
        raw_insertions = tuple(insertions)
        self.label = str(label)
        if self.genus < 0 or self.degree < 0:
            raise ValueError("genus and degree must be nonnegative")
        if any(not isinstance(item, TwistedInsertion) for item in raw_insertions):
            raise TypeError("twisted vertex insertions must be TwistedInsertion objects")
        # Stable-map invariants are symmetric in their marked insertions.  A
        # canonical order both makes sparse lift planning deterministic and
        # collapses the factorial family of edge-label permutations before
        # fixed-graph or admcycles work begins.
        self.insertions = tuple(sorted(
            raw_insertions, key=lambda item: item.localization_signature()
        ))

    @property
    def valence(self):
        return ZZ(len(self.insertions))

    @property
    def virtual_dimension(self):
        return PlaneCubicDimension.twisted_zero_virtual_dimension(
            self.genus, self.degree, self.valence
        )

    def signature(self):
        return (
            self.genus,
            self.degree,
            tuple(item.localization_signature() for item in self.insertions),
        )

    def normalized_scales(self):
        r"""Factor multilinear scalar coefficients out of the request."""
        factor = prod(item.scale for item in self.insertions)
        if not factor:
            return QQ.zero(), self
        if all(item.scale == 1 for item in self.insertions):
            return QQ.one(), self
        normalized = tuple(item.with_scale(1) for item in self.insertions)
        return factor, TwistedZeroVertexRequest(
            self.genus, self.degree, normalized, label=self.label
        )

    def to_record(self):
        return {
            "genus": int(self.genus),
            "degree": int(self.degree),
            "insertions": [item.to_record() for item in self.insertions],
            "label": self.label,
        }

    def localization_record(self):
        return {
            "genus": int(self.genus),
            "degree": int(self.degree),
            "insertions": [
                item.localization_record() for item in self.insertions
            ],
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["genus"],
            record["degree"],
            tuple(TwistedInsertion.from_record(item)
                  for item in record.get("insertions", ())),
            label=record.get("label", ""),
        )

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, TwistedZeroVertexRequest) and self.signature() == other.signature()

    def _repr_(self):
        return "TwistedZeroVertexRequest(g=%s,d=%s,insertions=%s)" % (
            self.genus, self.degree, self.insertions
        )


class TwistedZeroVertexBackend(SageObject):
    r"""Stable zero-vertex evaluator with explicit exact support boundaries."""

    def __init__(self, rings=None, hodge_backend=None):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.hodge = hodge_backend or HodgeIntegralBackend()
        self._registry = {}

    def register(self, request, value, provenance="external exact value"):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        self._registry[request] = (self.rings.full(value), str(provenance))

    def provenance(self, request):
        if request in self._registry:
            return self._registry[request][1]
        if request.genus == 0 and request.degree == 0:
            return "constant genus-zero map: Mbar_0,n x P2"
        return "unsupported"

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        if request in self._registry:
            return self._registry[request][0]
        if request.genus == 0 and request.degree == 0:
            return self._constant_genus_zero(request)
        raise UnsupportedGeometryError(
            "no exact individual O(3)-twisted zero-vertex backend for %s; "
            "implement base-torus localization/general Hodge integrals or register "
            "a proven value" % request
        )

    def _constant_genus_zero(self, request):
        if request.valence < 3:
            raise ValueError("a constant genus-zero stable vertex needs at least three flags")
        psi_powers = tuple(item.psi_power for item in request.insertions)
        psi_value = self.hodge.integral(0, psi_powers)
        if not psi_value:
            return self.rings.full(0)
        h_power = sum(item.h_power for item in request.insertions)
        scale = prod(item.scale for item in request.insertions)

        # e((R pi_*f^*O(3))^vee)=t-3H for a constant rational map.
        if h_power == 2:
            p2_value = self.rings.t
        elif h_power == 1:
            p2_value = -3
        else:
            p2_value = 0
        return self.rings.full(scale * psi_value * p2_value)


class CollectingTwistedZeroVertexBackend(TwistedZeroVertexBackend):
    r"""Record zero-vertex requests without performing their integrations.

    Returning zero is intentional: the graph compiler still visits every
    local zero-vertex option, but it does not spend time assembling infinity
    contributions that are irrelevant to a precomputation inventory.
    """

    def __init__(self, rings=None):
        super().__init__(rings=rings)
        self._requests = {}

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        self._requests.setdefault(request, request)
        return self.rings.full(0)

    @property
    def requests(self):
        return tuple(sorted(
            self._requests.values(),
            key=lambda request: (
                int(request.degree), int(request.genus), int(request.valence),
                repr(request.signature()),
            ),
        ))


class O3TwistedIBackend(SageObject):
    """Fixed-point restrictions of the equivariant genus-zero I-function."""

    def __init__(self):
        self.polynomial = PolynomialRing(
            QQ, names=("lambda0", "lambda1", "lambda2", "s", "z")
        )
        self.field = self.polynomial.fraction_field()
        self.lambda0, self.lambda1, self.lambda2, self.s, self.z = \
            self.polynomial.gens()
        self.lambdas = (self.lambda0, self.lambda1, self.lambda2)

    def restriction(self, fixed_point, degree):
        fixed_point = ZZ(fixed_point)
        degree = ZZ(degree)
        if fixed_point < 0 or fixed_point > 2 or degree < 0:
            raise ValueError("fixed point is 0,1,2 and degree is nonnegative")
        if degree == 0:
            return self.field(1)
        weight = self.lambdas[fixed_point]
        numerator = prod(
            self.s + 3 * weight + m * self.z
            for m in range(1, 3 * degree + 1)
        )
        denominator = prod(
            weight - other + m * self.z
            for other in self.lambdas
            for m in range(1, degree + 1)
        )
        return self.field(numerator / denominator)

    def section_weights_on_line(self, fixed_left, fixed_right, degree):
        fixed_left = ZZ(fixed_left)
        fixed_right = ZZ(fixed_right)
        degree = ZZ(degree)
        if fixed_left == fixed_right or fixed_left not in range(3) \
                or fixed_right not in range(3) or degree <= 0:
            raise ValueError("choose two distinct fixed points and positive degree")
        left = self.lambdas[fixed_left]
        right = self.lambdas[fixed_right]
        return tuple(
            self.s + QQ(3 * degree - k) / degree * left
            + QQ(k) / degree * right
            for k in range(3 * degree + 1)
        )


class ResummedO3TwistedTheory(SageObject):
    r"""Flat-coordinate primitive stationary blocks for the plane cubic."""

    def __init__(self, max_intrinsic_degree):
        self.max_degree = ZZ(max_intrinsic_degree)
        if self.max_degree < 0:
            raise ValueError("degree bound must be nonnegative")
        self.Q = PowerSeriesRing(QQ, "q", default_prec=self.max_degree + 1)
        self.q = self.Q.gen()

    def primitive_block(self, genus):
        genus = ZZ(genus)
        if genus <= 0:
            raise ValueError("the stationary primitive genus is positive")
        scale = QQ(2) / factorial(2 * genus)
        return self.Q(sum(
            scale * sum(QQ(k) ** (2 * genus - 1) for k in divisors(degree))
            * self.q ** degree
            for degree in range(1, self.max_degree + 1)
        )).add_bigoh(self.max_degree + 1)

    def genus_two_zero_level(self):
        a2 = self.primitive_block(1)
        a4 = self.primitive_block(2)
        return a4 + a2 ** 2 / 2
