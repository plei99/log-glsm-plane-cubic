r"""Chow-degree bookkeeping for numerical CJR localization coefficients.

The localization formula of CJR III is an identity in localized equivariant
Chow homology.  If a probe leaves complex dimension ``delta`` after its
insertions, the coefficient of ``t^k`` has ordinary Chow dimension

    delta + k.

Consequently a coefficient can be passed to the scalar infinity-vertex
solver exactly when this dimension is nonpositive: dimension zero is an
honest zero-cycle, while a negative-dimensional coefficient vanishes.  A
positive-dimensional coefficient is a Chow-class relation and needs a
tautological test class before it becomes numerical.

This module centralizes that distinction.  In particular, an
under-dimensioned probe of defect ``delta > 0`` first produces a scalar row at
``t^(-delta)``.  Treating the coefficients ``t^(-1), ..., t^(1-delta)`` as
numbers is invalid and was the source of the former genus-one string-equation
contradiction.
"""

load("log_glsm_conventions.sage")


class ChowCoefficientDegree(SageObject):
    """The ordinary Chow dimension of one Laurent coefficient."""

    def __init__(self, probe, t_power):
        if not isinstance(probe, ProbeSpec):
            raise TypeError("Chow coefficients require a ProbeSpec")
        self.probe = probe
        self.t_power = ZZ(t_power)
        self.probe_dimension = ZZ(probe.dimension_defect())
        self.residual_dimension = self.probe_dimension + self.t_power

    @property
    def is_zero_cycle(self):
        return self.residual_dimension == 0

    @property
    def is_forced_zero(self):
        return self.residual_dimension < 0

    @property
    def is_scalar(self):
        return self.residual_dimension <= 0

    @property
    def is_class_valued(self):
        return self.residual_dimension > 0

    def to_record(self):
        return {
            "probe_dimension": int(self.probe_dimension),
            "t_power": int(self.t_power),
            "residual_dimension": int(self.residual_dimension),
            "is_zero_cycle": bool(self.is_zero_cycle),
            "is_forced_zero": bool(self.is_forced_zero),
        }

    def _repr_(self):
        return "ChowCoefficientDegree(delta=%s,t^%s,dimension=%s)" % (
            self.probe_dimension, self.t_power, self.residual_dimension
        )


class CJRInfinityChowBackend(SageObject):
    r"""Project localized Chow identities to their numerical coefficients.

    The backend intentionally refuses positive-dimensional coefficients.  It
    does not discard them or silently integrate them: callers must first
    supply a complementary tautological class in a future class-pairing
    extension.  The zero-cycle and negative-dimensional projections are
    already sufficient for the first primary genus-three infinity rows.
    """

    def degree(self, probe, t_power):
        return ChowCoefficientDegree(probe, t_power)

    def require_scalar(self, probe, t_power):
        degree = self.degree(probe, t_power)
        if degree.is_class_valued:
            raise UnsupportedGeometryError(
                "the coefficient t^%s has residual Chow dimension %s; "
                "pair it with a codimension-%s tautological class before "
                "using the scalar infinity solver"
                % (degree.t_power, degree.residual_dimension,
                   degree.residual_dimension)
            )
        return degree

    def numerical_t_power(self, probe):
        """Return the highest Laurent power giving a numerical coefficient."""
        if not isinstance(probe, ProbeSpec):
            raise TypeError("expected a ProbeSpec")
        defect = ZZ(probe.dimension_defect())
        if defect < 0:
            raise UnsupportedGeometryError(
                "an over-dimensioned probe vanishes before Chow extraction"
            )
        return -defect

    def scalar_powers(self, probe, requested_powers):
        """Filter a Laurent schedule to its valid numerical projections."""
        return tuple(
            ZZ(power) for power in requested_powers
            if self.degree(probe, power).is_scalar
        )
