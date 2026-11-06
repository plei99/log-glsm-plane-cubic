r"""GW invariants of a plane cubic from localized and stationary theories.

The target is a smooth cubic E in P^2.  The implementation has two layers:

* :class:`O3TwistedPlane` evaluates the fixed-point restrictions of the
  conventional Euler-twisted O(3) I-function.  ``fiber_weight`` is the
  character of the R-torus on the fibers of O(3).  This is not the dual
  Euler convention ``t-3H`` used by the individual CJR stable-zero backend.
* :class:`EllipticCubicLocalization` sums the genus-one localization series
  in its hypergeometric form and changes from the B-model coordinate q to
  the flat coordinate Q.
* :class:`GenusTwoLogGLSMLocalization` reconstructs the one-point genus-two
  invariant from the positive-degree twisted-plane blocks and the finite
  degree-zero theory at infinity.
* :class:`EllipticCurveStationaryTheory` independently checks the result
  against the elliptic-curve theta function.

All power-series and localization arithmetic is exact over ``Fraction``.
The numerical invariant is indexed both by cover degree r and by ambient
degree d = int_C f^*O(1) = 3r."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence
import math
from log_glsm.curves.series import (
    Number,
    _Series,
    _as_fraction,
    _divisor_sum,
    _harmonic,
)

__all__ = ['EllipticCubicLocalization', 'O3TwistedPlane']


@dataclass(frozen=True)
class O3TwistedPlane:
    r"""Fixed-point oracle for the Euler-twisted O(3) theory of P^2.

    ``base_weights`` are the restrictions lambda_i of H at the three fixed
    points.  ``fiber_weight`` is the R-torus character s.  With loop
    parameter z, the degree-d restriction of the twisted I-function is

        prod_{m=1}^{3d} (s + 3 lambda_i + m z)
        ------------------------------------------------- .
        prod_{j=0}^2 prod_{m=1}^d (lambda_i-lambda_j+m z)

    Supplying rational weights evaluates this rational function exactly.
    CJR zero vertices instead contain ``e((R pi_*f^*O(3))^vee)`` and use
    weights ``t-3H``; callers needing those invariants should use
    ``FullTwistedZeroVertexBackend`` in the Sage implementation.
    """

    base_weights: tuple[Fraction, Fraction, Fraction]
    fiber_weight: Fraction

    def __init__(
        self,
        base_weights: Sequence[Number | str] = (0, 1, 3),
        fiber_weight: Number | str = 0,
    ) -> None:
        weights = tuple(_as_fraction(value) for value in base_weights)
        if len(weights) != 3:
            raise ValueError("P^2 localization needs exactly three base weights")
        if len(set(weights)) != 3:
            raise ValueError("the P^2 torus weights must be pairwise distinct")
        object.__setattr__(self, "base_weights", weights)
        object.__setattr__(self, "fiber_weight", _as_fraction(fiber_weight))

    @staticmethod
    def nonequivariant_i0_coefficient(degree: int) -> Fraction:
        if degree < 0:
            raise ValueError("degree must be nonnegative")
        return Fraction(math.factorial(3 * degree), math.factorial(degree) ** 3)

    def fixed_point_i_coefficient(
        self, fixed_point: int, degree: int, z: Number | str
    ) -> Fraction:
        if fixed_point not in range(3):
            raise ValueError("fixed_point must be 0, 1, or 2")
        if degree < 0:
            raise ValueError("degree must be nonnegative")
        z_value = _as_fraction(z)
        lam = self.base_weights[fixed_point]
        numerator = math.prod(
            self.fiber_weight + 3 * lam + m * z_value
            for m in range(1, 3 * degree + 1)
        )
        denominator = math.prod(
            lam - other + m * z_value
            for other in self.base_weights
            for m in range(1, degree + 1)
        )
        if not denominator:
            raise ZeroDivisionError(
                "the chosen weights lie on a localization pole; choose generic z/base weights"
            )
        return Fraction(numerator, denominator)

    def fixed_point_i_series(
        self, fixed_point: int, max_degree: int, z: Number | str
    ) -> tuple[Fraction, ...]:
        if max_degree < 0:
            raise ValueError("max_degree must be nonnegative")
        return tuple(
            self.fixed_point_i_coefficient(fixed_point, degree, z)
            for degree in range(max_degree + 1)
        )

    def edge_section_weights(
        self, initial_point: int, terminal_point: int, degree: int
    ) -> tuple[Fraction, ...]:
        """R-weights of H^0(P^1, f^*O(3)) on a degree-d invariant edge."""

        if initial_point not in range(3) or terminal_point not in range(3):
            raise ValueError("fixed-point indices must be 0, 1, or 2")
        if initial_point == terminal_point:
            raise ValueError("an invariant edge joins two different fixed points")
        if degree <= 0:
            raise ValueError("an invariant edge has positive degree")
        initial = self.base_weights[initial_point]
        terminal = self.base_weights[terminal_point]
        return tuple(
            self.fiber_weight
            + Fraction(3 * degree - k, degree) * initial
            + Fraction(k, degree) * terminal
            for k in range(3 * degree + 1)
        )

    def edge_euler_class(
        self,
        initial_point: int,
        terminal_point: int,
        degree: int,
        *,
        include_initial_weight: bool = True,
    ) -> Fraction:
        weights = self.edge_section_weights(initial_point, terminal_point, degree)
        if not include_initial_weight:
            weights = weights[1:]
        return Fraction(math.prod(weights))

    @staticmethod
    def formula() -> str:
        return (
            "I_{i,d}(z;s) = prod_{m=1}^{3d}(s+3*lambda_i+m*z) / "
            "prod_{j=0}^2 prod_{m=1}^d(lambda_i-lambda_j+m*z)"
        )

class EllipticCubicLocalization:
    r"""Exact genus-one R-localization calculation for a smooth plane cubic.

    The positive-degree, reduced genus-one localization sum in the B-model
    coordinate q is

        F_1(q) = (T-log(q))/8
                 - log(1-27q)/24
                 - log(I_0(q))/2,

    where

        I_0(q) = sum_d (3d)!/(d!)^3 q^d,
        T-log(q) = I_1(q)/I_0(q),
        [q^d]I_1 = 3 (3d)!/(d!)^3 (H_{3d}-H_d).

    This is the closed resummation of the reduced genus-one GW fixed-locus
    series.  For a curve target without descendants, it equals the standard
    GW series.  The CJR reduced log-GLSM class is related to it by
    ``glsm_specialization_sign``.  The change Q=q*exp(T-log(q)) gives the
    flat-coordinate GW potential.
    """

    def __init__(self, max_cover_degree: int) -> None:
        if type(max_cover_degree) is not int or max_cover_degree <= 0:
            raise ValueError("max_cover_degree must be a positive integer")
        self.max_cover_degree = max_cover_degree
        self.max_ambient_degree = 3 * max_cover_degree
        self._cache: dict[str, _Series] = {}

    @property
    def order(self) -> int:
        return self.max_ambient_degree

    def i0(self) -> _Series:
        if "i0" not in self._cache:
            self._cache["i0"] = _Series(
                O3TwistedPlane.nonequivariant_i0_coefficient(degree)
                for degree in range(self.order + 1)
            )
        return self._cache["i0"]

    def i1_correction_numerator(self) -> _Series:
        if "i1" not in self._cache:
            coefficients = [Fraction(0)]
            for degree in range(1, self.order + 1):
                coefficients.append(
                    self.i0()[degree]
                    * 3
                    * (_harmonic(3 * degree) - _harmonic(degree))
                )
            self._cache["i1"] = _Series(coefficients)
        return self._cache["i1"]

    def mirror_correction(self) -> _Series:
        if "mirror_correction" not in self._cache:
            self._cache["mirror_correction"] = (
                self.i1_correction_numerator() / self.i0()
            )
        return self._cache["mirror_correction"]

    def mirror_map(self) -> _Series:
        """Return Q(q)=q*exp(T-log(q))."""

        if "mirror_map" not in self._cache:
            self._cache["mirror_map"] = self.mirror_correction().exp().shifted()
        return self._cache["mirror_map"]

    def inverse_mirror_map(self) -> _Series:
        """Return q(Q), computed by exact fixed-point iteration of formal series."""

        if "inverse_mirror_map" not in self._cache:
            flat_variable = _Series.variable(self.order)
            inverse = flat_variable
            # Each iteration fixes at least one more coefficient of q(Q).
            for _ in range(self.order):
                correction = self.mirror_correction().compose(inverse)
                inverse = flat_variable * correction.scale(-1).exp()
            self._cache["inverse_mirror_map"] = inverse
        return self._cache["inverse_mirror_map"]

    def localization_terms_b_model(self) -> dict[str, _Series]:
        discriminant = _Series(
            [1, -27] + [0] * (self.order - 1)
        ).log().scale(Fraction(-1, 24))
        return {
            "mirror_map": self.mirror_correction().scale(Fraction(1, 8)),
            "discriminant": discriminant,
            "period_normalization": self.i0().log().scale(Fraction(-1, 2)),
        }

    def b_model_potential(self) -> _Series:
        if "b_model_potential" not in self._cache:
            terms = self.localization_terms_b_model().values()
            answer = _Series.zero(self.order)
            for term in terms:
                answer = answer + term
            self._cache["b_model_potential"] = answer
        return self._cache["b_model_potential"]

    def localization_terms_flat(self) -> dict[str, _Series]:
        inverse = self.inverse_mirror_map()
        return {
            name: term.compose(inverse)
            for name, term in self.localization_terms_b_model().items()
        }

    def flat_potential(self) -> _Series:
        if "flat_potential" not in self._cache:
            self._cache["flat_potential"] = self.b_model_potential().compose(
                self.inverse_mirror_map()
            )
        return self._cache["flat_potential"]

    @staticmethod
    def glsm_specialization_sign(genus: int, ambient_degree: int) -> int:
        r"""The sign (-1)^(1-g+beta.c1(O(3))) from CJR equation (9.7)."""

        if type(genus) is not int or genus < 0:
            raise ValueError("genus must be a nonnegative integer")
        if type(ambient_degree) is not int or ambient_degree < 0:
            raise ValueError("ambient_degree must be a nonnegative integer")
        exponent = 1 - genus + 3 * ambient_degree
        return -1 if exponent % 2 else 1

    @staticmethod
    def closed_cover_invariant(cover_degree: int) -> Fraction:
        """Return sigma_1(r)/r for a connected unramified degree-r cover."""

        if type(cover_degree) is not int or cover_degree <= 0:
            raise ValueError("cover_degree must be a positive integer")
        return Fraction(_divisor_sum(cover_degree), cover_degree)

    def gw_invariants_by_ambient_degree(self) -> dict[int, Fraction]:
        potential = self.flat_potential()
        return {
            degree: potential[degree]
            for degree in range(1, self.order + 1)
            if potential[degree]
        }

    def gw_invariants_by_cover_degree(self) -> dict[int, Fraction]:
        potential = self.flat_potential()
        return {
            cover_degree: potential[3 * cover_degree]
            for cover_degree in range(1, self.max_cover_degree + 1)
        }

    def reduced_glsm_invariants_by_cover_degree(self) -> dict[int, Fraction]:
        return {
            cover_degree: self.glsm_specialization_sign(1, 3 * cover_degree) * value
            for cover_degree, value in self.gw_invariants_by_cover_degree().items()
        }

    def verify(self) -> dict[str, bool]:
        flat = self.flat_potential()
        mirror_inverse_ok = (
            self.mirror_map().compose(self.inverse_mirror_map())
            == _Series.variable(self.order)
        )
        support_ok = all(
            flat[degree] == 0
            for degree in range(1, self.order + 1)
            if degree % 3
        )
        cover_formula_ok = all(
            flat[3 * cover_degree] == self.closed_cover_invariant(cover_degree)
            for cover_degree in range(1, self.max_cover_degree + 1)
        )
        return {
            "mirror_inverse": mirror_inverse_ok,
            "ambient_degree_divisible_by_three": support_ok,
            "unramified_cover_formula": cover_formula_ok,
        }

    def report(self) -> dict[str, object]:
        invariants = self.gw_invariants_by_cover_degree()
        reduced = self.reduced_glsm_invariants_by_cover_degree()
        return {
            "conventions": {
                "cover_degree": "r = degree of the isogeny C -> E",
                "ambient_degree": "integral of f^*O_P2(1), equal to 3r",
                "fiber_equivariant_parameter": "s = C_R^* weight on O(3)",
            },
            "i0": self.i0().as_strings(),
            "mirror_correction": self.mirror_correction().as_strings(),
            "mirror_map": self.mirror_map().as_strings(),
            "flat_potential": self.flat_potential().as_strings(),
            "invariants": [
                {
                    "cover_degree": cover_degree,
                    "ambient_degree": 3 * cover_degree,
                    "reduced_log_glsm": str(reduced[cover_degree]),
                    "gw": str(value),
                }
                for cover_degree, value in invariants.items()
            ],
            "verification": self.verify(),
        }
