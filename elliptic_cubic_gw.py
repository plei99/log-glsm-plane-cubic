#!/usr/bin/env python3
"""GW invariants of a plane cubic from localized and stationary theories.

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
degree d = int_C f^*O(1) = 3r.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence


Number = int | Fraction


def _as_fraction(value: Number | str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    return Fraction(value)


def _harmonic(number: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, number + 1)), Fraction(0))


def _divisor_power_sum(number: int, power: int) -> int:
    if number <= 0:
        raise ValueError("divisor power sums are defined here for positive integers")
    if power < 0:
        raise ValueError("the divisor power must be nonnegative")
    total = 0
    root = math.isqrt(number)
    for divisor in range(1, root + 1):
        if number % divisor == 0:
            total += divisor**power
            quotient = number // divisor
            if quotient != divisor:
                total += quotient**power
    return total


def _divisor_sum(number: int) -> int:
    return _divisor_power_sum(number, 1)


def _bernoulli_number(index: int) -> Fraction:
    """Return ``B_index`` with the convention ``B_1=-1/2``."""

    if type(index) is not int or index < 0:
        raise ValueError("the Bernoulli index must be a nonnegative integer")
    values = [Fraction(1)]
    for degree in range(1, index + 1):
        values.append(
            -sum(
                Fraction(math.comb(degree + 1, k)) * values[k]
                for k in range(degree)
            )
            / Fraction(degree + 1)
        )
    return values[index]


@dataclass(frozen=True)
class _Series:
    """A truncated power series over Q with a fixed inclusive order."""

    coefficients: tuple[Fraction, ...]

    def __init__(self, coefficients: Iterable[Number]) -> None:
        object.__setattr__(
            self, "coefficients", tuple(_as_fraction(value) for value in coefficients)
        )
        if not self.coefficients:
            raise ValueError("a truncated series needs at least its constant term")

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    @classmethod
    def zero(cls, order: int) -> "_Series":
        return cls([0] * (order + 1))

    @classmethod
    def one(cls, order: int) -> "_Series":
        return cls([1] + [0] * order)

    @classmethod
    def variable(cls, order: int) -> "_Series":
        if order < 1:
            return cls.zero(order)
        return cls([0, 1] + [0] * (order - 1))

    def __getitem__(self, degree: int) -> Fraction:
        return self.coefficients[degree]

    def _check_order(self, other: "_Series") -> None:
        if self.order != other.order:
            raise ValueError("power series must have the same truncation order")

    def __add__(self, other: "_Series") -> "_Series":
        self._check_order(other)
        return _Series(a + b for a, b in zip(self.coefficients, other.coefficients))

    def __sub__(self, other: "_Series") -> "_Series":
        self._check_order(other)
        return _Series(a - b for a, b in zip(self.coefficients, other.coefficients))

    def scale(self, scalar: Number) -> "_Series":
        scalar = _as_fraction(scalar)
        return _Series(scalar * value for value in self.coefficients)

    def __mul__(self, other: "_Series") -> "_Series":
        self._check_order(other)
        answer = [Fraction(0)] * (self.order + 1)
        for left_degree, left in enumerate(self.coefficients):
            if not left:
                continue
            for right_degree in range(self.order - left_degree + 1):
                right = other[right_degree]
                if right:
                    answer[left_degree + right_degree] += left * right
        return _Series(answer)

    def reciprocal(self) -> "_Series":
        if not self[0]:
            raise ZeroDivisionError("a series with zero constant term is not invertible")
        answer = [1 / self[0]] + [Fraction(0)] * self.order
        for degree in range(1, self.order + 1):
            convolution = sum(
                self[k] * answer[degree - k] for k in range(1, degree + 1)
            )
            answer[degree] = -convolution / self[0]
        return _Series(answer)

    def __truediv__(self, other: "_Series") -> "_Series":
        return self * other.reciprocal()

    def exp(self) -> "_Series":
        if self[0]:
            raise ValueError("exact formal exp is implemented for zero constant term")
        answer = [Fraction(1)] + [Fraction(0)] * self.order
        for degree in range(1, self.order + 1):
            answer[degree] = sum(
                k * self[k] * answer[degree - k] for k in range(1, degree + 1)
            ) / degree
        return _Series(answer)

    def log(self) -> "_Series":
        if self[0] != 1:
            raise ValueError("exact formal log is implemented for constant term one")
        derivative = _Series(
            [(degree + 1) * self[degree + 1] for degree in range(self.order)] + [0]
        )
        quotient = derivative / self
        answer = [Fraction(0)] * (self.order + 1)
        for degree in range(1, self.order + 1):
            answer[degree] = quotient[degree - 1] / degree
        return _Series(answer)

    def compose(self, inner: "_Series") -> "_Series":
        self._check_order(inner)
        if inner[0]:
            raise ValueError("the inner series in a formal composition must vanish at zero")
        answer = _Series.zero(self.order)
        power = _Series.one(self.order)
        for coefficient in self.coefficients:
            if coefficient:
                answer = answer + power.scale(coefficient)
            power = power * inner
        return answer

    def shifted(self, amount: int = 1) -> "_Series":
        if amount < 0:
            raise ValueError("only nonnegative shifts are supported")
        if amount > self.order:
            return _Series.zero(self.order)
        return _Series([0] * amount + list(self.coefficients[: self.order + 1 - amount]))

    def as_strings(self) -> list[str]:
        return [str(value) for value in self.coefficients]


def _format_series(
    series: _Series,
    variable: str,
    *,
    start_degree: int = 0,
    include_zero: bool = False,
) -> str:
    pieces: list[tuple[str, str]] = []
    for degree in range(start_degree, series.order + 1):
        coefficient = series[degree]
        if not coefficient and not include_zero:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = abs(coefficient)
        if degree == 0:
            monomial = str(magnitude)
        else:
            variable_part = variable if degree == 1 else f"{variable}^{degree}"
            monomial = variable_part if magnitude == 1 else f"{magnitude}*{variable_part}"
        pieces.append((sign, monomial))
    if not pieces:
        return "0"
    first_sign, first = pieces[0]
    answer = ("-" if first_sign == "-" else "") + first
    for sign, monomial in pieces[1:]:
        answer += f" {sign} {monomial}"
    return answer


def _integer_partitions(number: int, maximum: int | None = None):
    """Yield integer partitions in weakly decreasing order."""

    if number == 0:
        yield ()
        return
    if maximum is None or maximum > number:
        maximum = number
    for first in range(maximum, 0, -1):
        for rest in _integer_partitions(number - first, first):
            yield (first,) + rest


class EllipticCurveStationaryTheory:
    r"""One-point stationary GW theory of an elliptic curve.

    In the normalization of Okounkov-Pandharipande/Pixton,

        sum_g <tau_{2g-2}(pt)>_g z^(2g)
            = exp(sum_{k>=1} C_{2k}(q) z^(2k)).

    Hence ``<pt psi^2>_{2,1} = C_4 + C_2^2/2``.  Here ``degree`` is
    intrinsic degree on the elliptic curve; its degree after the plane-cubic
    embedding is three times as large.
    """

    def __init__(self, max_degree: int) -> None:
        if type(max_degree) is not int or max_degree < 0:
            raise ValueError("max_degree must be a nonnegative integer")
        self.max_degree = max_degree

    def c2(self) -> _Series:
        return _Series(
            [Fraction(-1, 24)]
            + [
                _divisor_power_sum(degree, 1)
                for degree in range(1, self.max_degree + 1)
            ]
        )

    def c4(self) -> _Series:
        return _Series(
            [Fraction(1, 2880)]
            + [
                Fraction(_divisor_power_sum(degree, 3), 12)
                for degree in range(1, self.max_degree + 1)
            ]
        )

    def genus_two_tau2_point_series(self) -> _Series:
        c2 = self.c2()
        return self.c4() + (c2 * c2).scale(Fraction(1, 2))

    @staticmethod
    def genus_two_tau2_point_closed(degree: int) -> Fraction:
        """Coefficient formula for ``<pt psi^2>_{2,1,d}``."""

        if type(degree) is not int or degree < 0:
            raise ValueError("degree must be a nonnegative integer")
        if degree == 0:
            return Fraction(7, 5760)
        convolution = sum(
            _divisor_power_sum(left_degree, 1)
            * _divisor_power_sum(degree - left_degree, 1)
            for left_degree in range(1, degree)
        )
        return (
            Fraction(_divisor_power_sum(degree, 3), 12)
            - Fraction(_divisor_power_sum(degree, 1), 24)
            + Fraction(convolution, 2)
        )

    def genus_two_tau2_point_invariants(self) -> dict[int, Fraction]:
        series = self.genus_two_tau2_point_series()
        return {degree: series[degree] for degree in range(self.max_degree + 1)}

    def genus_two_tau2_point_by_ambient_degree(self) -> dict[int, Fraction]:
        return {
            3 * degree: value
            for degree, value in self.genus_two_tau2_point_invariants().items()
        }

    def reduced_glsm_genus_two_tau2_point(self) -> dict[int, Fraction]:
        """CJR reduced-class coefficients, indexed by intrinsic degree."""

        return {
            degree: EllipticCubicLocalization.glsm_specialization_sign(
                genus=2, ambient_degree=3 * degree
            )
            * value
            for degree, value in self.genus_two_tau2_point_invariants().items()
        }

    def _partition_q_bracket_series(self, order: int | None = None) -> _Series:
        r"""Independent finite partition evaluation of ``[z^3] 1/Theta(z)``."""

        if order is None:
            order = self.max_degree
        if order < 0 or order > self.max_degree:
            raise ValueError("partition verification order is outside the series range")
        vacuum = Fraction(7, 5760)
        denominator = []
        numerator = []
        for degree in range(order + 1):
            partitions = tuple(_integer_partitions(degree))
            denominator.append(len(partitions))
            shifted_cubic_total = Fraction(0)
            for partition in partitions:
                for index, part in enumerate(partition, start=1):
                    occupied = Fraction(part - index, 1) + Fraction(1, 2)
                    vacuum_occupied = Fraction(-index, 1) + Fraction(1, 2)
                    shifted_cubic_total += (
                        occupied**3 - vacuum_occupied**3
                    ) / math.factorial(3)
            numerator.append(len(partitions) * vacuum + shifted_cubic_total)
        return _Series(numerator) / _Series(denominator)

    def verify(self) -> dict[str, bool]:
        series = self.genus_two_tau2_point_series()
        closed_formula = all(
            series[degree] == self.genus_two_tau2_point_closed(degree)
            for degree in range(self.max_degree + 1)
        )
        partition_order = min(self.max_degree, 12)
        partition_series = self._partition_q_bracket_series(partition_order)
        partition_q_bracket = (
            series.coefficients[: partition_order + 1]
            == partition_series.coefficients
        )
        return {
            "divisor_sum_formula": closed_formula,
            "partition_q_bracket": partition_q_bracket,
        }

    def report(self) -> dict[str, object]:
        invariants = self.genus_two_tau2_point_invariants()
        reduced = self.reduced_glsm_genus_two_tau2_point()
        return {
            "invariant": "<pt psi^2>_{2,1,d}",
            "generating_series": "C4 + C2^2/2",
            "values": [
                {
                    "intrinsic_degree": degree,
                    "ambient_plane_degree": 3 * degree,
                    "reduced_log_glsm": str(reduced[degree]),
                    "gw": str(value),
                }
                for degree, value in invariants.items()
            ],
            "verification": self.verify(),
        }


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


class DegreeZeroInfinityTheory:
    r"""The infinity-vertex input needed by genus-two cubic localization.

    For an infinity vertex of genus ``h`` and ambient plane degree ``D``,
    CJR I, equation (7.2), specializes for ``O_P2(3)`` to

        3*D - (2*h-2) = sum_E (c(E_infinity)+1).

    Every infinity contact has ``c(E_infinity) <= -1``.  In total genus two,
    the left side is positive whenever ``D > 0``; hence all infinity vertices
    have degree zero.  Their primitive Hodge weights are the coefficients of

        log((z/2)/sinh(z/2))
          = sum_{h>=1} -B_(2h)/(2h*(2h)!) z^(2h).

    Only the genus-one and genus-two values are needed below.
    """

    @staticmethod
    def balance_left(genus: int, ambient_degree: int) -> int:
        if type(genus) is not int or genus < 0:
            raise ValueError("genus must be a nonnegative integer")
        if type(ambient_degree) is not int or ambient_degree < 0:
            raise ValueError("ambient_degree must be a nonnegative integer")
        return 3 * ambient_degree - (2 * genus - 2)

    @classmethod
    def positive_degree_can_satisfy_balance(
        cls, genus: int, ambient_degree: int
    ) -> bool:
        """Return a necessary balance test for an infinity vertex.

        The right side of the balance equation is nonpositive.  Thus a
        positive left side proves that the proposed vertex cannot occur.
        """

        return cls.balance_left(genus, ambient_degree) <= 0

    @staticmethod
    def primitive_vertex(genus: int) -> Fraction:
        if type(genus) is not int or genus <= 0:
            raise ValueError("a primitive infinity Hodge vertex has positive genus")
        return -_bernoulli_number(2 * genus) / (
            2 * genus * math.factorial(2 * genus)
        )

    @classmethod
    def genus_two_terms(cls) -> dict[str, Fraction]:
        genus_one = cls.primitive_vertex(1)
        genus_two = cls.primitive_vertex(2)
        return {
            "genus_one_primitive": genus_one,
            "genus_two_primitive": genus_two,
            "two_genus_one_vertices": genus_one * genus_one / 2,
            "assembled_degree_zero": genus_two + genus_one * genus_one / 2,
        }

    @classmethod
    def assembled_degree_zero(cls, genus: int) -> Fraction:
        """Assemble degree-zero infinity graphs through the requested genus."""

        if type(genus) is not int or genus < 0:
            raise ValueError("genus must be a nonnegative integer")
        primitive = _Series(
            [0] + [cls.primitive_vertex(h) for h in range(1, genus + 1)]
        )
        return primitive.exp()[genus]


class GenusTwoLogGLSMLocalization:
    r"""Specialized CJR graph reconstruction of ``<pt psi^2>`` for a cubic.

    This is a numerical genus-two specialization, not a general-purpose
    implementation of punctured R-map moduli.  The positive-degree zero-side
    of the CJR graph sum is packaged into two resummed O(3)-twisted blocks

        A_2(q) = sum_{d>0} sigma_1(d) q^d,
        A_4(q) = (1/12) sum_{d>0} sigma_3(d) q^d.

    ``A_2`` is the divisor insertion in the genus-one twisted-plane theory:
    ``[q^d]A_2 = d * N_{1,d}``.  ``A_4`` is the primitive genus-two edge
    block, equivalently ``sum_{k|d} k^3/12``.  All new infinity data are the
    degree-zero Hodge weights ``b_2=-1/24`` and ``b_4=1/2880``.

    The four localization graph types sum to

        A_4 + A_2^2/2 + b_2*A_2 + (b_4+b_2^2/2).

    The stationary theta formula is used only in :meth:`verify`, as an
    independent check; it is not used to construct the answer.
    """

    def __init__(self, max_degree: int) -> None:
        if type(max_degree) is not int or max_degree < 0:
            raise ValueError("max_degree must be a nonnegative integer")
        self.max_degree = max_degree
        self.infinity = DegreeZeroInfinityTheory()

    def zero_level_genus_one_block(self) -> _Series:
        coefficients = [Fraction(0)]
        for degree in range(1, self.max_degree + 1):
            genus_one = EllipticCubicLocalization.closed_cover_invariant(degree)
            coefficients.append(degree * genus_one)
        return _Series(coefficients)

    def zero_level_genus_two_primitive(self) -> _Series:
        return _Series(
            [Fraction(0)]
            + [
                Fraction(_divisor_power_sum(degree, 3), 12)
                for degree in range(1, self.max_degree + 1)
            ]
        )

    def localization_graph_terms(self) -> dict[str, _Series]:
        a2 = self.zero_level_genus_one_block()
        infinity = self.infinity.genus_two_terms()
        pure_infinity = _Series(
            [infinity["assembled_degree_zero"]] + [0] * self.max_degree
        )
        return {
            "zero_level_genus_two_primitive": self.zero_level_genus_two_primitive(),
            "zero_level_split_genus_one": (a2 * a2).scale(Fraction(1, 2)),
            "mixed_degree_zero_infinity_tail": a2.scale(
                infinity["genus_one_primitive"]
            ),
            "pure_degree_zero_infinity": pure_infinity,
        }

    def gw_series(self) -> _Series:
        answer = _Series.zero(self.max_degree)
        for contribution in self.localization_graph_terms().values():
            answer = answer + contribution
        return answer

    def gw_invariants(self) -> dict[int, Fraction]:
        series = self.gw_series()
        return {degree: series[degree] for degree in range(self.max_degree + 1)}

    def reduced_log_glsm_invariants(self) -> dict[int, Fraction]:
        return {
            degree: EllipticCubicLocalization.glsm_specialization_sign(
                genus=2, ambient_degree=3 * degree
            )
            * value
            for degree, value in self.gw_invariants().items()
        }

    def contributions_by_degree(self) -> dict[int, dict[str, Fraction]]:
        terms = self.localization_graph_terms()
        total = self.gw_series()
        return {
            degree: {
                **{name: series[degree] for name, series in terms.items()},
                "gw_total": total[degree],
            }
            for degree in range(self.max_degree + 1)
        }

    def verify(self) -> dict[str, bool]:
        terms = self.localization_graph_terms()
        reconstructed = _Series.zero(self.max_degree)
        for term in terms.values():
            reconstructed = reconstructed + term
        stationary = EllipticCurveStationaryTheory(
            self.max_degree
        ).genus_two_tau2_point_series()
        infinity_vanishing = all(
            not self.infinity.positive_degree_can_satisfy_balance(genus, degree)
            for genus in range(3)
            for degree in range(1, self.max_degree + 1)
        )
        a2 = self.zero_level_genus_one_block()
        divisor_equation = all(
            a2[degree]
            == degree * EllipticCubicLocalization.closed_cover_invariant(degree)
            for degree in range(1, self.max_degree + 1)
        )
        return {
            "positive_degree_infinity_vanishes": infinity_vanishing,
            "degree_zero_hodge_assembly": (
                self.infinity.assembled_degree_zero(2) == Fraction(7, 5760)
            ),
            "genus_one_divisor_block": divisor_equation,
            "graph_sum_is_stationary_series": reconstructed == stationary,
        }

    def report(self) -> dict[str, object]:
        infinity = self.infinity.genus_two_terms()
        reduced = self.reduced_log_glsm_invariants()
        contributions = self.contributions_by_degree()
        return {
            "invariant": "<pt psi^2>_{2,1,d}",
            "method": "CJR log-GLSM localization with degree-zero infinity vertices",
            "infinity_balance": (
                "3*D-(2*g-2)=sum_E(c_E+1), "
                "D=ambient degree and c_E<=-1"
            ),
            "infinity_vertices": {
                name: str(value) for name, value in infinity.items()
            },
            "values": [
                {
                    "intrinsic_degree": degree,
                    "ambient_plane_degree": 3 * degree,
                    "contributions": {
                        name: str(value)
                        for name, value in contributions[degree].items()
                        if name != "gw_total"
                    },
                    "reduced_log_glsm": str(reduced[degree]),
                    "gw": str(contributions[degree]["gw_total"]),
                }
                for degree in range(self.max_degree + 1)
            ],
            "verification": self.verify(),
        }


def _fraction_argument(value: str) -> Fraction:
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(f"not a rational number: {value}") from error


def _print_invariants(calculator: EllipticCubicLocalization) -> None:
    print("O(3)-twisted period I0(q):")
    print("  " + _format_series(calculator.i0(), "q"))
    print("Mirror map Q(q):")
    print("  " + _format_series(calculator.mirror_map(), "q"))
    print("Genus-one flat potential F1(Q):")
    print("  " + _format_series(calculator.flat_potential(), "Q", start_degree=1))
    print("\ncover r | ambient d | reduced log-GLSM | GW N_{1,d}")
    print("--------+-----------+------------------+-----------")
    gw = calculator.gw_invariants_by_cover_degree()
    glsm = calculator.reduced_glsm_invariants_by_cover_degree()
    for cover_degree in range(1, calculator.max_cover_degree + 1):
        print(
            f"{cover_degree:7d} | {3 * cover_degree:9d} | "
            f"{str(glsm[cover_degree]):>16} | {gw[cover_degree]}"
        )
    checks = calculator.verify()
    print("\nchecks:", ", ".join(f"{name}={value}" for name, value in checks.items()))


def _print_genus_two_stationary(theory: EllipticCurveStationaryTheory) -> None:
    print("Genus-two stationary potential:")
    print("  <pt psi^2>_{2,1}(q) = C4(q) + C2(q)^2/2")
    print("  " + _format_series(theory.genus_two_tau2_point_series(), "q"))
    print("\nintrinsic d | plane degree 3d | reduced log-GLSM | GW")
    print("------------+-----------------+------------------+-----------")
    invariants = theory.genus_two_tau2_point_invariants()
    reduced = theory.reduced_glsm_genus_two_tau2_point()
    for degree, value in invariants.items():
        print(
            f"{degree:11d} | {3 * degree:15d} | "
            f"{str(reduced[degree]):>16} | {value}"
        )
    checks = theory.verify()
    print("\nchecks:", ", ".join(f"{name}={value}" for name, value in checks.items()))


def _print_genus_two_log_glsm(calculator: GenusTwoLogGLSMLocalization) -> None:
    infinity = calculator.infinity.genus_two_terms()
    print("Degree-zero infinity vertex data:")
    print(
        "  b2 = {genus_one_primitive}, b4 = {genus_two_primitive}, "
        "[degree 0] = {assembled_degree_zero}".format(**infinity)
    )
    print("Genus-two CJR localization graph sum:")
    print("  A4 + A2^2/2 + b2*A2 + (b4+b2^2/2)")
    print("  " + _format_series(calculator.gw_series(), "q"))
    print("\nintrinsic d | plane degree 3d | reduced log-GLSM | GW")
    print("------------+-----------------+------------------+-----------")
    invariants = calculator.gw_invariants()
    reduced = calculator.reduced_log_glsm_invariants()
    for degree, value in invariants.items():
        print(
            f"{degree:11d} | {3 * degree:15d} | "
            f"{str(reduced[degree]):>16} | {value}"
        )
    checks = calculator.verify()
    print("\nchecks:", ", ".join(f"{name}={value}" for name, value in checks.items()))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")

    invariants = subparsers.add_parser(
        "invariants", help="compute genus-one plane-cubic GW invariants"
    )
    invariants.add_argument("--max-cover-degree", type=int, default=5)
    invariants.add_argument("--json", action="store_true")

    stationary = subparsers.add_parser(
        "stationary-g2", help="compute <pt psi^2>_{2,1,d} for the plane cubic"
    )
    stationary.add_argument("--max-degree", type=int, default=5)
    stationary.add_argument("--json", action="store_true")

    log_glsm_g2 = subparsers.add_parser(
        "log-glsm-g2",
        help="compute <pt psi^2>_{2,1,d} by the specialized CJR graph sum",
    )
    log_glsm_g2.add_argument("--max-degree", type=int, default=5)
    log_glsm_g2.add_argument("--json", action="store_true")

    twisted = subparsers.add_parser(
        "twisted-i", help="evaluate the R-equivariant O(3)-twisted fixed-point I-series"
    )
    twisted.add_argument("--max-degree", type=int, default=3)
    twisted.add_argument("--fixed-point", type=int, choices=(0, 1, 2), default=0)
    twisted.add_argument(
        "--base-weights", type=_fraction_argument, nargs=3, default=(0, 1, 3)
    )
    twisted.add_argument("--fiber-weight", type=_fraction_argument, default=0)
    twisted.add_argument("--z", type=_fraction_argument, default=5)
    twisted.add_argument("--json", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    arguments = parser.parse_args(argv)
    command = arguments.command or "invariants"
    if command == "invariants":
        max_cover_degree = getattr(arguments, "max_cover_degree", 5)
        as_json = getattr(arguments, "json", False)
        calculator = EllipticCubicLocalization(max_cover_degree)
        if as_json:
            print(json.dumps(calculator.report(), indent=2))
        else:
            _print_invariants(calculator)
        return 0

    if command == "stationary-g2":
        theory = EllipticCurveStationaryTheory(arguments.max_degree)
        if arguments.json:
            print(json.dumps(theory.report(), indent=2))
        else:
            _print_genus_two_stationary(theory)
        return 0

    if command == "log-glsm-g2":
        calculator = GenusTwoLogGLSMLocalization(arguments.max_degree)
        if arguments.json:
            print(json.dumps(calculator.report(), indent=2))
        else:
            _print_genus_two_log_glsm(calculator)
        return 0

    theory = O3TwistedPlane(arguments.base_weights, arguments.fiber_weight)
    coefficients = theory.fixed_point_i_series(
        arguments.fixed_point, arguments.max_degree, arguments.z
    )
    payload = {
        "formula": theory.formula(),
        "fixed_point": arguments.fixed_point,
        "base_weights": [str(value) for value in theory.base_weights],
        "fiber_weight": str(theory.fiber_weight),
        "z": str(arguments.z),
        "coefficients": [str(value) for value in coefficients],
    }
    if arguments.json:
        print(json.dumps(payload, indent=2))
    else:
        print(payload["formula"])
        print("I_i(q) = " + " + ".join(
            f"({coefficient})*q^{degree}"
            for degree, coefficient in enumerate(coefficients)
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
