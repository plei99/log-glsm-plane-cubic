r"""curves / stationary."""

from __future__ import annotations

from fractions import Fraction
import math
from log_glsm.curves.elliptic import (
    EllipticCubicLocalization,
)
from log_glsm.curves.series import (
    _Series,
    _divisor_power_sum,
    _integer_partitions,
)

__all__ = ['EllipticCurveStationaryTheory']


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
