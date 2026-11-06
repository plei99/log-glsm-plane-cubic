r"""curves / genus_two."""

from __future__ import annotations

from fractions import Fraction
import math
from log_glsm.curves.elliptic import (
    EllipticCubicLocalization,
)
from log_glsm.curves.series import (
    _Series,
    _bernoulli_number,
    _divisor_power_sum,
)
from log_glsm.curves.stationary import (
    EllipticCurveStationaryTheory,
)

__all__ = ['DegreeZeroInfinityTheory', 'GenusTwoLogGLSMLocalization']


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
