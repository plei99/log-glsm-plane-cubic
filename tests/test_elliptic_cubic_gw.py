"""Tests for the exact plane-cubic R-localization calculator."""

from fractions import Fraction
import json
import subprocess
import sys
import unittest

from elliptic_cubic_gw import (
    DegreeZeroInfinityTheory,
    EllipticCubicLocalization,
    EllipticCurveStationaryTheory,
    GenusTwoLogGLSMLocalization,
    O3TwistedPlane,
)


class TwistedPlaneTests(unittest.TestCase):
    def test_nonequivariant_hypergeometric_coefficients(self) -> None:
        self.assertEqual(
            [O3TwistedPlane.nonequivariant_i0_coefficient(d) for d in range(5)],
            [1, 6, 90, 1680, 34650],
        )

    def test_equivariant_fixed_point_coefficient(self) -> None:
        theory = O3TwistedPlane(base_weights=(0, 1, 3), fiber_weight=5)
        self.assertEqual(
            theory.fixed_point_i_coefficient(fixed_point=0, degree=1, z=2),
            Fraction(-693, 2),
        )

    def test_edge_euler_weights_include_fiber_character(self) -> None:
        theory = O3TwistedPlane(base_weights=(0, 1, 3), fiber_weight=2)
        self.assertEqual(
            theory.edge_section_weights(0, 1, 1),
            (Fraction(2), Fraction(3), Fraction(4), Fraction(5)),
        )
        self.assertEqual(theory.edge_euler_class(0, 1, 1), 120)
        self.assertEqual(
            theory.edge_euler_class(0, 1, 1, include_initial_weight=False), 60
        )


class EllipticCubicTests(unittest.TestCase):
    def setUp(self) -> None:
        self.calculator = EllipticCubicLocalization(max_cover_degree=5)

    def test_twisted_period_and_mirror_map(self) -> None:
        self.assertEqual(
            self.calculator.i0().coefficients[:5],
            (Fraction(1), Fraction(6), Fraction(90), Fraction(1680), Fraction(34650)),
        )
        self.assertEqual(
            self.calculator.mirror_correction().coefficients[:5],
            (Fraction(0), Fraction(15), Fraction(333, 2), Fraction(2669), Fraction(199269, 4)),
        )
        self.assertEqual(
            self.calculator.mirror_map().coefficients[:5],
            (Fraction(0), Fraction(1), Fraction(15), Fraction(279), Fraction(5729)),
        )

    def test_flat_potential_and_cover_formula(self) -> None:
        self.assertEqual(
            self.calculator.gw_invariants_by_cover_degree(),
            {
                1: Fraction(1),
                2: Fraction(3, 2),
                3: Fraction(4, 3),
                4: Fraction(7, 4),
                5: Fraction(6, 5),
            },
        )
        self.assertEqual(
            self.calculator.gw_invariants_by_ambient_degree(),
            {
                3: Fraction(1),
                6: Fraction(3, 2),
                9: Fraction(4, 3),
                12: Fraction(7, 4),
                15: Fraction(6, 5),
            },
        )
        self.assertTrue(all(self.calculator.verify().values()))

    def test_log_glsm_specialization_sign(self) -> None:
        self.assertEqual(
            self.calculator.reduced_glsm_invariants_by_cover_degree(),
            {
                1: Fraction(-1),
                2: Fraction(3, 2),
                3: Fraction(-4, 3),
                4: Fraction(7, 4),
                5: Fraction(-6, 5),
            },
        )

    def test_json_cli(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                "elliptic_cubic_gw.py",
                "invariants",
                "--max-cover-degree",
                "2",
                "--json",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["invariants"][1]["gw"], "3/2")
        self.assertTrue(all(payload["verification"].values()))


class EllipticStationaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.theory = EllipticCurveStationaryTheory(max_degree=6)

    def test_normalized_eisenstein_series(self) -> None:
        self.assertEqual(
            self.theory.c2().coefficients[:5],
            (Fraction(-1, 24), Fraction(1), Fraction(3), Fraction(4), Fraction(7)),
        )
        self.assertEqual(
            self.theory.c4().coefficients[:5],
            (
                Fraction(1, 2880),
                Fraction(1, 12),
                Fraction(3, 4),
                Fraction(7, 3),
                Fraction(73, 12),
            ),
        )

    def test_genus_two_tau2_point_values(self) -> None:
        self.assertEqual(
            self.theory.genus_two_tau2_point_invariants(),
            {
                0: Fraction(7, 5760),
                1: Fraction(1, 24),
                2: Fraction(9, 8),
                3: Fraction(31, 6),
                4: Fraction(343, 24),
                5: Fraction(117, 4),
                6: Fraction(111, 2),
            },
        )
        self.assertTrue(all(self.theory.verify().values()))

    def test_plane_degree_and_reduced_sign(self) -> None:
        self.assertEqual(
            self.theory.genus_two_tau2_point_by_ambient_degree()[6],
            Fraction(9, 8),
        )
        self.assertEqual(
            self.theory.reduced_glsm_genus_two_tau2_point(),
            {
                0: Fraction(-7, 5760),
                1: Fraction(1, 24),
                2: Fraction(-9, 8),
                3: Fraction(31, 6),
                4: Fraction(-343, 24),
                5: Fraction(117, 4),
                6: Fraction(-111, 2),
            },
        )

    def test_stationary_json_cli(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                "elliptic_cubic_gw.py",
                "stationary-g2",
                "--max-degree",
                "2",
                "--json",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["values"][2]["gw"], "9/8")
        self.assertTrue(all(payload["verification"].values()))


class GenusTwoLogGLSMTests(unittest.TestCase):
    def setUp(self) -> None:
        self.calculator = GenusTwoLogGLSMLocalization(max_degree=6)

    def test_degree_zero_infinity_hodge_weights(self) -> None:
        infinity = DegreeZeroInfinityTheory.genus_two_terms()
        self.assertEqual(infinity["genus_one_primitive"], Fraction(-1, 24))
        self.assertEqual(infinity["genus_two_primitive"], Fraction(1, 2880))
        self.assertEqual(infinity["two_genus_one_vertices"], Fraction(1, 1152))
        self.assertEqual(infinity["assembled_degree_zero"], Fraction(7, 5760))
        self.assertEqual(
            DegreeZeroInfinityTheory.assembled_degree_zero(2),
            Fraction(7, 5760),
        )

    def test_positive_degree_infinity_is_excluded_by_balance(self) -> None:
        for genus in range(3):
            for degree in range(1, 7):
                self.assertFalse(
                    DegreeZeroInfinityTheory.positive_degree_can_satisfy_balance(
                        genus, degree
                    )
                )

    def test_localization_graph_terms(self) -> None:
        terms = self.calculator.localization_graph_terms()
        self.assertEqual(
            terms["zero_level_genus_two_primitive"].coefficients[:4],
            (Fraction(0), Fraction(1, 12), Fraction(3, 4), Fraction(7, 3)),
        )
        self.assertEqual(
            terms["zero_level_split_genus_one"].coefficients[:4],
            (Fraction(0), Fraction(0), Fraction(1, 2), Fraction(3)),
        )
        self.assertEqual(
            terms["mixed_degree_zero_infinity_tail"].coefficients[:4],
            (Fraction(0), Fraction(-1, 24), Fraction(-1, 8), Fraction(-1, 6)),
        )
        self.assertEqual(
            terms["pure_degree_zero_infinity"].coefficients[:4],
            (Fraction(7, 5760), Fraction(0), Fraction(0), Fraction(0)),
        )

    def test_reconstruction_matches_requested_values(self) -> None:
        self.assertEqual(
            self.calculator.gw_invariants(),
            {
                0: Fraction(7, 5760),
                1: Fraction(1, 24),
                2: Fraction(9, 8),
                3: Fraction(31, 6),
                4: Fraction(343, 24),
                5: Fraction(117, 4),
                6: Fraction(111, 2),
            },
        )
        self.assertTrue(all(self.calculator.verify().values()))

    def test_log_glsm_json_cli(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                "elliptic_cubic_gw.py",
                "log-glsm-g2",
                "--max-degree",
                "2",
                "--json",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(
            payload["infinity_vertices"]["genus_two_primitive"], "1/2880"
        )
        self.assertEqual(payload["values"][2]["gw"], "9/8")
        self.assertTrue(all(payload["verification"].values()))


if __name__ == "__main__":
    unittest.main()
