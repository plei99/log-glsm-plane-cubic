r"""cli / elliptic_cubic_gw."""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence
import argparse
import json
from log_glsm.curves.elliptic import (
    EllipticCubicLocalization,
    O3TwistedPlane,
)
from log_glsm.curves.genus_two import (
    GenusTwoLogGLSMLocalization,
)
from log_glsm.curves.series import (
    _format_series,
)
from log_glsm.curves.stationary import (
    EllipticCurveStationaryTheory,
)

__all__ = ['main']


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
