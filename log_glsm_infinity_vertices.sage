r"""
Reconstruct the degree-zero log-GLSM infinity vertices for a plane cubic.

The known Gromov-Witten input is supplied by ``bo_coefficient.sage``.  For
the one-point genus-two correlator it gives

    C_[2] = 1/2 E2^2 + 1/12 E4.

The positive-degree zero-level contribution in the CJR localization formula
is computed from the resummed O(3)-twisted plane blocks

    A2 = sum_{d>0} sigma_1(d) q^d,
    A4 = 1/12 sum_{d>0} sigma_3(d) q^d.

We first compare the known genus-one series ``C_[0]=E2`` with ``A2`` to
solve for ``b2``.  We then solve, rather than assume, ``b4`` in

    C_[2] = A4 + A2^2/2 + b2*A2 + (b4 + b2^2/2).

Here q records intrinsic degree on the elliptic curve.  Its ambient P^2
degree is three times as large.  The answer is b2=-1/24 and b4=1/2880.

Run from the repository directory with

    sage log_glsm_infinity_vertices.sage --max-degree 8

or use the functions after

    load("log_glsm_infinity_vertices.sage")
"""

import argparse
import json
import os
import sys


load("bo_coefficient.sage")


def _validate_max_degree(max_degree):
    max_degree = ZZ(max_degree)
    if max_degree < 1:
        raise ValueError("max_degree must be at least one for reconstruction")
    return max_degree


def paper_eisenstein_qexp(weight, max_degree):
    r"""Return the paper-normalized Eisenstein series through q^max_degree.

    The convention in ``bo_coefficient.sage`` is

        E_k = zeta(1-k)/2 + sum_{d>=1} sigma_(k-1)(d) q^d.
    """
    weight = ZZ(weight)
    max_degree = ZZ(max_degree)
    if weight < 2 or weight % 2:
        raise ValueError("the Eisenstein weight must be a positive even integer")
    if max_degree < 0:
        raise ValueError("max_degree must be nonnegative")

    Q = PowerSeriesRing(QQ, "q", default_prec=max_degree + 1)
    q = Q.gen()
    constant = -bernoulli(weight) / (2 * weight)
    return Q(constant + sum(
        sum(QQ(k)^(weight - 1) for k in divisors(degree)) * q^degree
        for degree in range(1, max_degree + 1)
    )).add_bigoh(max_degree + 1)


def connected_stationary_qseries(descendants, max_degree):
    r"""Expand a connected BO stationary series in intrinsic degree q."""
    max_degree = ZZ(max_degree)
    if max_degree < 0:
        raise ValueError("max_degree must be nonnegative")
    polynomial = connected_coefficient(descendants)
    e2 = paper_eisenstein_qexp(2, max_degree)
    e4 = paper_eisenstein_qexp(4, max_degree)
    e6 = paper_eisenstein_qexp(6, max_degree)
    return polynomial, polynomial(e2, e4, e6)


def o3_twisted_genus_two_blocks(max_degree):
    r"""Return the positive-degree zero-level O(3)-twisted blocks.

    ``A2[d]=d*N_(1,d)=sigma_1(d)`` is the divisor derivative of the
    genus-one O(3)-twisted plane potential.  The primitive genus-two edge
    block is ``A4[d]=sigma_3(d)/12``.  These are the zero-level inputs in
    flat (intrinsic-degree) coordinate after the twisted-plane graph sum.
    """
    max_degree = _validate_max_degree(max_degree)
    Q = PowerSeriesRing(QQ, "q", default_prec=max_degree + 1)
    q = Q.gen()

    a2 = Q(sum(
        sum(QQ(k) for k in divisors(degree)) * q^degree
        for degree in range(1, max_degree + 1)
    )).add_bigoh(max_degree + 1)
    a4 = Q(sum(
        sum(QQ(k)^3 for k in divisors(degree)) * q^degree / 12
        for degree in range(1, max_degree + 1)
    )).add_bigoh(max_degree + 1)
    return {
        "A2_genus_one_divisor": a2,
        "A4_genus_two_primitive": a4,
        "zero_level_total": a4 + a2^2 / 2,
    }


def infinity_contact_profiles(genus, valence):
    r"""List degree-zero infinity contact profiles allowed by balancing.

    For beta_infinity=0 the CJR hypersurface balance equation is

        2-2g = sum_i (1-k_i),  k_i=-c_i >= 1,

    or ``sum_i(k_i-1)=2g-2``.  Returned entries are the negative contact
    orders c_i, sorted increasingly.  For example genus two gives ``(-3,)``
    at valence one and ``(-3,-1), (-2,-2)`` at valence two.
    """
    genus = ZZ(genus)
    valence = ZZ(valence)
    if genus < 0:
        raise ValueError("genus must be nonnegative")
    if valence <= 0:
        raise ValueError("valence must be positive")
    excess = 2 * genus - 2
    if excess < 0:
        return tuple()

    profiles = set()
    for allocation in IntegerVectors(excess, valence):
        positive_orders = sorted((ZZ(x) + 1 for x in allocation), reverse=True)
        profiles.add(tuple(-order for order in positive_orders))
    return tuple(sorted(profiles))


def reconstruct_genus_two_infinity_vertices(max_degree=8):
    r"""Solve the genus-two CJR localization equation for infinity data.

    The genus-one localization equation ``C_[0]=A2+b2`` determines the basic
    contact-minus-one infinity vertex.  In genus two, after subtracting
    ``b2*A2``, the residual must be supported in degree zero.  Its constant
    term determines the primitive genus-two combination b4 after removing
    the graph with two genus-one infinity vertices.

    The one-point equation determines the aggregate b4.  It does not by
    itself separate the individual genus-two profiles ``(-3,)`` and
    ``(-2,-2)`` (or versions with additional contact ``-1`` legs).
    """
    max_degree = _validate_max_degree(max_degree)
    genus_one_polynomial, known_genus_one = connected_stationary_qseries(
        [0], max_degree
    )
    gw_polynomial, known_gw = connected_stationary_qseries([2], max_degree)
    twisted = o3_twisted_genus_two_blocks(max_degree)
    a2 = twisted["A2_genus_one_divisor"]
    zero_level = twisted["zero_level_total"]

    genus_one_residual = known_genus_one - a2
    if any(genus_one_residual[d] for d in range(1, max_degree + 1)):
        raise ArithmeticError(
            "the genus-one O(3)-twisted block leaves positive-degree "
            "infinity data"
        )
    b2 = genus_one_residual[0]

    localization_residual = known_gw - zero_level
    pure_infinity = localization_residual - b2 * a2

    positive_residual = tuple(pure_infinity[d]
                              for d in range(1, max_degree + 1))
    if any(positive_residual):
        raise ArithmeticError(
            "the proposed O(3)-twisted blocks leave positive-degree "
            "infinity data: %s" % (positive_residual,)
        )

    two_genus_one_gluing = b2^2 / 2
    b4 = pure_infinity[0] - two_genus_one_gluing
    reconstructed = zero_level + b2 * a2 + b4 + two_genus_one_gluing

    return {
        "max_degree": max_degree,
        "known_genus_one_polynomial": genus_one_polynomial,
        "known_genus_one_series": known_genus_one,
        "genus_one_localization_residual": genus_one_residual,
        "known_gw_polynomial": gw_polynomial,
        "known_gw_series": known_gw,
        "twisted_blocks": twisted,
        "localization_residual": localization_residual,
        "genus_one_basic_contact_minus_one": b2,
        "genus_two_primitive_profile_combination": b4,
        "two_genus_one_gluing": two_genus_one_gluing,
        "assembled_degree_zero_infinity": b4 + two_genus_one_gluing,
        "reconstructed_gw_series": reconstructed,
        "checks": {
            "positive_degree_infinity_vanishes": not any(positive_residual),
            "reconstructs_known_genus_one": a2 + b2 == known_genus_one,
            "reconstructs_known_gw": reconstructed == known_gw,
            "genus_one_profile_is_contact_minus_one": (
                infinity_contact_profiles(1, 1) == ((-1,),)
            ),
            "genus_two_profiles_have_excess_two": (
                infinity_contact_profiles(2, 1) == ((-3,),)
                and infinity_contact_profiles(2, 2)
                    == ((-3, -1), (-2, -2))
            ),
        },
    }


def _series_coefficients(series, max_degree):
    return [str(series[d]) for d in range(max_degree + 1)]


def infinity_vertex_report(max_degree=8):
    """Return a JSON-serializable report of the reconstruction."""
    result = reconstruct_genus_two_infinity_vertices(max_degree)
    max_degree = result["max_degree"]
    twisted = result["twisted_blocks"]
    return {
        "invariant": "<tau_2(pt)>_{g=2,d}",
        "degree_coordinate": "intrinsic q-degree; ambient P2 degree is 3d",
        "known_genus_one_polynomial": str(
            result["known_genus_one_polynomial"]
        ),
        "known_gw_polynomial": str(result["known_gw_polynomial"]),
        "known_gw_coefficients": _series_coefficients(
            result["known_gw_series"], max_degree
        ),
        "o3_twisted_zero_level": {
            name: _series_coefficients(series, max_degree)
            for name, series in twisted.items()
        },
        "infinity_vertices": {
            "genus_one_basic_contact_minus_one": str(
                result["genus_one_basic_contact_minus_one"]
            ),
            "genus_two_primitive_profile_combination": str(
                result["genus_two_primitive_profile_combination"]
            ),
            "two_genus_one_gluing": str(result["two_genus_one_gluing"]),
            "assembled_degree_zero": str(
                result["assembled_degree_zero_infinity"]
            ),
        },
        "allowed_degree_zero_profiles": {
            "genus_1_valence_1": [
                [int(c) for c in p] for p in infinity_contact_profiles(1, 1)
            ],
            "genus_2_valence_1": [
                [int(c) for c in p] for p in infinity_contact_profiles(2, 1)
            ],
            "genus_2_valence_2": [
                [int(c) for c in p] for p in infinity_contact_profiles(2, 2)
            ],
        },
        "localization_residual_coefficients": _series_coefficients(
            result["localization_residual"], max_degree
        ),
        "checks": result["checks"],
        "identifiability": (
            "The one-point equation determines the aggregate genus-two "
            "primitive combination, not the individual (-3) and (-2,-2) "
            "punctured-vertex integrals."
        ),
    }


def print_infinity_vertex_report(max_degree=8, as_json=False):
    report = infinity_vertex_report(max_degree)
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return

    vertices = report["infinity_vertices"]
    print("Known connected stationary series:")
    print("  genus 1: %s" % report["known_genus_one_polynomial"])
    print("  genus 2: %s" % report["known_gw_polynomial"])
    print("Reconstructed degree-zero infinity vertices:")
    print("  genus 1, basic contact -1: %s" %
          vertices["genus_one_basic_contact_minus_one"])
    print("  genus 2, primitive profile combination: %s" %
          vertices["genus_two_primitive_profile_combination"])
    print("  two genus-1 gluing: %s" % vertices["two_genus_one_gluing"])
    print("  assembled degree zero: %s" % vertices["assembled_degree_zero"])
    print("GW coefficients q^0 through q^%s:" % max_degree)
    print("  %s" % ", ".join(report["known_gw_coefficients"]))
    print("Checks:")
    for name, passed in sorted(report["checks"].items()):
        print("  %s=%s" % (name, passed))
    print("Note: %s" % report["identifiability"])


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    print_infinity_vertex_report(arguments.max_degree, arguments.json)


_entrypoint = os.path.basename(sys.argv[0])
if _entrypoint in (
        "log_glsm_infinity_vertices.sage",
        "log_glsm_infinity_vertices.sage.py",
):
    _main()
