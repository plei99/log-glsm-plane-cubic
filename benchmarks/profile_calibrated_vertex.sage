r"""Compare calibrated Givental--Teleman with virtual localization.

The default request is the first nontrivial descendant check
``<H psi>_(1,1,1)``. Run from the repository root.
"""

import argparse
import json
import time

load("o3_givental_teleman.sage")


def profile_calibrated_vertex(genus=1, degree=1, h_power=1, psi_power=1,
                              lift_strategy="sparse"):
    request = TwistedZeroVertexRequest(
        genus, degree, (TwistedInsertion(h_power, psi_power),)
    )
    localization = FullTwistedZeroVertexBackend(
        lift_strategy=lift_strategy
    )
    givental = O3CalibratedGiventalBackend(
        rings=localization.rings, lift_strategy=lift_strategy
    )

    started = time.time()
    givental_value = givental.evaluate(request)
    givental_seconds = time.time() - started

    started = time.time()
    localization_value = localization.evaluate(request)
    localization_seconds = time.time() - started

    return {
        "request": str(request),
        "lift_strategy": str(lift_strategy),
        "equal": bool(givental_value == localization_value),
        "givental_value": str(givental_value),
        "localization_value": str(localization_value),
        "givental_seconds": float(givental_seconds),
        "localization_seconds": float(localization_seconds),
        "givental_over_localization": (
            None if localization_seconds == 0
            else float(givental_seconds / localization_seconds)
        ),
    }


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, default=1)
    parser.add_argument("--degree", type=int, default=1)
    parser.add_argument("--h-power", type=int, choices=(0, 1, 2), default=1)
    parser.add_argument("--psi-power", type=int, default=1)
    parser.add_argument("--standard-lifts", action="store_true")
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    report = profile_calibrated_vertex(
        genus=arguments.genus,
        degree=arguments.degree,
        h_power=arguments.h_power,
        psi_power=arguments.psi_power,
        lift_strategy=(
            "standard" if arguments.standard_lifts else "sparse"
        ),
    )
    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print(report["request"])
    print("lift strategy:", report["lift_strategy"])
    print("exact values agree:", report["equal"])
    print("Givental seconds: %.6f" % report["givental_seconds"])
    print("localization seconds: %.6f" % report["localization_seconds"])
    print("Givental/localization: %.3fx" %
          report["givental_over_localization"])


if __name__ == "__main__":
    _main()
