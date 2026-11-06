r"""cli / cjr_full_equation_provider."""

from sage.all import *

import argparse
import json
from log_glsm.core.conventions import (
    UnsupportedGeometryError,
)
from log_glsm.infinity.equations import (
    genus_two_end_to_end_report,
    genus_two_profile_split_report,
)

__all__ = []


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=8,
                        help="maximum intrinsic elliptic-curve degree")
    parser.add_argument("--json", action="store_true")
    parser.add_argument(
        "--profile-split-rank", action="store_true",
        help="compile the slower string/dilaton contact-splitting rank witness",
    )
    arguments = parser.parse_args()
    if arguments.profile_split_rank:
        try:
            report = genus_two_profile_split_report()
        except UnsupportedGeometryError as error:
            parser.error(str(error))
        if arguments.json:
            print(json.dumps(report, indent=2, sort_keys=True))
            return
        print("Genus-two contact-profile rank witness")
        print("matrix:")
        for row in report["rows"]:
            print("  ", row)
        print("determinant =", report["determinant"])
        print("full rank =", report["full_rank"])
        return
    report = genus_two_end_to_end_report(arguments.max_degree)
    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print("Plane-cubic log-GLSM genus-two end-to-end aggregate")
    print("intrinsic degree | ambient degree | <pt psi^2>_(2,1,d)")
    print("-----------------+----------------+--------------------")
    for item in report["values"]:
        print("%16s | %14s | %s" % (
            item["intrinsic_degree"], item["ambient_degree"], item["value"]
        ))
    print("b2 =", report["b2"])
    print("b4 (aggregate) =", report["b4_aggregate"])
    print("reconstruction check =", report["reconstruction_check"])
