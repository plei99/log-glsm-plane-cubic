r"""cli / cjr_plane_cubic_equation_provider."""

from sage.all import *

import argparse
from log_glsm.infinity.one_point import (
    print_plane_cubic_provider_report,
)

__all__ = []


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-genus", type=int, default=4)
    parser.add_argument("--max-degree", type=int, default=6)
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    print_plane_cubic_provider_report(
        arguments.max_genus, arguments.max_degree, arguments.json
    )
