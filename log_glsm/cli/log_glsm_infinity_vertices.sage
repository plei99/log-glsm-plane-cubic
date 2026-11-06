r"""cli / log_glsm_infinity_vertices."""

from sage.all import *

import argparse
from log_glsm.infinity.stationary import (
    print_infinity_vertex_report,
)

__all__ = []


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    print_infinity_vertex_report(arguments.max_degree, arguments.json)
