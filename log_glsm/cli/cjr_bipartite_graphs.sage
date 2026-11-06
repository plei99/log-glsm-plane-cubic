r"""cli / cjr_bipartite_graphs."""

from sage.all import *


from log_glsm.infinity.graphs import (
    PlaneCubicGraphEnumerator,
)

__all__ = []


def _command_line(argv):
    if len(argv) != 4:
        print("usage: sage cjr_bipartite_graphs.sage GENUS MARKINGS DEGREE")
        return Integer(2)
    enumerator = PlaneCubicGraphEnumerator(
        Integer(argv[1]), Integer(argv[2]), Integer(argv[3])
    )
    print(enumerator.report())
    return Integer(0)
