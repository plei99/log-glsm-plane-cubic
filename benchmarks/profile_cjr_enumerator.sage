r"""Profile one finite CJR bipartite graph family."""

import argparse
import cProfile
import io
import pstats
import time

load("cjr_bipartite_graphs.sage")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--genus", type=int, default=3)
parser.add_argument("--markings", type=int, default=0)
parser.add_argument("--degree", type=int, default=3)
parser.add_argument("--limit", type=int, default=30)
arguments = parser.parse_args()

profiler = cProfile.Profile()
profiler.enable()
started = time.time()
graphs = PlaneCubicGraphEnumerator(
    arguments.genus, arguments.markings, arguments.degree
).graphs()
elapsed = time.time() - started
profiler.disable()

stream = io.StringIO()
pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats(
    "cumulative"
).print_stats(arguments.limit)
print("graphs=%s elapsed=%.6f" % (len(graphs), elapsed))
print(stream.getvalue())
