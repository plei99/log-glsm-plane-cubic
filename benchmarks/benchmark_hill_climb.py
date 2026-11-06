"""Cold-object exact workloads; run with the Sage environment's Python.

From the repository root: python benchmarks/benchmark_hill_climb.py --case all
Every case constructs fresh evaluators, without persistent result caches.
Module loading and output fingerprinting are outside the timed region.
"""
import argparse
import cProfile
import hashlib
import json
import platform
import time
from sage.all import *
from sage.env import SAGE_VERSION


from log_glsm.core.conventions import (
    ProbeSpec,
)
from log_glsm.infinity.equations import (
    PlaneCubicFullEquationProvider,
)
from log_glsm.infinity.graphs import (
    PlaneCubicGraphEnumerator,
)
from log_glsm.localization.plane import (
    P2FixedLocusEvaluator,
)
from log_glsm.localization.twisted_plane import (
    TwistedInsertion,
    TwistedZeroVertexRequest,
)

CASES = ('graphs-3-0-3', 'graphs-3-1-3', 'o3-3-0', 'o3-3-psi4',
         'o3-2-1', 'relation-3-1-0')


def workload(case):
    if case.startswith('graphs-'):
        g, n, d = map(int, case.split('-')[1:])
        return PlaneCubicGraphEnumerator(g, n, d).graphs()
    if case == 'relation-3-1-0':
        provider = PlaneCubicFullEquationProvider(laurent_precision=8)
        return provider.relation(ProbeSpec.stationary(3, 0, (0,)), -4)
    backend = P2FixedLocusEvaluator()
    if case == 'o3-3-0':
        request = TwistedZeroVertexRequest(3, 0, ())
    elif case == 'o3-3-psi4':
        request = TwistedZeroVertexRequest(3, 0, (TwistedInsertion(1, 4),))
    elif case == 'o3-2-1':
        request = TwistedZeroVertexRequest(2, 1, (TwistedInsertion(1, 3),))
    else:
        raise ValueError(case)
    return backend.evaluate(request)


def exact_record(case, value):
    if case.startswith('graphs-'):
        return sorted((str(g.canonical_key()), str(g.automorphism_order()))
                      for g in value)
    if case.startswith('relation-'):
        return {'known': str(value.known_gw),
                'zero': str(value.twisted_zero_level),
                'terms': sorted(((str(c), [v.to_record() for v in factors])
                                 for c, factors in value.terms),
                                key=lambda term: json.dumps(term, sort_keys=True))}
    return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=CASES + ('all',), default='all')
    parser.add_argument('--output')
    parser.add_argument('--profile')
    args = parser.parse_args()
    report = {'sage': SAGE_VERSION, 'python': platform.python_version(),
              'platform': platform.platform(), 'cases': []}
    profiler = cProfile.Profile() if args.profile else None
    for case in CASES if args.case == 'all' else (args.case,):
        if profiler:
            profiler.enable()
        started = time.perf_counter()
        value = workload(case)
        elapsed = time.perf_counter() - started
        if profiler:
            profiler.disable()
        record = exact_record(case, value)
        digest = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
        row = {'case': case, 'seconds': elapsed, 'sha256': digest,
               'exact_result': record}
        report['cases'].append(row)
        print(json.dumps({k: v for k, v in row.items() if k != 'exact_result'}), flush=True)
        if args.output:
            with open(args.output, 'w') as stream:
                json.dump(report, stream, indent=2, sort_keys=True)
    if profiler:
        profiler.dump_stats(args.profile)


if __name__ == '__main__':
    main()
