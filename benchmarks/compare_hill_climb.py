"""Compare exact outputs and repeated timings, using Sage rational equality."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
from sage.all import *


from log_glsm.core.conventions import (
    PlaneCubicCoefficientRing,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('measurement_directory', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    before = [json.loads((args.measurement_directory / f'before-{i}.json').read_text())
              for i in range(1, 4)]
    after = [json.loads((args.measurement_directory / f'after-{i}.json').read_text())
             for i in range(1, 4)]
    field = PlaneCubicCoefficientRing().full_field
    cases = [r['case'] for r in before[0]['cases']]
    records = []
    for index, name in enumerate(cases):
        reference = before[0]['cases'][index]['exact_result']
        for run in before + after:
            row = run['cases'][index]
            assert row['case'] == name
            result = row['exact_result']
            if name.startswith('o3-'):
                assert field(reference) == field(result), name
            else:
                assert reference == result, name
        old = [run['cases'][index]['seconds'] for run in before]
        new = [run['cases'][index]['seconds'] for run in after]
        records.append({'case': name, 'before_seconds': old, 'after_seconds': new,
                        'before_median': statistics.median(old), 'after_median': statistics.median(new),
                        'speedup': statistics.median(old)/statistics.median(new), 'exact_equal': True,
                        'graph_count': len(reference) if name.startswith('graphs-') else None})
    old_totals = [sum(row['seconds'] for row in run['cases']) for run in before]
    new_totals = [sum(row['seconds'] for row in run['cases']) for run in after]
    summary = {'environment': {k: v for k, v in before[0].items() if k != 'cases'},
               'cases': records, 'total_before_seconds': old_totals,
               'total_after_seconds': new_totals,
               'total_before_median': statistics.median(old_totals),
               'total_after_median': statistics.median(new_totals),
               'total_speedup': statistics.median(old_totals)/statistics.median(new_totals),
               'all_exact_equal': True}
    args.output_directory.mkdir(parents=True, exist_ok=True)
    (args.output_directory/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    # Keep all exact results, including graph representatives and symmetry
    # orders, so the equality audit can be repeated independently.
    for label, runs in (('before', before), ('after', after)):
        for i, run in enumerate(runs, 1):
            (args.output_directory/f'{label}-{i}.json').write_text(json.dumps(run, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
