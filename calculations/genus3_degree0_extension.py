"""Extend the saved first stage with corrected degree-zero two-marking rows.
Run with the Sage environment's Python from the repository root.
"""
import json
from pathlib import Path
import time
from sage.all import *


from log_glsm.infinity.genus_three import (
    GenusThreeVertexComputation,
    print_computation_report,
)

OUTPUT = Path('results/genus3/corrected-2026-11-05/d0-extension')
OUTPUT.mkdir(parents=True, exist_ok=True)


def save(name, value):
    (OUTPUT/name).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def progress(message):
    print(message, flush=True)


computation = GenusThreeVertexComputation(
    genus=3, max_valence=1, max_infinity_degree=1,
    max_markings=2, t_powers=(0,-1,-2,-3,-4,-5,-6),
    include_unit_relatives=False, max_unit_insertions=0,
    include_chow_relations=True, max_chow_unit_insertions=0,
    chow_primary_only=True, effective_basis_only=True,
    additional_probe_ambient_degrees=(1,2),
    zero_vertex_cache='results/genus3/zero-vertices.sqlite',
    hodge_cache='results/genus3/hodge.sqlite',
    progress=progress,
)
progress('Targets: '+str(computation.roots))
computation.load_checkpoint('results/genus3/corrected-2026-11-05/stage1-checkpoint.json')
computation.restrict_future_probe_degree(0)
progress('Keeping all first-stage and imported primary rows; adding only degree-zero probes')
report = computation.run(checkpoint_path=str(OUTPUT/'checkpoint.json'), time_budget=120)
save('report.json', report)
save('values.json', [dict(vertex=v.to_record(), value=str(value))
     for v,value in sorted(computation.solved_values().items(), key=lambda item:item[0].order_key())])
print_computation_report(report)
