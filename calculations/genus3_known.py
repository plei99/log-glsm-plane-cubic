"""Exact genus-three stationary input and corrected degree-zero CJR equations."""
import json
from pathlib import Path
from sage.all import *

from log_glsm.core.conventions import (
    ProbeSpec,
)
from log_glsm.infinity.equations import (
    PlaneCubicFullEquationProvider,
)
from log_glsm.infinity.stationary import (
    connected_stationary_qseries,
)

out=Path('results/genus3/corrected-2026-11-05');out.mkdir(parents=True,exist_ok=True)
def save(name,data):
    (out/name).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

known=[]
for descendants in ((4,), (3,1), (2,2), (2,1,1)):
    polynomial,series=connected_stationary_qseries(descendants,6)
    known.append({'descendants':list(descendants),'polynomial':str(polynomial),
                  'coefficients':[str(series[d]) for d in range(7)]})
    print(descendants,polynomial,series,flush=True)
    save('stationary-gw.json',known)
save('stationary-gw.json',known)

provider=PlaneCubicFullEquationProvider(laurent_precision=8,
    zero_vertex_cache='results/genus3/zero-vertices.sqlite',
    autosave_zero_vertices=True,hodge_cache='results/genus3/hodge.sqlite')
relations=[]
for descendants,power in (((4,),0),((0,),-4),((3,1),0),((2,2),0),((4,0),0)):
    probe=ProbeSpec.stationary(3,0,descendants)
    relation=provider.boundary_compared_relation(probe,power)
    record=relation.to_record()
    record['graph_count']=len(relation.compilation.contributions)
    relations.append(record)
    print('CJR',descendants,'t power',power,'known',relation.known_gw,
          'twisted',relation.twisted_zero_level,'terms',len(relation.terms),flush=True)
save('selected-relations.json',relations)
