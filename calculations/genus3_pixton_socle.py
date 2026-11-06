"""Check the Pixton quotient in codimension six on Mbar_3.

This checks the ambient tautological relation space, not the unavailable
decorated-strata lifts of the log-GLSM effective cycles.
"""
import json
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path('vendor').resolve()))
from sage.all import *
from admcycles import DR, kappaclass
from admcycles.moduli import get_moduli

started = time.time()
relations = DR.rels_matrix(3, 6, 0, symm=0,
    moduli_type=get_moduli('st', DRpy=True), usespin=True)
rank = relations.rank()
degree = kappaclass(6, 3, 0).evaluate()
assert relations.ncols()-rank == 1
assert degree == QQ(1)/82944
result = {
    'genus': 3, 'markings': 0, 'codimension': 6,
    'method': 'Pixton 3-spin formula, usespin=True',
    'rows': int(relations.nrows()), 'columns': int(relations.ncols()),
    'rank': int(rank), 'quotient_dimension': int(relations.ncols()-rank),
    'integral_kappa_6': str(degree),
    'boxed_rhs_as_kappa_6_coefficient': str((QQ(113)/144)/degree),
    'effective_cycle_strata_lifts_computed': False,
    'elapsed_seconds': time.time()-started,
}
out = Path('results/genus3/corrected-2026-11-05/d1-extension/pixton-socle.json')
out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2), flush=True)
