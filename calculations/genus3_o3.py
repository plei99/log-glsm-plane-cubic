"""Independent genus-three ordinary Euler ancestors and CJR dual descendants."""
import json,sys,time
from pathlib import Path
from types import SimpleNamespace
from sage.all import *

from log_glsm.localization.o3_givental import (
    O3CalibratedGiventalBackend,
)
from log_glsm.localization.plane import (
    P2FixedLocusEvaluator,
)
from log_glsm.localization.twisted_plane import (
    TwistedInsertion,
    TwistedZeroVertexRequest,
)
sys.path.insert(0,str(Path.cwd().parent))
from t_twisted_cubic import TTwistedCubic
from formal_cubic_modular import modular_series

out=Path('results/genus3/corrected-2026-11-05')
P=PolynomialRing(QQ,'epsilon');epsilon=P.gen();F=P.fraction_field()
rings=SimpleNamespace(full_field=F,t=F.one(),
    base_weights=lambda:(F.zero(),F(epsilon),F(3*epsilon)),full=F)
calibrated=O3CalibratedGiventalBackend(rings=rings,lift_strategy='standard')
localization=P2FixedLocusEvaluator(rings=rings,lift_strategy='standard')
direct=TTwistedCubic()
records=[]
for insertion in ((1,4),(2,3)):
    series=modular_series(direct.ancestor(3,(insertion,)),7)
    for degree in (1,2):
        request=TwistedZeroVertexRequest(3,degree,(TwistedInsertion(*insertion),))
        sign=ZZ(-1)**(3*degree+1-3)
        ancestor=QQ(F(sign*calibrated.evaluate_ancestor(request))(0))
        assert ancestor==QQ(series[degree])
        row={'genus':3,'ambient_degree':degree,'insertion':list(insertion),
             'ordinary_ancestor':str(ancestor),'dual_ancestor':str(sign*ancestor),
             'direct_series':[str(series[d]) for d in range(7)]}
        if degree==1:
            descendant=F(calibrated.evaluate(request))
            assert descendant==localization.evaluate(request)
            row['dual_descendant']=str(QQ(descendant(0)))
            row['ordinary_descendant']=str(sign*QQ(descendant(0)))
        records.append(row)
        print(row,flush=True)
        (out/'o3-genus3.json').write_text(json.dumps(records,indent=2)+'\n')
print('All genus-three O(3) comparisons passed',flush=True)
