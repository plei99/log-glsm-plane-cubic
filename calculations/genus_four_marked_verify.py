"""Replay saved certificates, numerical identities, and the displayed relation."""
import gzip,json
from pathlib import Path
from sage.all import *

from log_glsm.cycles.marked import (
    genus_four_backend,
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
from log_glsm.targets.grassmannian import (
    grassmannian_seven_hyperplane_cy3,
)
OUT=Path('results/genus-four-marked/2026-11-06')
summary=json.loads((OUT/'class-summary.json').read_text())
certificate_rows=0
for row in summary['degrees']:
    n,r=row['markings'],row['codimension']
    M=load(str(OUT/('pixton-n%s-r%s.sobj'%(n,r))))
    U=load(str(OUT/('certificate-n%s-r%s.sobj'%(n,r))))
    V=load(str(OUT/('relations-n%s-r%s.sobj'%(n,r))))
    assert U*M==V and V.rank()==row['combined_rank']
    certificate_rows+=V.nrows()
sample=json.loads((OUT/'sample-relation.json').read_text())
v=vector(QQ,sample['coefficients']);M=load(str(OUT/'pixton-n0-r2.sobj'))
assert M.transpose().solve_right(v)*M==v
with gzip.open(OUT/'p2-3-classes-r3.json.gz','rt') as f:d=json.load(f)
row=next(x for x in d['relations'] if x['markings']==0 and x['codimension']==2)['rows'][0]
w=vector(QQ,17,{int(i):QQ(x) for i,x in row.items()})
assert v==QQ(13824)/169*w
sectors=0;profiles=0;unknown=0
for target in genus_four_targets():
    d=json.loads((OUT/(target.slug+'-numbers.json')).read_text())
    for row in d['rows']:
        z=QQ(row['twisted_unmarked']);ordinary=row['ordinary_unmarked'];sectors+=1;unknown+=ordinary is None
        for record in row['profiles']:
            factor=marked_multiplier(4,row['H_degree'],record['insertions'])
            assert QQ(record['twisted'])==factor*z
            if ordinary is not None:assert QQ(record['ordinary'])==factor*QQ(ordinary)
            elif factor==0:assert record['ordinary']=='0'
            else:assert record['ordinary'] is None
            profiles+=1
# Input validation must occur even when the grading shortcut applies.
for target,bad in ((genus_four_targets()[0],(-1,)),(grassmannian_seven_hyperplane_cy3(),(0,1))):
    try:genus_four_backend(target).marked_coefficients(bad,('1',),3)
    except ValueError:pass
    else:raise AssertionError('invalid curve class accepted')
result={'exact_certificate_rows_replayed':int(certificate_rows),'sample_relation_normalization_verified':True,
    'numerical_sectors_replayed':sectors,'ordered_profiles_replayed':profiles,'unknown_ordinary_sectors':unknown,
    'invalid_curve_classes_rejected':True}
(OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
