"""Replay numerical values, old seeds, and every saved Pixton certificate."""
import json,gzip
from pathlib import Path
from sage.all import *

from log_glsm.cycles.marked import (
    marked_multiplier,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersection,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
)
OUT=Path('results/higher-marked/2026-11-06')


def main():
    profile_count=0;baseline_count=0;unknown=0
    for target in genus_four_targets():
        olddir=('complete-intersections' if isinstance(target,ProjectiveCompleteIntersection)
                else 'grassmannian-27' if target.slug=='G15'
                else 'grassmannian-complete-intersections' if isinstance(target,GrassmannianCompleteIntersection)
                else 'toric-complete-intersections')
        oldpath=Path('results')/olddir/'2026-11-06'/target.slug/'summary.json'
        old=json.loads(oldpath.read_text()) if oldpath.exists() else None
        for g in (2,3):
            d=json.loads((OUT/(target.slug+'-g'+str(g)+'-numbers.json')).read_text())
            for row in d['rows']:
                z=QQ(row['twisted_unmarked']);ordinary=row['ordinary_unmarked'];unknown+=ordinary is None
                if old:
                    seed=next(x for x in old['genera'] if x['genus']==g)['top_twisted_integrals'][','.join(map(str,row['beta']))]
                    assert z==QQ(seed),(target.slug,g,row['beta'],z,seed);baseline_count+=1
                for rec in row['profiles']:
                    m=marked_multiplier(g,row['H_degree'],rec['insertions'])
                    assert QQ(rec['twisted'])==m*z
                    if ordinary is not None:assert QQ(rec['ordinary'])==m*QQ(ordinary)
                    elif not m:assert rec['ordinary']=='0'
                    else:assert rec['ordinary'] is None
                    profile_count+=1
    summary=json.loads((OUT/'class-summary.json').read_text());total=0;target_rows=0
    for item in summary['degrees']:
        g,n,r=item['genus'],item['markings'],item['codimension'];suffix='g%s-n%s-r%s'%(g,n,r)
        M=load(str(OUT/('pixton-'+suffix+'.sobj')));U=load(str(OUT/('certificate-'+suffix+'.sobj')));V=load(str(OUT/('relations-'+suffix+'.sobj')))
        assert U*M==V and V.rank()==item['found_rank'];total+=V.rank()
        W=V.row_space()
        for target in genus_four_targets():
            with gzip.open(OUT/(target.slug+'-g'+str(g)+'-classes.json.gz'),'rt') as f:d=json.load(f)
            rows=next(x for x in d['relations'] if x['markings']==n and x['codimension']==r)['rows']
            for row in rows:
                assert vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in row.items()}) in W
                target_rows+=1
    result={'ordered_profiles_replayed':profile_count,'previous_unmarked_seeds_matched':baseline_count,
        'unknown_ordinary_sectors':unknown,'combined_certificate_rows_replayed':int(total),'per_target_rows_replayed':target_rows}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result)


if __name__=='__main__':main()
