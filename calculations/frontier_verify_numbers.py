"""Replay numerical fractions and provenance independently of Sage localization."""
import gzip,json,math
from fractions import Fraction as F
from pathlib import Path
OUT=Path('results/marked-frontier/2026-11-06');OLD=Path('results/higher-marked/2026-11-06')
sectors=profiles=unknown=old_seeds=0;targets=set()
for path in sorted((OUT/'numbers').glob('*-n5.json.gz')):
    with gzip.open(path,'rt') as f:d=json.load(f)
    slug=d['target'].get('slug',d['target'].get('id'));g=d['genus'];targets.add(slug)
    prior=OLD/(slug+'-g%s-numbers.json'%g)
    old=json.loads(prior.read_text()) if prior.exists() else None
    for row in d['rows']:
        sectors+=1;unknown+=row['ordinary_unmarked'] is None
        z=F(row['twisted_unmarked']);N=None if row['ordinary_unmarked'] is None else F(row['ordinary_unmarked'])
        if old:
            before=next(x for x in old['rows'] if x['beta']==row['beta'])
            assert F(before['twisted_unmarked'])==z and before['ordinary_unmarked']==row['ordinary_unmarked']
            old_seeds+=1
        assert len(row['profiles'])==364 and len({tuple(x['insertions']) for x in row['profiles']})==364
        for v in row['profiles']:
            p=v['insertions'];k=p.count('psi');h=p.count('H')
            m=0 if '1' in p else row['H_degree']**h*math.prod(2*g-2+h+j for j in range(k))
            assert F(v['multiplier'])==m and F(v['twisted'])==m*z
            compact=F(0) if not m else None if N is None else m*N
            assert (None if v['ordinary'] is None else F(v['ordinary']))==compact
            power=row['alpha']+k+h-(3*g-3+len(p));assert v['top_t_power']==power
            B=-m*z if power<0 else row['sign']*compact-m*z if power==0 and compact is not None else None
            assert (None if v['infinity_aggregate'] is None else F(v['infinity_aggregate']))==B
            profiles+=1
assert len(targets)==54 and sectors==534 and old_seeds==362
summary={'presentations':54,'genera':[2,3],'max_markings':5,'sectors':sectors,
         'ordered_profiles':profiles,'unknown_ordinary_sectors':unknown,'previous_seeds_matched':old_seeds}
(OUT/'numerical-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
