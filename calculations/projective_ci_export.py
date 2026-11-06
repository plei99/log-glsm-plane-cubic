"""Replay exact certificates, remove cross-target duplicates, and export all examples."""
import gzip,hashlib,json
from pathlib import Path
from sage.all import *

from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.targets.projective import (
    CompleteIntersectionCycles,
    cy3_line_count_grassmannian,
    cy3_ordinary_gw_number,
    requested_complete_intersections,
)
OUT=Path('results/complete-intersections/2026-11-06')

def read(path):
    path=Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_text())

def write(path,data,compressed=False):
    raw=(json.dumps(data,indent=None if compressed else 2,default=int)+'\n').encode()
    Path(str(path)+('.gz' if compressed else '')).write_bytes(gzip.compress(raw,mtime=0) if compressed else raw)

targets=requested_complete_intersections()
catalog=[read(OUT/t.slug/'summary.json') for t in targets]
assert len(catalog)==14 and all(len(s['genera'])==2 for s in catalog)
if any(s['max_degree']!=3 or s['weight_bases']!=[4,5] for s in catalog):
    raise ValueError('this report exports the full degree-three scan with weight bases 4 and 5')
write(OUT/'catalog.json',catalog)
matrices={(g,r):load(str(OUT/('pixton-g%s-r%s.sobj'%(g,r))))
          for g in (2,3) for r in range(3*g-2)}
rows={key:[] for key in matrices}
union_bases={}
certificates=0
for t,summary in zip(targets,catalog):
    for gs in summary['genera']:
        g=gs['genus']
        finite=read(OUT/t.slug/('finite-classes-g%s.json.gz'%g))
        assert len(finite)==4*(3*g-2)
        assert gs['negative_cjr_identities']==sum(t.t_power(g,D,r)<0 for D in range(4) for r in range(3*g-2))
        for ds in gs['degrees']:
            r=ds['codimension'];M=matrices[g,r]
            records=read(OUT/t.slug/('relations-g%s-r%s.json.gz'%(g,r)))
            assert len(records)==ds['relation_rank']
            if not records:continue
            V=matrix(QQ,[rec['strata_vector'] for rec in records])
            U=matrix(QQ,len(records),M.nrows(),sparse=True)
            for i,rec in enumerate(records):
                for j,co in rec['pixton_weights']:U[i,j]=QQ(co)
            assert U*M==V and V.rank()==len(records),(t.slug,g,r)
            rows[g,r].extend(V.rows());certificates+=len(records)
    print('Replayed certificates',t.slug,flush=True)

union=[]
readable=['Every expression below equals zero in rational Chow. Coordinates S_i refer to the dictionary for the same (g,r). These bases remove repetitions across all 14 targets.\n']
for (g,r),M in matrices.items():
    V=matrix(QQ,rows[g,r]) if rows[g,r] else matrix(QQ,0,M.ncols())
    B=V.row_space().basis_matrix()
    union_bases[g,r]=B
    records=[]
    if B.nrows():
        U=M.transpose().solve_right(B.transpose()).transpose()
        assert U*M==B
        for v,u in zip(B.rows(),U.rows()):
            records.append({'strata_vector':[str(c) for c in v],
                'pixton_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
    write(OUT/('union-relations-g%s-r%s.json'%(g,r)),records,True)
    union.append({'genus':g,'codimension':r,'rank':B.nrows(),'pixton_rank':M.rank(),'strata_columns':M.ncols()})
    readable.append('\nGenus %s, codimension %s: %s independent expressions\n'%(g,r,B.nrows()))
    for i,s in enumerate(read(OUT/('strata-g%s-r%s.json'%(g,r)))):readable.append('S_%s = %s\n'%(i,s))
    for i,v in enumerate(B.rows()):readable.append('R_%s: '%i+' + '.join('(%s)*S_%s'%(c,j) for j,c in enumerate(v) if c)+' = 0\n')
    print('Union',g,r,B.nrows(),'/',M.rank(),flush=True)
(OUT/'independent-relations.txt').write_text(''.join(readable))
(OUT/'independent-relations.txt.gz').write_bytes(gzip.compress(''.join(readable).encode(),mtime=0))
write(OUT/'union-summary.json',{'individual_certificates':certificates,'degrees':union,
    'total_union_rank':sum(x['rank'] for x in union)})

ordinary=[]
for t in targets[-5:]:
    values={str(g):{str(D):str(cy3_ordinary_gw_number(t,g,D)) for D in (0,1)} for g in (2,3)}
    infinity={}
    summary=next(s for s in catalog if s['target']['slug']==t.slug)
    for gs in summary['genera']:
        g=gs['genus'];infinity[str(g)]={}
        for D in (0,1):
            infinity[str(g)][str(D)]=str(t.comparison_sign(g,D)*QQ(values[str(g)][str(D)])-QQ(gs['top_twisted_integrals'][str(D)]))
    ordinary.append({'target':t.record(),'line_count':str(cy3_line_count_grassmannian(t)),
        'ordinary_gw_numbers':values,'infinity_aggregate_integrals':infinity})
write(OUT/'ordinary-cy3-gw.json',ordinary)

# Short geometric expressions in the span actually found in this scan.
adm=_load_vendored_admcycles()
R2=adm.TautologicalRing(2,0);R11=adm.TautologicalRing(1,1)
xi11=adm.StableGraph([1,1],[[1],[2]],[(1,2)])
short2=R2.kappa(1)-2*R2.lambdaclass(1)-xi11.boundary_pushforward([R11.one(),R11.one()])/2
R3=adm.TautologicalRing(3,0);R21=adm.TautologicalRing(2,1)
xi12=adm.StableGraph([1,2],[[1],[2]],[(1,2)])
short3=R3.kappa(2)-R3.lambdaclass(1)*R3.kappa(1)+4*R3.lambdaclass(2)-xi12.boundary_pushforward([R11.one(),R21.psi(1)-R21.lambdaclass(1)])
short=[]
for g,r,cl in ((2,1,short2),(3,2,short3)):
    v=vector(QQ,cl.vector(r));B=union_bases[g,r]
    u=B.transpose().solve_right(v)
    assert u*B==v
    short.append({'genus':g,'codimension':r,'strata_vector':[str(c) for c in v],
        'union_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
write(OUT/'short-geometric-relations.json',short)


# Independent direct Hodge integration, without the Pixton-defined integration
# functional used by the batch scan. Cover every ambient dimension and both genera.
checks=[]
for N in range(3,8):
    t=next(t for t in targets if t.N==N)
    s=next(s for s in catalog if s['target']['slug']==t.slug)
    for gs in s['genera']:
        g=gs['genus'];actual=CompleteIntersectionCycles(t,g).top_integral(1)
        assert actual==QQ(gs['top_twisted_integrals']['1']),(t.slug,g,actual)
        checks.append({'target':t.slug,'genus':g,'degree':1,'integral':str(actual)})
    print('Direct Hodge integration checked',t.slug,flush=True)
write(OUT/'verification.json',{'individual_certificates_replayed':certificates,
    'union_certificates_replayed':sum(x['rank'] for x in union),'direct_hodge_checks':checks})

# Tables are included from this generated appendix in the explanatory report.
lines=['All results use unmarked genera two and three, ambient degrees D=0,1,2,3, and two nonresonant weight rays. Each rank counts independent strata expressions before quotienting by tautological relations.\n\n',
'| Target | Genus 2 ranks r=1,2,3 | Genus 3 ranks r=2,3,4,5,6 | Negative-t coefficient slots (g2,g3) |\n|---|---|---|---|\n']
for t,s in zip(targets,catalog):
    a,b=s['genera']
    rank=lambda gs,start: ', '.join(str(x['relation_rank']) for x in gs['degrees'][start:])
    lines.append('| %s | %s | %s | %s, %s |\n'%(t.label,rank(a,1),rank(b,2),a['negative_cjr_identities'],b['negative_cjr_identities']))
lines+=['\nAcross all targets, the independent ranks are:\n\n| Genus | Codimension | Found | Generated Pixton rank |\n|---|---|---:|---:|\n']
for x in union:
    if x['pixton_rank']:lines.append('| %(genus)s | %(codimension)s | %(rank)s | %(pixton_rank)s |\n'%x)
for g in (2,3):
    lines+=['\nExact ambient twisted top integrals, genus %s:\n\n'%g,'| Target | D=0 | D=1 | D=2 | D=3 |\n|---|---|---|---|---|\n']
    for t,s in zip(targets,catalog):
        gs=next(gs for gs in s['genera'] if gs['genus']==g)
        lines.append('| %s | %s |\n'%(t.label,' | '.join(gs['top_twisted_integrals'][str(D)] for D in range(4))))
lines+=['\nOrdinary CY3 GW numbers (distinct from the ambient twisted values above):\n\n',
'| CY3 | Lines | N₂,₀ | N₃,₀ | N₂,₁ | N₃,₁ |\n|---|---:|---|---|---|---|\n']
for rec in ordinary:
    v=rec['ordinary_gw_numbers'];lines.append('| %s | %s | %s | %s | %s | %s |\n'%(rec['target']['name'],rec['line_count'],v['2']['0'],v['3']['0'],v['2']['1'],v['3']['1']))
lines+=['\nIntegrated CY3 infinity aggregates at t⁰, obtained from the signed compact comparison:\n\n',
'| CY3 | g2,D0 | g3,D0 | g2,D1 | g3,D1 |\n|---|---|---|---|---|\n']
for rec in ordinary:
    v=rec['infinity_aggregate_integrals'];lines.append('| %s | %s | %s | %s | %s |\n'%(rec['target']['name'],v['2']['0'],v['3']['0'],v['2']['1'],v['3']['1']))
(OUT/'tables.md').write_text(''.join(lines))
print('VERIFIED AND EXPORTED',certificates,'individual certificates;',sum(x['rank'] for x in union),'independent union relations',flush=True)
