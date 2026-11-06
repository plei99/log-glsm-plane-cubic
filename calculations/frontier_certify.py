"""Certificates and inherited spans at each completed marking level."""
import argparse,json,gzip,time
from pathlib import Path
from itertools import combinations
from sage.all import *
from calculations import marked_frontier as run
from log_glsm.core.hodge import _load_vendored_admcycles
_load_vendored_admcycles()
from admcycles import DR
from admcycles.moduli import get_moduli
OUT=run.OUT;OLD=Path('results/higher-marked/2026-11-06')


def suffix(g,n,r):return 'g%s-n%s-r%s'%(g,n,r)


def exact_span(A):
    """Select a small row set modulo a prime, then prove its span over QQ."""
    if not A.nrows():return A
    for prime in (2147483647,2147483629,2147483587):
        try:
            selected=A.change_ring(GF(prime)).transpose().pivots();break
        except ZeroDivisionError:continue
    else:raise ArithmeticError('all modular row-selection primes divide denominators')
    E=A.matrix_from_rows(selected).echelon_form(algorithm='multimodular',height_guess=1000,proof=True)
    pivots=E.pivots();V=E[:len(pivots)]
    # This exact identity makes the modular selection a heuristic only. If
    # a prime lost rank, fall back to a proved full rational echelon form.
    if A.matrix_from_columns(pivots)*V!=A:
        E=A.echelon_form(algorithm='multimodular',height_guess=1000,proof=True)
        V=E[:len(E.pivots())]
        assert A.matrix_from_columns(V.pivots())*V==A
    return V


def pixton(g,n,r):
    tag=suffix(g,n,r);path=OUT/('pixton-'+tag+'.sobj')
    if path.exists():return load(str(path))
    R=run.common.converter(g).adm.TautologicalRing(g,n);names=[str(v) for v in R.generators(r)]
    raw=DR.rels_matrix(g,r,n,symm=0,moduli_type=get_moduli('st',DRpy=True),usespin=True)
    M=matrix(QQ,[DR.convert_vector_to_monomial_basis(v,g,r,tuple(range(1,n+1)),
                moduli_type=get_moduli('st',DRpy=True)) for v in raw.rows()]) if raw.nrows() else matrix(QQ,0,len(names))
    save(M,str(path));run.dump(OUT/('strata-'+tag+'.json'),names)
    print('PIXTON',tag,M.dimensions(),M.rank(),flush=True);return M


def previous(g,n,r):
    path=OUT/('relations-'+suffix(g,n,r)+'.sobj')
    if n<=3 and not path.exists():path=OLD/path.name
    return load(str(path)) if path.exists() else None


def pullback_matrix(g,n,r):
    """Linear on raw strata; compute each basis image once, before dense sums."""
    path=OUT/('forgetful-'+suffix(g,n,r)+'.sobj')
    if path.exists():return load(str(path))
    R=run.common.converter(g).adm.TautologicalRing(g,n-1)
    T=matrix(QQ,[vector(QQ,b.forgetful_pullback([n]).vector(r)) for b in R.generators(r)])
    save(T,str(path));return T


def certify(n,max_r,allow_partial=False):
    missing=[];records=[]
    for i in range(54):
        for g in (2,3):
            path=run.result_path(i,g,n,max_r)
            if not path.exists():missing.append((i,g));continue
            d=run.read(path)
            if n>=4 and d['profiles']!=[list(p) for p in run.profiles(n)]:
                missing.append((i,g));continue
            records.append(d)
    if missing and not allow_partial:raise ValueError('missing target/genus checkpoints: '+str(missing))
    summary=[]
    for g in (2,3):
        vectors={}
        for r in range(max_r+1):
            tag=suffix(g,n,r);M=pixton(g,n,r);rows=[]
            for d in records:
                if d['genus']!=g:continue
                item=next(s for s in d['relations'] if s['codimension']==r)
                rows.extend(vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in row.items()}) for row in item['rows'])
            direct_count=len(rows)
            if n>=4:
                from calculations.frontier_forgetting import section_matrix
                P=previous(g,n-1,r)
                if P is None:raise ValueError('certify the preceding marking level first')
                for row in (P*pullback_matrix(g,n,r)).rows():
                    for i in range(1,n+1):
                        mapping=list(range(1,n+1));mapping[i-1],mapping[n-1]=n,i
                        rows.append(run.common.permute_vector(row,g,n,r,tuple(mapping)))
                if r:
                    lower=previous(g,n-1,r-1)
                    for row in (lower*section_matrix(g,n,r,1)).rows():
                        for a,b in combinations(range(1,n+1),2):
                            middle=[j for j in range(1,n+1) if j not in (a,b)]
                            mapping=(a,)+tuple(middle)+(b,)
                            rows.append(run.common.permute_vector(row,g,n,r,mapping))
            print('SPAN',tag,len(rows),'rows',flush=True)
            V=exact_span(matrix(QQ,rows)) if rows else matrix(QQ,0,M.ncols());vectors[r]=V
            print('CERTIFY',tag,V.nrows(),flush=True)
            U=M.transpose().solve_right(V.transpose()).transpose() if V.nrows() else matrix(QQ,0,M.nrows())
            assert U*M==V
            # exact_span already proved that every input row lies in V.
            if not missing:
                save(V,str(OUT/('relations-'+tag+'.sobj')));save(U,str(OUT/('certificate-'+tag+'.sobj')))
            run.dump(OUT/(('partial-' if missing else '')+'certificate-'+tag+'.json'),
                     {'equation':'U M = V','V':[run.common.sparse(v) for v in V.rows()],
                      'U':[run.common.sparse(v) for v in U.rows()],'complete':not missing})
            inherited=[];P=previous(g,n-1,r) if n else None
            if P is not None:
                for row in (P*pullback_matrix(g,n,r)).rows():
                    for i in range(1,n+1):
                        mapping=list(range(1,n+1));mapping[i-1],mapping[n-1]=n,i
                        inherited.append(run.common.permute_vector(row,g,n,r,tuple(mapping)))
            B=exact_span(matrix(QQ,inherited)) if inherited else matrix(QQ,0,M.ncols())
            if B.nrows():assert B.row_space().is_subspace(M.row_space())
            ideal=B
            if r==2:
                R=run.common.converter(g).adm.TautologicalRing(g,n);basis=R.generators(1);products=[]
                for row in vectors[1].rows():
                    relation=sum((c*basis[j] for j,c in enumerate(row) if c),R.zero())
                    products.extend(vector(QQ,(relation*divisor).vector(2)) for divisor in basis)
                if products:ideal=B.stack(matrix(QQ,products)).row_space().basis_matrix()
                assert ideal.row_space().is_subspace(M.row_space())
            item={'g':g,'n':n,'r':r,'found':V.nrows(),'pixton':M.rank(),'strata':M.ncols(),
                  'per_target_rows':direct_count,'span_input_rows':len(rows),'unit_pullback':B.rank(),
                  'direct_localization_rows':direct_count,
                  'beyond_unit_pullback':V.stack(B).rank()-B.rank(),
                  'unit_and_divisor_products':ideal.rank(),
                  'beyond_unit_and_divisor_products':V.stack(ideal).rank()-ideal.rank(),
                  'certificate_verified':True}
            summary.append(item);print(item,flush=True)
    run.dump(OUT/(('partial-' if missing else '')+'summary-n%s-r%s.json'%(n,max_r)),
             {'markings':n,'max_codimension':max_r,'complete':not missing,'completed':len(records),
              'missing':missing,'spaces':summary})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--markings',type=int,required=True)
    p.add_argument('--codimension',type=int,default=2);p.add_argument('--partial',action='store_true')
    p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    if a.prepare_only:
        for g in (2,3):
            for r in range(a.codimension+1):pixton(g,a.markings,r)
    else:certify(a.markings,a.codimension,a.partial)
