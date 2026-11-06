"""Exact Pixton certificates and comparison with unit pullbacks."""
import argparse,gzip,json
from pathlib import Path
from sage.all import *
from itertools import permutations
from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
adm=_load_vendored_admcycles()
from admcycles import DR
from admcycles.moduli import get_moduli
OUT=Path('results/higher-marked/2026-11-06')


def dump(path,data):path.write_text(json.dumps(data,indent=2)+'\n')
def suffix(g,n,r):return 'g%s-n%s-r%s'%(g,n,r)
def sparse(v):return {str(i):str(x) for i,x in enumerate(v) if x}


def prepare():
    for g in (2,3):
        for n in range(4):
            R=adm.TautologicalRing(g,n)
            for r in range(3):
                s=suffix(g,n,r);path=OUT/('pixton-'+s+'.sobj')
                if path.exists():continue
                names=[str(x) for x in R.generators(r)]
                raw=DR.rels_matrix(g,r,n,symm=0,moduli_type=get_moduli('st',DRpy=True),usespin=True)
                M=matrix(QQ,[DR.convert_vector_to_monomial_basis(row,g,r,tuple(range(1,n+1)),
                    moduli_type=get_moduli('st',DRpy=True)) for row in raw.rows()]) if raw.nrows() else matrix(QQ,0,len(names))
                assert M.ncols()==len(names)
                save(M,str(path));dump(OUT/('strata-'+s+'.json'),names)
                print('Pixton',s,M.dimensions(),'rank',M.rank(),flush=True)


def certify():
    sectors={(g,n,r):[] for g in (2,3) for n in range(4) for r in range(3)};candidates=0;individual=0;vectors={}
    for target in genus_four_targets():
        for g in (2,3):
            with gzip.open(OUT/(target.slug+'-g'+str(g)+'-classes.json.gz'),'rt') as f:data=json.load(f)
            for item in data['relations']:
                n,r=item['markings'],item['codimension'];candidates+=item['candidates']
                for i,row in enumerate(item['rows']):sectors[g,n,r].append((target.slug,i,row))
    summary=[]
    for (g,n,r),records in sectors.items():
        s=suffix(g,n,r);M=load(str(OUT/('pixton-'+s+'.sobj')));echelon={};rows=[];origins=[]
        for slug,i,record in records:
            v=vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in record.items()});rem=copy(v)
            for pivot,row in sorted(echelon.items()):
                if rem[pivot]:rem-=rem[pivot]*row
            if rem:
                pivot=next(j for j,x in enumerate(rem) if x);echelon[pivot]=rem/rem[pivot]
                rows.append(v);origins.append({'target':slug,'row':i})
        V=matrix(QQ,rows) if rows else matrix(QQ,0,M.ncols());vectors[g,n,r]=V
        U=M.transpose().solve_right(V.transpose()).transpose() if rows else matrix(QQ,0,M.nrows())
        assert U*M==V
        W=V.row_space()
        assert all(vector(QQ,M.ncols(),{int(j):QQ(c) for j,c in rec.items()}) in W for _,_,rec in records)
        save(V,str(OUT/('relations-'+s+'.sobj')));save(U,str(OUT/('certificate-'+s+'.sobj')))
        dump(OUT/('certificate-'+s+'.json'),{'genus':g,'markings':n,'codimension':r,'equation':'U M = V',
            'origins':origins,'U':[sparse(v) for v in U.rows()],'V':[sparse(v) for v in V.rows()]})
        individual+=len(records)
        item={'genus':g,'markings':n,'codimension':r,'strata':M.ncols(),'pixton_rank':M.rank(),
            'found_rank':V.nrows(),'per_target_rows':len(records),'certificate_verified':True}
        summary.append(item);print(item,flush=True)
    # Compare with the span obtained only by adding unit markings to the
    # preceding marking level, including every choice of forgotten label.
    inheritance=[]
    for g in (2,3):
        for n in range(1,4):
            for r in (1,2):
                R=adm.TautologicalRing(g,n-1);basis=R.generators(r);inherited=[]
                for row in vectors[g,n-1,r].rows():
                    cycle=sum((c*basis[i] for i,c in enumerate(row) if c),R.zero()).forgetful_pullback([n])
                    for order in permutations(range(1,n+1)):
                        inherited.append(vector(QQ,cycle.rename_legs({i+1:j for i,j in enumerate(order)},inplace=False).vector(r)))
                V=vectors[g,n,r];B=matrix(QQ,inherited).row_space().basis_matrix() if inherited else matrix(QQ,0,V.ncols())
                M=load(str(OUT/('pixton-'+suffix(g,n,r)+'.sobj')))
                assert all(row in M.row_space() for row in B.rows())
                combined=V.stack(B).rank()
                inheritance.append({'genus':g,'markings':n,'codimension':r,'found_rank':V.rank(),
                    'unit_pullback_rank':B.rank(),'union_rank':combined,
                    'found_directions_beyond_unit_pullbacks':combined-B.rank()})
                save(B,str(OUT/('unit-pullbacks-'+suffix(g,n,r)+'.sobj')))
    dump(OUT/'class-summary.json',{'presentations':43,'genera':[2,3],'max_markings':3,'curve_bound':2,'max_codimension':2,
        'candidate_vectors':candidates,'per_target_independent_rows':individual,'combined_rank_sum':sum(x['found_rank'] for x in summary),
        'degrees':summary,'unit_pullback_comparison':inheritance})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare-only',action='store_true');args=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    prepare()
    if not args.prepare_only:certify()
