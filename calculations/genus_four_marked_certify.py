"""Build and replay exact Pixton certificates for the bounded marked search."""
import argparse,gzip,json,time
from pathlib import Path
from sage.all import *

from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.targets.catalogue import (
    genus_four_targets,
)
_load_vendored_admcycles()
from admcycles import DR
from admcycles.moduli import get_moduli
OUT=Path('results/genus-four-marked/2026-11-06')


def dump(path,data):path.write_text(json.dumps(data,indent=2)+'\n')


def prepare(max_r):
    adm=_load_vendored_admcycles()
    for n in range(3):
        R=adm.TautologicalRing(4,n)
        for r in range(max_r+1):
            path=OUT/('pixton-n%s-r%s.sobj'%(n,r))
            if path.exists():continue
            names=[str(s) for s in R.generators(r)]
            raw=DR.rels_matrix(4,r,n,symm=0,moduli_type=get_moduli('st',DRpy=True),usespin=True)
            M=matrix(QQ,[DR.convert_vector_to_monomial_basis(row,4,r,tuple(range(1,n+1)),
                moduli_type=get_moduli('st',DRpy=True)) for row in raw.rows()]) if raw.nrows() else matrix(QQ,0,len(names))
            assert M.ncols()==len(names)
            save(M,str(path));dump(OUT/('strata-n%s-r%s.json'%(n,r)),names)
            print('Pixton',n,r,M.dimensions(),'rank',M.rank(),flush=True)


def certify(max_r):
    sectors={(n,r):[] for n in range(3) for r in range(max_r+1)};individual=0;candidates=0
    for target in genus_four_targets():
        path=OUT/(target.slug+'-classes-r'+str(max_r)+'.json.gz')
        with gzip.open(path,'rt') as f:data=json.load(f)
        for item in data['relations']:
            n,r=item['markings'],item['codimension'];candidates+=item['candidates']
            for i,row in enumerate(item['rows']):sectors[n,r].append((target.slug,i,row))
    summary=[]
    for (n,r),records in sectors.items():
        M=load(str(OUT/('pixton-n%s-r%s.sobj'%(n,r))));vectors=[];echelon={};origins=[]
        for slug,i,record in records:
            v=vector(QQ,M.ncols(),{int(j):QQ(co) for j,co in record.items()});remainder=copy(v)
            for pivot,row in sorted(echelon.items()):
                if remainder[pivot]:remainder-=remainder[pivot]*row
            if remainder:
                pivot=next(j for j,x in enumerate(remainder) if x);echelon[pivot]=remainder/remainder[pivot]
                vectors.append(v);origins.append({'target':slug,'row':i})
        V=matrix(QQ,vectors) if vectors else matrix(QQ,0,M.ncols())
        # A certificate uses the original generated rows, not an assumed rank.
        U=M.transpose().solve_right(V.transpose()).transpose() if V.nrows() else matrix(QQ,0,M.nrows())
        assert U*M==V
        save(U,str(OUT/('certificate-n%s-r%s.sobj'%(n,r))))
        save(V,str(OUT/('relations-n%s-r%s.sobj'%(n,r))))
        # Every target row must lie in the certified combined span.
        W=V.row_space()
        assert all(vector(QQ,M.ncols(),{int(j):QQ(co) for j,co in row.items()}) in W for _,_,row in records)
        dump(OUT/('certificate-n%s-r%s.json'%(n,r)),{'genus':4,'markings':n,'codimension':r,
            'equation':'U * M = V','origins':origins,
            'U':[{str(j):str(x) for j,x in enumerate(row) if x} for row in U.rows()],
            'V':[{str(j):str(x) for j,x in enumerate(row) if x} for row in V.rows()]})
        individual+=len(records)
        summary.append({'markings':n,'codimension':r,'strata':M.ncols(),'pixton_rank':M.rank(),
            'individual_target_rows':len(records),'combined_rank':V.nrows(),'exact_certificate_replayed':True})
        print(summary[-1],flush=True)
    result={'presentations':43,'genus':4,'max_codimension':max_r,'candidate_vectors':candidates,
        'individual_target_rows':individual,'combined_rank_sum':sum(x['combined_rank'] for x in summary),'degrees':summary}
    dump(OUT/'class-summary.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-r',type=int,default=3);p.add_argument('--prepare-only',action='store_true');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);prepare(args.max_r)
    if not args.prepare_only:certify(args.max_r)
