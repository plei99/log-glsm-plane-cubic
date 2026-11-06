"""Exact unit insertion with stable-map descendant boundary corrections."""
from log_glsm.cycles.higher_marked import marked_profile_permutations
from functools import lru_cache
from sage.all import *
from calculations import marked_frontier as run
from calculations.frontier_certify import pullback_matrix,suffix,OLD


@lru_cache(maxsize=None)
def section_matrix(g,n,r,i):
    """Push Mbar(g,n-1) into D_{i,n}; r is the resulting codimension."""
    path=run.OUT/('section-'+suffix(g,n,r)+'-i%s.sobj'%i)
    if path.exists():return load(str(path))
    adm=run.common.converter(g).adm
    legs=[n+1 if j==i else j for j in range(1,n)]
    graph=adm.StableGraph([g,0],[legs,[i,n,n+2]],[(n+1,n+2)])
    unit=adm.TautologicalRing(0,3).fundamental_class()
    R=adm.TautologicalRing(g,n-1)
    D=matrix(QQ,[vector(QQ,graph.boundary_pushforward([b,unit]).vector(r)) for b in R.generators(r-1)])
    save(D,str(path));return D


def prior_data(index,g,n,r):
    if n<=3 and index<43:
        slug=run.targets()[index].slug;return run.read(OLD/(slug+'-g%s-classes.json.gz'%g))
    return run.read(run.result_path(index,g,n,r))


def unit_coefficients(index,g,n,max_r):
    old=prior_data(index,g,n-1,max_r);stored={}
    for item in old.get('finite_canonical_classes',old.get('finite_classes',[])):
        p=tuple(item.get('profile',item.get('insertions')))
        if len(p)==n-1:stored[tuple(item['beta']),p,item['codimension']]=item['vector']
    def get(beta,p,r):
        p=tuple(p);canonical=tuple(sorted(p,key=('1','psi','H').index))
        size=len(run.common.converter(g).adm.TautologicalRing(g,n-1).generators(r))
        row=stored.get((beta,canonical,r),{})
        v=vector(QQ,size,{int(j):QQ(c) for j,c in row.items()})
        if p==canonical:return v
        mapping=next(m for pp,m in marked_profile_permutations(canonical) if pp==p)
        return run.common.permute_vector(v,g,n-1,r,mapping)
    added=[]
    for beta in run.common.genus_four_degrees(run.targets()[index]):
        for p in run.profiles(n):
            if '1' not in p:continue
            q=p[1:]
            for r in range(max_r+1):
                v=get(beta,q,r)*pullback_matrix(g,n,r)
                if r:
                    for i,x in enumerate(q,1):
                        if x=='psi':
                            lower=list(q);lower[i-1]='1'
                            v+=get(beta,tuple(lower),r-1)*section_matrix(g,n,r,i)
                # q followed by the new unit becomes the canonical profile p.
                v=run.common.permute_vector(v,g,n,r,tuple(range(2,n+1))+(1,))
                if v:added.append({'beta':list(map(int,beta)),'profile':p,'codimension':r,
                                  'vector':run.common.sparse(v),'method':'unit forgetting with descendant tails'})
    return added


def complete(index,g,n,max_r):
    path=run.result_path(index,g,n,max_r);data=run.read(path)
    if data.get('unit_profiles_reconstructed'):return
    if data['profiles']==[list(p) for p in run.profiles(n)]:return
    added=unit_coefficients(index,g,n,max_r)
    data['direct_profiles']=data['profiles'];data['profiles']=run.profiles(n)
    data['finite_canonical_classes']+=added;data['unit_profiles_reconstructed']=True
    data['unit_formula']='C(I,1)=pi^*C(I)+sum_{I_i=psi} (section_i)_* C(I_i replaced by 1)'
    run.dump(path,data)
    print('FORGOTTEN',run.targets()[index].slug,g,n,len(added),'nonzero canonical coefficients',flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--markings',type=int,default=5);a=p.parse_args()
    for i in range(54):
        for g in (2,3):complete(i,g,a.markings,2)
