r"""Convex complete intersections in non-Fano toric products and type-A flags.

The Hirzebruch factors use the effective basis (S,F), S^2=-e. Partial flags
use the degrees of det(S_k)^vee. GKM normal summands are retained individually;
in particular negative normal degrees on F_e still contribute H^1 factors."""

from sage.all import *

from functools import lru_cache
from itertools import permutations
from math import factorial
from log_glsm.core.hodge import (
    HodgeIntegralBackend,
)
from log_glsm.targets.grassmannian import (
    GrassmannianAmbient,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersection,
    ToricProductAmbient,
    toric_genus_zero_number,
)

__all__ = ['ConvexGKMCompleteIntersection', 'PartialFlagAmbient', 'convex_gkm_ordinary_seed', 'convex_gkm_targets', 'partial_flag_fibre_lines', 'partial_flag_primitive_counts']


class PartialFlagAmbient:
    def __init__(self,n,steps):
        self.n=ZZ(n);self.steps=tuple(ZZ(k) for k in steps)
        if not self.steps or tuple(sorted(set(self.steps)))!=self.steps or not 0<self.steps[0]<=self.steps[-1]<n:
            raise ValueError('strictly increasing proper flag dimensions required')
        cuts=(0,)+self.steps+(self.n,)
        self.blocks=tuple(tuple(range(a,b)) for a,b in zip(cuts,cuts[1:]))
        self.sizes=tuple(len(b) for b in self.blocks)
        self.block_of={i:k for k,b in enumerate(self.blocks) for i in b}
        self.dimension=sum(a*b for i,a in enumerate(self.sizes) for b in self.sizes[i+1:])
        self.picard_rank=len(self.steps);self.polarization=(ZZ(1),)*self.picard_rank
        self.anticanonical=tuple(a+b for a,b in zip(self.sizes,self.sizes[1:]))
        self.slug='flag%s-%s'%(n,'-'.join(map(str,steps)))
        self.points=tuple(p for p in permutations(range(n)) if all(tuple(sorted(p[i] for i in b))==tuple(p[i] for i in b) for b in self.blocks))
        self.point_index={p:i for i,p in enumerate(self.points)}
        roots=tuple((i,j) for i in range(n) for j in range(n) if self.block_of[i]<self.block_of[j])
        basis=identity_matrix(ZZ,n);self.edges={};self.tangent=[]
        for pi,p in enumerate(self.points):
            weight=lambda i,j:tuple(basis[p[j]]-basis[p[i]])
            self.tangent.append(tuple(weight(i,j) for i,j in roots))
            for a,b in roots:
                A,B=self.block_of[a],self.block_of[b];q=list(p);q[a],q[b]=q[b],q[a]
                for block in self.blocks:
                    values=sorted(q[k] for k in block)
                    for k,v in zip(block,values):q[k]=v
                normals=[]
                for i,j in roots:
                    if (i,j)==(a,b):continue
                    degree=ZZ((i==a and self.block_of[j]<=B) or (j==b and self.block_of[i]>=A))
                    normals.append((degree,weight(i,j)))
                beta=tuple(ZZ(A<=k<B) for k in range(self.picard_rank))
                assert 2+sum(d for d,_ in normals)==sum(x*y for x,y in zip(self.anticanonical,beta))
                self.edges[pi,self.point_index[tuple(q)]]={'weight':weight(a,b),'normals':tuple(normals),'class':beta}
        self.P=PolynomialRing(QQ,n,['x%s'%i for i in range(n)]);x=self.P.gens();self.chow=self.P
        self.total_chern=prod(1+x[j]-x[i] for i,j in roots)

    def divisor(self,L):
        if len(L)!=self.picard_rank:raise ValueError('wrong flag Picard rank')
        x=self.P.gens();return -sum(a*sum(x[:k]) for a,k in zip(L,self.steps))

    def homogeneous(self,p,d):return self.P({m:c for m,c in self.P(p).dict().items() if sum(m)==d})

    def integral(self,p):
        p=self.homogeneous(p,self.dimension);weights=tuple(QQ(3)**i for i in range(self.n))
        return sum(p(*[weights[i] for i in point])/prod(self.character_value(w,weights) for w in tangent)
                   for point,tangent in zip(self.points,self.tangent))

    def pairing(self,L,beta):
        if len(L)!=self.picard_rank or len(beta)!=self.picard_rank:raise ValueError('wrong flag Picard rank')
        return sum(a*b for a,b in zip(L,beta))

    def fibre_character(self,L,point):
        p=self.points[point]
        return tuple(-sum(a*ZZ(i in p[:k]) for a,k in zip(L,self.steps)) for i in range(self.n))

    character_value=ToricProductAmbient.character_value
    ray=GrassmannianAmbient.ray

    def record(self):
        return {'slug':self.slug,'flag_steps':list(map(int,self.steps)),'n':int(self.n),
                'dimension':int(self.dimension),'picard_rank':self.picard_rank,
                'fixed_points':len(self.points),'invariant_curves':len(self.edges)//2,
                'anticanonical':list(map(int,self.anticanonical))}

class ConvexGKMCompleteIntersection(ToricCompleteIntersection):
    """The geometric kind is proved for each sample, not inferred from ambient ampleness."""
    def __init__(self,ambient,bundles,identifier,kind,geometry):
        self.ambient=ambient;self.bundles=tuple(tuple(ZZ(c) for c in L) for L in bundles)
        self.rank=len(self.bundles);self.dimension=ambient.dimension-self.rank
        if self.dimension not in (2,3):raise ValueError('surface/threefold sample required')
        self.slug=identifier;self.kind=kind;self.geometry=geometry
        for L in self.bundles:
            if not any(L) or len(L)!=ambient.picard_rank:raise ValueError('nontrivial defining bundle required')
            for (i,j),edge in ambient.edges.items():
                if ambient.pairing(L,edge['class'])<0:raise ValueError('defining bundle is not nef/convex')
                assert vector(ambient.fibre_character(L,j))==vector(ambient.fibre_character(L,i))-ambient.pairing(L,edge['class'])*vector(edge['weight'])
        self.c1=tuple(a-sum(L[i] for L in self.bundles) for i,a in enumerate(ambient.anticanonical))
        if kind in ('K3','CY3'):assert not any(self.c1)
        ct=ambient.total_chern;self.fundamental=ambient.chow.one()
        for L in self.bundles:
            line=ambient.divisor(L);self.fundamental*=line
            ct*=sum((-line)**k for k in range(self.dimension+1))
            ct=sum(ambient.homogeneous(ct,k) for k in range(self.dimension+1))
        self.chern=tuple(ambient.homogeneous(ct,k) for k in range(self.dimension+1))

    def record(self):
        result=super().record();result['geometry']=self.geometry;result['convexity']='nef globally generated line summands'
        return result

def convex_gkm_targets():
    A=ToricProductAmbient(('F2',1));B=ToricProductAmbient(('F2',2))
    C=ToricProductAmbient(('F3',1));D=ToricProductAmbient(('F3',2))
    F=PartialFlagAmbient(3,(1,2));G=PartialFlagAmbient(4,(1,2));J=PartialFlagAmbient(4,(1,3))
    data=(
        (A,((1,2,0),),'NF1','Fano','P1 x P1, a section S+2F times P1'),
        (A,((2,4,2),),'NF2','K3','smooth anticanonical K3 in F2 x P1'),
        (B,((2,4,3),),'NF3','CY3','smooth anticanonical CY3 in F2 x P2'),
        (C,((1,3,0),),'NF4','Fano','P1 x P1, a section S+3F times P1'),
        (D,((1,3,0),),'NF5','Fano','P1 x P2, a section S+3F times P2'),
        (F,((1,1),),'PF1','Fano','degree-six del Pezzo surface in Fl(1,2;3)'),
        (F,((2,2),),'PF2','K3','anticanonical K3 in Fl(1,2;3)'),
        (G,((1,0),(0,1)),'PF3','Fano','Fano threefold in Fl(1,2;4)'),
        (G,((1,1),(1,2)),'PF4','CY3','CY3 in Fl(1,2;4)'),
        (J,((1,1),(1,1)),'PF5','Fano','Fano threefold in Fl(1,3;4)'),
        (J,((1,1),(2,2)),'PF6','CY3','CY3 in Fl(1,3;4)'))
    return tuple(ConvexGKMCompleteIntersection(*args) for args in data)

def convex_gkm_ordinary_seed(target,g,beta):
    if target.kind=='K3':return QQ.zero()
    if (1-g)*(target.dimension-3)+target.ambient.pairing(target.c1,beta):return QQ.zero()
    if not any(beta):return (-1)**g*target.chern_numbers()['C']*HodgeIntegralBackend.top_hodge_triple_constant(g)
    if target.slug=='PF4':
        counts=partial_flag_primitive_counts();factor=abs(bernoulli(2*g))/(2*g*factorial(2*g-2))
        if beta[1]==0:return factor*counts[1,0]*QQ(beta[0])**(2*g-3)
        if beta[0]==0:
            k=ZZ(beta[1]);answer=counts[0,1]*k**(2*g-3)
            if k%2==0:answer+=counts[0,2]*(k/2)**(2*g-3)
            return factor*answer
        if tuple(beta)==(1,1):return factor*counts[1,1]
    # A section S+eF avoids S and has numerical class (1,e). Below our
    # Mori-sum bound it contributes no nonconstant degree in the F_e factor.
    if target.slug in ('NF1','NF4','NF5') and (beta[0] or beta[1]):return QQ.zero()
    return None

def partial_flag_fibre_lines():
    """Top Chern count on Gr(2,4): c2(S* tensor O(1)) c2(S* tensor O(2))."""
    A=GrassmannianAmbient(4);x,y=A.P.gens();H=x+y
    return A.integral(A.chow((x+H)*(y+H)*(x+2*H)*(y+2*H)))

@lru_cache(maxsize=None)
def partial_flag_primitive_counts():
    """Parameter-space counts for pure fibers; mixed count by convex localization.

    Isolation and exclusion of intersecting components are proved in
    docs/CONVEX_GKM_GEOMETRY_2026-11-06.md. Genus-zero values alone would
    not justify the higher-genus use of these counts.
    """
    A=PartialFlagAmbient(4,(1,3));weights=A.ray(4);lines=QQ.zero()
    for point,tangent in zip(A.points,A.tangent):
        a,b,c,d=[weights[i] for i in point]
        numerator=(-b-2*a)*(-c-2*a)*(-2*b-3*a)*(-b-c-3*a)*(-2*c-3*a)
        lines+=numerator/prod(A.character_value(w,weights) for w in tangent)
    # On P3, c(Q*)=1-h+h²-h³ and c3(Q* tensor O(2))=(-1+2-4+8)h³.
    conics=ZZ(-1+2-4+8)
    target=convex_gkm_targets()[8]
    mixed=toric_genus_zero_number(target,(1,1),4)
    return {(1,0):partial_flag_fibre_lines(),(0,1):lines,(0,2):conics,(1,1):mixed}
