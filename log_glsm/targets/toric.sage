r"""Bounded toric/GKM dual-Euler localization with vector-valued curve classes.

Supported ambients are products of projective spaces and Hirzebruch surfaces,
of dimension <=4. The original sample has twelve Fano/K3/CY3 targets. Graphs retain Mori tuples;
negative normal degrees use H1 in the numerator. See the toric source audit."""

from sage.all import *

from functools import lru_cache
from itertools import product
from math import factorial
from log_glsm.core.hodge import (
    AdmcyclesHodgeIntegralBackend,
    _load_vendored_admcycles,
)
from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.localization.plane_graphs import (
    P2FixedLocusGraph,
    P2FixedLocusGraphEnumerator,
    _fl_weak_compositions,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
)

__all__ = ['ToricCompleteIntersection', 'ToricCompleteIntersectionCycles', 'ToricProductAmbient', 'hirzebruch_index', 'toric_degree_tuples', 'toric_fixed_graphs', 'toric_genus_zero_number', 'toric_known_ordinary_gw', 'toric_primitive_curve_counts', 'toric_sample']


def _toric_embed(values,offset,length):
    return tuple([ZZ(0)]*offset+list(values)+[ZZ(0)]*(length-offset-len(values)))

def hirzebruch_index(factor):
    if isinstance(factor,str) and factor.startswith('F') and factor[1:].isdigit():return ZZ(factor[1:])
    return None

class ToricProductAmbient:
    """Explicit GKM data and a separate ordinary Chow-ring calculation."""
    def __init__(self,factors):
        self.factors=tuple(factors)  # positive integers denote P^n; 'Fe' denotes F_e
        if not self.factors or any(hirzebruch_index(f) is None and (f not in ZZ or f<1) for f in self.factors):
            raise ValueError('factors must be positive projective dimensions or Fe, e >= 0')
        self.dimension=sum(2 if hirzebruch_index(f) is not None else f for f in self.factors)
        if self.dimension>4:raise ValueError('this bounded module supports ambient dimension at most four')
        self.picard_rank=sum(2 if hirzebruch_index(f) is not None else 1 for f in self.factors)
        self.slug='x'.join(f.lower() if hirzebruch_index(f) is not None else 'p%s'%f for f in self.factors)
        self.points=tuple(product(*(range(4 if hirzebruch_index(f) is not None else f+1) for f in self.factors)))
        self.point_index={p:i for i,p in enumerate(self.points)}
        self.edges={};self.tangent=[[] for p in self.points]
        self._local=[];toffset=0;poffset=0;antican=[]
        for factor in self.factors:
            if hirzebruch_index(factor) is not None:
                e=hirzebruch_index(factor)
                tau={(0,1):(1,0),(1,0):(-1,0),(0,2):(0,1),(2,0):(0,-1),
                     (1,3):(e,1),(3,1):(-e,-1),(2,3):(1,0),(3,2):(-1,0)}
                classes={(0,1):(1,0),(0,2):(0,1),(1,3):(0,1),(2,3):(1,e)}
                normal={(0,1):(-e,(0,1)),(1,0):(-e,(e,1)),
                        (0,2):(0,(1,0)),(2,0):(0,(1,0)),
                        (1,3):(0,(-1,0)),(3,1):(0,(-1,0)),
                        (2,3):(e,(0,-1)),(3,2):(e,(-e,-1))}
                local={'tau':tau,'classes':classes,'normal':{e:(n,) for e,n in normal.items()},
                       'tdim':2,'pdim':2};antican.extend((2,e+2))
            else:
                aa=[tuple([ZZ(0)]*factor)]+[tuple(ZZ(i==j) for j in range(factor)) for i in range(factor)]
                tau={(i,j):tuple(a-b for a,b in zip(aa[i],aa[j])) for i in range(factor+1) for j in range(factor+1) if i!=j}
                local={'tau':tau,'classes':{(i,j):(1,) for i in range(factor+1) for j in range(i+1,factor+1)},
                    'normal':{(i,j):tuple((1,tau[i,k]) for k in range(factor+1) if k not in (i,j)) for i,j in tau},
                    'aa':aa,'tdim':factor,'pdim':1};antican.append(factor+1)
            local.update(toffset=toffset,poffset=poffset)
            self._local.append(local);toffset+=local['tdim'];poffset+=local['pdim']
        self.anticanonical=tuple(antican)
        for i,p in enumerate(self.points):
            for k,local in enumerate(self._local):
                for (left,right),w in local['tau'].items():
                    if left!=p[k]:continue
                    q=list(p);q[k]=right;j=self.point_index[tuple(q)]
                    weight=_toric_embed(w,local['toffset'],self.dimension)
                    beta=_toric_embed(local['classes'][tuple(sorted((left,right)))],local['poffset'],self.picard_rank)
                    normals=[(a,_toric_embed(ell,local['toffset'],self.dimension)) for a,ell in local['normal'][left,right]]
                    for kk,other in enumerate(self._local):
                        if kk==k:continue
                        normals.extend((0,_toric_embed(ww,other['toffset'],self.dimension))
                                       for (ll,rr),ww in other['tau'].items() if ll==p[kk])
                    self.edges[i,j]={'weight':weight,'normals':tuple(normals),'class':beta}
                    self.tangent[i].append(weight)
        assert all(len(t)==self.dimension for t in self.tangent)
        self._build_chow()

    def _build_chow(self):
        self.P=PolynomialRing(QQ,self.picard_rank,['h%s'%i for i in range(self.picard_rank)])
        variables=self.P.gens();ideal=[];top=self.P.one();ct=self.P.one();offset=0
        for f in self.factors:
            if hirzebruch_index(f) is not None:
                e=hirzebruch_index(f)
                s,b=variables[offset:offset+2];ideal.extend((s*s+e*s*b,b*b));top*=s*b
                ct*=(1+s)*(1+s+e*b)*(1+b)**2;offset+=2
            else:
                h=variables[offset];ideal.append(h**(f+1));top*=h**f;ct*=(1+h)**(f+1);offset+=1
        self.chow=self.P.quotient(self.P.ideal(ideal));self.top_monomial=top
        self.total_chern=self.chow(ct)

    def divisor(self,coefficients):
        if len(coefficients)!=self.picard_rank:raise ValueError('wrong divisor rank')
        return self.chow(sum(a*h for a,h in zip(coefficients,self.P.gens())))

    def homogeneous(self,poly,degree):
        p=self.chow(poly).lift()
        return self.chow(self.P({powers:co for powers,co in p.dict().items() if sum(powers)==degree}))

    def integral(self,poly):
        return QQ(self.chow(poly).lift().monomial_coefficient(self.top_monomial))

    def pairing(self,divisor,beta):
        if len(divisor)!=self.picard_rank or len(beta)!=self.picard_rank:raise ValueError('wrong divisor or curve rank')
        ans=ZZ(0);offset=0
        for f in self.factors:
            if hirzebruch_index(f) is not None:
                a,b=divisor[offset:offset+2];s,v=beta[offset:offset+2]
                ans+=(b-hirzebruch_index(f)*a)*s+a*v;offset+=2
            else:ans+=divisor[offset]*beta[offset];offset+=1
        return ans

    def fibre_character(self,divisor,point):
        answer=vector(ZZ,self.dimension)
        for f,loc,p in zip(self.factors,self._local,self.points[point]):
            off=loc['poffset']
            if hirzebruch_index(f) is not None:
                a,b=divisor[off:off+2];chars=((0,0),(hirzebruch_index(f)*a-b,0),(0,-a),(-b,-a));char=chars[p]
            else:char=tuple(divisor[off]*x for x in loc['aa'][p])
            answer+=vector(ZZ,_toric_embed(char,loc['toffset'],self.dimension))
        return tuple(answer)

    def character_value(self,char,weights):
        return sum(c*w for c,w in zip(char,weights))

    def ray(self,base=4):
        weights=tuple(QQ(base)**i for i in range(self.dimension))
        for (i,j),edge in self.edges.items():
            w=self.character_value(edge['weight'],weights)
            if not w:raise ValueError('resonant tangent weight')
            for e in (1,2):
                for a,ell in edge['normals']:
                    if a>=0 and any(self.character_value(ell,weights)-QQ(r)*w/e==0 for r in range(a*e+1)):
                        raise ValueError('resonant normal H0 weight')
        for tangent in self.tangent:
            values=[self.character_value(w,weights) for w in tangent]
            if any(u/e+v/f==0 for u in values for v in values for e in (1,2) for f in (1,2)):
                raise ValueError('resonant node flag sum')
        return weights

    def record(self):
        return {'slug':self.slug,'factors':list(self.factors),'dimension':int(self.dimension),
                'picard_rank':int(self.picard_rank),'fixed_points':len(self.points),
                'invariant_curves':len(self.edges)//2,'anticanonical':list(self.anticanonical)}

class ToricCompleteIntersection:
    def __init__(self,ambient,bundles,identifier=None):
        self.ambient=ambient;self.bundles=tuple(tuple(ZZ(x) for x in L) for L in bundles)
        self.rank=len(self.bundles);self.dimension=ambient.dimension-self.rank
        if not self.bundles or self.dimension not in (2,3):raise ValueError('this sample supports surface and threefold complete intersections')
        for L in self.bundles:
            if len(L)!=ambient.picard_rank:raise ValueError('wrong bundle Picard rank')
            if not any(L):raise ValueError('defining line bundles must be nontrivial')
            if any(ambient.pairing(L,e['class'])<0 for e in ambient.edges.values()):raise ValueError('twisting bundles must be nef')
            for (i,j),edge in ambient.edges.items():
                assert vector(ambient.fibre_character(L,j))==vector(ambient.fibre_character(L,i))-ambient.pairing(L,edge['class'])*vector(edge['weight'])
        self.c1=tuple(a-sum(L[i] for L in self.bundles) for i,a in enumerate(ambient.anticanonical))
        if all(x==0 for x in self.c1):
            if self.dimension==2 and not all(ambient.pairing(L,e['class'])>0 for L in self.bundles for e in ambient.edges.values()):
                raise NotImplementedError('CY surface classification here requires ample summands; c1=0 alone does not imply K3')
            self.kind='K3' if self.dimension==2 else 'CY3'
        elif all(ambient.pairing(self.c1,e['class'])>0 for e in ambient.edges.values()):self.kind='Fano'
        else:raise ValueError('the residual anticanonical divisor must be ample or zero')
        self.slug=identifier or ambient.slug+'-'+'.'.join('-'.join(map(str,L)) for L in self.bundles)
        ct=ambient.total_chern;self.fundamental=ambient.chow.one()
        for L in self.bundles:
            line=ambient.divisor(L);self.fundamental*=line
            ct*=sum((-line)**k for k in range(ambient.dimension+1))
        self.chern=tuple(ambient.homogeneous(ct,k) for k in range(self.dimension+1))

    def chern_numbers(self):
        integrate=lambda p:self.ambient.integral(self.fundamental*p)
        if self.dimension==2:return {'c1_squared':integrate(self.chern[1]**2),'euler':integrate(self.chern[2])}
        c12=integrate(self.chern[1]*self.chern[2]);c3=integrate(self.chern[3])
        return {'c1_cubed':integrate(self.chern[1]**3),'c1_c2':c12,'euler':c3,'C':c3-c12}

    def alpha(self,g,beta):return self.dimension*(g-1)-self.ambient.pairing(self.c1,beta)
    def t_power(self,g,beta,r):return self.alpha(g,beta)-r
    def comparison_sign(self,g,beta):return QQ(-1)**(self.rank*(1-g)+sum(self.ambient.pairing(L,beta) for L in self.bundles))

    def compact_degree_zero_terms(self,g):
        if g not in (2,3):raise NotImplementedError('implemented genera are two and three')
        c=self.chern_numbers()
        if self.dimension==2:
            return {('interior',None,(2,) if g==2 else (1,3)):c['c1_squared']} if c['c1_squared'] else {}
        return {('interior',None,(1,2) if g==2 else (1,2,3)):(-1)**g*c['C']}

    def record(self):
        return {'id':self.slug,'ambient':self.ambient.record(),'bundles':[list(L) for L in self.bundles],
                'dimension':int(self.dimension),'kind':self.kind,'c1':list(self.c1),
                'chern_numbers':{k:str(v) for k,v in self.chern_numbers().items()}}

def toric_sample():
    A=ToricProductAmbient((1,2));B=ToricProductAmbient((1,3));C=ToricProductAmbient((2,2))
    D=ToricProductAmbient(('F1',1));E=ToricProductAmbient(('F1',2))
    data=((A,((1,1),)),(A,((2,3),)),(B,((1,2),)),(B,((2,4),)),
          (C,((1,2),)),(C,((3,3),)),(C,((1,1),(1,1))),(C,((1,1),(2,2))),
          (D,((1,1,1),)),(D,((2,3,2),)),(E,((1,1,1),)),(E,((2,3,3),)))
    return tuple(ToricCompleteIntersection(a,b,'T%s'%(i+1)) for i,(a,b) in enumerate(data))

def toric_degree_tuples(ambient,bound=2):
    return tuple(beta for d in range(bound+1) for beta in _fl_weak_compositions(d,ambient.picard_rank))

_TORIC_GRAPH_CACHE={}

def toric_fixed_graphs(ambient,genus,bound=2):
    if bound not in (1,2):raise NotImplementedError('the bounded toric search supports Mori coefficient sum at most two')
    key=(ambient.slug,int(genus),bound)
    if key in _TORIC_GRAPH_CACHE:return _TORIC_GRAPH_CACHE[key]
    class ToricGraph(P2FixedLocusGraph):
        fixed_point_count=len(ambient.points)
        def flag_weight(self,edge_index,vertex,lambdas):
            other=self.other_endpoint(edge_index,vertex);e=self.edges[edge_index][2]
            char=ambient.edges[self.vertex_colours[vertex],self.vertex_colours[other]]['weight']
            return ambient.character_value(char,lambdas)/e
        def mori_degree(self):
            beta=vector(ZZ,ambient.picard_rank)
            for i,j,e in self.edges:beta+=e*vector(ambient.edges[self.vertex_colours[i],self.vertex_colours[j]]['class'])
            return tuple(beta)
    class ToricEnumerator(P2FixedLocusGraphEnumerator):
        fixed_point_count=len(ambient.points);graph_type=ToricGraph
        @classmethod
        def _proper_colours(cls,vertex_count,endpoints):
            for colours in product(range(cls.fixed_point_count),repeat=vertex_count):
                if all((colours[i],colours[j]) in ambient.edges for i,j in endpoints):yield colours
    result={b:[] for b in toric_degree_tuples(ambient,bound)}
    for cover_degree in range(bound+1):
        for graph in ToricEnumerator(genus,0,cover_degree).graphs():
            beta=graph.mori_degree()
            if beta in result:result[beta].append(graph)
    result={b:tuple(gs) for b,gs in result.items()};_TORIC_GRAPH_CACHE[key]=result
    return result

def _toric_edge_factor(ambient,bundles,i,j,e,weights,q):
    edge=ambient.edges[i,j];w=ambient.character_value(edge['weight'],weights)
    factor=(-1)**e*e**(2*e)/(factorial(e)**2*w**(2*e))
    for a,ell in edge['normals']:
        ell=ambient.character_value(ell,weights);ae=a*e
        if ae>=0:factor/=prod(ell-QQ(r)*w/e for r in range(ae+1))
        else:factor*=prod(ell-QQ(r)*w/e for r in range(ae+1,0))
    for k,L in enumerate(bundles):
        factor*=prod(q[i][k]+QQ(r)*w/e for r in range(e*ambient.pairing(L,edge['class'])+1))
    return factor

class ToricCompleteIntersectionCycles(QuarticK3GraphCycles):
    def __init__(self,target,genus=3,base=4):
        if genus not in (2,3):raise NotImplementedError('stabilized cycles support genera two and three')
        self.target=target;self.ambient=target.ambient;self.genus=ZZ(genus);self.dimension=3*genus-3
        self.weights=self.ambient.ray(base)
        self.P=PolynomialRing(QQ,'epsilon');self.epsilon=self.P.gen();self.F=self.P.fraction_field()
        self.a=tuple(self.F(w*self.epsilon) for w in self.weights)
        self.tangent=tuple(tuple(self.ambient.character_value(w,self.a) for w in ws) for ws in self.ambient.tangent)
        self.q=tuple(tuple(1-self.ambient.character_value(self.ambient.fibre_character(L,i),self.a)
                           for L in target.bundles) for i in range(len(self.ambient.points)))
        self.adm=_load_vendored_admcycles();self.R=self.adm.TautologicalRing(genus,0)
        self.hodge=AdmcyclesHodgeIntegralBackend()

    def graphs(self,beta):
        beta=tuple(ZZ(x) for x in beta)
        if len(beta)!=self.ambient.picard_rank or any(x<0 for x in beta) or sum(beta)>2:
            raise ValueError('expected a nonnegative Mori tuple of coefficient sum at most two')
        return toric_fixed_graphs(self.ambient,self.genus)[beta]

    @lru_cache(maxsize=None)
    def _hodge_vertex(self,genus,colour,flags):
        dimension=3*genus-3+flags;answer={():self.F.one()}
        for w in self.tangent[colour]+self.q[colour]:
            factor={(() if ell==0 else (ell,)):(-1)**ell*w**(genus-ell) for ell in range(genus+1)}
            answer=self._hodge_product(genus,answer,factor,dimension)
        scalar=prod(self.tangent[colour])**(flags-1)*prod(q**(1-flags-2*genus) for q in self.q[colour])
        return {indices:self.F(co*scalar) for indices,co in answer.items()}

    @lru_cache(maxsize=None)
    def _general_edge(self,i,j,degree):
        return self.F(_toric_edge_factor(self.ambient,self.target.bundles,i,j,degree,self.a,self.q))

    def _nodal_factor(self,colour,weights):
        if not sum(weights):raise ValueError('resonant node')
        return prod(self.tangent[colour])/(sum(weights)*prod(self.q[colour]))

    @lru_cache(maxsize=None)
    def rational_terms(self,beta):
        self.graphs(beta)  # validate and populate all Mori sectors once
        if any(beta):return super().rational_terms(beta)
        result={}
        for i in range(len(self.ambient.points)):
            for indices,co in self._hodge_vertex(self.genus,i,0).items():
                key=('interior',None,indices);result[key]=result.get(key,self.F.zero())+co
        return {key:self.F(co) for key,co in result.items() if co}

    def top_integral(self,beta):
        if not any(beta):return QuarticK3Cycles.top_integral(self,beta)
        return super().top_integral(beta)

    def aggregate_relation(self,beta,codimension):
        power=self.target.t_power(self.genus,beta,codimension)
        if codimension not in range(self.dimension+1) or power>0:raise ValueError('expected a nonpositive t power in the Chow range')
        zero=self.cycle(beta,codimension);known=power<0 or self.target.kind=='K3' or not any(beta)
        compact=(self.R.zero() if power<0 or self.target.kind=='K3' else
                 self.class_from_terms(self.target.compact_degree_zero_terms(self.genus)) if not any(beta) else None)
        sign=self.target.comparison_sign(self.genus,beta)
        return {'genus':self.genus,'beta':tuple(beta),'codimension':codimension,'t_power':power,
                'zero_vertex':zero,'compact_known':known,'compact_class':compact,'comparison_sign':sign,
                'infinity_aggregate':sign*compact-zero if known else None}

def toric_genus_zero_number(target,beta,base=4):
    """Ordinary CY3 genus-zero numbers by convex Lefschetz, for the bounded tuples."""
    if target.kind!='CY3' or not any(beta):raise ValueError('expected a positive CY3 class')
    X=target.ambient;weights=X.ray(base)
    # Generic character shifts avoid zero node fibres at t=0. The ordinary
    # invariant is independent of these choices of auxiliary linearization.
    q=tuple(tuple(-X.character_value(X.fibre_character(L,i),weights)-101*(k+1)
                  for k,L in enumerate(target.bundles)) for i in range(len(X.points)))
    total=QQ.zero()
    for graph in toric_fixed_graphs(X,0)[tuple(beta)]:
        co=QQ(graph.localization_weight())
        for i,j,e in graph.edges:co*=_toric_edge_factor(X,target.bundles,graph.vertex_colours[i],graph.vertex_colours[j],e,weights,q)
        for v in range(graph.vertex_count):
            colour=graph.vertex_colours[v]
            ws=[graph.flag_weight(e,v,weights) for e in graph.incident_edges(v)]
            if len(ws)==1:co*=ws[0]
            elif len(ws)==2:
                co*=prod(X.character_value(w,weights) for w in X.tangent[colour])/(sum(ws)*prod(q[colour]))
            else:raise NotImplementedError('genus-zero vertices of valence three require psi integration')
        total+=co
    return target.comparison_sign(0,beta)*total

def toric_primitive_curve_counts(target):
    """Independent parameter-space top-Chern counts; see the results report."""
    geometry=(target.ambient.slug,target.bundles)
    if geometry==('p1xp3',((2,4),)):
        P=PolynomialRing(QQ,names=('x','y','h'));x,y,h=P.gens()
        lines=((x-y)*prod(2*h+(4-r)*x+r*y for r in range(5))).monomial_coefficient(x**3*y**2*h)
        return {(1,0):ZZ(4**3),(0,1):ZZ(lines)}
    if geometry==('p2xp2',((3,3),)):
        return {(1,0):ZZ(21*9),(0,1):ZZ(21*9)}
    if geometry==('f1xp2',((2,3,3),)):
        return {(1,0,0):ZZ(3**2),(0,1,0):ZZ(9*(1+2+3)),(0,0,1):ZZ(21*8)}
    return {}

def toric_known_ordinary_gw(target,g,beta):
    """Proved CY3 scalar values only; other positive degrees return None."""
    if target.kind!='CY3' or g not in (2,3):raise ValueError('expected a CY3 and genus two or three')
    if len(beta)!=target.ambient.picard_rank or any(x<0 or x not in ZZ for x in beta):raise ValueError('invalid Mori tuple')
    if not any(beta):return target.chern_numbers()['euler']*(QQ(1)/5760 if g==2 else -QQ(1)/1451520)
    counts=toric_primitive_curve_counts(target)
    if tuple(beta) in counts:return QQ(counts[tuple(beta)])/(240 if g==2 else 6048)
    geometry=(target.ambient.slug,target.bundles)
    pure_double=(geometry==('p1xp3',((2,4),)) and tuple(beta)==(2,0)) or (
        geometry==('f1xp2',((2,3,3),)) and tuple(beta) in ((2,0,0),(0,2,0)))
    if pure_double:
        primitive=tuple(ZZ(x/2) for x in beta)
        return QQ(counts[primitive])*2**(2*g-3)/(240 if g==2 else 6048)
    return None
