r"""Dual-Euler localization and CJR aggregates for projective complete intersections.

All summands of E = direct_sum O(d_j) have the same CJR parameter t.
Compact positive-degree GW classes remain unknown; K3 vanishing is not
inherited. Genus-two/three curve stabilization is shared with the quartic
backend. See docs/PROJECTIVE_COMPLETE_INTERSECTION_SOURCES.md."""

from sage.all import *

from functools import lru_cache
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
)

__all__ = ['CompleteIntersectionCycles', 'ProjectiveCompleteIntersection', 'complete_intersection_weights', 'cy3_line_count_grassmannian', 'cy3_line_count_localization', 'cy3_ordinary_gw_number', 'requested_complete_intersections']


class ProjectiveCompleteIntersection:
    def __init__(self,ambient_dimension,degrees,label=None):
        self.N=ZZ(ambient_dimension)
        self.degrees=tuple(sorted(ZZ(d) for d in degrees))
        if self.N<1 or not self.degrees or any(d<1 for d in self.degrees):
            raise ValueError('positive ambient dimension and positive defining degrees are required')
        self.rank=len(self.degrees)
        self.dimension=self.N-self.rank
        if self.dimension<1:raise ValueError('the complete intersection must have positive dimension')
        self.index=self.N+1-sum(self.degrees)
        self.degree=prod(self.degrees)
        self.label=label or 'P%s-(%s)'%(self.N,','.join(map(str,self.degrees)))
        self.slug='p%s-'%self.N+'-'.join(map(str,self.degrees))

    def virtual_dimension(self,genus,degree):
        return (1-genus)*(self.dimension-3)+self.index*degree

    def alpha(self,genus,degree):
        return self.dimension*(genus-1)-self.index*degree

    def t_power(self,genus,degree,codimension):
        return self.alpha(genus,degree)-codimension

    def comparison_sign(self,genus,degree):
        return QQ(-1)**(self.rank*(1-genus)+sum(self.degrees)*degree)

    def chern_coefficients(self):
        ring=PowerSeriesRing(QQ,'H',default_prec=self.dimension+1)
        H=ring.gen()
        series=(1+H)**(self.N+1)/prod(1+d*H for d in self.degrees)
        return tuple(ZZ(series[k]) for k in range(self.dimension+1))

    def euler_characteristic(self):
        return self.degree*self.chern_coefficients()[-1]

    def compact_degree_zero_terms(self,genus):
        if genus not in (2,3):raise NotImplementedError('constant-map class formulas are implemented in genera 2 and 3')
        ch=self.chern_coefficients()
        if self.dimension==2:
            indices=(2,) if genus==2 else (1,3)
            return {('interior',None,indices):self.degree*ch[1]**2}
        if self.dimension==3:
            indices=(1,2) if genus==2 else (1,2,3)
            return {('interior',None,indices):(-1)**genus*self.degree*(ch[3]-ch[1]*ch[2])}
        if self.dimension>=4:return {}
        raise NotImplementedError('curve-target constants are outside this module')

    def record(self):
        return {'name':self.label,'slug':self.slug,'ambient_dimension':int(self.N),
                'degrees':[int(d) for d in self.degrees],'dimension':int(self.dimension),
                'index':int(self.index),'degree':int(self.degree),
                'chern_coefficients':[int(x) for x in self.chern_coefficients()],
                'euler_characteristic':int(self.euler_characteristic())}

def requested_complete_intersections():
    names={2:'quadric',3:'cubic',4:'quartic',5:'quintic'}
    dimensions={2:'surface',3:'threefold',4:'fourfold'}
    targets=[ProjectiveCompleteIntersection(m+1,(d,),names[d]+' '+dimensions[m])
             for m in (2,3,4) for d in range(2,m+2)]
    targets += [ProjectiveCompleteIntersection(N,degrees,'CY3 '+str(degrees))
                for N,degrees in ((4,(5,)),(5,(2,4)),(5,(3,3)),
                                  (6,(2,2,3)),(7,(2,2,2,2)))]
    return tuple(targets)

def complete_intersection_weights(number,base=4):
    """Small nonresonant weights for degree at most three, using base-b digits 0,1."""
    return tuple(ZZ(sum(((i>>k)&1)*base**k for k in range(i.bit_length()))) for i in range(number))

class CompleteIntersectionCycles(QuarticK3GraphCycles):
    _graph_cache={}

    def __init__(self,target,genus=3,base_weights=None):
        if genus not in (2,3):raise NotImplementedError('cycle stabilization supports genera two and three')
        self.target=target
        self.genus=ZZ(genus);self.dimension=3*self.genus-3
        self.weights=tuple(QQ(a) for a in (base_weights if base_weights is not None
            else complete_intersection_weights(target.N+1)))
        if len(self.weights)!=target.N+1 or len(set(self.weights))!=target.N+1:
            raise ValueError('one distinct rational weight is required for every projective fixed point')
        self.P=PolynomialRing(QQ,'epsilon');self.epsilon=self.P.gen();self.F=self.P.fraction_field()
        self.a=tuple(self.F(a*self.epsilon) for a in self.weights)
        self.adm=_load_vendored_admcycles();self.R=self.adm.TautologicalRing(genus,0)
        self.hodge=AdmcyclesHodgeIntegralBackend()

    @lru_cache(maxsize=None)
    def _hodge_vertex(self,genus,colour,flags):
        dimension=3*genus-3+flags
        answer={():self.F.one()}
        tangent=[self.a[colour]-self.a[j] for j in range(self.target.N+1) if j!=colour]
        fibres=[1-d*self.a[colour] for d in self.target.degrees]
        for w in tangent+fibres:
            factor={(() if ell==0 else (ell,)):(-1)**ell*w**(genus-ell)
                    for ell in range(genus+1)}
            answer=self._hodge_product(genus,answer,factor,dimension)
        scalar=prod(tangent)**(flags-1)*prod(q**(1-flags-2*genus) for q in fibres)
        return {indices:self.F(co*scalar) for indices,co in answer.items()}

    def graphs(self,degree):
        if degree not in (0,1,2,3):raise NotImplementedError('the scan supports ambient degrees zero through three')
        key=(int(self.target.N),int(self.genus),int(degree))
        if key not in self._graph_cache:
            graph_type=type('ProjectiveFixedLocusGraph',(P2FixedLocusGraph,),
                            {'fixed_point_count':int(self.target.N+1)})
            enumerator=type('ProjectiveFixedLocusGraphEnumerator',(P2FixedLocusGraphEnumerator,),
                           {'fixed_point_count':int(self.target.N+1),'graph_type':graph_type})
            self._graph_cache[key]=enumerator(self.genus,0,degree).graphs()
        return self._graph_cache[key]

    @lru_cache(maxsize=None)
    def _general_edge(self,i,j,degree):
        delta=self.a[i]-self.a[j]
        denominator=prod(QQ(m)/degree*self.a[i]+QQ(degree-m)/degree*self.a[j]-self.a[k]
                         for k in range(self.target.N+1) if k not in (i,j)
                         for m in range(degree+1))
        if not denominator:raise ValueError('resonant auxiliary weights on an edge')
        factor=(-1)**degree*degree**(2*degree)/(factorial(degree)**2*delta**(2*degree)*denominator)
        factor*=prod(1-QQ(d*degree-m)/degree*self.a[i]-QQ(m)/degree*self.a[j]
                     for d in self.target.degrees for m in range(d*degree+1))
        return self.F(factor)

    def _nodal_factor(self,colour,weights):
        if not sum(weights):raise ValueError('resonant auxiliary weights at a node')
        return prod(self.a[colour]-self.a[k] for k in range(self.target.N+1) if k!=colour)/(
            sum(weights)*prod(1-d*self.a[colour] for d in self.target.degrees))

    @lru_cache(maxsize=None)
    def rational_terms(self,degree):
        if degree:
            return super().rational_terms(degree)
        result={}
        for colour in range(self.target.N+1):
            for indices,co in self._hodge_vertex(self.genus,colour,0).items():
                key=('interior',None,indices)
                result[key]=result.get(key,self.F.zero())+co
        return {k:self.F(co) for k,co in result.items() if co}

    def aggregate_relation(self,degree,codimension):
        if codimension not in range(self.dimension+1):raise ValueError('codimension outside Mbar_g')
        power=self.target.t_power(self.genus,degree,codimension)
        if power>0:raise ValueError('positive t powers do not give automatic CJR relations')
        record={'genus':int(self.genus),'ambient_degree':int(degree),
                'codimension':int(codimension),'t_power':int(power),
                'zero_vertex':self.cycle(degree,codimension)}
        if power<0:
            record.update(compact_known=True,compact_class=self.R.zero(),comparison_sign=1)
        else:
            record.update(compact_known=(degree==0),
                compact_class=(self.class_from_terms(self.target.compact_degree_zero_terms(self.genus))
                               if degree==0 else None),
                comparison_sign=self.target.comparison_sign(self.genus,degree))
        record['infinity_aggregate']=(record['comparison_sign']*record['compact_class']-record['zero_vertex']
                                      if record['compact_known'] else None)
        return record

def cy3_line_count_localization(target,weights=None):
    """Degree-one genus-zero dual-Euler fixed edges, independently of Grassmannian integration."""
    if target.dimension!=3 or target.index!=0:
        raise ValueError('the unmarked line count here is a number only for CY3 targets')
    a=tuple(QQ(x) for x in (weights or range(target.N+1)))
    if len(a)!=target.N+1 or len(set(a))!=len(a):raise ValueError('distinct fixed-point weights are required')
    answer=QQ.zero()
    for i in range(target.N+1):
        for j in range(i+1,target.N+1):
            numerator=prod(-(d-u)*a[i]-u*a[j] for d in target.degrees for u in range(d+1))
            denominator=prod((a[i]-a[k])*(a[j]-a[k])
                             for k in range(target.N+1) if k not in (i,j))
            # The two unstable flag factors cancel the -Delta^2 tangent factor.
            answer+=numerator/denominator
    return answer

def cy3_line_count_grassmannian(target):
    if target.dimension!=3 or target.index!=0:raise ValueError('expected a CY3 complete intersection')
    ring=PolynomialRing(QQ,names=('x','y'));x,y=ring.gens()
    polynomial=(x-y)*prod(u*x+(d-u)*y for d in target.degrees for u in range(d+1))
    return polynomial.monomial_coefficient(x**target.N*y**(target.N-1))

def cy3_ordinary_gw_number(target,genus,degree):
    """Ordinary connected GW numbers in degrees 0 and 1, via FP Theorems 3/4.

    A general member has transverse isolated lines, with normal bundle O(-1)^2.
    Deformation invariance extends their local contributions to every smooth
    member. This function is separate from ambient dual-Euler localization.
    """
    if target.dimension!=3 or target.index!=0:
        raise ValueError('expected a Calabi-Yau threefold complete intersection')
    if genus not in (2,3):raise NotImplementedError('implemented genera are two and three')
    if degree==0:
        return target.euler_characteristic()*(QQ(1)/5760 if genus==2 else -QQ(1)/1451520)
    if degree==1:
        return cy3_line_count_grassmannian(target)/(240 if genus==2 else 6048)
    raise NotImplementedError('ordinary positive-degree invariants beyond degree one are not supplied')
