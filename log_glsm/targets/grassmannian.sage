r"""Gr(2,n) Pluecker complete intersections with exact GKM localization.

The bounded samples use n=4,5,7 and positive line-bundle summands. The GKM
backend and stabilization are shared with the toric module; no toric fan or
higher-genus quantum Lefschetz identification is assumed."""

from sage.all import *

from functools import lru_cache
from itertools import combinations
from itertools import product
from log_glsm.targets.toric import (
    ToricCompleteIntersection,
    ToricCompleteIntersectionCycles,
    ToricProductAmbient,
    toric_genus_zero_number,
)

__all__ = ['GrassmannianAmbient', 'GrassmannianCompleteIntersection', 'GrassmannianCompleteIntersectionCycles', 'grassmannian_known_ordinary_gw', 'grassmannian_line_count', 'grassmannian_sample', 'grassmannian_seven_hyperplane_cy3']


class GrassmannianAmbient:
    def __init__(self,n):
        if n not in (3,4,5,7):raise ValueError('implemented ambients are Gr(2,n), n in {3,4,5,7}')
        n=ZZ(n)
        self.n=ZZ(n);self.dimension=2*(n-2);self.picard_rank=1
        self.slug='gr2-%s'%n;self.anticanonical=(ZZ(n),)
        self.points=tuple(combinations(range(n),2));self.point_index={p:i for i,p in enumerate(self.points)}
        basis=identity_matrix(ZZ,n)
        weight=lambda i,j:tuple(basis[j]-basis[i])
        self.tangent=[];self.edges={}
        for p,I in enumerate(self.points):
            outside=tuple(j for j in range(n) if j not in I)
            self.tangent.append(tuple(weight(i,j) for i in I for j in outside))
            for a in I:
                for b in outside:
                    J=tuple(sorted((set(I)-{a})|{b}));q=self.point_index[J]
                    normals=[]
                    for i in I:
                        for j in outside:
                            if (i,j)==(a,b):continue
                            normals.append((ZZ(i==a or j==b),weight(i,j)))
                    self.edges[p,q]={'weight':weight(a,b),'normals':tuple(normals),'class':(ZZ(1),)}
        self._build_chow()

    def _build_chow(self):
        # Chern-root splitting algebra. Only symmetric classes are used.
        # Integration is the Grassmannian coefficient formula, not the
        # degree of the top monomial in this larger splitting algebra.
        self.P=PolynomialRing(QQ,2,('x','y'));x,y=self.P.gens()
        h=lambda d:sum(x**i*y**(d-i) for i in range(d+1))
        self.chow=self.P.quotient(self.P.ideal([h(self.n-1),h(self.n)]))
        self.total_chern=self.chow((1+x)**self.n*(1+y)**self.n)*sum(
            self.chow((x-y)**(2*j)) for j in range(self.dimension//2+1))

    def divisor(self,coefficients):
        if len(coefficients)!=1:raise ValueError('a Pluecker divisor has one coefficient')
        return self.chow(coefficients[0]*sum(self.P.gens()))

    homogeneous=ToricProductAmbient.homogeneous

    def integral(self,poly):
        x,y=self.P.gens();p=self.homogeneous(poly,self.dimension).lift()
        return QQ(((x-y)*p).monomial_coefficient(x**(self.n-1)*y**(self.n-2)))

    def pairing(self,divisor,beta):
        if len(divisor)!=1 or len(beta)!=1:raise ValueError('expected a single Pluecker degree')
        return divisor[0]*beta[0]

    def fibre_character(self,divisor,point):
        if len(divisor)!=1:raise ValueError('expected a single Pluecker coefficient')
        return tuple(-divisor[0]*ZZ(i in self.points[point]) for i in range(self.n))

    character_value=ToricProductAmbient.character_value

    def ray(self,base=4):
        if base<4:raise ValueError('use a base of at least four for the degree-two bound')
        weights=tuple(QQ(base)**i for i in range(self.n))
        for edge in self.edges.values():
            w=self.character_value(edge['weight'],weights)
            if not w:raise ValueError('resonant tangent weight')
            for e in (1,2):
                for a,ell in edge['normals']:
                    if any(self.character_value(ell,weights)-QQ(r)*w/e==0 for r in range(a*e+1)):
                        raise ValueError('resonant normal weight')
        for tangent in self.tangent:
            values=[self.character_value(w,weights) for w in tangent]
            if any(u/e+v/f==0 for u in values for v in values for e in (1,2) for f in (1,2)):
                raise ValueError('resonant node flag sum')
        return weights

    def record(self):
        return {'slug':self.slug,'k':2,'n':int(self.n),'dimension':int(self.dimension),
                'picard_rank':1,'fixed_points':len(self.points),'invariant_curves':len(self.edges)//2,
                'anticanonical':list(self.anticanonical)}

class GrassmannianCompleteIntersection(ToricCompleteIntersection):
    def __init__(self,n,degrees,identifier=None):
        self.degrees=tuple(sorted(ZZ(d) for d in degrees))
        if not self.degrees or any(d<1 for d in self.degrees):raise ValueError('positive Pluecker degrees are required')
        super().__init__(GrassmannianAmbient(n),tuple((d,) for d in self.degrees),identifier)

    def record(self):
        record=super().record();record['degrees']=list(self.degrees)
        record['pluecker_degree']=int(self.ambient.integral(self.fundamental*self.ambient.divisor((1,))**self.dimension))
        return record

def grassmannian_sample():
    data=[(4,d) for d in ((1,),(2,),(3,),(4,),(1,1),(1,2),(1,3),(2,2))]
    data += [(5,d) for d in ((1,1,1),(1,1,2),(1,1,3),(1,2,2),(1,1,1,1),(1,1,1,2))]
    return tuple(GrassmannianCompleteIntersection(n,d,'G%s'%(i+1)) for i,(n,d) in enumerate(data))

def grassmannian_seven_hyperplane_cy3():
    """The degree-42 linear section Gr(2,7) intersect P^13."""
    return GrassmannianCompleteIntersection(7,(1,)*7,'G15')

class GrassmannianCompleteIntersectionCycles(ToricCompleteIntersectionCycles):
    """The GKM edge data, rather than a toric hypothesis, drive this backend."""

    def __init__(self,target,genus=3,base=4,*,graphwise_series=None):
        super().__init__(target,genus,base)
        self.graphwise_series=(target.ambient.n==7 if graphwise_series is None else bool(graphwise_series))

    @lru_cache(maxsize=None)
    def _bounded_coefficients(self,beta):
        """Sum exact Laurent coefficients before combining rational functions.

        Only powers <= 0 are requested by the bounded relation search.
        Expanding a single graph summand avoids a common denominator across
        the 21 fixed points. Every quotient has unit constant denominator
        after removing its valuation, so its needed finite jet is exact.
        """
        result={r:{} for r in range(self.dimension+1)}
        S=PowerSeriesRing(QQ,'z');maximum_jet=0
        for graph in self.graphs(beta):
            scalar=self.F(graph.localization_weight())
            for left,right,d in graph.edges:
                scalar*=self._general_edge(graph.vertex_colours[left],graph.vertex_colours[right],d)
            stable,data=self._stabilized_graph(graph);local=[]
            for v in range(graph.vertex_count):
                colour=graph.vertex_colours[v]
                weights=tuple(graph.flag_weight(e,v,self.a) for e in graph.incident_edges(v))
                kind=graph.vertex_type(v)
                if kind=='stable':local.append(self._general_vertex(graph.vertex_genera[v],colour,weights))
                elif kind=='univalent':scalar*=weights[0]
                elif kind=='nodal':scalar*=self._nodal_factor(colour,weights)
                else:raise ValueError('unsupported unstable vertex')
            for entries in product(*(p.items() for p in local)):
                powers=tuple(item[0][0] for item in entries)
                lambdas=tuple(item[0][1] for item in entries);canonical=[]
                for flags,ps in zip(data[1],powers):
                    ps=list(ps);slots=[i for i,flag in enumerate(flags) if flag<=data[3]]
                    for i,value in zip(slots,sorted(ps[i] for i in slots)):ps[i]=value
                    canonical.append(tuple(ps))
                key=('graph',data,tuple(canonical),lambdas);r=self.codimension(key)
                if not 0<=r<=self.dimension:continue
                co=scalar*prod(item[1] for item in entries)
                if not co:continue
                numerator,denominator=co.numerator(),co.denominator()
                nv,dv=numerator.valuation(),denominator.valuation();valuation=nv-dv
                if valuation>0:continue
                precision=1-valuation;maximum_jet=max(maximum_jet,precision)
                a=S(numerator.shift(-nv)).add_bigoh(precision)
                b=S(denominator.shift(-dv)).add_bigoh(precision)
                quotient=a/b
                if quotient.prec()<precision:raise ArithmeticError('insufficient graphwise Laurent precision')
                for j in range(precision):
                    value=quotient[j]
                    if value:
                        terms=result[r].setdefault(int(valuation+j),{})
                        terms[key]=terms.get(key,QQ.zero())+value
        self.maximum_laurent_jet=maximum_jet
        return {r:{power:{key:co for key,co in terms.items() if co}
                   for power,terms in powers.items() if any(terms.values())}
                for r,powers in result.items()}

    def coefficients(self,beta,codimension):
        if not self.graphwise_series or not any(beta):
            return super().coefficients(beta,codimension)
        if codimension not in range(self.dimension+1):raise ValueError('codimension outside Mbar_g')
        return self._bounded_coefficients(tuple(beta))[codimension]

def grassmannian_line_count(target,weights=None):
    """Independent integration on the flag variety F(1,3;n) of ambient lines.

    For a flag A subset B, the line is P(B/A), and its Pluecker O(1) is
    A^* tensor O_P(B/A)(1). The restriction bundle is
    (A^*)^d tensor Sym^d(B/A)^*. Its top Chern class is localized here
    directly on the flag variety, with no stable-map graph factors.
    """
    if target.kind!='CY3':raise ValueError('a numerical line count requires a CY threefold')
    n=target.ambient.n;a=tuple(QQ(x) for x in (weights if weights is not None else range(n)))
    if len(a)!=n or len(set(a))!=n:raise ValueError('one distinct weight per coordinate is required')
    total=QQ.zero()
    for i in range(n):
        for pair in combinations([j for j in range(n) if j!=i],2):
            j,k=pair;outside=[h for h in range(n) if h not in (i,j,k)]
            numerator=prod(-d*a[i]-r*a[j]-(d-r)*a[k] for d in target.degrees for r in range(d+1))
            denominator=(a[j]-a[i])*(a[k]-a[i])*prod(a[h]-a[v] for h in outside for v in (i,j,k))
            total+=numerator/denominator
    return total

def grassmannian_known_ordinary_gw(target,g,beta):
    """Constants and CY3 lines; also conics for the audited Gr(2,7)[1^7]."""
    if g not in (2,3):raise ValueError('implemented genera are two and three')
    if tuple(beta)==(0,) and target.kind=='CY3':
        chi=target.chern_numbers()['euler']
        return chi/5760 if g==2 else -chi/1451520
    if tuple(beta)==(1,) and target.kind=='CY3':
        return grassmannian_line_count(target)/(240 if g==2 else 6048)
    if (tuple(beta)==(2,) and target.kind=='CY3' and target.ambient.n==7
            and target.degrees==(1,)*7):
        # The incidence and normal-obstruction proof is recorded in
        # docs/GRASSMANNIAN_27_SOURCES.md. It excludes distinct meeting
        # lines and proves that the conics are isolated (-1,-1) curves.
        lines=grassmannian_line_count(target)
        conics=toric_genus_zero_number(target,(2,))-lines/8
        return (conics+2**(2*g-3)*lines)/(240 if g==2 else 6048)
    if target.kind=='K3':return QQ.zero()
    return None
