r"""Ordinary quartic-K3 CJR aggregates and class-valued O(4)/P3 localization.

The cycle backend supports unmarked genera two and three, ambient degrees
zero and one. It does not identify individual punctured effective cycles.
The ordinary K3 compact class is zero in this range, including degree zero;
the reduced K3 theory is a different theory and is deliberately rejected.

Conventions and independent formulas: docs/QUARTIC_K3_SOURCES.md.
Only the target-independent Hodge/admcycles backend is shared with the cubic."""

from sage.all import *

from functools import lru_cache
from log_glsm.core.hodge import (
    AdmcyclesHodgeIntegralBackend,
    HodgeIntegralRequest,
    _load_vendored_admcycles,
)

__all__ = ['QuarticK3Cycles', 'QuarticK3Theory']


class QuarticK3Theory:
    ambient_dimension = 3
    hypersurface_degree = 4

    @staticmethod
    def compact_virtual_class(genus, degree, reduced=False):
        if reduced:
            raise NotImplementedError('reduced K3 theory is not ordinary K3 theory')
        if genus < 2 or degree < 0:
            raise ValueError('this unmarked module requires genus >= 2 and degree >= 0')
        return QQ.zero()

    @staticmethod
    def t_power(genus, codimension):
        return ZZ(2*genus-2-codimension)

    @staticmethod
    def relation_codimensions(genus):
        if genus < 2:
            raise ValueError('unmarked stable curves require genus >= 2')
        return tuple(range(2*genus-2,3*genus-2))

    @staticmethod
    def balanced(genus, degree, contacts):
        contacts=tuple(contacts)
        return (genus >= 0 and degree >= 0
                and all(c <= -1 and c in ZZ for c in contacts)
                and sum(c+1 for c in contacts) == 4*degree-(2*genus-2))

class QuarticK3Cycles:
    """Exact stabilized dual-Euler classes; auxiliary poles are retained."""

    def __init__(self, genus=3, base_weights=(0,1,3,7)):
        if genus not in (2,3):
            raise NotImplementedError('cycle stabilization currently supports genera 2 and 3')
        weights=tuple(QQ(a) for a in base_weights)
        if len(weights)!=4 or len(set(weights))!=4:
            raise ValueError('P3 requires four distinct rational auxiliary weights')
        self.genus=ZZ(genus)
        self.dimension=3*self.genus-3
        self.weights=weights
        self.P=PolynomialRing(QQ,'epsilon')
        self.epsilon=self.P.gen()
        self.F=self.P.fraction_field()
        self.a=tuple(self.F(a*self.epsilon) for a in weights)
        self.adm=_load_vendored_admcycles()
        self.R=self.adm.TautologicalRing(genus,0)
        self.hodge=AdmcyclesHodgeIntegralBackend()

    @staticmethod
    def _normalize(genus,indices):
        scale,request=HodgeIntegralRequest(genus,(),indices).low_genus_normal_form()
        return scale,request.lambda_indices

    def _hodge_product(self, genus, left, right, dimension):
        result={}
        for l,cl in left.items():
            for r,cr in right.items():
                if sum(l)+sum(r)>dimension:continue
                scale,indices=self._normalize(genus,l+r)
                if scale:
                    result[indices]=result.get(indices,self.F.zero())+scale*cl*cr
        return {k:c for k,c in result.items() if c}

    @lru_cache(maxsize=None)
    def _hodge_vertex(self, genus, colour, flags):
        dimension=3*genus-3+flags
        answer={():self.F.one()}
        q=1-4*self.a[colour]
        weights=[self.a[colour]-self.a[j] for j in range(4) if j!=colour]
        # Mumford: 1/Lambda_g^+(q) = Lambda_g^vee(q)/q^(2g).
        # This is an exact Chow identity, not a finite numerical approximation.
        for w in weights+[q]:
            factor={(() if ell==0 else (ell,)):(-1)**ell*w**(genus-ell)
                    for ell in range(genus+1)}
            answer=self._hodge_product(genus,answer,factor,dimension)
        scalar=q**(1-flags-2*genus)*prod(weights)**(flags-1)
        return {indices:self.F(co*scalar) for indices,co in answer.items()}

    @lru_cache(maxsize=None)
    def _vertex(self,genus,colour,other):
        w=self.a[colour]-self.a[other]
        dimension=3*genus-2
        return {(power,indices):co/w**(power+1)
                for indices,co in self._hodge_vertex(genus,colour,1).items()
                for power in range(dimension-sum(indices)+1)}

    @lru_cache(maxsize=None)
    def _edge(self,i,j):
        numerator=prod(1-(4-m)*self.a[i]-m*self.a[j] for m in range(5))
        denominator=(self.a[i]-self.a[j])**2*prod(
            (self.a[i]-self.a[k])*(self.a[j]-self.a[k])
            for k in range(4) if k not in (i,j))
        return -numerator/denominator

    @staticmethod
    def codimension(key):
        if key[0]=='interior':
            return (0 if key[1] is None else key[1])+sum(key[2])
        return 1+key[2]+sum(key[3])+key[4]+sum(key[5])

    @lru_cache(maxsize=None)
    def rational_terms(self,degree):
        if degree not in (0,1):
            raise NotImplementedError('ambient degrees greater than one require further fixed graphs')
        result={}
        def add(key,co):
            if self.codimension(key)<=self.dimension and co:
                result[key]=result.get(key,self.F.zero())+co
        if degree==0:
            for colour in range(4):
                for indices,co in self._hodge_vertex(self.genus,colour,0).items():
                    add(('interior',None,indices),co)
        else:
            for i in range(4):
                for j in range(4):
                    if i==j:continue
                    edge=self._edge(i,j)
                    for (power,indices),co in self._vertex(self.genus,i,j).items():
                        if power:
                            add(('interior',power-1,indices),
                                edge*(self.a[j]-self.a[i])*co)
                    # Ordered colours and all genus splits double-index bridges.
                    for left_genus in range(1,self.genus):
                        right_genus=self.genus-left_genus
                        for (p,l),cl in self._vertex(left_genus,i,j).items():
                            for (q,m),cr in self._vertex(right_genus,j,i).items():
                                if left_genus<=right_genus:
                                    key=('boundary',left_genus,p,l,q,m)
                                else:
                                    key=('boundary',right_genus,q,m,p,l)
                                add(key,edge*cl*cr/2)
        return {k:self.F(c) for k,c in result.items() if c}

    @lru_cache(maxsize=None)
    def coefficients(self,degree,codimension):
        if codimension not in range(self.dimension+1):
            raise ValueError('codimension outside the dimension of Mbar_g')
        result={}
        series_ring=LaurentSeriesRing(QQ,'epsilon',default_prec=32)
        for key,co in self.rational_terms(degree).items():
            if self.codimension(key)!=codimension:continue
            series=series_ring(co.numerator())/series_ring(co.denominator())
            if not series:continue
            if series.prec()<=0:raise ArithmeticError('insufficient Laurent precision')
            for power in range(int(series.valuation()),1):
                if series[power]:result.setdefault(int(power),{})[key]=QQ(series[power])
        return result

    @lru_cache(maxsize=None)
    def hodge_class(self,genus,markings,indices):
        ring=self.adm.TautologicalRing(genus,markings)
        return prod((ring.lambdaclass(i) for i in indices),ring.one())

    @lru_cache(maxsize=None)
    def monomial_class(self,key):
        if key[0]=='interior':
            result=self.hodge_class(self.genus,0,key[2])
            if key[1] is not None:result*=self.R.kappa(key[1])
            return result
        h=key[1]; k=self.genus-h
        left=self.adm.TautologicalRing(h,1)
        right=self.adm.TautologicalRing(k,1)
        graph=self.adm.StableGraph([h,k],[[1],[2]],[(1,2)])
        return graph.boundary_pushforward([
            left.psi(1)**key[2]*self.hodge_class(h,1,key[3]),
            right.psi(1)**key[4]*self.hodge_class(k,1,key[5])])

    def class_from_terms(self,terms):
        return sum((co*self.monomial_class(key) for key,co in terms.items()),self.R.zero())

    def cycle(self,degree,codimension,epsilon_power=0):
        return self.class_from_terms(self.coefficients(degree,codimension).get(epsilon_power,{}))

    def aggregate_relation(self,degree,codimension):
        if codimension not in QuarticK3Theory.relation_codimensions(self.genus):
            raise ValueError('positive t coefficients are not forced to vanish by ordinary K3 vanishing')
        return {'genus':int(self.genus),'ambient_degree':int(degree),
                'codimension':int(codimension),
                't_power':int(QuarticK3Theory.t_power(self.genus,codimension)),
                'compact_class':'0','infinity_aggregate':-self.cycle(degree,codimension)}

    def top_integral(self,degree):
        result=QQ.zero()
        for key,co in self.coefficients(degree,self.dimension).get(0,{}).items():
            if key[0]=='interior':
                psi=() if key[1] is None else (key[1]+1,)
                result+=co*self.hodge.integral(self.genus,psi,key[2])
            else:
                result+=co*self.hodge.integral(key[1],(key[2],),key[3])*self.hodge.integral(
                    self.genus-key[1],(key[4],),key[5])
        return result
