r"""Cycle-valued dual-Euler O(3) theory on Mbar_3, in degrees zero and one.

The torus-fixed vertex polynomials are retained as classes. Degree-one
rational tails are forgotten (psi^a pushes to kappa_(a-1)); rational
bridges are contracted and the genus-one and genus-two factors are glued.
The auxiliary torus limit is taken only after summing all fixed graphs.
Negative powers of its parameter are returned for independent verification
in the tautological ring, never silently dropped as zero classes."""

from sage.all import *

from functools import lru_cache
from types import SimpleNamespace
from log_glsm.core.hodge import (
    HodgeIntegralRequest,
    _load_vendored_admcycles,
)
from log_glsm.localization.plane import (
    P2FixedLocusEvaluator,
)
from log_glsm.localization.plane_graphs import (
    P2FixedLocusGraphEnumerator,
)

__all__ = ['O3UnmarkedCycleBackend']


class O3UnmarkedCycleBackend:
    def __init__(self, base_weights=(0,1,3)):
        if len(base_weights)!=3 or len(set(base_weights))!=3:
            raise ValueError('three distinct auxiliary base weights are required')
        self.P=PolynomialRing(QQ,'epsilon')
        self.epsilon=self.P.gen()
        self.F=self.P.fraction_field()
        self.rings=SimpleNamespace(full_field=self.F,full=self.F,t=self.F.one(),
            base_weights=lambda:tuple(self.F(a*self.epsilon) for a in base_weights))
        self.evaluator=P2FixedLocusEvaluator(self.rings,lift_strategy='standard')
        self.adm=_load_vendored_admcycles()
        self.R=self.adm.TautologicalRing(3,0)
        self._cache={}

    @staticmethod
    def _normal(genus,psi,lambdas):
        factor,request=HodgeIntegralRequest(genus,psi,lambdas).low_genus_normal_form()
        return factor,request.lambda_indices

    @staticmethod
    def codimension(key):
        if key[0]=='interior':
            return (0 if key[1] is None else key[1])+sum(key[2])
        return 1+key[1]+sum(key[2])+key[3]+sum(key[4])

    def rational_terms(self,degree):
        if degree not in (0,1):
            raise ValueError('cycle stabilization is implemented only for degrees zero and one')
        if degree in self._cache:
            return self._cache[degree]
        terms={}
        def add(key,value):
            if value:
                terms[key]=terms.get(key,self.F.zero())+value
        graphs=P2FixedLocusGraphEnumerator(3,0,degree).graphs()
        for graph in graphs:
            scalar=self.F(graph.localization_weight())
            stable=[]
            for e in range(graph.edge_count):
                scalar*=self.evaluator._edge_factor(graph,e)
            for v in range(graph.vertex_count):
                if graph.vertex_type(v)=='stable':
                    weights=tuple(graph.flag_weight(e,v,self.evaluator.lambdas)
                                  for e in graph.incident_edges(v))
                    poly=self.evaluator._stable_vertex_polynomial(
                        graph.vertex_genera[v],graph.vertex_colours[v],weights,())
                    stable.append((graph.vertex_genera[v],poly))
                else:
                    scalar*=self.evaluator._unstable_vertex_factor(graph,v,())
            if len(stable)==1:
                genus,poly=stable[0]
                assert genus==3
                for (psi,lambdas),co in poly.terms.items():
                    factor,lambdas=self._normal(genus,psi,lambdas)
                    if not factor:continue
                    if degree==0:
                        assert not psi
                        add(('interior',None,lambdas),scalar*co*factor)
                    else:
                        assert len(psi)==1
                        if psi[0]:
                            # Hodge classes pull back through the forgetful map.
                            add(('interior',int(psi[0]-1),lambdas),scalar*co*factor)
            else:
                assert degree==1 and sorted(g for g,_ in stable)==[1,2]
                stable.sort(key=lambda pair:pair[0])
                for (psi1,lam1),co1 in stable[0][1].terms.items():
                    f1,lam1=self._normal(1,psi1,lam1)
                    if not f1:continue
                    for (psi2,lam2),co2 in stable[1][1].terms.items():
                        f2,lam2=self._normal(2,psi2,lam2)
                        if f2:
                            add(('boundary12',int(psi1[0]),lam1,int(psi2[0]),lam2),
                                scalar*co1*co2*f1*f2)
        self._cache[degree]={k:self.F(v) for k,v in terms.items() if v}
        return self._cache[degree]

    @lru_cache(maxsize=None)
    def hodge_class(self,genus,markings,indices):
        ring=self.adm.TautologicalRing(genus,markings)
        if not indices:return ring.one()
        return self.hodge_class(genus,markings,indices[:-1])*ring.lambdaclass(indices[-1])

    @lru_cache(maxsize=None)
    def monomial_class(self,key):
        if key[0]=='interior':
            answer=self.hodge_class(3,0,key[2])
            if key[1] is not None:answer*=self.R.kappa(key[1])
            return answer
        left=self.adm.TautologicalRing(1,1)
        right=self.adm.TautologicalRing(2,1)
        a=left.psi(1)**key[1]*self.hodge_class(1,1,key[2])
        b=right.psi(1)**key[3]*self.hodge_class(2,1,key[4])
        graph=self.adm.StableGraph([1,2],[[1],[2]],[(1,2)])
        return graph.boundary_pushforward([a,b])

    def top_integral(self,degree):
        """Evaluate the top class through Hodge integrals, without expanding graphs."""
        hodge=self.evaluator.hodge
        answer=QQ.zero()
        for key,co in self.coefficients(degree,6).get(0,{}).items():
            if key[0]=='interior':
                psi=() if key[1] is None else (key[1]+1,)
                answer+=co*hodge.integral(3,psi,key[2])
            else:
                answer+=co*hodge.integral(1,(key[1],),key[2])*hodge.integral(2,(key[3],),key[4])
        return answer

    def coefficients(self,degree,codimension):
        """Return epsilon powers <= 0, each with formal monomial coefficients."""
        if codimension<0 or codimension>6:
            raise ValueError('codimension must be between zero and six on Mbar_3')
        result={}
        series_ring=LaurentSeriesRing(QQ,'epsilon',default_prec=32)
        for key,co in self.rational_terms(degree).items():
            if self.codimension(key)!=codimension:continue
            series=series_ring(co.numerator())/series_ring(co.denominator())
            if not series:continue
            if series.prec()<=0:
                raise ArithmeticError('insufficient precision for the auxiliary torus limit')
            for power in range(int(series.valuation()),1):
                value=QQ(series[power])
                if value:result.setdefault(power,{})[key]=value
        return result

    def cycle(self,degree,codimension,epsilon_power=0):
        terms=self.coefficients(degree,codimension).get(epsilon_power,{})
        return sum((value*self.monomial_class(key) for key,value in terms.items()),self.R.zero())

    def record(self,degree,codimension):
        coefficients=self.coefficients(degree,codimension)
        return {'genus':3,'ambient_degree':int(degree),'codimension':int(codimension),
            't_power':int(2-codimension),'twist':'dual Euler',
            'terms':[{'key':key,'coefficient':str(value)}
                     for key,value in coefficients.get(0,{}).items()],
            'auxiliary_pole_powers':sorted(int(k) for k in coefficients if k<0)}
