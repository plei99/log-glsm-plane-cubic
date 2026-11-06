r"""Quartic cycle localization through ambient degree three.

Retain contracted stable components. Replace rational tails by temporary
markings, collapse bivalent rational vertices, glue, then forget temporary
markings using admcycles. This includes all multiple-tail correction terms."""

from sage.all import *

from functools import lru_cache
from itertools import product
from math import factorial
from log_glsm.localization.plane_graphs import (
    P2FixedLocusGraph,
    P2FixedLocusGraphEnumerator,
    _fl_weak_compositions,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
)

__all__ = ['P3FixedLocusGraph', 'P3FixedLocusGraphEnumerator', 'QuarticK3GraphCycles']


class P3FixedLocusGraph(P2FixedLocusGraph):
    fixed_point_count=4

    def _repr_(self):
        return 'P3FixedLocusGraph%s'%(self.signature(),)

class P3FixedLocusGraphEnumerator(P2FixedLocusGraphEnumerator):
    fixed_point_count=4
    graph_type=P3FixedLocusGraph

class QuarticK3GraphCycles(QuarticK3Cycles):
    @lru_cache(maxsize=None)
    def graphs(self,degree):
        if degree not in (0,1,2,3):
            raise NotImplementedError('the graph calculation covers degrees zero through three')
        return P3FixedLocusGraphEnumerator(self.genus,0,degree).graphs()

    @lru_cache(maxsize=None)
    def _general_edge(self,i,j,degree):
        difference=self.a[i]-self.a[j]
        factor=(-1)**degree*degree**(2*degree)/(factorial(degree)**2*difference**(2*degree))
        denominator=prod((QQ(m)/degree*self.a[i]+QQ(degree-m)/degree*self.a[j]-self.a[k])
                         for k in range(4) if k not in (i,j) for m in range(degree+1))
        if not denominator:raise ValueError('resonant auxiliary weights on a multiple-cover edge')
        factor/=denominator
        factor*=prod(1-QQ(4*degree-m)/degree*self.a[i]-QQ(m)/degree*self.a[j]
                     for m in range(4*degree+1))
        return self.F(factor)

    @lru_cache(maxsize=None)
    def _general_vertex(self,genus,colour,weights):
        dimension=3*genus-3+len(weights)
        answer={}
        for indices,co in self._hodge_vertex(genus,colour,len(weights)).items():
            for total in range(dimension-sum(indices)+1):
                for powers in _fl_weak_compositions(total,len(weights)):
                    answer[powers,indices]=co/prod(w**(p+1) for w,p in zip(weights,powers))
        return answer

    def _stabilized_graph(self,graph):
        stable=tuple(v for v in range(graph.vertex_count) if graph.vertex_type(v)=='stable')
        assert stable
        vmap={v:i for i,v in enumerate(stable)}
        tails=[]; pairs=[]; visited=set()
        for v in stable:
            for e in graph.incident_edges(v):
                stub=(v,e)
                if stub in visited:continue
                visited.add(stub)
                next_vertex=graph.other_endpoint(e,v); next_edge=e
                while graph.vertex_type(next_vertex)=='nodal':
                    next_edge=next(x for x in graph.incident_edges(next_vertex) if x!=next_edge)
                    next_vertex=graph.other_endpoint(next_edge,next_vertex)
                if graph.vertex_type(next_vertex)=='univalent':
                    tails.append(stub)
                else:
                    assert next_vertex in vmap
                    other=(next_vertex,next_edge)
                    visited.add(other)
                    pairs.append((stub,other))
        labels={stub:i+1 for i,stub in enumerate(tails)}
        edges=[]; next_label=len(tails)+1
        for left,right in pairs:
            labels[left]=next_label; labels[right]=next_label+1
            edges.append((next_label,next_label+1));next_label+=2
        legs=tuple(tuple(labels[v,e] for e in graph.incident_edges(v)) for v in stable)
        data=(tuple(graph.vertex_genera[v] for v in stable),legs,tuple(edges),len(tails))
        return stable,data

    @staticmethod
    def codimension(key):
        if key[0]!='graph':return QuarticK3Cycles.codimension(key)
        data,powers,lambdas=key[1:]
        return len(data[2])-data[3]+sum(sum(p) for p in powers)+sum(sum(l) for l in lambdas)

    def _nodal_factor(self,colour,weights):
        if not sum(weights):raise ValueError('resonant weights at a bivalent vertex')
        return prod(self.a[colour]-self.a[k] for k in range(4) if k!=colour)/(
            sum(weights)*(1-4*self.a[colour]))

    @lru_cache(maxsize=None)
    def rational_terms(self,degree):
        if degree==0:return super().rational_terms(degree)
        result={}
        for graph in self.graphs(degree):
            scalar=self.F(graph.localization_weight())
            for left,right,d in graph.edges:
                scalar*=self._general_edge(graph.vertex_colours[left],graph.vertex_colours[right],d)
            stable,data=self._stabilized_graph(graph)
            local=[]
            for v in range(graph.vertex_count):
                colour=graph.vertex_colours[v]
                weights=tuple(graph.flag_weight(e,v,self.a) for e in graph.incident_edges(v))
                kind=graph.vertex_type(v)
                if kind=='stable':
                    local.append(self._general_vertex(graph.vertex_genera[v],colour,weights))
                elif kind=='univalent':scalar*=weights[0]
                elif kind=='nodal':
                    scalar*=self._nodal_factor(colour,weights)
                else:raise ValueError('unsupported unstable vertex')
            for entries in product(*(p.items() for p in local)):
                powers=tuple(item[0][0] for item in entries)
                lambdas=tuple(item[0][1] for item in entries)
                # Forgotten markings on a fixed component are interchangeable.
                # Canonicalize their exponents before collecting graph terms.
                canonical=[]
                for flags,ps in zip(data[1],powers):
                    ps=list(ps)
                    slots=[i for i,flag in enumerate(flags) if flag<=data[3]]
                    values=sorted(ps[i] for i in slots)
                    for i,value in zip(slots,values):ps[i]=value
                    canonical.append(tuple(ps))
                powers=tuple(canonical)
                key=('graph',data,powers,lambdas)
                if not 0<=self.codimension(key)<=self.dimension:continue
                co=scalar*prod(item[1] for item in entries)
                result[key]=result.get(key,self.F.zero())+co
        return {k:self.F(c) for k,c in result.items() if c}

    @lru_cache(maxsize=None)
    def monomial_class(self,key):
        if key[0]!='graph':return super().monomial_class(key)
        data,powers,lambdas=key[1:]
        genera,legs,edges,markings=data
        remaining_legs=[[flag for flag in flags if flag>markings] for flags in legs]
        if any(2*g-2+len(flags)<=0 for g,flags in zip(genera,remaining_legs)):
            # A degree-three trivalent contracted rational component can
            # destabilize when its temporary tail markings are forgotten.
            # Let admcycles perform the global stabilization in that case.
            graph=self.adm.StableGraph(list(genera),[list(l) for l in legs],list(edges))
            local=[]
            positive=[v for v,g in enumerate(genera) if g]
            global_hodge=(len(positive)==1 and genera[positive[0]]==self.genus)
            for genus,flags,psi,indices in zip(genera,legs,powers,lambdas):
                ring=self.adm.TautologicalRing(genus,len(flags))
                part=prod((ring.psi(i+1)**p for i,p in enumerate(psi)),ring.one())
                if not global_hodge:part*=self.hodge_class(genus,len(flags),indices)
                local.append(part)
            result=graph.boundary_pushforward(local).forgetful_pushforward(list(range(1,markings+1)))
            if global_hodge:result*=self.hodge_class(self.genus,0,lambdas[positive[0]])
            return result
        graph=self.adm.StableGraph(list(genera),remaining_legs,list(edges))
        local=[]
        for genus,flags,psi,indices in zip(genera,legs,powers,lambdas):
            ring=self.adm.TautologicalRing(genus,len(flags))
            forget=[i+1 for i,flag in enumerate(flags) if flag<=markings]
            # Hodge classes pull back under forgetting. Push the psi monomial
            # first to avoid a large, unnecessary lambda expansion on Mbar_(g,n).
            part=prod((ring.psi(i+1)**p for i,p in enumerate(psi)),ring.one())
            if forget:part=part.forgetful_pushforward(forget)
            # Renumber retained local markings to match boundary_pushforward.
            if forget and len(flags)>len(forget):
                part=part.rename_legs({old:new for new,old in enumerate(
                    [i for i in range(1,len(flags)+1) if i not in forget],1)},inplace=False)
            part*=self.hodge_class(genus,len(flags)-len(forget),indices)
            local.append(part)
        if len(genera)==1 and not edges:
            return local[0]
        return graph.boundary_pushforward(local)

    def top_integral(self,degree):
        if degree==0:return super().top_integral(degree)
        result=QQ.zero()
        for key,co in self.coefficients(degree,self.dimension).get(0,{}).items():
            data,powers,lambdas=key[1:]
            result+=co*prod(self.hodge.integral(g,p,l)
                           for g,p,l in zip(data[0],powers,lambdas))
        return result
