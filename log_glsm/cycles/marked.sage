r"""Genus-four dual-Euler calculations and true descendant reductions.

The recorded genus-four pass has n <= 2. The shared localization backend
accepts labeled profiles of arbitrary length; batch drivers impose their bounds.
An unknown ordinary invariant is None, never zero. Numerical reductions
are not certificates of equality in the strata algebra."""

from sage.all import *

from functools import lru_cache
from itertools import product
from log_glsm.core.hodge import (
    _load_vendored_admcycles,
)
from log_glsm.core.hodge_grr import (
    GenusFourHodgeIntegrals,
    mumford_monomial,
)
from log_glsm.localization.plane_graphs import (
    P2FixedLocusGraph,
    P2FixedLocusGraphEnumerator,
    _fl_weak_compositions,
)
from log_glsm.targets.grassmannian import (
    GrassmannianAmbient,
)
from log_glsm.targets.projective import (
    CompleteIntersectionCycles,
    ProjectiveCompleteIntersection,
    complete_intersection_weights,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersectionCycles,
    hirzebruch_index,
    toric_degree_tuples,
)

__all__ = ['GenusFourGKMCycles', 'GenusFourMixin', 'GenusFourProjectiveCycles', 'ambient_polarization', 'genus_four_backend', 'marked_gkm_graphs', 'marked_multiplier', 'nonpositive_laurent']


def nonpositive_laurent(co):
    """Exact finite jet through epsilon^0; no fixed precision guess."""
    if not co:return {}
    numerator,denominator=co.numerator(),co.denominator()
    nv,dv=numerator.valuation(),denominator.valuation();v=nv-dv
    if v>0:return {}
    precision=1-v;S=PowerSeriesRing(QQ,'z')
    series=S(numerator.shift(-nv)).add_bigoh(precision)/S(denominator.shift(-dv)).add_bigoh(precision)
    if series.prec()<precision:raise ArithmeticError('insufficient Laurent precision')
    return {int(v+j):QQ(series[j]) for j in range(precision) if series[j]}

class GenusFourMixin:
    def marked_stabilization(self,graph):
        """Retain real labels, collapse map edges, forget only empty tails."""
        n=graph.marking_count
        stable=tuple(v for v in range(graph.vertex_count) if graph.vertex_type(v)=='stable')
        labels={};tails=[];pairs=[];visited=set()
        for v in stable:
            for e in graph.incident_edges(v):
                stub=(v,e)
                if stub in visited:continue
                visited.add(stub);w=graph.other_endpoint(e,v);f=e
                while graph.vertex_type(w)=='nodal':
                    f=next(x for x in graph.incident_edges(w) if x!=f);w=graph.other_endpoint(f,w)
                kind=graph.vertex_type(w)
                if kind=='univalent':tails.append(stub)
                elif kind=='marked':labels[stub]=graph.markings_at_vertex(w)[0]+1
                else:
                    if w not in stable:raise ValueError('unrecognized stable-map endpoint')
                    other=(w,f);visited.add(other);pairs.append((stub,other))
        labels.update({stub:n+i+1 for i,stub in enumerate(tails)})
        edges=[];next_label=n+len(tails)+1
        for left,right in pairs:
            labels[left]=next_label;labels[right]=next_label+1
            edges.append((next_label,next_label+1));next_label+=2
        legs=tuple(tuple(labels[v,e] for e in graph.incident_edges(v))+
                   tuple(m+1 for m in graph.markings_at_vertex(v)) for v in stable)
        return stable,(tuple(graph.vertex_genera[v] for v in stable),legs,tuple(edges),n,len(tails))

    @lru_cache(maxsize=None)
    def marked_vertex_terms(self,genus,colour,weights,marking_psi,local_bound):
        dimension=3*genus-3+len(weights)+len(marking_psi)
        result={}
        for indices,co in self._marked_hodge_vertex(genus,colour,len(weights),len(marking_psi)).items():
            remaining=min(dimension,local_bound)-sum(indices)-sum(marking_psi)
            for total in range(remaining+1):
                for powers in _fl_weak_compositions(total,len(weights)):
                    result[powers+marking_psi,indices]=co/prod(w**(p+1) for w,p in zip(weights,powers))
        return result

    def marked_coefficients(self,beta,insertions=(),max_codimension=3):
        """Exact decorated-graph classes in a declared codimension window.

        The graph data contain genuine marked moduli factors. Descendants at
        unstable marked endpoints use -w, not the stabilized curve psi.
        """
        insertions=tuple(insertions);n=len(insertions)
        if any(x not in ('1','psi','H') for x in insertions):raise ValueError('invalid insertion profile')
        if not 0<=max_codimension<=3*self.genus-3+n:raise ValueError('invalid codimension bound')
        result={r:{} for r in range(max_codimension+1)}
        if isinstance(self.target,ProjectiveCompleteIntersection):
            if isinstance(beta,(tuple,list)) and len(beta)!=1:raise ValueError('projective degree must have one component')
            D=beta[0] if isinstance(beta,(tuple,list)) else beta
            if D not in (0,1,2):raise ValueError('degree must be zero, one, or two')
            ambient_shift=self.target.N*(self.genus-1)-(self.target.N+1)*D
        else:
            if len(beta)!=self.ambient.picard_rank or any(d not in ZZ or d<0 for d in beta) or sum(beta)>2:
                raise ValueError('expected a nonnegative Mori tuple of coefficient sum at most two')
            ambient_shift=self.ambient.dimension*(self.genus-1)-self.ambient.pairing(self.ambient.anticanonical,beta)
        # Theta(t)=t^rank sum_{j>=0} theta_j/t^j, deg(theta_j)=j.
        # A coefficient of curve codimension r has epsilon degree
        # ambient_shift + insertion_degree + j - r. Thus every individual
        # graph term has strictly positive epsilon valuation in this case.
        if ambient_shift+sum(x!='1' for x in insertions)>max_codimension:return result
        for graph in self.marked_graphs(beta,n):
            stable,data=self.marked_stabilization(graph);offset=len(data[2])-data[4]
            # Every Hodge or flag-psi decoration has nonnegative degree.
            # Descendants already on stable vertices give a lower bound
            # before any rational-function edge or vertex work is needed.
            minimum=offset+sum(insertions[i]=='psi' for v in stable for i in graph.markings_at_vertex(v))
            if minimum>max_codimension:continue
            bound=max_codimension-offset;scalar=self.F(graph.localization_weight());local=[]
            for left,right,d in graph.edges:scalar*=self._general_edge(graph.vertex_colours[left],graph.vertex_colours[right],d)
            for v in range(graph.vertex_count):
                colour=graph.vertex_colours[v]
                weights=tuple(graph.flag_weight(e,v,self.a) for e in graph.incident_edges(v))
                marks=tuple(insertions[i] for i in graph.markings_at_vertex(v))
                scalar*=self.hyperplane_character(colour)**marks.count('H')
                kind=graph.vertex_type(v)
                if kind=='stable':local.append(self.marked_vertex_terms(graph.vertex_genera[v],colour,weights,tuple(ZZ(x=='psi') for x in marks),bound))
                elif kind=='univalent':scalar*=weights[0]
                elif kind=='marked':scalar*=(-weights[0])**marks.count('psi')
                elif kind=='nodal':scalar*=self._nodal_factor(colour,weights)
                else:raise ValueError('unsupported fixed vertex')
            for entries in product(*(p.items() for p in local)):
                powers=tuple(item[0][0] for item in entries);lambdas=tuple(item[0][1] for item in entries)
                r=offset+sum(map(sum,powers))+sum(map(sum,lambdas))
                if not 0<=r<=max_codimension:continue
                # Hodge classes pull back through forgetting. A psi class of
                # degree below the relative dimension pushes to zero before
                # multiplying any Hodge factor; do not expand such factors.
                forgotten=[sum(n<f<=n+data[4] for f in flags) for flags in data[1]]
                if any(2*g-2+len(flags)-count>0 and sum(ps)<count
                       for g,flags,ps,count in zip(data[0],data[1],powers,forgotten)):continue
                key=(data,powers,lambdas)
                for power,value in nonpositive_laurent(scalar*prod(item[1] for item in entries)).items():
                    terms=result[r].setdefault(power,{})
                    terms[key]=terms.get(key,QQ.zero())+value
        return {r:{p:{k:v for k,v in terms.items() if v} for p,terms in powers.items() if any(terms.values())} for r,powers in result.items()}

    @lru_cache(maxsize=None)
    def marked_monomial_class(self,key):
        data,powers,lambdas=key;genera,legs,edges,n,t=data
        forget_labels=tuple(range(n+1,n+t+1))
        remaining=[[f for f in flags if f not in forget_labels] for flags in legs]
        global_forget=any(2*g-2+len(flags)<=0 for g,flags in zip(genera,remaining))
        if global_forget:
            graph=self.adm.StableGraph(list(genera),[list(flags) for flags in legs],list(edges));local=[]
            positive=[v for v,g in enumerate(genera) if g]
            global_hodge=len(positive)==1 and genera[positive[0]]==self.genus
            for g,flags,ps,ls in zip(genera,legs,powers,lambdas):
                ring=self.adm.TautologicalRing(g,len(flags))
                part=prod((ring.psi(i+1)**p for i,p in enumerate(ps)),ring.one())
                if not global_hodge:part*=self.hodge_class(g,len(flags),ls)
                local.append(part)
            result=graph.boundary_pushforward(local)
            if t:result=result.forgetful_pushforward(list(forget_labels))
            if global_hodge:
                if all(not result.vector(d) for d in result.degree_list()):return self.adm.TautologicalRing(self.genus,n).zero()
                result*=self.hodge_class(self.genus,n,lambdas[positive[0]])
            return result
        graph=self.adm.StableGraph(list(genera),remaining,list(edges));local=[]
        for g,flags,ps,ls in zip(genera,legs,powers,lambdas):
            ring=self.adm.TautologicalRing(g,len(flags))
            part=prod((ring.psi(i+1)**p for i,p in enumerate(ps)),ring.one())
            forget=[i+1 for i,f in enumerate(flags) if f in forget_labels]
            if sum(ps)<len(forget):return self.adm.TautologicalRing(self.genus,n).zero()
            if forget:
                part=part.forgetful_pushforward(forget)
                part=part.rename_legs({old:new for new,old in enumerate([i for i in range(1,len(flags)+1) if i not in forget],1)},inplace=False)
            if all(not part.vector(d) for d in part.degree_list()):return self.adm.TautologicalRing(self.genus,n).zero()
            part*=self.hodge_class(g,len(flags)-len(forget),ls);local.append(part)
        if len(genera)==1 and not edges:
            # Local order follows edge flags, then original vertex markings.
            return local[0].rename_legs({i+1:f for i,f in enumerate(remaining[0])},inplace=False)
        return graph.boundary_pushforward(local)

    def _hodge_product(self,genus,left,right,dimension):
        result={}
        for l,cl in left.items():
            for r,cr in right.items():
                if sum(l)+sum(r)>dimension:continue
                for indices,scale in mumford_monomial(genus,tuple(sorted(l+r))):
                    result[indices]=result.get(indices,self.F.zero())+scale*cl*cr
        return {key:co for key,co in result.items() if co}

    @lru_cache(maxsize=None)
    def integrated_vertex(self,genus,colour,weights,marking_psi=()):
        flags=len(weights);dimension=3*genus-3+flags+len(marking_psi)
        # Hodge factor's flag count is the edge valence, NOT number of marks.
        # Its old dimension truncation must be enlarged for actual marks.
        answer=self.F.zero()
        for indices,co in self._marked_hodge_vertex(genus,colour,flags,len(marking_psi)).items():
            remaining=dimension-sum(indices)-sum(marking_psi)
            if remaining<0:continue
            for powers in _fl_weak_compositions(remaining,flags):
                integral=self.hodge.integral(genus,powers+marking_psi,indices)
                if integral:answer+=co*integral/prod(w**(p+1) for w,p in zip(weights,powers))
        return self.F(answer)

    @lru_cache(maxsize=None)
    def _marked_hodge_vertex(self,genus,colour,flags,marks):
        dimension=3*genus-3+flags+marks;answer={():self.F.one()}
        tangent,fibres=self.vertex_weights(colour)
        for w in tangent+fibres:
            factor={(() if ell==0 else (ell,)):(-1)**ell*w**(genus-ell) for ell in range(genus+1)}
            answer=self._hodge_product(genus,answer,factor,dimension)
        scalar=prod(tangent)**(flags-1)*prod(q**(1-flags-2*genus) for q in fibres)
        return {indices:self.F(co*scalar) for indices,co in answer.items()}

    def numerical_localization(self,beta,insertions=()):
        """Integrate each fixed locus, then sum finite Laurent jets.

        Direct marked localization is available to check the universal
        string/dilaton/divisor reductions. psi means psi on stable maps.
        """
        insertions=tuple(insertions)
        if any(x not in ('1','psi','H') for x in insertions):
            raise ValueError('each insertion must be 1, psi, or H')
        result={};graphs=self.marked_graphs(beta,len(insertions))
        for graph in graphs:
            co=self.F(graph.localization_weight())
            for left,right,d in graph.edges:
                co*=self._general_edge(graph.vertex_colours[left],graph.vertex_colours[right],d)
            for v in range(graph.vertex_count):
                colour=graph.vertex_colours[v]
                weights=tuple(sorted(graph.flag_weight(e,v,self.a) for e in graph.incident_edges(v)))
                marks=tuple(insertions[i] for i in graph.markings_at_vertex(v))
                co*=self.hyperplane_character(colour)**marks.count('H')
                kind=graph.vertex_type(v)
                if kind=='stable':co*=self.integrated_vertex(graph.vertex_genera[v],colour,weights,tuple(ZZ(x=='psi') for x in marks))
                elif kind=='univalent':co*=weights[0]
                elif kind=='marked':co*=(-weights[0])**marks.count('psi')
                elif kind=='nodal':co*=self._nodal_factor(colour,weights)
                else:raise ValueError('unsupported fixed vertex')
            for power,value in nonpositive_laurent(co).items():result[power]=result.get(power,QQ.zero())+value
        return {'finite':result.get(0,QQ.zero()),'poles':{k:v for k,v in result.items() if k<0 and v},'graphs':len(graphs)}

class GenusFourProjectiveCycles(GenusFourMixin,CompleteIntersectionCycles):
    def __init__(self,target,genus=4,base=4,hodge=None):
        if genus not in (2,3,4):raise ValueError('supported genera are 2, 3, 4')
        self.target=target;self.genus=ZZ(genus);self.dimension=3*genus-3
        self.weights=complete_intersection_weights(target.N+1,base)
        self.P=PolynomialRing(QQ,'epsilon');self.epsilon=self.P.gen();self.F=self.P.fraction_field()
        self.a=tuple(self.F(a*self.epsilon) for a in self.weights)
        self.adm=_load_vendored_admcycles();self.R=self.adm.TautologicalRing(genus,0)
        self.hodge=hodge or GenusFourHodgeIntegrals()

    def vertex_weights(self,colour):
        return (tuple(self.a[colour]-self.a[j] for j in range(self.target.N+1) if j!=colour),
                tuple(1-d*self.a[colour] for d in self.target.degrees))

    def hyperplane_character(self,colour):return self.a[colour]

    def marked_graphs(self,beta,markings):
        degree=ZZ(beta[0]) if isinstance(beta,(list,tuple)) else ZZ(beta)
        if degree not in (0,1,2):raise ValueError('degree must be zero, one, or two')
        key=(int(self.target.N),int(self.genus),int(degree),int(markings))
        if key not in self._graph_cache:
            graph_type=type('MarkedProjectiveGraph',(P2FixedLocusGraph,),{'fixed_point_count':int(self.target.N+1)})
            enumerator=type('MarkedProjectiveEnumerator',(P2FixedLocusGraphEnumerator,),{'fixed_point_count':int(self.target.N+1),'graph_type':graph_type})
            self._graph_cache[key]=enumerator(self.genus,markings,degree).graphs()
        return self._graph_cache[key]

class GenusFourGKMCycles(GenusFourMixin,ToricCompleteIntersectionCycles):
    def __init__(self,target,genus=4,base=4,hodge=None):
        if genus not in (2,3,4):raise ValueError('supported genera are 2, 3, 4')
        self.target=target;self.ambient=target.ambient;self.genus=ZZ(genus);self.dimension=3*genus-3
        self.weights=self.ambient.ray(base)
        self.P=PolynomialRing(QQ,'epsilon');self.epsilon=self.P.gen();self.F=self.P.fraction_field()
        self.a=tuple(self.F(w*self.epsilon) for w in self.weights)
        self.tangent=tuple(tuple(self.ambient.character_value(w,self.a) for w in ws) for ws in self.ambient.tangent)
        self.q=tuple(tuple(1-self.ambient.character_value(self.ambient.fibre_character(L,i),self.a) for L in target.bundles) for i in range(len(self.ambient.points)))
        self.adm=_load_vendored_admcycles();self.R=self.adm.TautologicalRing(genus,0)
        self.hodge=hodge or GenusFourHodgeIntegrals()
        self.polarization=ambient_polarization(target)

    def vertex_weights(self,colour):return self.tangent[colour],self.q[colour]

    def hyperplane_character(self,colour):
        return self.ambient.character_value(self.ambient.fibre_character(self.polarization,colour),self.a)

    def marked_graphs(self,beta,markings):
        if not markings:return self.graphs(beta)
        return marked_gkm_graphs(self.ambient,self.genus,markings)[tuple(beta)]

_MARKED_GKM_GRAPHS={}

def marked_gkm_graphs(ambient,genus,markings):
    key=(ambient.slug,int(genus),int(markings))
    if key in _MARKED_GKM_GRAPHS:return _MARKED_GKM_GRAPHS[key]
    class Graph(P2FixedLocusGraph):
        fixed_point_count=len(ambient.points)
        def flag_weight(self,edge_index,vertex,lambdas):
            other=self.other_endpoint(edge_index,vertex);e=self.edges[edge_index][2]
            return ambient.character_value(ambient.edges[self.vertex_colours[vertex],self.vertex_colours[other]]['weight'],lambdas)/e
        def mori_degree(self):
            beta=vector(ZZ,ambient.picard_rank)
            for i,j,e in self.edges:beta+=e*vector(ambient.edges[self.vertex_colours[i],self.vertex_colours[j]]['class'])
            return tuple(beta)
    class Enumerator(P2FixedLocusGraphEnumerator):
        fixed_point_count=len(ambient.points);graph_type=Graph
        @classmethod
        def _proper_colours(cls,vertex_count,endpoints):
            for colours in product(range(cls.fixed_point_count),repeat=vertex_count):
                if all((colours[i],colours[j]) in ambient.edges for i,j in endpoints):yield colours
    result={b:[] for b in toric_degree_tuples(ambient)}
    for d in range(3):
        for graph in Enumerator(genus,markings,d).graphs():
            beta=graph.mori_degree()
            if beta in result:result[beta].append(graph)
    _MARKED_GKM_GRAPHS[key]={b:tuple(gs) for b,gs in result.items()}
    return _MARKED_GKM_GRAPHS[key]

def ambient_polarization(target):
    if isinstance(target,ProjectiveCompleteIntersection):return (ZZ(1),)
    if isinstance(target.ambient,GrassmannianAmbient):return (ZZ(1),)
    if hasattr(target.ambient,'polarization'):return target.ambient.polarization
    return tuple(x for f in target.ambient.factors for x in ((1,hirzebruch_index(f)+1) if hirzebruch_index(f) is not None else (1,)))

def marked_multiplier(genus,hyperplane_degree,insertions):
    insertions=tuple(insertions)
    if genus<2 or any(x not in ('1','psi','H') for x in insertions):
        raise ValueError('genus >= 2 and insertions 1/psi/H required')
    if '1' in insertions:return ZZ.zero()
    k=insertions.count('psi');h=insertions.count('H')
    return ZZ(hyperplane_degree)**h*prod(2*genus-2+h+i for i in range(k))

def genus_four_backend(target,base=4,hodge=None,genus=4):
    cls=GenusFourProjectiveCycles if isinstance(target,ProjectiveCompleteIntersection) else GenusFourGKMCycles
    return cls(target,genus,base,hodge)
