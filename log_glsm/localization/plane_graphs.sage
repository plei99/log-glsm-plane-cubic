r"""localization / plane_graphs."""

from sage.all import *

from collections import defaultdict
from collections import deque
from itertools import product
from log_glsm.localization.twisted_plane import (
    TwistedInsertion,
)

__all__ = ['P2FixedLocusGraph', 'P2FixedLocusGraphEnumerator', 'SparseEquivariantLiftPlanner']


def _fl_weak_compositions(total, length):
    total = ZZ(total)
    length = ZZ(length)
    if total < 0 or length < 0:
        return
    if length == 0:
        if total == 0:
            yield tuple()
        return
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in _fl_weak_compositions(total - first, length - 1):
            yield (ZZ(first),) + tail

def _fl_positive_compositions(total, length):
    total = ZZ(total)
    length = ZZ(length)
    if length == 0:
        if total == 0:
            yield tuple()
        return
    if total < length:
        return
    for shifted in _fl_weak_compositions(total - length, length):
        yield tuple(value + 1 for value in shifted)

def _fl_connected(vertex_count, endpoint_pairs):
    vertex_count = ZZ(vertex_count)
    if vertex_count == 1:
        return len(endpoint_pairs) == 0
    adjacency = [[] for _ in range(vertex_count)]
    for left, right in endpoint_pairs:
        adjacency[left].append(right)
        adjacency[right].append(left)
    seen = {0}
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        for other in adjacency[vertex]:
            if other not in seen:
                seen.add(other)
                queue.append(other)
    return len(seen) == vertex_count

class SparseEquivariantLiftPlanner(SageObject):
    r"""Choose deterministic sparse lifts of ``H`` and ``H^2``.

    Point insertions are supported cyclically at ``p_0,p_1,p_2``.  Divisor
    insertions use ``H-lambda_i``, again cycling through the fixed points.
    Explicitly chosen lifts are preserved.  The cyclic rule is inexpensive,
    symmetric, and forces all three colours as soon as three point insertions
    occur; for five points in degree two this reduces 1557 graph components
    to the three all-colour paths.
    """

    def plan(self, insertions):
        insertions = tuple(insertions)
        counters = {1: ZZ.zero(), 2: ZZ.zero()}
        planned = []
        for insertion in insertions:
            if not isinstance(insertion, TwistedInsertion):
                raise TypeError("sparse lifts require TwistedInsertion objects")
            if insertion.h_power in (1, 2) \
                    and insertion.lift_kind == "standard":
                index = counters[insertion.h_power] % 3
                counters[insertion.h_power] += 1
                planned.append(insertion.with_sparse_lift(index))
            else:
                planned.append(insertion)
        return tuple(planned)

class P2FixedLocusGraph(SageObject):
    """One decorated fixed graph for stable maps to ``P^2``."""

    fixed_point_count = 3

    def __init__(self, vertex_genera, vertex_colours,
                 marking_vertices=(), edges=()):
        self.vertex_genera = tuple(ZZ(value) for value in vertex_genera)
        self.vertex_colours = tuple(ZZ(value) for value in vertex_colours)
        self.marking_vertices = tuple(ZZ(value) for value in marking_vertices)
        self.edges = tuple(sorted(
            (min(ZZ(left), ZZ(right)), max(ZZ(left), ZZ(right)), ZZ(degree))
            for left, right, degree in edges
        ))
        self._canonical_key_cache = None
        self._automorphism_order_cache = None

    @property
    def vertex_count(self):
        return ZZ(len(self.vertex_genera))

    @property
    def edge_count(self):
        return ZZ(len(self.edges))

    @property
    def marking_count(self):
        return ZZ(len(self.marking_vertices))

    def first_betti_number(self):
        return self.edge_count - self.vertex_count + 1

    def total_genus(self):
        return ZZ(sum(self.vertex_genera) + self.first_betti_number())

    def total_degree(self):
        return ZZ(sum(degree for _, _, degree in self.edges))

    def incident_edges(self, vertex):
        vertex = ZZ(vertex)
        return tuple(
            edge_index for edge_index, (left, right, _) in enumerate(self.edges)
            if left == vertex or right == vertex
        )

    def markings_at_vertex(self, vertex):
        vertex = ZZ(vertex)
        return tuple(
            marking for marking, support in enumerate(self.marking_vertices)
            if support == vertex
        )

    def edge_valence(self, vertex):
        return ZZ(len(self.incident_edges(vertex)))

    def marking_valence(self, vertex):
        return ZZ(len(self.markings_at_vertex(vertex)))

    def vertex_type(self, vertex):
        genus = self.vertex_genera[vertex]
        edge_valence = self.edge_valence(vertex)
        marking_valence = self.marking_valence(vertex)
        if 2 * genus - 2 + edge_valence + marking_valence > 0:
            return "stable"
        if genus == 0 and edge_valence == 1 and marking_valence == 0:
            return "univalent"
        if genus == 0 and edge_valence == 1 and marking_valence == 1:
            return "marked"
        if genus == 0 and edge_valence == 2 and marking_valence == 0:
            return "nodal"
        return "invalid"

    def other_endpoint(self, edge_index, vertex):
        left, right, _ = self.edges[edge_index]
        if vertex == left:
            return right
        if vertex == right:
            return left
        raise ValueError("the edge is not incident to the vertex")

    def flag_weight(self, edge_index, vertex, lambdas):
        other = self.other_endpoint(edge_index, vertex)
        degree = self.edges[edge_index][2]
        return (lambdas[self.vertex_colours[vertex]]
                - lambdas[self.vertex_colours[other]]) / degree

    def validation_errors(self):
        errors = []
        if len(self.vertex_genera) != len(self.vertex_colours):
            errors.append("vertex genus/colour arrays have different lengths")
            return errors
        if self.vertex_count == 0:
            errors.append("the graph has no vertices")
            return errors
        if any(genus < 0 for genus in self.vertex_genera):
            errors.append("vertex genera must be nonnegative")
        if any(colour not in range(self.fixed_point_count) for colour in self.vertex_colours):
            errors.append("vertex colour is outside the target's fixed points")
        for support in self.marking_vertices:
            if support < 0 or support >= self.vertex_count:
                errors.append("a marking has an invalid supporting vertex")
        endpoint_pairs = []
        for left, right, degree in self.edges:
            if left < 0 or right >= self.vertex_count or left >= right:
                errors.append("an edge has invalid endpoints")
                continue
            endpoint_pairs.append((left, right))
            if degree <= 0:
                errors.append("edge degrees must be positive")
            if self.vertex_colours[left] == self.vertex_colours[right]:
                errors.append("an edge joins equal fixed-point colours")
        if not errors and not _fl_connected(self.vertex_count, endpoint_pairs):
            errors.append("the fixed graph is disconnected")
        for vertex in range(self.vertex_count):
            if self.vertex_type(vertex) == "invalid":
                errors.append("vertex %s has an invalid unstable type" % vertex)
        return errors

    def is_valid(self):
        return not self.validation_errors()

    def _coloured_incidence_graph(self):
        graph = Graph(multiedges=False, loops=False)
        buckets = defaultdict(list)
        for vertex in range(self.vertex_count):
            node = ("v", vertex)
            graph.add_vertex(node)
            decoration = (
                "vertex", self.vertex_colours[vertex],
                self.vertex_genera[vertex], self.markings_at_vertex(vertex),
            )
            buckets[decoration].append(node)
        for edge_index, (left, right, degree) in enumerate(self.edges):
            node = ("e", edge_index)
            graph.add_vertex(node)
            graph.add_edge(node, ("v", left))
            graph.add_edge(node, ("v", right))
            buckets[("edge", degree)].append(node)
        colours = tuple(sorted(buckets))
        partition = [buckets[colour] for colour in colours]
        return graph, colours, partition

    def canonical_key(self):
        if self._canonical_key_cache is None:
            graph, colours, partition = self._coloured_incidence_graph()
            canonical = graph.canonical_label(partition=partition)
            colour_profile = tuple(
                (colour, len(block))
                for colour, block in zip(colours, partition)
            )
            self._canonical_key_cache = (
                colour_profile, tuple(canonical.edges(labels=False, sort=True))
            )
        return self._canonical_key_cache

    def automorphism_order(self):
        if self._automorphism_order_cache is None:
            graph, _, partition = self._coloured_incidence_graph()
            self._automorphism_order_cache = ZZ(
                graph.automorphism_group(partition=partition).order()
            )
        return self._automorphism_order_cache

    def deck_order(self):
        return ZZ(prod(degree for _, _, degree in self.edges))

    def localization_weight(self):
        return QQ(1) / (self.automorphism_order() * self.deck_order())

    def signature(self):
        vertices = tuple(
            (self.vertex_colours[vertex], self.vertex_genera[vertex],
             self.vertex_type(vertex), self.markings_at_vertex(vertex))
            for vertex in range(self.vertex_count)
        )
        return vertices, self.edges, self.automorphism_order(), self.deck_order()

    def __hash__(self):
        return hash(self.canonical_key())

    def __eq__(self, other):
        return (isinstance(other, P2FixedLocusGraph)
                and self.canonical_key() == other.canonical_key())

    def _repr_(self):
        return "P2FixedLocusGraph%s" % (self.signature(),)

class P2FixedLocusGraphEnumerator(SageObject):
    """Enumerate all fixed graphs for ``Mbar_(g,n)(P^2,d)``."""

    fixed_point_count = 3
    graph_type = P2FixedLocusGraph

    def __init__(self, genus, markings, degree,
                 allowed_marking_colours=None):
        self.genus = ZZ(genus)
        self.markings = ZZ(markings)
        self.degree = ZZ(degree)
        if self.genus < 0 or self.markings < 0 or self.degree < 0:
            raise ValueError("genus, markings, and degree must be nonnegative")
        if allowed_marking_colours is None:
            self.allowed_marking_colours = tuple(
                tuple(ZZ(i) for i in range(self.fixed_point_count)) for _ in range(self.markings)
            )
        else:
            if len(allowed_marking_colours) != self.markings:
                raise ValueError("one fixed-point choice is needed per marking")
            normalized = []
            for choices in allowed_marking_colours:
                choices = tuple(sorted(set(ZZ(value) for value in choices)))
                if not choices or any(value not in range(self.fixed_point_count)
                                      for value in choices):
                    raise ValueError(
                        "allowed marking colours must be nonempty subsets of the target's fixed points"
                    )
                normalized.append(choices)
            self.allowed_marking_colours = tuple(normalized)
        self._cache = None

    @staticmethod
    def _endpoint_multigraphs(vertex_count, edge_count):
        cells = tuple(
            (left, right)
            for left in range(vertex_count)
            for right in range(left + 1, vertex_count)
        )
        for multiplicities in _fl_weak_compositions(edge_count, len(cells)):
            endpoints = []
            for pair, multiplicity in zip(cells, multiplicities):
                endpoints.extend([pair] * multiplicity)
            if _fl_connected(vertex_count, endpoints):
                yield tuple(endpoints)

    @classmethod
    def _proper_colours(cls, vertex_count, endpoints):
        for colours in product(range(cls.fixed_point_count), repeat=vertex_count):
            if all(colours[left] != colours[right] for left, right in endpoints):
                yield tuple(ZZ(colour) for colour in colours)

    def _degree_zero_graphs(self):
        for colour in range(self.fixed_point_count):
            if any(colour not in choices
                   for choices in self.allowed_marking_colours):
                continue
            graph = self.graph_type(
                (self.genus,), (colour,),
                tuple(0 for _ in range(self.markings)), (),
            )
            if graph.is_valid():
                yield graph

    def _positive_degree_graphs(self):
        for edge_count in range(1, self.degree + 1):
            for vertex_count in range(2, edge_count + 2):
                first_betti = edge_count - vertex_count + 1
                residual_genus = self.genus - first_betti
                if residual_genus < 0:
                    continue
                for endpoints in self._endpoint_multigraphs(
                        vertex_count, edge_count):
                    for colours in self._proper_colours(vertex_count, endpoints):
                        for genera in _fl_weak_compositions(
                                residual_genus, vertex_count):
                            marking_choices = tuple(
                                tuple(
                                    vertex for vertex in range(vertex_count)
                                    if colours[vertex] in allowed_colours
                                )
                                for allowed_colours
                                in self.allowed_marking_colours
                            )
                            if any(not choices for choices in marking_choices):
                                continue
                            for marking_vertices in product(*marking_choices):
                                provisional = self.graph_type(
                                    genera, colours, marking_vertices,
                                    tuple((left, right, 1)
                                          for left, right in endpoints),
                                )
                                if any(provisional.vertex_type(vertex) == "invalid"
                                       for vertex in range(vertex_count)):
                                    continue
                                for degrees in _fl_positive_compositions(
                                        self.degree, edge_count):
                                    graph = self.graph_type(
                                        genera, colours, marking_vertices,
                                        tuple(
                                            (left, right, degree)
                                            for (left, right), degree
                                            in zip(endpoints, degrees)
                                        ),
                                    )
                                    if graph.is_valid():
                                        yield graph

    def graphs(self):
        if self._cache is None:
            representatives = {}
            source = (self._degree_zero_graphs() if self.degree == 0
                      else self._positive_degree_graphs())
            for graph in source:
                representatives.setdefault(graph.canonical_key(), graph)
            self._cache = tuple(
                representatives[key] for key in sorted(representatives)
            )
        return self._cache

    def graph_count(self):
        return ZZ(len(self.graphs()))
