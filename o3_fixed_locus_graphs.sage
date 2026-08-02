r"""Base-torus localization for full ``O(3)``-twisted ``P^2`` vertices.

The fixed loci of ``Mbar_(g,n)(P^2,d)`` are decorated connected graphs:

* a vertex is labelled by a genus and one of the three torus fixed points;
* an edge joins vertices of different colours and has a positive degree;
* ordinary markings are labelled and lie on vertices.

This module enumerates those graphs, includes graph and deck automorphisms,
and integrates every stable-vertex Hodge/psi expression with ``admcycles``.
The implementation can also omit the ``O(3)`` Euler class, which is useful
for testing the virtual-normal conventions against ordinary ``P^2`` GW
invariants.

Conventions agree with ``log_glsm_conventions.sage``: ``H|p_i=lambda_i``
and the fibre weight of ``(O(3))^vee`` at ``p_i`` is ``t-3*lambda_i``.
"""

load("o3_twisted_plane_vertices.sage")

from collections import defaultdict, deque
from itertools import product
from math import factorial


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
        if any(colour not in (0, 1, 2) for colour in self.vertex_colours):
            errors.append("vertex colours must be 0, 1, or 2")
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

    def __init__(self, genus, markings, degree,
                 allowed_marking_colours=None):
        self.genus = ZZ(genus)
        self.markings = ZZ(markings)
        self.degree = ZZ(degree)
        if self.genus < 0 or self.markings < 0 or self.degree < 0:
            raise ValueError("genus, markings, and degree must be nonnegative")
        if allowed_marking_colours is None:
            self.allowed_marking_colours = tuple(
                (ZZ(0), ZZ(1), ZZ(2)) for _ in range(self.markings)
            )
        else:
            if len(allowed_marking_colours) != self.markings:
                raise ValueError("one fixed-point choice is needed per marking")
            normalized = []
            for choices in allowed_marking_colours:
                choices = tuple(sorted(set(ZZ(value) for value in choices)))
                if not choices or any(value not in (0, 1, 2)
                                      for value in choices):
                    raise ValueError(
                        "allowed marking colours are nonempty subsets of 0,1,2"
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

    @staticmethod
    def _proper_colours(vertex_count, endpoints):
        for colours in product(range(3), repeat=vertex_count):
            if all(colours[left] != colours[right] for left, right in endpoints):
                yield tuple(ZZ(colour) for colour in colours)

    def _degree_zero_graphs(self):
        for colour in range(3):
            if any(colour not in choices
                   for choices in self.allowed_marking_colours):
                continue
            graph = P2FixedLocusGraph(
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
                                provisional = P2FixedLocusGraph(
                                    genera, colours, marking_vertices,
                                    tuple((left, right, 1)
                                          for left, right in endpoints),
                                )
                                if any(provisional.vertex_type(vertex) == "invalid"
                                       for vertex in range(vertex_count)):
                                    continue
                                for degrees in _fl_positive_compositions(
                                        self.degree, edge_count):
                                    graph = P2FixedLocusGraph(
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


class TautologicalPolynomial(SageObject):
    r"""A truncated polynomial in vertex psi and lambda classes.

    A term key is ``(psi_powers, lambda_indices)``.  Coefficients lie in the
    full equivariant fraction field, while tautological codimension is kept
    separately so rational functions in torus weights never reach admcycles.
    """

    def __init__(self, field, marking_count, dimension, terms=None):
        self.field = field
        self.marking_count = ZZ(marking_count)
        self.dimension = ZZ(dimension)
        self.terms = {}
        for key, coefficient in (terms or {}).items():
            psi_powers, lambda_indices = key
            psi_powers = tuple(ZZ(value) for value in psi_powers)
            lambda_indices = tuple(sorted(
                ZZ(value) for value in lambda_indices if value
            ))
            if len(psi_powers) != self.marking_count:
                raise ValueError("a tautological monomial has the wrong arity")
            if sum(psi_powers) + sum(lambda_indices) > self.dimension:
                continue
            coefficient = self.field(coefficient)
            if coefficient:
                self.terms[(psi_powers, lambda_indices)] = coefficient

    @classmethod
    def one(cls, field, marking_count, dimension):
        return cls(field, marking_count, dimension, {
            ((ZZ.zero(),) * marking_count, ()): field.one(),
        })

    @classmethod
    def monomial(cls, field, marking_count, dimension,
                 psi_powers=None, lambda_indices=(), coefficient=1):
        psi_powers = ((ZZ.zero(),) * marking_count if psi_powers is None
                      else tuple(psi_powers))
        return cls(field, marking_count, dimension, {
            (psi_powers, tuple(lambda_indices)): coefficient,
        })

    def _compatible(self, other):
        return (isinstance(other, TautologicalPolynomial)
                and self.field == other.field
                and self.marking_count == other.marking_count
                and self.dimension == other.dimension)

    def __add__(self, other):
        if not self._compatible(other):
            raise TypeError("incompatible tautological polynomials")
        terms = dict(self.terms)
        for key, coefficient in other.terms.items():
            terms[key] = terms.get(key, self.field.zero()) + coefficient
            if not terms[key]:
                del terms[key]
        return TautologicalPolynomial(
            self.field, self.marking_count, self.dimension, terms
        )

    def __neg__(self):
        return self.scaled(-1)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not self._compatible(other):
            raise TypeError("incompatible tautological polynomials")
        terms = {}
        for (left_psi, left_lambda), left_coefficient in self.terms.items():
            for (right_psi, right_lambda), right_coefficient in other.terms.items():
                psi_powers = tuple(
                    left + right for left, right in zip(left_psi, right_psi)
                )
                lambda_indices = tuple(sorted(left_lambda + right_lambda))
                if sum(psi_powers) + sum(lambda_indices) > self.dimension:
                    continue
                key = psi_powers, lambda_indices
                terms[key] = (terms.get(key, self.field.zero())
                              + left_coefficient * right_coefficient)
        return TautologicalPolynomial(
            self.field, self.marking_count, self.dimension, terms
        )

    def __pow__(self, exponent):
        exponent = ZZ(exponent)
        if exponent < 0:
            raise ValueError("tautological powers must be nonnegative")
        answer = TautologicalPolynomial.one(
            self.field, self.marking_count, self.dimension
        )
        factor = self
        while exponent:
            if exponent % 2:
                answer = answer * factor
            exponent //= 2
            if exponent:
                factor = factor * factor
        return answer

    def scaled(self, scalar):
        scalar = self.field(scalar)
        return TautologicalPolynomial(
            self.field, self.marking_count, self.dimension,
            {key: scalar * coefficient
             for key, coefficient in self.terms.items()},
        )

    def integrate(self, genus, hodge_backend):
        answer = self.field.zero()
        for (psi_powers, lambda_indices), coefficient in self.terms.items():
            answer += coefficient * hodge_backend.integral(
                genus, psi_powers, lambda_indices
            )
        return self.field(answer)


class P2FixedLocusEvaluator(SageObject):
    r"""Evaluate ordinary or ``O(3)``-twisted ``P^2`` fixed graphs."""

    def __init__(self, rings=None, hodge_backend=None, include_twist=True,
                 lift_strategy="sparse"):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.field = self.rings.full_field
        self.lambdas = self.rings.base_weights()
        self.hodge = hodge_backend or AdmcyclesHodgeIntegralBackend()
        self.include_twist = bool(include_twist)
        self.lift_strategy = str(lift_strategy)
        if self.lift_strategy not in ("sparse", "standard"):
            raise ValueError("lift strategy must be sparse or standard")
        self.lift_planner = SparseEquivariantLiftPlanner()
        self._cache = {}

    def planned_insertions(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        if self.lift_strategy == "sparse":
            return self.lift_planner.plan(request.insertions)
        return request.insertions

    def fixed_graphs(self, request):
        insertions = self.planned_insertions(request)
        allowed = tuple(
            insertion.allowed_fixed_points() for insertion in insertions
        )
        return P2FixedLocusGraphEnumerator(
            request.genus, request.valence, request.degree,
            allowed_marking_colours=allowed,
        ).graphs()

    def graph_count(self, request):
        return ZZ(len(self.fixed_graphs(request)))

    def _lambda_dual(self, genus, weight, marking_count, dimension):
        answer = TautologicalPolynomial(
            self.field, marking_count, dimension
        )
        for index in range(genus + 1):
            term = TautologicalPolynomial.monomial(
                self.field, marking_count, dimension,
                lambda_indices=(() if index == 0 else (index,)),
                coefficient=(-1) ** index * weight ** (genus - index),
            )
            answer = answer + term
        return answer

    def _lambda_plus_inverse(self, genus, weight,
                             marking_count, dimension):
        one = TautologicalPolynomial.one(
            self.field, marking_count, dimension
        )
        positive = TautologicalPolynomial(
            self.field, marking_count, dimension
        )
        for index in range(1, genus + 1):
            positive = positive + TautologicalPolynomial.monomial(
                self.field, marking_count, dimension,
                lambda_indices=(index,), coefficient=weight ** (-index),
            )
        answer = one
        power = one
        for exponent in range(1, dimension + 1):
            power = power * positive
            answer = answer + power.scaled((-1) ** exponent)
        return answer.scaled(weight ** (-genus))

    def _flag_denominator(self, local_index, weight,
                          marking_count, dimension):
        answer = TautologicalPolynomial(
            self.field, marking_count, dimension
        )
        for power in range(dimension + 1):
            psi_powers = [ZZ.zero()] * marking_count
            psi_powers[local_index] = ZZ(power)
            answer = answer + TautologicalPolynomial.monomial(
                self.field, marking_count, dimension,
                psi_powers=tuple(psi_powers),
                coefficient=weight ** (-power - 1),
            )
        return answer

    def _edge_factor(self, graph, edge_index):
        left, right, degree = graph.edges[edge_index]
        left_colour = graph.vertex_colours[left]
        right_colour = graph.vertex_colours[right]
        left_weight = self.lambdas[left_colour]
        right_weight = self.lambdas[right_colour]
        difference = left_weight - right_weight
        factor = ((-1) ** degree * degree ** (2 * degree)
                  / (factorial(degree) ** 2
                     * difference ** (2 * degree)))
        third_colour = next(
            colour for colour in range(3)
            if colour not in (left_colour, right_colour)
        )
        third_weight = self.lambdas[third_colour]
        factor /= prod(
            (QQ(step) / degree * left_weight
             + QQ(degree - step) / degree * right_weight
             - third_weight)
            for step in range(degree + 1)
        )
        if self.include_twist:
            factor *= prod(
                self.rings.t
                - QQ(3 * degree - step) / degree * left_weight
                - QQ(step) / degree * right_weight
                for step in range(3 * degree + 1)
            )
        return self.field(factor)

    def _stable_vertex_factor(self, graph, vertex, insertions):
        genus = graph.vertex_genera[vertex]
        colour = graph.vertex_colours[vertex]
        edge_indices = graph.incident_edges(vertex)
        marking_indices = graph.markings_at_vertex(vertex)
        edge_valence = len(edge_indices)
        local_count = edge_valence + len(marking_indices)
        dimension = 3 * genus - 3 + local_count
        answer = TautologicalPolynomial.one(
            self.field, local_count, dimension
        )

        tangent_product = self.field.one()
        for other_colour in range(3):
            if other_colour == colour:
                continue
            tangent_weight = self.lambdas[colour] - self.lambdas[other_colour]
            tangent_product *= tangent_weight
            answer = answer * self._lambda_dual(
                genus, tangent_weight, local_count, dimension
            )
        answer = answer.scaled(tangent_product ** (edge_valence - 1))

        for local_index, edge_index in enumerate(edge_indices):
            weight = graph.flag_weight(edge_index, vertex, self.lambdas)
            answer = answer * self._flag_denominator(
                local_index, weight, local_count, dimension
            )

        if self.include_twist:
            fibre_weight = self.rings.t - 3 * self.lambdas[colour]
            answer = answer * self._lambda_plus_inverse(
                genus, fibre_weight, local_count, dimension
            )
            answer = answer.scaled(fibre_weight ** (1 - edge_valence))

        for offset, marking_index in enumerate(marking_indices):
            insertion = insertions[marking_index]
            psi_powers = [ZZ.zero()] * local_count
            psi_powers[edge_valence + offset] = insertion.psi_power
            answer = answer * TautologicalPolynomial.monomial(
                self.field, local_count, dimension,
                psi_powers=tuple(psi_powers),
                coefficient=insertion.restriction(colour, self.lambdas),
            )
        return answer.integrate(genus, self.hodge)

    def _unstable_vertex_factor(self, graph, vertex, insertions):
        vertex_type = graph.vertex_type(vertex)
        colour = graph.vertex_colours[vertex]
        edge_indices = graph.incident_edges(vertex)
        if vertex_type == "univalent":
            return self.field(graph.flag_weight(
                edge_indices[0], vertex, self.lambdas
            ))
        if vertex_type == "marked":
            marking_index = graph.markings_at_vertex(vertex)[0]
            insertion = insertions[marking_index]
            flag_weight = graph.flag_weight(
                edge_indices[0], vertex, self.lambdas
            )
            return self.field(
                insertion.restriction(colour, self.lambdas)
                * (-flag_weight) ** insertion.psi_power
            )
        if vertex_type == "nodal":
            tangent_product = prod(
                self.lambdas[colour] - self.lambdas[other_colour]
                for other_colour in range(3) if other_colour != colour
            )
            flag_sum = sum(
                graph.flag_weight(edge_index, vertex, self.lambdas)
                for edge_index in edge_indices
            )
            factor = tangent_product / flag_sum
            if self.include_twist:
                factor /= self.rings.t - 3 * self.lambdas[colour]
            return self.field(factor)
        raise ValueError("expected an unstable fixed-graph vertex")

    def graph_contribution(self, graph, insertions):
        if not graph.is_valid():
            raise ValueError("invalid fixed graph: %s" % graph.validation_errors())
        if len(insertions) != graph.marking_count:
            raise ValueError("the insertion count does not match the graph")
        answer = self.field(graph.localization_weight())
        for edge_index in range(graph.edge_count):
            answer *= self._edge_factor(graph, edge_index)
        for vertex in range(graph.vertex_count):
            if graph.vertex_type(vertex) == "stable":
                answer *= self._stable_vertex_factor(graph, vertex, insertions)
            else:
                answer *= self._unstable_vertex_factor(graph, vertex, insertions)
        return self.field(answer)

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        key = request, self.include_twist, self.lift_strategy
        if key not in self._cache:
            insertions = self.planned_insertions(request)
            graphs = self.fixed_graphs(request)
            if not graphs:
                unrestricted = P2FixedLocusGraphEnumerator(
                    request.genus, request.valence, request.degree
                ).graphs()
                if unrestricted:
                    self._cache[key] = self.field.zero()
                    return self._cache[key]
                raise ValueError(
                    "the stable-map type (g,n,d)=(%s,%s,%s) is unstable"
                    % (request.genus, request.valence, request.degree)
                )
            self._cache[key] = self.field(sum(
                self.graph_contribution(graph, insertions)
                for graph in graphs
            ))
        return self._cache[key]


class FullTwistedZeroVertexBackend(TwistedZeroVertexBackend):
    """Registry-compatible full zero-vertex backend using fixed loci."""

    def __init__(self, rings=None, hodge_backend=None,
                 lift_strategy="sparse"):
        hodge_backend = hodge_backend or AdmcyclesHodgeIntegralBackend()
        super().__init__(rings=rings, hodge_backend=hodge_backend)
        self.fixed_loci = P2FixedLocusEvaluator(
            rings=self.rings, hodge_backend=self.hodge, include_twist=True,
            lift_strategy=lift_strategy,
        )

    def provenance(self, request):
        if request in self._registry:
            return self._registry[request][1]
        return "P2 torus localization with admcycles vertex integration"

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        if request in self._registry:
            return self._registry[request][0]
        return self.fixed_loci.evaluate(request)
