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
import json
import os
import sqlite3


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


class PersistentZeroVertexStore(SageObject):
    """Concurrent SQLite store for normalized exact zero-vertex values."""

    FORMAT = "log-glsm-o3-zero-vertices-sqlite"
    VERSION = 1

    def __init__(self, path, include_twist, lift_strategy,
                 base_weight_specialization=None):
        self.path = os.path.abspath(path)
        directory = os.path.dirname(self.path)
        if directory and not os.path.isdir(directory):
            os.makedirs(directory)
        self.connection = sqlite3.connect(self.path, timeout=float(60))
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("PRAGMA busy_timeout=60000")
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS metadata ("
            "key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS zero_vertices ("
            "request_key TEXT PRIMARY KEY, request_json TEXT NOT NULL, "
            "value TEXT NOT NULL)"
        )
        expected = {
            "format": self.FORMAT,
            "version": str(int(self.VERSION)),
            "include_twist": str(bool(include_twist)),
            "lift_strategy": str(lift_strategy),
        }
        if base_weight_specialization is not None:
            if isinstance(base_weight_specialization, str):
                expected["base_weight_specialization"] = json.dumps(
                    base_weight_specialization
                )
                if base_weight_specialization == "nonequivariant":
                    expected["nonequivariant_path"] = json.dumps(
                        NONEQUIVARIANT_BASE_PATH
                    )
            else:
                expected["base_weight_specialization"] = json.dumps([
                    str(QQ(value)) for value in base_weight_specialization
                ])
        metadata = dict(self.connection.execute(
            "SELECT key, value FROM metadata"
        ).fetchall())
        if metadata:
            if any(metadata.get(key) != value
                   for key, value in expected.items()):
                raise ValueError("incompatible O(3) zero-vertex SQLite cache")
        else:
            self.connection.executemany(
                "INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)",
                tuple(expected.items()),
            )
            self.connection.commit()

    @staticmethod
    def request_key(request):
        return json.dumps(
            request.localization_record(),
            sort_keys=True, separators=(",", ":")
        )

    def get(self, request):
        row = self.connection.execute(
            "SELECT value FROM zero_vertices WHERE request_key = ?",
            (self.request_key(request),),
        ).fetchone()
        return None if row is None else row[0]

    def put(self, request, value):
        self.connection.execute(
            "INSERT OR IGNORE INTO zero_vertices"
            "(request_key, request_json, value) VALUES (?, ?, ?)",
            (
                self.request_key(request),
                json.dumps(request.to_record(), sort_keys=True),
                str(value),
            ),
        )
        self.connection.commit()

    def __len__(self):
        return int(self.connection.execute(
            "SELECT COUNT(*) FROM zero_vertices"
        ).fetchone()[0])

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None


class P2FixedLocusEvaluator(SageObject):
    r"""Evaluate ordinary or ``O(3)``-twisted ``P^2`` fixed graphs."""

    CACHE_FORMAT = "log-glsm-o3-zero-vertices"
    CACHE_VERSION = 1

    def __init__(self, rings=None, hodge_backend=None, include_twist=True,
                 lift_strategy="sparse", cache_path=None, autosave=False,
                 base_weight_specialization=None):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.field = self.rings.full_field
        if base_weight_specialization is None:
            self.base_weight_specialization = None
            self.lambdas = self.rings.base_weights()
        elif base_weight_specialization == "nonequivariant":
            # Keep one scale epsilon while evaluating individual fixed loci,
            # then let epsilon tend to zero only after their complete sum.
            # A linear ray such as epsilon*(0,1,3) is resonant: at a nodal
            # unstable fixed vertex, two flag weights can cancel exactly
            # (already in degree three).  The curved path below approaches
            # the same nonequivariant point without imposing any constant
            # rational relation among the three weights.
            self.base_weight_specialization = "nonequivariant"
            epsilon = self.rings.full_lambda0
            self.lambdas = (
                self.field.zero(), epsilon, epsilon ** 2
            )
        elif isinstance(base_weight_specialization, str):
            raise ValueError(
                "base weights must be 'nonequivariant', None, or a triple"
            )
        else:
            if len(base_weight_specialization) != 3 \
                    or len(set(base_weight_specialization)) != 3:
                raise ValueError("base weights must be three distinct values")
            self.base_weight_specialization = tuple(
                QQ(value) for value in base_weight_specialization
            )
            self.lambdas = tuple(
                self.field(value) for value in self.base_weight_specialization
            )
        self.hodge = hodge_backend or AdmcyclesHodgeIntegralBackend()
        self.include_twist = bool(include_twist)
        self.lift_strategy = str(lift_strategy)
        if self.lift_strategy not in ("sparse", "standard"):
            raise ValueError("lift strategy must be sparse or standard")
        self.lift_planner = SparseEquivariantLiftPlanner()
        self._cache = {}
        self._fixed_graph_cache = {}
        self._edge_factor_cache = {}
        self._stable_vertex_value_cache = {}
        self._cache_hits = ZZ.zero()
        self._cache_misses = ZZ.zero()
        self.cache_path = (
            None if cache_path is None else os.path.abspath(cache_path)
        )
        self._sqlite_store = None
        self.autosave = bool(autosave)
        if self.autosave and self.cache_path is None:
            raise ValueError("autosave requires a zero-vertex cache path")
        if self.cache_path is not None:
            if self.cache_path.endswith((".sqlite", ".sqlite3", ".db")):
                self._sqlite_store = PersistentZeroVertexStore(
                    self.cache_path, self.include_twist, self.lift_strategy,
                    self.base_weight_specialization,
                )
            elif os.path.exists(self.cache_path):
                self.load_cache(self.cache_path)

    def _cache_key(self, request):
        return request, self.include_twist, self.lift_strategy

    def is_cached(self, request):
        scale, normalized = request.normalized_scales()
        if not scale or self._cache_key(normalized) in self._cache:
            return True
        if self._sqlite_store is not None:
            stored = self._sqlite_store.get(normalized)
            if stored is not None:
                self._cache[self._cache_key(normalized)] = self.field(stored)
                return True
        return False

    def cache_info(self):
        return {
            "entries": len(self._cache),
            "hits": int(self._cache_hits),
            "misses": int(self._cache_misses),
            "fixed_graph_families": len(self._fixed_graph_cache),
            "edge_factors": len(self._edge_factor_cache),
            "stable_vertex_values": len(self._stable_vertex_value_cache),
            "hodge_integrals": len(getattr(self.hodge, "_cache", {})),
            "hodge_cache": (
                self.hodge.cache_info()
                if hasattr(self.hodge, "cache_info") else None
            ),
            "path": self.cache_path,
            "persistent_entries": (
                int(0) if self._sqlite_store is None
                else int(len(self._sqlite_store))
            ),
            "cache_backend": (
                "sqlite" if self._sqlite_store is not None else "json"
            ),
        }

    def _remember_value(self, request, value):
        key = self._cache_key(request)
        self._cache[key] = self.field(value)
        if self._sqlite_store is not None:
            self._sqlite_store.put(request, self._cache[key])
        elif self.autosave:
            self.save_cache()
        return self._cache[key]

    def save_cache(self, path=None):
        r"""Atomically persist every exact zero-vertex value computed so far."""
        path = self.cache_path if path is None else os.path.abspath(path)
        if path is None:
            raise ValueError("a zero-vertex cache path is required")
        if path.endswith((".sqlite", ".sqlite3", ".db")):
            if self._sqlite_store is None or self._sqlite_store.path != path:
                self._sqlite_store = PersistentZeroVertexStore(
                    path, self.include_twist, self.lift_strategy,
                    self.base_weight_specialization,
                )
            for (request, include_twist, lift_strategy), value \
                    in self._cache.items():
                if include_twist == self.include_twist \
                        and lift_strategy == self.lift_strategy:
                    self._sqlite_store.put(request, value)
            self.cache_path = path
            return path
        directory = os.path.dirname(path)
        if directory and not os.path.isdir(directory):
            os.makedirs(directory)
        records = []
        for (request, include_twist, lift_strategy), value in self._cache.items():
            if include_twist != self.include_twist \
                    or lift_strategy != self.lift_strategy:
                continue
            records.append({
                "request": request.to_record(),
                "value": str(value),
            })
        records.sort(key=lambda record: json.dumps(
            record["request"], sort_keys=True
        ))
        payload = {
            "format": self.CACHE_FORMAT,
            "version": int(self.CACHE_VERSION),
            "include_twist": self.include_twist,
            "lift_strategy": self.lift_strategy,
            "base_weight_specialization": (
                None if self.base_weight_specialization is None
                else (
                    self.base_weight_specialization
                    if isinstance(self.base_weight_specialization, str)
                    else [str(value) for value
                          in self.base_weight_specialization]
                )
            ),
            "nonequivariant_path": (
                NONEQUIVARIANT_BASE_PATH
                if self.base_weight_specialization == "nonequivariant"
                else None
            ),
            "values": records,
        }
        temporary = path + ".tmp"
        with open(temporary, "w") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary, path)
        self.cache_path = path
        return path

    def load_cache(self, path=None):
        r"""Load exact values, rejecting caches for different conventions."""
        path = self.cache_path if path is None else os.path.abspath(path)
        if path is None:
            raise ValueError("a zero-vertex cache path is required")
        if path.endswith((".sqlite", ".sqlite3", ".db")):
            self._sqlite_store = PersistentZeroVertexStore(
                path, self.include_twist, self.lift_strategy,
                self.base_weight_specialization,
            )
            self.cache_path = path
            return len(self._sqlite_store)
        with open(path) as stream:
            payload = json.load(stream)
        if payload.get("format") != self.CACHE_FORMAT \
                or payload.get("version") != int(self.CACHE_VERSION):
            raise ValueError("unsupported O(3) zero-vertex cache format")
        if bool(payload.get("include_twist")) != self.include_twist:
            raise ValueError("zero-vertex cache has the wrong twist convention")
        if payload.get("lift_strategy") != self.lift_strategy:
            raise ValueError("zero-vertex cache has the wrong lift strategy")
        expected_weights = (
            None if self.base_weight_specialization is None
            else (
                self.base_weight_specialization
                if isinstance(self.base_weight_specialization, str)
                else [str(value) for value
                      in self.base_weight_specialization]
            )
        )
        if payload.get("base_weight_specialization") != expected_weights:
            raise ValueError("zero-vertex cache has the wrong base weights")
        expected_path = (
            NONEQUIVARIANT_BASE_PATH
            if self.base_weight_specialization == "nonequivariant" else None
        )
        if payload.get("nonequivariant_path") != expected_path:
            raise ValueError(
                "zero-vertex cache has the wrong nonequivariant path"
            )
        for record in payload.get("values", ()):
            request = TwistedZeroVertexRequest.from_record(record["request"])
            scale, normalized = request.normalized_scales()
            if scale:
                self._cache[self._cache_key(normalized)] = self.field(
                    self.field(record["value"]) / scale
                )
        self.cache_path = path
        return len(payload.get("values", ()))

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
        key = (
            ZZ(request.genus), ZZ(request.valence), ZZ(request.degree),
            allowed,
        )
        if key not in self._fixed_graph_cache:
            self._fixed_graph_cache[key] = P2FixedLocusGraphEnumerator(
                request.genus, request.valence, request.degree,
                allowed_marking_colours=allowed,
            ).graphs()
        return self._fixed_graph_cache[key]

    def graph_count(self, request):
        return ZZ(len(self.fixed_graphs(request)))

    def _constant_map_value(self, request):
        r"""Evaluate any stable degree-zero vertex without fixed graphs.

        ``Mbar_(g,n)(P2,0)`` is ``Mbar_(g,n) x P2``.  We form the Hodge
        polynomial at the three auxiliary fixed points, add those
        polynomials first, and only then evaluate their distinct monomials.
        """
        if 2 * request.genus - 2 + request.valence <= 0:
            raise ValueError("the degree-zero stable-map type is unstable")
        insertions = self.planned_insertions(request)
        dimension = 3 * request.genus - 3 + request.valence
        allowed_colours = set(range(3))
        for insertion in insertions:
            allowed_colours.intersection_update(
                int(colour) for colour in insertion.allowed_fixed_points()
            )
        polynomials = tuple(
            self._stable_vertex_polynomial(
                request.genus, colour, (), insertions
            )
            for colour in sorted(allowed_colours)
        )
        if len(polynomials) < 3:
            # Sparse lifts already removed at least one fixed point.  In that
            # case adding rational-function coefficients before integration
            # costs more than the saved Hodge-cache lookups.
            return self.field(sum(
                polynomial.integrate(request.genus, self.hodge)
                for polynomial in polynomials
            ))
        aggregate = TautologicalPolynomial(
            self.field, request.valence, dimension
        )
        for polynomial in polynomials:
            aggregate = aggregate + polynomial
        return aggregate.integrate(request.genus, self.hodge)

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
        cache_key = left_colour, right_colour, ZZ(degree)
        if cache_key in self._edge_factor_cache:
            return self._edge_factor_cache[cache_key]
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
        self._edge_factor_cache[cache_key] = self.field(factor)
        return self._edge_factor_cache[cache_key]

    def _stable_vertex_polynomial(self, genus, colour, flag_weights,
                                  insertions):
        flag_weights = tuple(self.field(weight) for weight in flag_weights)
        insertions = tuple(insertions)
        edge_valence = len(flag_weights)
        local_count = edge_valence + len(insertions)
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

        for local_index, weight in enumerate(flag_weights):
            answer = answer * self._flag_denominator(
                local_index, weight, local_count, dimension
            )

        if self.include_twist:
            fibre_weight = self.rings.t - 3 * self.lambdas[colour]
            answer = answer * self._lambda_plus_inverse(
                genus, fibre_weight, local_count, dimension
            )
            answer = answer.scaled(fibre_weight ** (1 - edge_valence))

        for offset, insertion in enumerate(insertions):
            psi_powers = [ZZ.zero()] * local_count
            psi_powers[edge_valence + offset] = insertion.psi_power
            answer = answer * TautologicalPolynomial.monomial(
                self.field, local_count, dimension,
                psi_powers=tuple(psi_powers),
                coefficient=insertion.restriction(colour, self.lambdas),
            )
        return answer

    def _stable_vertex_factor(self, graph, vertex, insertions):
        genus = graph.vertex_genera[vertex]
        colour = graph.vertex_colours[vertex]
        edge_indices = graph.incident_edges(vertex)
        flag_weights = tuple(
            graph.flag_weight(edge_index, vertex, self.lambdas)
            for edge_index in edge_indices
        )
        local_insertions = tuple(
            insertions[index] for index in graph.markings_at_vertex(vertex)
        )
        cache_key = (
            ZZ(genus), ZZ(colour), flag_weights,
            tuple(item.localization_signature() for item in local_insertions),
        )
        if cache_key not in self._stable_vertex_value_cache:
            polynomial = self._stable_vertex_polynomial(
                genus, colour, flag_weights, local_insertions
            )
            self._stable_vertex_value_cache[cache_key] = polynomial.integrate(
                genus, self.hodge
            )
        return self._stable_vertex_value_cache[cache_key]

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
        scale, request = request.normalized_scales()
        if not scale:
            return self.field.zero()
        key = self._cache_key(request)
        if key in self._cache:
            self._cache_hits += 1
            return self.field(scale * self._cache[key])
        if self._sqlite_store is not None:
            stored = self._sqlite_store.get(request)
            if stored is not None:
                self._cache_hits += 1
                self._cache[key] = self.field(stored)
                return self.field(scale * self._cache[key])
        self._cache_misses += 1
        if request.degree == 0:
            value = self._remember_value(
                request, self._constant_map_value(request)
            )
            return self.field(scale * value)
        if key not in self._cache:
            insertions = self.planned_insertions(request)
            graphs = self.fixed_graphs(request)
            if not graphs:
                unrestricted = P2FixedLocusGraphEnumerator(
                    request.genus, request.valence, request.degree
                ).graphs()
                if unrestricted:
                    value = self._remember_value(request, self.field.zero())
                    return self.field(scale * value)
                raise ValueError(
                    "the stable-map type (g,n,d)=(%s,%s,%s) is unstable"
                    % (request.genus, request.valence, request.degree)
                )
            value = self.field(sum(
                self.graph_contribution(graph, insertions)
                for graph in graphs
            ))
            value = self._remember_value(request, value)
        return self.field(scale * value)


class FullTwistedZeroVertexBackend(TwistedZeroVertexBackend):
    """Registry-compatible full zero-vertex backend using fixed loci."""

    def __init__(self, rings=None, hodge_backend=None,
                 lift_strategy="sparse", cache_path=None, autosave=False,
                 hodge_cache_path=None, base_weight_specialization=None):
        hodge_backend = hodge_backend or AdmcyclesHodgeIntegralBackend(
            cache_path=hodge_cache_path
        )
        super().__init__(rings=rings, hodge_backend=hodge_backend)
        self.fixed_loci = P2FixedLocusEvaluator(
            rings=self.rings, hodge_backend=self.hodge, include_twist=True,
            lift_strategy=lift_strategy, cache_path=cache_path,
            autosave=autosave,
            base_weight_specialization=base_weight_specialization,
        )

    def is_cached(self, request):
        return (
            request in self._registry
            or (request.genus == 0 and request.degree == 0)
            or self.fixed_loci.is_cached(request)
        )

    def cache_info(self):
        return self.fixed_loci.cache_info()

    def save_cache(self, path=None):
        return self.fixed_loci.save_cache(path)

    def load_cache(self, path=None):
        return self.fixed_loci.load_cache(path)

    def provenance(self, request):
        if request in self._registry:
            return self._registry[request][1]
        if request.degree == 0:
            return "direct constant-map Hodge polynomial on Mbar_g,n x P2"
        return "P2 torus localization with admcycles vertex integration"

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        if request in self._registry:
            return self._registry[request][0]
        return self.fixed_loci.evaluate(request)
