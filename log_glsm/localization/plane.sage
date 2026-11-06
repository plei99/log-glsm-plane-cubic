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
and the fibre weight of ``(O(3))^vee`` at ``p_i`` is ``t-3*lambda_i``."""

from sage.all import *

from math import factorial
import json
import os
from log_glsm.core.conventions import (
    NONEQUIVARIANT_BASE_PATH,
    PlaneCubicCoefficientRing,
)
from log_glsm.core.hodge import (
    AdmcyclesHodgeIntegralBackend,
)
from log_glsm.core.tautological_polynomial import (
    TautologicalPolynomial,
)
from log_glsm.localization.plane_graphs import (
    P2FixedLocusGraphEnumerator,
    SparseEquivariantLiftPlanner,
)
from log_glsm.localization.twisted_plane import (
    TwistedZeroVertexBackend,
    TwistedZeroVertexRequest,
)
from log_glsm.localization.zero_store import (
    PersistentZeroVertexStore,
)

__all__ = ['FullTwistedZeroVertexBackend', 'P2FixedLocusEvaluator']


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
        r"""Invert the Hodge Euler polynomial using Mumford's relation.

        For Lambda_+(w)=sum lambda_i*w^(g-i), the identity
        c(E)c(E^vee)=1 implies Lambda_+(w)*Lambda_-(w)=w^(2g).
        Thus 1/Lambda_+(w)=Lambda_-(w)/w^(2g) in the Chow ring.
        This finite representative avoids generating and then integrating
        the vanishing higher Hodge products of a geometric series.
        """
        return self._lambda_dual(
            genus, weight, marking_count, dimension
        ).scaled(weight ** (-2 * genus))

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
        insertion_powers = (ZZ.zero(),) * edge_valence + tuple(
            insertion.psi_power for insertion in insertions
        )
        residual_dimension = dimension - sum(insertion_powers)
        insertion_coefficient = prod(
            insertion.restriction(colour, self.lambdas)
            for insertion in insertions
        )
        if residual_dimension < 0 or not insertion_coefficient:
            return TautologicalPolynomial(self.field, local_count, dimension)
        answer = TautologicalPolynomial.one(
            self.field, local_count, residual_dimension
        )

        tangent_product = self.field.one()
        for other_colour in range(3):
            if other_colour == colour:
                continue
            tangent_weight = self.lambdas[colour] - self.lambdas[other_colour]
            tangent_product *= tangent_weight
            answer = answer * self._lambda_dual(
                genus, tangent_weight, local_count, residual_dimension
            )
        answer = answer.scaled(tangent_product ** (edge_valence - 1))

        for local_index, weight in enumerate(flag_weights):
            answer = answer * self._flag_denominator(
                local_index, weight, local_count, residual_dimension
            )

        if self.include_twist:
            fibre_weight = self.rings.t - 3 * self.lambdas[colour]
            answer = answer * self._lambda_plus_inverse(
                genus, fibre_weight, local_count, residual_dimension
            )
            answer = answer.scaled(fibre_weight ** (1 - edge_valence))

        # Ordinary psi insertions already consume this much Chow degree.
        # Shift them back only after expanding the remaining characteristic
        # classes; no discarded monomial can reach the top degree.
        return TautologicalPolynomial(
            self.field, local_count, dimension,
            {(tuple(a + b for a, b in zip(psi, insertion_powers)), lambdas):
             coefficient * insertion_coefficient
             for (psi, lambdas), coefficient in answer.terms.items()},
        )

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
