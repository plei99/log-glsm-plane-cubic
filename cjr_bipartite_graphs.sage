r"""
Finite enumeration of compact-type CJR localization graphs for a plane cubic.

This module implements the decorated bipartite graphs of CJR III, Section 3.4,
in the specialization

    X = P^2,    E = O(3),    r = 1,

with all ordinary markings of compact type.  Thus every marking lies on a
zero-level vertex and all inertia sectors are trivial.  An infinity vertex V
satisfies the specialized balancing equation

    sum_{e incident to V} c(e) = 2*g(V) - 2 + val_E(V) - 3*d(V).

The implementation is intended to be the combinatorial front end of the
localization equation provider.  It enumerates isomorphism classes, computes
automorphism orders, and retains the full contact-order decoration.

The public entry point is ``PlaneCubicGraphEnumerator(g, n, d).graphs()``.
"""

from collections import Counter, defaultdict, deque
from itertools import product
from math import factorial


def _weak_compositions(total, length):
    """Yield tuples of ``length`` nonnegative integers summing to ``total``."""
    total = Integer(total)
    length = Integer(length)
    if length < 0 or total < 0:
        return
    if length == 0:
        if total == 0:
            yield tuple()
        return
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in _weak_compositions(total - first, length - 1):
            yield (Integer(first),) + tail


def _positive_compositions(total, length):
    """Yield tuples of ``length`` positive integers summing to ``total``."""
    total = Integer(total)
    length = Integer(length)
    if length == 0:
        if total == 0:
            yield tuple()
        return
    if total < length:
        return
    for shifted in _weak_compositions(total - length, length):
        yield tuple(x + 1 for x in shifted)


def _inverse_permutation(order):
    """For a new-to-old ordering, return the old-to-new index map."""
    answer = [None] * len(order)
    for new, old in enumerate(order):
        answer[old] = new
    return tuple(answer)


def _connected_bipartite(zero_count, infinity_count, pairs):
    """Test connectedness of a bipartite multigraph from its edge pairs."""
    vertex_count = zero_count + infinity_count
    if vertex_count == 1:
        return len(pairs) == 0
    if vertex_count == 0 or len(pairs) == 0:
        return False

    adjacency = [[] for _ in range(vertex_count)]
    for zero, infinity in pairs:
        left = zero
        right = zero_count + infinity
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


class BipartiteLocalizationGraph(SageObject):
    r"""A connected decorated CJR graph in the compact plane-cubic sector.

    ``edges`` consists of triples ``(zero_index, infinity_index, contact)``.
    Markings are labeled; ``marking_vertices[j]`` is the zero vertex carrying
    marking ``j``.  Edge instances are unlabeled.
    """

    def __init__(self, zero_genera, zero_degrees,
                 infinity_genera, infinity_degrees,
                 marking_vertices=(), edges=()):
        self.zero_genera = tuple(Integer(x) for x in zero_genera)
        self.zero_degrees = tuple(Integer(x) for x in zero_degrees)
        self.infinity_genera = tuple(Integer(x) for x in infinity_genera)
        self.infinity_degrees = tuple(Integer(x) for x in infinity_degrees)
        self.marking_vertices = tuple(Integer(x) for x in marking_vertices)
        self.edges = tuple(sorted(
            (Integer(z), Integer(i), Integer(c)) for z, i, c in edges
        ))
        self._canonical_key_cache = None
        self._automorphism_order_cache = None

    @property
    def zero_count(self):
        return len(self.zero_genera)

    @property
    def infinity_count(self):
        return len(self.infinity_genera)

    @property
    def edge_count(self):
        return len(self.edges)

    @property
    def vertex_count(self):
        return self.zero_count + self.infinity_count

    def zero_edge_valences(self):
        answer = [0] * self.zero_count
        for zero, _, _ in self.edges:
            if 0 <= zero < self.zero_count:
                answer[zero] += 1
        return tuple(answer)

    def infinity_edge_valences(self):
        answer = [0] * self.infinity_count
        for _, infinity, _ in self.edges:
            if 0 <= infinity < self.infinity_count:
                answer[infinity] += 1
        return tuple(answer)

    def zero_marking_valences(self):
        answer = [0] * self.zero_count
        for zero in self.marking_vertices:
            if 0 <= zero < self.zero_count:
                answer[zero] += 1
        return tuple(answer)

    def first_betti_number(self):
        if self.vertex_count == 0:
            return Integer(-1)
        return Integer(self.edge_count - self.vertex_count + 1)

    def total_genus(self):
        return Integer(sum(self.zero_genera) + sum(self.infinity_genera)
                       + self.first_betti_number())

    def total_degree(self):
        return Integer(sum(self.zero_degrees) + sum(self.infinity_degrees))

    def markings_at_zero_vertex(self, zero):
        return tuple(j for j, vertex in enumerate(self.marking_vertices)
                     if vertex == zero)

    def zero_vertex_type(self, zero):
        r"""Return stable/nonspecial/marked/nodal, or ``invalid``.

        The three unstable types are exactly Definition 3.11 of CJR III.
        """
        genus = self.zero_genera[zero]
        degree = self.zero_degrees[zero]
        edge_valence = self.zero_edge_valences()[zero]
        marking_valence = self.zero_marking_valences()[zero]
        valence = edge_valence + marking_valence
        if genus != 0 or degree != 0 or valence > 2:
            return "stable"
        if edge_valence == 1 and marking_valence == 0:
            return "nonspecial"
        if edge_valence == 1 and marking_valence == 1:
            return "marked"
        if edge_valence == 2 and marking_valence == 0:
            return "nodal"
        return "invalid"

    def infinity_balance(self, infinity):
        contact_sum = sum(c for _, i, c in self.edges if i == infinity)
        valence = self.infinity_edge_valences()[infinity]
        expected = (2 * self.infinity_genera[infinity] - 2 + valence
                    - 3 * self.infinity_degrees[infinity])
        return Integer(contact_sum), Integer(expected)

    def infinity_contact_profile(self, infinity):
        """Return the CJR infinity half-edge contacts, hence negative values."""
        if infinity < 0 or infinity >= self.infinity_count:
            raise IndexError("infinity vertex index is out of range")
        return tuple(sorted(-c for _, i, c in self.edges if i == infinity))

    def infinity_vertex_data(self, infinity):
        """Return ``(genus, ambient_degree, negative_contact_profile)``."""
        return (
            self.infinity_genera[infinity],
            self.infinity_degrees[infinity],
            self.infinity_contact_profile(infinity),
        )

    def all_infinity_vertex_data(self):
        """Return the undecorated effective-vertex data needed by this graph."""
        return tuple(self.infinity_vertex_data(i)
                     for i in range(self.infinity_count))

    def validation_errors(self):
        errors = []
        if len(self.zero_genera) != len(self.zero_degrees):
            errors.append("zero genus/degree arrays have different lengths")
        if len(self.infinity_genera) != len(self.infinity_degrees):
            errors.append("infinity genus/degree arrays have different lengths")
        if errors:
            return errors
        if self.vertex_count == 0:
            errors.append("the graph has no vertices")
            return errors

        decorations = (self.zero_genera + self.zero_degrees
                       + self.infinity_genera + self.infinity_degrees)
        if any(x < 0 for x in decorations):
            errors.append("genera and degrees must be nonnegative")
        incidence_error = False
        for zero in self.marking_vertices:
            if zero < 0 or zero >= self.zero_count:
                errors.append("a marking is not attached to a zero vertex")
                incidence_error = True
        for zero, infinity, contact in self.edges:
            if zero < 0 or zero >= self.zero_count:
                errors.append("an edge has an invalid zero endpoint")
                incidence_error = True
            if infinity < 0 or infinity >= self.infinity_count:
                errors.append("an edge has an invalid infinity endpoint")
                incidence_error = True
            if contact <= 0:
                errors.append("contact orders must be positive")
        if incidence_error:
            return errors

        endpoint_pairs = [(z, i) for z, i, _ in self.edges]
        if not _connected_bipartite(
                self.zero_count, self.infinity_count, endpoint_pairs):
            errors.append("the graph is not connected")

        for zero in range(self.zero_count):
            vertex_type = self.zero_vertex_type(zero)
            if vertex_type == "invalid":
                errors.append("zero vertex %s has no allowed stability type" % zero)
            if vertex_type == "nonspecial":
                contacts = [c for z, _, c in self.edges if z == zero]
                if len(contacts) != 1 or contacts[0] <= 1:
                    errors.append(
                        "nonspecial zero vertex %s must have contact > 1" % zero
                    )

        # With compact markings, Lemma 3.12 leaves no unstable infinity
        # vertices: the only possible unstable infinity type is marked.
        for infinity in range(self.infinity_count):
            genus = self.infinity_genera[infinity]
            degree = self.infinity_degrees[infinity]
            valence = self.infinity_edge_valences()[infinity]
            if genus == 0 and degree == 0 and valence <= 2:
                errors.append("infinity vertex %s is unstable" % infinity)
            actual, expected = self.infinity_balance(infinity)
            if actual != expected:
                errors.append(
                    "infinity vertex %s violates balance: %s != %s"
                    % (infinity, actual, expected)
                )
        return errors

    def is_valid(self):
        return len(self.validation_errors()) == 0

    def _transformed_key(self, zero_order, infinity_order):
        zero_inverse = _inverse_permutation(zero_order)
        infinity_inverse = _inverse_permutation(infinity_order)
        zero_data = tuple(
            (self.zero_genera[old], self.zero_degrees[old],
             self.markings_at_zero_vertex(old))
            for old in zero_order
        )
        infinity_data = tuple(
            (self.infinity_genera[old], self.infinity_degrees[old])
            for old in infinity_order
        )
        edges = tuple(sorted(
            (zero_inverse[z], infinity_inverse[i], c)
            for z, i, c in self.edges
        ))
        return zero_data, infinity_data, edges

    def _colored_incidence_graph(self):
        r"""Encode vertices, edge instances, and decorations as a colored graph.

        Turning every localization edge into an incidence vertex handles
        parallel edges and their permutations without any special cases.
        """
        graph = Graph(multiedges=False, loops=False)
        buckets = defaultdict(list)
        for zero in range(self.zero_count):
            node = ("z", zero)
            graph.add_vertex(node)
            color = ("zero", self.zero_genera[zero], self.zero_degrees[zero],
                     self.markings_at_zero_vertex(zero))
            buckets[color].append(node)
        for infinity in range(self.infinity_count):
            node = ("i", infinity)
            graph.add_vertex(node)
            color = ("infinity", self.infinity_genera[infinity],
                     self.infinity_degrees[infinity])
            buckets[color].append(node)
        for edge_index, (zero, infinity, contact) in enumerate(self.edges):
            node = ("e", edge_index)
            graph.add_vertex(node)
            graph.add_edge(node, ("z", zero))
            graph.add_edge(node, ("i", infinity))
            buckets[("edge", contact)].append(node)
        colors = tuple(sorted(buckets))
        partition = [buckets[color] for color in colors]
        return graph, colors, partition

    def identity_key(self):
        return self._transformed_key(
            tuple(range(self.zero_count)),
            tuple(range(self.infinity_count)),
        )

    def canonical_key(self):
        if self._canonical_key_cache is None:
            graph, colors, partition = self._colored_incidence_graph()
            canonical = graph.canonical_label(partition=partition)
            color_profile = tuple(
                (color, len(block))
                for color, block in zip(colors, partition)
            )
            adjacency = tuple(canonical.edges(labels=False, sort=True))
            self._canonical_key_cache = color_profile, adjacency
        return self._canonical_key_cache

    def vertex_automorphism_order(self):
        return self.automorphism_order() // self.edge_automorphism_order()

    def edge_automorphism_order(self):
        multiplicities = Counter(self.edges)
        return Integer(prod(factorial(m) for m in multiplicities.values()))

    def automorphism_order(self):
        if self._automorphism_order_cache is None:
            graph, _, partition = self._colored_incidence_graph()
            self._automorphism_order_cache = Integer(
                graph.automorphism_group(partition=partition).order()
            )
        return self._automorphism_order_cache

    def localization_weight(self):
        """Return the graph factor ``1 / |Aut(Gamma)|``."""
        return QQ(1) / self.automorphism_order()

    def __eq__(self, other):
        return (isinstance(other, BipartiteLocalizationGraph)
                and self.canonical_key() == other.canonical_key())

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.canonical_key())

    def signature(self):
        """A compact, deterministic description useful in reports and tests."""
        zero = []
        for vertex in range(self.zero_count):
            zero.append((
                self.zero_genera[vertex], self.zero_degrees[vertex],
                self.zero_vertex_type(vertex),
                self.markings_at_zero_vertex(vertex),
            ))
        infinity = []
        for vertex in range(self.infinity_count):
            infinity.append((
                self.infinity_genera[vertex], self.infinity_degrees[vertex],
            ))
        return (tuple(zero), tuple(infinity), self.edges,
                self.automorphism_order())

    def _repr_(self):
        return ("BipartiteLocalizationGraph(zero=%s, infinity=%s, "
                "edges=%s, aut=%s)"
                % (self.signature()[0], self.signature()[1], self.edges,
                   self.automorphism_order()))


class PlaneCubicGraphEnumerator(SageObject):
    r"""Enumerate all compact-type CJR graphs with fixed ``(g,n,d)``.

    The search is finite.  For graphs using both levels we use the bound

        |V| <= 4*g + 2*n + d - 1.

    Here is a useful bookkeeping proof.  Write I for the number of infinity
    vertices; A,B,C for nonspecial, marked, and nodal unstable zero vertices;
    P for stable zero vertices of positive genus or degree; and T for stable
    rational degree-zero zero vertices.  Balance gives I <= g and
    A <= 2*g-2*I.  Also B <= n and P <= g+d.  Stability and the Euler formula
    give C+2*T <= g+I-1+n.  Adding these inequalities gives the displayed
    bound.  One-level exceptional graphs are added separately.  Canonical
    deduplication removes isomorphic copies generated by labeled vertices and
    edge instances.
    """

    def __init__(self, genus, markings, degree):
        self.genus = Integer(genus)
        self.markings = Integer(markings)
        self.degree = Integer(degree)
        if self.genus < 0 or self.markings < 0 or self.degree < 0:
            raise ValueError("genus, markings, and degree must be nonnegative")
        self._cache = None
        self._row_cache = {}
        self._degree_cache = {}

    def mixed_vertex_bound(self):
        if self.genus == 0:
            return Integer(0)
        return 4 * self.genus + 2 * self.markings + self.degree - 1

    def _one_level_graphs(self):
        # The CJR vertex classification is made inside globally X-stable
        # data.  Reject the exceptional globally unstable degree-zero types
        # before applying it; otherwise (g,n,d)=(1,0,0) produces spurious
        # zero-only and infinity-only graphs even though the connected
        # stable-map probe itself does not exist.
        if self.degree == 0 and 2 * self.genus - 2 + self.markings <= 0:
            return tuple()
        candidates = []
        zero_only = BipartiteLocalizationGraph(
            (self.genus,), (self.degree,), (), (),
            tuple(0 for _ in range(self.markings)), (),
        )
        if zero_only.is_valid():
            candidates.append(zero_only)

        if self.markings == 0:
            infinity_only = BipartiteLocalizationGraph(
                (), (), (self.genus,), (self.degree,), (), (),
            )
            if infinity_only.is_valid():
                candidates.append(infinity_only)
        return candidates

    def _rows_with_sum(self, infinity_count, row_sum):
        """Cache the incidence rows of one zero vertex with a fixed valence."""
        key = (infinity_count, row_sum)
        if key not in self._row_cache:
            self._row_cache[key] = tuple(sorted(
                _weak_compositions(row_sum, infinity_count), reverse=True
            ))
        return self._row_cache[key]

    def _canonical_incidence_matrices(self, zero_count, infinity_count,
                                      edge_count):
        r"""Yield incidence matrices up to relabeling of the zero vertices.

        Row ``z`` records how many edges join zero vertex ``z`` to each
        infinity vertex.  Rows are emitted in non-increasing ``(sum, row)``
        order, so exactly one representative of each orbit under the symmetric
        group on zero vertices is produced.

        This loses no isomorphism class.  Given any decorated graph, let
        ``pi`` be a permutation sorting the rows of its incidence matrix.
        Applying ``pi`` to the genus, degree, and marking decorations as well
        gives an isomorphic decorated graph whose matrix is canonical, and
        ``_mixed_graphs`` enumerates every decoration of every emitted matrix.

        Each row sum is positive: a zero vertex with no incident edge would
        disconnect a mixed-level graph, and disconnected graphs are rejected.
        """
        def extend(rows_left, remaining, previous):
            if rows_left == 0:
                if remaining == 0:
                    yield tuple()
                return
            # Every remaining row still needs at least one edge.
            if remaining < rows_left:
                return
            largest = remaining - (rows_left - 1)
            if previous is not None:
                largest = min(largest, previous[0])
            for row_sum in range(largest, 0, -1):
                # Later rows have sum at most row_sum, so this bounds the
                # whole tail; smaller sums only make it worse.
                if remaining > row_sum * rows_left:
                    break
                for row in self._rows_with_sum(infinity_count, row_sum):
                    key = (row_sum, row)
                    if previous is not None and key > previous:
                        continue
                    for tail in extend(rows_left - 1, remaining - row_sum,
                                       key):
                        yield (row,) + tail

        return extend(zero_count, edge_count, None)

    def _endpoint_multigraphs(self, zero_count, infinity_count, edge_count):
        for matrix in self._canonical_incidence_matrices(
                zero_count, infinity_count, edge_count):
            # An infinity vertex with no incident edge disconnects the graph.
            if any(not any(row[infinity] for row in matrix)
                   for infinity in range(infinity_count)):
                continue
            pairs = tuple(
                (zero, infinity)
                for zero, row in enumerate(matrix)
                for infinity, multiplicity in enumerate(row)
                for _ in range(multiplicity)
            )
            if _connected_bipartite(zero_count, infinity_count, pairs):
                yield pairs

    def _genus_decorations(self, vertex_count, infinity_count, first_betti):
        # Infinity stability plus balance forces g_infinity >= 1: for g=0,
        # balance would give sum(c_e) <= valence-2, impossible for c_e >= 1.
        residual = self.genus - first_betti - infinity_count
        if residual < 0:
            return
        for excess in _weak_compositions(residual, vertex_count):
            zero_count = vertex_count - infinity_count
            zero_genera = excess[:zero_count]
            infinity_genera = tuple(
                1 + x for x in excess[zero_count:]
            )
            yield zero_genera, infinity_genera

    def _degree_decorations(self, vertex_count):
        """Cache degree distributions shared by all incidence matrices."""
        if vertex_count not in self._degree_cache:
            self._degree_cache[vertex_count] = tuple(
                _weak_compositions(self.degree, vertex_count)
            )
        return self._degree_cache[vertex_count]

    @staticmethod
    def _canonical_marking_decorations(marking_count, zero_genera,
                                       zero_degrees, pairs):
        r"""Assign labeled markings modulo interchangeable zero vertices.

        Zero vertices with the same genus, degree, and incidence row can be
        permuted while fixing every infinity vertex.  Generating all
        ``zero_count^marking_count`` assignments and deduplicating only after
        graph construction is the dominant cost of the four-mark genus-three
        family.  Within each interchangeable group we use restricted-growth
        labels: a new vertex may be used only after every earlier vertex in
        that group has appeared.  This selects exactly one representative of
        every orbit under that evident product of symmetric groups.
        """
        marking_count = Integer(marking_count)
        zero_count = len(zero_genera)
        if marking_count < 0:
            raise ValueError("marking_count must be nonnegative")
        if marking_count == 0:
            yield tuple()
            return

        infinity_count = 1 + max(infinity for _, infinity in pairs)
        rows = [[0] * infinity_count for _ in range(zero_count)]
        for zero, infinity in pairs:
            rows[zero][infinity] += 1

        grouped = defaultdict(list)
        for zero in range(zero_count):
            grouped[(
                zero_genera[zero], zero_degrees[zero], tuple(rows[zero])
            )].append(zero)
        groups = tuple(tuple(vertices) for _, vertices in sorted(
            grouped.items(), key=lambda item: item[1][0]
        ))
        vertex_data = {}
        for group_index, vertices in enumerate(groups):
            for position, zero in enumerate(vertices):
                vertex_data[zero] = (group_index, position)

        def extend(marking, assignment, largest_used):
            if marking == marking_count:
                yield tuple(assignment)
                return
            for zero in range(zero_count):
                group, position = vertex_data[zero]
                if position > largest_used[group] + 1:
                    continue
                previous = largest_used[group]
                if position > previous:
                    largest_used[group] = position
                assignment.append(zero)
                for result in extend(
                        marking + 1, assignment, largest_used):
                    yield result
                assignment.pop()
                largest_used[group] = previous

        for result in extend(
                0, [], [-1] * len(groups)):
            yield result

    @staticmethod
    def _contact_budget_available(pairs, zero_types, infinity_genera,
                                  infinity_degrees):
        r"""Necessary condition for any positive contact decoration to exist.

        The contacts at infinity vertex ``v`` sum to
        ``2g(v)-2+val(v)-3d(v)``.  Every edge contributes at least one, and an
        edge meeting a *nonspecial* zero vertex contributes at least two,
        because ``validation_errors`` rejects a nonspecial vertex whose
        contact equals one.  Testing this before ``_contact_decorations``
        discards configurations whose contact budget can never balance,
        instead of building each of their contact assignments first.
        """
        for infinity in range(len(infinity_genera)):
            valence = 0
            nonspecial = 0
            for zero, other in pairs:
                if other != infinity:
                    continue
                valence += 1
                if zero_types[zero] == "nonspecial":
                    nonspecial += 1
            required = (2 * infinity_genera[infinity] - 2 + valence
                        - 3 * infinity_degrees[infinity])
            if required < valence + nonspecial:
                return False
        return True

    @staticmethod
    def _zero_vertex_types(zero_genera, zero_degrees, pairs,
                           marking_vertices):
        r"""Classify all zero vertices without constructing a graph object.

        This is the innermost rejection test in ``_mixed_graphs``.  Calling
        ``BipartiteLocalizationGraph.zero_vertex_type`` once per vertex used
        to rebuild the complete edge- and marking-valence arrays on every
        call.  Computing both arrays once applies the same cases in time
        linear in the size of the decoration.
        """
        zero_count = len(zero_genera)
        edge_valences = [0] * zero_count
        marking_valences = [0] * zero_count
        for zero, _ in pairs:
            edge_valences[zero] += 1
        for zero in marking_vertices:
            marking_valences[zero] += 1

        types = []
        for zero in range(zero_count):
            genus = zero_genera[zero]
            degree = zero_degrees[zero]
            edge_valence = edge_valences[zero]
            marking_valence = marking_valences[zero]
            valence = edge_valence + marking_valence
            if genus != 0 or degree != 0 or valence > 2:
                kind = "stable"
            elif edge_valence == 1 and marking_valence == 0:
                kind = "nonspecial"
            elif edge_valence == 1 and marking_valence == 1:
                kind = "marked"
            elif edge_valence == 2 and marking_valence == 0:
                kind = "nodal"
            else:
                kind = "invalid"
            types.append(kind)
        return tuple(types)

    def _contact_decorations(self, pairs, infinity_genera,
                             infinity_degrees):
        by_infinity = defaultdict(list)
        for edge_index, (_, infinity) in enumerate(pairs):
            by_infinity[infinity].append(edge_index)

        choices = []
        for infinity in range(len(infinity_genera)):
            incident = by_infinity[infinity]
            valence = len(incident)
            required = (2 * infinity_genera[infinity] - 2 + valence
                        - 3 * infinity_degrees[infinity])
            compositions = tuple(_positive_compositions(required, valence))
            if not compositions:
                return
            choices.append((incident, compositions))

        for selected in product(*(entry[1] for entry in choices)):
            contacts = [None] * len(pairs)
            for (incident, _), local_contacts in zip(choices, selected):
                for edge_index, contact in zip(incident, local_contacts):
                    contacts[edge_index] = contact
            yield tuple(contacts)

    def _mixed_graphs(self):
        vertex_bound = self.mixed_vertex_bound()
        if vertex_bound < 2 or self.genus < 1:
            return

        # Each infinity vertex has genus at least one.
        for infinity_count in range(1, self.genus + 1):
            for zero_count in range(1, vertex_bound - infinity_count + 1):
                vertex_count = zero_count + infinity_count
                for edge_count in range(vertex_count - 1,
                                        self.genus - infinity_count
                                        + vertex_count):
                    first_betti = edge_count - vertex_count + 1
                    if first_betti > self.genus:
                        continue
                    for pairs in self._endpoint_multigraphs(
                            zero_count, infinity_count, edge_count):
                        for zero_genera, infinity_genera in \
                                self._genus_decorations(
                                    vertex_count, infinity_count, first_betti):
                            for degrees in self._degree_decorations(
                                    vertex_count):
                                zero_degrees = degrees[:zero_count]
                                infinity_degrees = degrees[zero_count:]
                                if any(2 * genus - 2 - 3 * degree < 0
                                       for genus, degree in zip(
                                           infinity_genera,
                                           infinity_degrees)):
                                    continue
                                for marking_vertices in \
                                        self._canonical_marking_decorations(
                                            self.markings,
                                            zero_genera,
                                            zero_degrees,
                                            pairs,
                                        ):
                                    # Reject forbidden zero unstable types before
                                    # the usually larger contact-order loop.
                                    zero_types = self._zero_vertex_types(
                                        zero_genera, zero_degrees,
                                        pairs, marking_vertices,
                                    )
                                    if any(kind == "invalid"
                                           for kind in zero_types):
                                        continue
                                    if not self._contact_budget_available(
                                            pairs, zero_types,
                                            infinity_genera,
                                            infinity_degrees):
                                        continue
                                    for contacts in self._contact_decorations(
                                            pairs, infinity_genera,
                                            infinity_degrees):
                                        edges = tuple(
                                            (z, i, contacts[j])
                                            for j, (z, i) in enumerate(pairs)
                                        )
                                        graph = BipartiteLocalizationGraph(
                                            zero_genera, zero_degrees,
                                            infinity_genera, infinity_degrees,
                                            marking_vertices, edges,
                                        )
                                        if graph.is_valid():
                                            yield graph

    def graphs(self):
        """Return all isomorphism classes, sorted by canonical key."""
        if self._cache is not None:
            return self._cache
        if self.degree == 0 and 2 * self.genus - 2 + self.markings <= 0:
            self._cache = tuple()
            return self._cache
        representatives = {}
        for graph in self._one_level_graphs():
            representatives.setdefault(graph.canonical_key(), graph)
        for graph in self._mixed_graphs():
            representatives.setdefault(graph.canonical_key(), graph)
        self._cache = tuple(representatives[key]
                            for key in sorted(representatives))
        return self._cache

    def graph_count(self):
        return Integer(len(self.graphs()))

    def report(self):
        lines = [
            "Plane-cubic CJR localization graphs for (g,n,d)=(%s,%s,%s)"
            % (self.genus, self.markings, self.degree),
            "mixed-level vertex bound: %s" % self.mixed_vertex_bound(),
            "isomorphism classes: %s" % self.graph_count(),
        ]
        for index, graph in enumerate(self.graphs(), start=1):
            lines.append("%s. %s" % (index, graph.signature()))
        return "\n".join(lines)


def enumerate_plane_cubic_graphs(genus, markings, degree):
    """Convenience wrapper returning all graph isomorphism classes."""
    return PlaneCubicGraphEnumerator(genus, markings, degree).graphs()


def _command_line(argv):
    if len(argv) != 4:
        print("usage: sage cjr_bipartite_graphs.sage GENUS MARKINGS DEGREE")
        return Integer(2)
    enumerator = PlaneCubicGraphEnumerator(
        Integer(argv[1]), Integer(argv[2]), Integer(argv[3])
    )
    print(enumerator.report())
    return Integer(0)


if __name__ == "__main__":
    import os
    import sys
    if os.path.basename(sys.argv[0]) in (
            "cjr_bipartite_graphs.sage",
            "cjr_bipartite_graphs.sage.py"):
        sys.exit(int(_command_line(sys.argv)))
