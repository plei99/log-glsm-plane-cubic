r"""Exact truncated Givental--Teleman reconstruction utilities.

The implementation uses the unit-preserving stable-graph action

``leg = R^(-1)(psi)``,
``edge = (eta^(-1)-R^(-1)(psi')eta^(-1)R^(-1)(psi'')^T)/(psi'+psi'')``,
``T(z) = z(1-R^(-1)(z)1)``.

The underlying topological field theory is supplied in an idempotent basis,
so its metric is diagonal and every stable vertex has a single colour.  All
integrals on the boundary factors are pure psi intersections and are handled
by the DVV backend; admcycles is used only to enumerate stable graph types.
"""

import itertools

if "PsiIntersectionBackend" not in globals():
    load("hodge_integrals.sage")


class TruncatedRMatrix(SageObject):
    """A finite matrix series ``R(z)=I+R_1 z+...+R_N z^N``."""

    def __init__(self, coefficient_field, coefficients):
        self.field = coefficient_field
        raw = tuple(coefficients)
        if not raw:
            raise ValueError("an R-matrix needs its constant coefficient")
        self.rank = ZZ(raw[0].nrows())
        if self.rank < 1 or any(
                matrix.nrows() != self.rank or matrix.ncols() != self.rank
                for matrix in raw):
            raise ValueError("all R-matrix coefficients must be square")
        self.coefficients = tuple(
            Matrix(self.field, matrix) for matrix in raw
        )
        if self.coefficients[0] != identity_matrix(self.field, self.rank):
            raise ValueError("R_0 must be the identity")
        self.max_power = ZZ(len(self.coefficients) - 1)
        self._inverse = None

    def coefficient(self, power):
        power = ZZ(power)
        if power < 0 or power > self.max_power:
            return zero_matrix(self.field, self.rank, self.rank)
        return self.coefficients[int(power)]

    def inverse_coefficients(self):
        """Return the truncated coefficients of ``R(z)^(-1)``."""
        if self._inverse is None:
            inverse = [identity_matrix(self.field, self.rank)]
            for degree in range(1, self.max_power + 1):
                inverse.append(-sum(
                    self.coefficients[index] * inverse[degree - index]
                    for index in range(1, degree + 1)
                ))
            self._inverse = tuple(inverse)
        return self._inverse

    def symplectic_defects(self, metric):
        r"""Coefficients of ``R(z) eta^-1 R(-z)^T-eta^-1``."""
        metric = tuple(self.field(value) for value in metric)
        if len(metric) != self.rank or any(not value for value in metric):
            raise ValueError("the idempotent metric must be nonzero")
        inverse_metric = diagonal_matrix(
            self.field, tuple(value ** (-1) for value in metric)
        )
        defects = []
        # With coefficients known only through z^N, the symplectic identity
        # is determined only through z^N. Higher apparent defects merely use
        # missing R_(N+1),...,R_(2N) coefficients.
        for total in range(self.max_power + 1):
            coefficient = zero_matrix(self.field, self.rank, self.rank)
            for left in range(total + 1):
                right = total - left
                if left <= self.max_power and right <= self.max_power:
                    coefficient += (
                        self.coefficient(left) * inverse_metric
                        * ((-1) ** right * self.coefficient(right).transpose())
                    )
            if total == 0:
                coefficient -= inverse_metric
            defects.append(coefficient)
        return tuple(defects)


class SemisimpleTFTData(SageObject):
    r"""Diagonal metric, unit, and truncated R-matrix at one semisimple point."""

    def __init__(self, coefficient_field, metric, r_matrix, unit=None,
                 label=""):
        if not isinstance(r_matrix, TruncatedRMatrix):
            raise TypeError("r_matrix must be a TruncatedRMatrix")
        if r_matrix.field != coefficient_field:
            raise ValueError("the metric and R-matrix need the same field")
        self.field = coefficient_field
        self.metric = tuple(self.field(value) for value in metric)
        self.rank = ZZ(len(self.metric))
        if self.rank != r_matrix.rank or any(not value for value in self.metric):
            raise ValueError("the idempotent metric has the wrong rank")
        self.r_matrix = r_matrix
        self.unit = tuple(
            self.field.one() for _ in range(self.rank)
        ) if unit is None else tuple(self.field(value) for value in unit)
        if len(self.unit) != self.rank:
            raise ValueError("the unit has the wrong rank")
        self.label = str(label)

    def vertex_weight(self, colour, genus):
        return self.metric[int(colour)] ** (1 - ZZ(genus))


class GiventalTelemanCorrelator(SageObject):
    """Integrate the truncated R-action over every stable graph."""

    def __init__(self, data, psi_backend=None):
        if not isinstance(data, SemisimpleTFTData):
            raise TypeError("data must be SemisimpleTFTData")
        self.data = data
        self.field = data.field
        self.psi = psi_backend or PsiIntersectionBackend()
        self.admcycles = _load_vendored_admcycles()
        self._graph_cache = {}
        self._edge_cache = {}
        self._translation_cache = {}
        self._inverse_r = self.data.r_matrix.inverse_coefficients()

    def stable_graphs(self, genus, markings):
        genus = ZZ(genus)
        markings = ZZ(markings)
        if 2 * genus - 2 + markings <= 0:
            raise ValueError("Givental reconstruction needs a stable (g,n)")
        key = genus, markings
        if key not in self._graph_cache:
            maximum_edges = 3 * genus - 3 + markings
            self._graph_cache[key] = tuple(
                graph
                for edge_count in range(maximum_edges + 1)
                for graph in self.admcycles.list_strata(
                    int(genus), int(markings), int(edge_count)
                )
            )
        return self._graph_cache[key]

    def _leg_polynomial(self, colour, vector, descendant_power):
        vector = tuple(self.field(entry) for entry in vector)
        if len(vector) != self.data.rank:
            raise ValueError("an insertion vector has the wrong rank")
        answer = {}
        for power, coefficient in enumerate(self._inverse_r):
            value = sum(
                coefficient[int(colour), column] * vector[column]
                for column in range(self.data.rank)
            )
            exponent = ZZ(power) + ZZ(descendant_power)
            if value:
                answer[exponent] = answer.get(
                    exponent, self.field.zero()
                ) + value
        return answer

    def _edge_polynomial(self, left_colour, right_colour):
        key = ZZ(left_colour), ZZ(right_colour)
        if key in self._edge_cache:
            return self._edge_cache[key]
        rank = int(self.data.rank)
        maximum = int(self.data.r_matrix.max_power)
        inverse_metric = tuple(value ** (-1) for value in self.data.metric)

        def numerator(left_power, right_power):
            if left_power > maximum or right_power > maximum:
                return self.field.zero()
            value = -sum(
                self._inverse_r[left_power][int(left_colour), middle]
                * inverse_metric[middle]
                * self._inverse_r[right_power][int(right_colour), middle]
                for middle in range(rank)
            )
            if left_power == 0 and right_power == 0 \
                    and left_colour == right_colour:
                value += inverse_metric[int(left_colour)]
            return self.field(value)

        # Divide the numerator by x+y coefficient-by-coefficient.  Terms of
        # total degree d need numerator coefficients of total degree d+1.
        answer = {}
        for total in range(maximum):
            previous = self.field.zero()
            for left_power in range(total + 1):
                right_power = total - left_power
                value = numerator(left_power, right_power + 1) - previous
                if value:
                    answer[(ZZ(left_power), ZZ(right_power))] = value
                previous = value
            boundary = numerator(total + 1, 0)
            if boundary != previous:
                raise ValueError(
                    "R-matrix edge numerator is not divisible by psi'+psi''; "
                    "increase the truncation or correct the symplectic data"
                )
        self._edge_cache[key] = answer
        return answer

    def _translation_polynomial(self, colour):
        colour = ZZ(colour)
        if colour in self._translation_cache:
            return self._translation_cache[colour]
        answer = {}
        for inverse_power in range(1, len(self._inverse_r)):
            value = -sum(
                self._inverse_r[inverse_power][int(colour), column]
                * self.data.unit[column]
                for column in range(self.data.rank)
            )
            if value:
                answer[ZZ(inverse_power + 1)] = self.field(value)
        self._translation_cache[colour] = answer
        return answer

    @staticmethod
    def _multiply_slot_polynomial(terms, slot, polynomial):
        answer = {}
        for exponents, coefficient in terms.items():
            for power, scalar in polynomial.items():
                changed = list(exponents)
                changed[slot] += power
                key = tuple(changed)
                answer[key] = answer.get(key, 0) + coefficient * scalar
        return {key: value for key, value in answer.items() if value}

    def _vertex_integral(self, genus, base_exponents, colour):
        genus = ZZ(genus)
        base_exponents = tuple(ZZ(value) for value in base_exponents)
        translation = self._translation_polynomial(colour)
        base_deficit = (
            3 * genus - 3 + len(base_exponents) - sum(base_exponents)
        )
        if base_deficit < 0:
            return self.field.zero()
        answer = self.field(self.psi.integral(genus, base_exponents))
        if not translation:
            return answer

        translation_terms = tuple(sorted(translation.items()))
        for count in range(1, int(base_deficit) + 1):
            for choices in itertools.product(translation_terms, repeat=count):
                powers = tuple(choice[0] for choice in choices)
                if sum(power - 1 for power in powers) != base_deficit:
                    continue
                coefficient = prod(choice[1] for choice in choices) \
                    / factorial(count)
                answer += coefficient * self.psi.integral(
                    genus, base_exponents + powers
                )
        return self.field(answer)

    def _graph_contribution(self, graph, colours, insertions):
        vertex_count = graph.num_verts()
        legs = tuple(tuple(graph.legs(vertex))
                     for vertex in range(vertex_count))
        positions = {
            leg: (vertex, offset)
            for vertex, vertex_legs in enumerate(legs)
            for offset, leg in enumerate(vertex_legs)
        }
        exponent_shape = tuple(
            tuple(ZZ.zero() for _ in vertex_legs) for vertex_legs in legs
        )
        terms = {exponent_shape: self.field.one()}

        # Ordinary markings are labelled 1,...,n by admcycles.
        for marking, (vector, descendant) in enumerate(insertions, start=1):
            vertex, slot = positions[marking]
            polynomial = self._leg_polynomial(
                colours[vertex], vector, descendant
            )
            updated = {}
            for exponents, coefficient in terms.items():
                for power, scalar in polynomial.items():
                    changed = [list(local) for local in exponents]
                    changed[vertex][slot] += power
                    key = tuple(tuple(local) for local in changed)
                    updated[key] = updated.get(
                        key, self.field.zero()
                    ) + coefficient * scalar
            terms = {key: value for key, value in updated.items() if value}
            if not terms:
                return self.field.zero()

        for left_leg, right_leg in graph.edges():
            left_vertex, left_slot = positions[left_leg]
            right_vertex, right_slot = positions[right_leg]
            polynomial = self._edge_polynomial(
                colours[left_vertex], colours[right_vertex]
            )
            updated = {}
            for exponents, coefficient in terms.items():
                for (left_power, right_power), scalar in polynomial.items():
                    changed = [list(local) for local in exponents]
                    changed[left_vertex][left_slot] += left_power
                    changed[right_vertex][right_slot] += right_power
                    key = tuple(tuple(local) for local in changed)
                    updated[key] = updated.get(
                        key, self.field.zero()
                    ) + coefficient * scalar
            terms = {key: value for key, value in updated.items() if value}
            if not terms:
                return self.field.zero()

        answer = self.field.zero()
        vertex_prefactor = prod(
            self.data.vertex_weight(colours[vertex], graph.genera()[vertex])
            for vertex in range(vertex_count)
        )
        for exponents, coefficient in terms.items():
            local_value = prod(
                self._vertex_integral(
                    graph.genera()[vertex], exponents[vertex], colours[vertex]
                )
                for vertex in range(vertex_count)
            )
            answer += coefficient * vertex_prefactor * local_value
        return self.field(answer / graph.automorphism_number())

    def evaluate(self, genus, insertions):
        r"""Evaluate descendants supplied as ``(idempotent_vector, psi)``."""
        genus = ZZ(genus)
        insertions = tuple(
            (tuple(self.field(value) for value in vector), ZZ(descendant))
            for vector, descendant in insertions
        )
        if any(descendant < 0 for _, descendant in insertions):
            raise ValueError("descendant powers must be nonnegative")
        graphs = self.stable_graphs(genus, len(insertions))
        answer = self.field.zero()
        for graph in graphs:
            for colours in itertools.product(
                    range(int(self.data.rank)), repeat=graph.num_verts()):
                answer += self._graph_contribution(graph, colours, insertions)
        return self.field(answer)


def exponential_series_coefficients(field, logarithm, maximum):
    """Exponentiate a scalar series represented by a coefficient dictionary."""
    maximum = ZZ(maximum)
    answer = [field.one()]
    for degree in range(1, maximum + 1):
        answer.append(field(sum(
            index * field(logarithm.get(ZZ(index), 0))
            * answer[degree - index]
            for index in range(1, degree + 1)
        ) / degree))
    return tuple(answer)
