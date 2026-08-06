r"""
Dynamic-programming infrastructure for reconstructing log-GLSM infinity data.

For a chosen probe correlator, write the CJR localization equation as

    known_GW = twisted_zero_level + sum_m coefficient_m * product(V in m) V.

An equation provider is responsible for enumerating the localization graphs
and returning this finite expression.  ``InfinityVertexDP`` recursively
substitutes lower effective vertices, isolates the requested diagonal term,
and memoizes the answer.  If several vertices of the same rank mix linearly,
``solve_block`` solves their finite matrix simultaneously.

The engine is exact and agnostic about how twisted invariants are computed.
For a plane cubic, an equation provider should use O(3)-twisted P^2 theory for
zero-level vertices and Bloch--Okounkov stationary series for ``known_GW``.
The function ``genus_two_bo_demo`` instantiates the already-known genus-one
and genus-two degree-zero equations.

Important: infinite computing power handles enumeration, but not a singular
diagonal.  A proposed probe must still have a nonzero target coefficient, or
a full-rank same-level block of probes must be supplied.
"""

load("log_glsm_infinity_vertices.sage")


class EffectiveVertex:
    """Hashable discrete data for a punctured infinity invariant."""

    def __init__(self, genus, ambient_degree, contacts, psi_min=0,
                 insertions=(), label=""):
        self.genus = ZZ(genus)
        self.ambient_degree = ZZ(ambient_degree)
        self.contacts = tuple(sorted(ZZ(c) for c in contacts))
        self.psi_min = ZZ(psi_min)
        self.insertions = tuple(ZZ(insertion) for insertion in insertions)
        self.label = str(label)

        if self.genus < 0:
            raise ValueError("genus must be nonnegative")
        if self.ambient_degree < 0:
            raise ValueError("ambient_degree must be nonnegative")
        if self.psi_min < 0:
            raise ValueError("psi_min must be nonnegative")
        if any(c >= 0 for c in self.contacts):
            raise ValueError("infinity contacts must be strictly negative")

    def signature(self):
        return (
            self.genus,
            self.ambient_degree,
            self.contacts,
            self.psi_min,
            self.insertions,
            self.label,
        )

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, EffectiveVertex) and self.signature() == other.signature()

    def __repr__(self):
        suffix = ", label=%r" % self.label if self.label else ""
        insertion_suffix = ", insertions=%r" % (self.insertions,) \
            if self.insertions else ""
        return (
            "EffectiveVertex(g=%s, D=%s, contacts=%s, psi_min=%s%s%s)"
            % (self.genus, self.ambient_degree, self.contacts, self.psi_min,
               insertion_suffix, suffix)
        )

    def to_record(self):
        """Return a versioned-checkpoint-friendly record."""
        return {
            "genus": int(self.genus),
            "ambient_degree": int(self.ambient_degree),
            "contacts": [int(value) for value in self.contacts],
            "psi_min": int(self.psi_min),
            "insertions": [int(value) for value in self.insertions],
            "label": self.label,
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["genus"], record["ambient_degree"], record["contacts"],
            psi_min=record.get("psi_min", 0),
            insertions=tuple(record.get("insertions", ())),
            label=record.get("label", ""),
        )

    def balance_defect(self, hypersurface_degree=3):
        r"""Return the defect in the hypersurface infinity balance equation.

        A balanced vertex satisfies

            degree(E)*D - (2g-2) = sum_i(c_i+1).
        """
        hypersurface_degree = ZZ(hypersurface_degree)
        return (
            hypersurface_degree * self.ambient_degree
            - (2 * self.genus - 2)
            - sum(c + 1 for c in self.contacts)
        )

    def is_balanced(self, hypersurface_degree=3):
        return self.balance_defect(hypersurface_degree) == 0

    def order_key(self):
        r"""Default well-founded candidate order for triangular recursion.

        The graph provider must verify triangularity for its chosen probes;
        this order alone is not a mathematical proof that a probe works.
        """
        contact_excess = sum(-c - 1 for c in self.contacts)
        return (
            self.genus,
            self.ambient_degree,
            contact_excess,
            len(self.contacts),
            self.psi_min,
            self.insertions,
            self.contacts,
            self.label,
        )


class LocalizationEquation:
    r"""A finite equation ``known = twisted + sum coefficient*monomial``."""

    def __init__(self, target, known_gw, twisted_zero_level=0,
                 terms=(), probe_label="", coefficient_field=QQ):
        if not isinstance(target, EffectiveVertex):
            raise TypeError("target must be an EffectiveVertex")
        self.coefficient_field = coefficient_field
        self.target = target
        self.known_gw = coefficient_field(known_gw)
        self.twisted_zero_level = coefficient_field(twisted_zero_level)
        self.probe_label = str(probe_label)
        normalized = []
        for coefficient, factors in terms:
            factors = tuple(factors)
            if any(not isinstance(v, EffectiveVertex) for v in factors):
                raise TypeError("every monomial factor must be an EffectiveVertex")
            normalized.append((coefficient_field(coefficient), factors))
        self.terms = tuple(normalized)


class NonTriangularLocalizationError(RuntimeError):
    pass


class SingularProbeError(RuntimeError):
    pass


class InfinityVertexDP:
    r"""Memoized scalar/block solver for triangular localization equations."""

    def __init__(self, equation_provider, initial_values=None, order_key=None,
                 coefficient_field=QQ):
        self.equation_provider = equation_provider
        self.order_key = order_key or (lambda vertex: vertex.order_key())
        self.coefficient_field = coefficient_field
        self.values = {
            vertex: coefficient_field(value)
            for vertex, value in dict(initial_values or {}).items()
        }
        self.equations = {}

    def equation(self, vertex):
        if vertex not in self.equations:
            equation = self.equation_provider(vertex)
            if equation.target != vertex:
                raise ValueError("the equation provider returned the wrong target")
            self.equations[vertex] = equation
        return self.equations[vertex]

    def _lower_value(self, factor, target, block):
        if factor in self.values:
            return self.values[factor]
        if factor in block:
            raise AssertionError("block factors are handled before lower substitution")
        if not self.order_key(factor) < self.order_key(target):
            raise NonTriangularLocalizationError(
                "probe for %r depends on non-lower vertex %r; include same-rank "
                "vertices in solve_block or choose another probe" % (target, factor)
            )
        return self.solve(factor)

    def solve(self, vertex):
        if vertex in self.values:
            return self.values[vertex]
        self.solve_block((vertex,))
        return self.values[vertex]

    def solve_block(self, vertices):
        r"""Solve a finite linearly-coupled block after lower substitution."""
        vertices = tuple(vertices)
        if not vertices:
            return {}
        if len(set(vertices)) != len(vertices):
            raise ValueError("a solve block cannot contain duplicate vertices")

        unsolved = tuple(vertex for vertex in vertices if vertex not in self.values)
        if not unsolved:
            return {vertex: self.values[vertex] for vertex in vertices}
        block = set(unsolved)
        index = {vertex: i for i, vertex in enumerate(unsolved)}
        rows = []
        right_hand_sides = []

        for target in unsolved:
            equation = self.equation(target)
            if equation.coefficient_field is not self.coefficient_field:
                raise TypeError("equation and solver coefficient fields must agree")
            row = [self.coefficient_field.zero()] * len(unsolved)
            off_diagonal = self.coefficient_field.zero()

            for coefficient, factors in equation.terms:
                block_occurrences = [factor for factor in factors if factor in block]
                if len(block_occurrences) > 1:
                    raise NonTriangularLocalizationError(
                        "probe %r contains a nonlinear monomial in unsolved block "
                        "vertices: %r" % (equation.probe_label, factors)
                    )

                multiplier = coefficient
                skipped_block_factor = False
                for factor in factors:
                    if factor in block and not skipped_block_factor:
                        skipped_block_factor = True
                        continue
                    multiplier *= self._lower_value(factor, target, block)

                if block_occurrences:
                    row[index[block_occurrences[0]]] += multiplier
                else:
                    off_diagonal += multiplier

            rows.append(row)
            right_hand_sides.append(
                equation.known_gw
                - equation.twisted_zero_level
                - off_diagonal
            )

        matrix = Matrix(self.coefficient_field, rows)
        rhs = vector(self.coefficient_field, right_hand_sides)
        if matrix.rank() < len(unsolved):
            raise SingularProbeError(
                "the chosen probes have singular diagonal matrix %s" % matrix
            )
        solution = matrix.solve_right(rhs)
        for vertex, value in zip(unsolved, solution):
            self.values[vertex] = value
        return {vertex: self.values[vertex] for vertex in vertices}

    def solve_up_to(self, vertices):
        """Solve a caller-supplied finite truncation in the configured order."""
        for vertex in sorted(set(vertices), key=self.order_key):
            self.solve(vertex)
        return {vertex: self.values[vertex] for vertex in vertices}

    def dependency_graph(self, roots):
        """Return the finite dependency DAG exposed by the equation provider."""
        graph = {}
        pending = list(roots)
        while pending:
            vertex = pending.pop()
            if vertex in graph or vertex in self.values:
                continue
            equation = self.equation(vertex)
            dependencies = set(
                factor for coefficient, factors in equation.terms
                for factor in factors if factor != vertex
            )
            graph[vertex] = tuple(sorted(dependencies, key=self.order_key))
            pending.extend(dependencies)
        return graph

    def checkpoint(self):
        """Return solved values in a deterministic versioned exact format."""
        entries = sorted(
            self.values.items(), key=lambda item: item[0].order_key()
        )
        return {
            "format": "log-glsm-infinity-dp",
            "version": int(1),
            "coefficient_field": str(self.coefficient_field),
            "values": [
                {"vertex": vertex.to_record(), "value": str(value)}
                for vertex, value in entries
            ],
        }

    def restore_checkpoint(self, checkpoint):
        """Merge a checkpoint after validating its version and duplicates."""
        if checkpoint.get("format") != "log-glsm-infinity-dp" \
                or checkpoint.get("version") != 1:
            raise ValueError("unsupported infinity-DP checkpoint format")
        restored = {}
        for entry in checkpoint.get("values", ()):
            vertex = EffectiveVertex.from_record(entry["vertex"])
            value = self.coefficient_field(entry["value"])
            if vertex in restored and restored[vertex] != value:
                raise ValueError("checkpoint contains conflicting duplicate values")
            if vertex in self.values and self.values[vertex] != value:
                raise ValueError("checkpoint conflicts with an existing solved value")
            restored[vertex] = value
        self.values.update(restored)
        return restored


def genus_two_bo_demo(max_degree=4):
    r"""Solve the two primitive degree-zero cumulants with the DP engine.

    This is the scalar base case already reconstructed by
    ``log_glsm_infinity_vertices.sage``.  The genus-two key is explicitly
    labeled as the aggregate contact-profile combination; it is not claimed
    to equal either individual ``(-3)`` or ``(-2,-2)`` vertex integral.
    """
    max_degree = _validate_max_degree(max_degree)
    genus_one_polynomial, genus_one = connected_stationary_qseries([0], max_degree)
    genus_two_polynomial, genus_two = connected_stationary_qseries([2], max_degree)

    b2 = EffectiveVertex(1, 0, (-1,), label="basic")
    b4 = EffectiveVertex(
        2, 0, (-3,), label="aggregate_genus_two_profile_combination"
    )

    def provider(vertex):
        if vertex == b2:
            return LocalizationEquation(
                b2,
                known_gw=genus_one[0],
                terms=((1, (b2,)),),
                probe_label="<tau_0(pt)>_(1,1,0)",
            )
        if vertex == b4:
            return LocalizationEquation(
                b4,
                known_gw=genus_two[0],
                terms=((1, (b4,)), (QQ(1) / 2, (b2, b2))),
                probe_label="<tau_2(pt)>_(2,1,0)",
            )
        raise KeyError("no demonstration probe for %r" % vertex)

    solver = InfinityVertexDP(provider)
    dependency_graph = solver.dependency_graph((b4,))
    values = solver.solve_up_to((b2, b4))
    return {
        "genus_one_polynomial": genus_one_polynomial,
        "genus_two_polynomial": genus_two_polynomial,
        "genus_one_vertex": b2,
        "genus_two_vertex": b4,
        "values": values,
        "dependency_graph": dependency_graph,
        "solver": solver,
    }
