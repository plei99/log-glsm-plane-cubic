r"""Compile enumerated CJR graphs into effective-vertex polynomials.

The compiler implements every universal edge, unstable-zero, diagonal, and
infinity-descendant factor.  Stable zero vertices are delegated to the full
base-torus localization backend.  A caller may still inject a deliberately
partial backend; in that case the structured unsupported-geometry diagnostic
is retained and no missing vertex is silently assigned zero.
"""

from itertools import product

load("cjr_bipartite_graphs.sage")
load("log_glsm_infinity_dp.sage")
load("cjr_graph_factors.sage")
load("o3_fixed_locus_graphs.sage")


class GraphContribution(SageObject):
    """Exact polynomial or explicit unsupported result for one graph."""

    def __init__(self, graph, polynomial=None, unsupported=(), provenance=()):
        self.graph = graph
        self.polynomial = polynomial
        self.unsupported = tuple(str(item) for item in unsupported)
        self.provenance = tuple(str(item) for item in provenance)

    @property
    def supported(self):
        return self.polynomial is not None and not self.unsupported

    def _repr_(self):
        if self.supported:
            return "GraphContribution(%s)" % self.polynomial
        return "UnsupportedGraphContribution(%s)" % (self.unsupported,)

    def report(self):
        terms = []
        if self.polynomial is not None:
            terms = [
                {
                    "coefficient": str(coefficient),
                    "factors": [factor.to_record() for factor in factors],
                }
                for factors, coefficient in sorted(
                    self.polynomial.terms.items(),
                    key=lambda item: tuple(factor.order_key() for factor in item[0]),
                )
            ]
        return {
            "graph": str(self.graph.signature()),
            "canonical_key": str(self.graph.canonical_key()),
            "supported": self.supported,
            "automorphism_order": int(self.graph.automorphism_order()),
            "terms": terms,
            "unsupported": list(self.unsupported),
            "provenance": list(self.provenance),
        }


class ProbeCompilation(SageObject):
    """Sum of all supported graph terms plus any exact diagnostics."""

    def __init__(self, probe, contributions, coefficient_field):
        self.probe = probe
        self.contributions = tuple(contributions)
        self.coefficient_field = coefficient_field
        self.polynomial = EffectivePolynomial(coefficient_field)
        for contribution in self.contributions:
            if contribution.supported:
                self.polynomial = self.polynomial + contribution.polynomial

    @property
    def unsupported(self):
        return tuple(
            (contribution.graph, contribution.unsupported)
            for contribution in self.contributions if not contribution.supported
        )

    @property
    def is_complete(self):
        return not self.unsupported

    def require_complete(self):
        if self.unsupported:
            raise UnsupportedGeometryError(
                "probe compilation has %s unsupported graph(s): %s"
                % (len(self.unsupported), self.unsupported)
            )
        return self

    def report(self):
        return {
            "probe": str(self.probe),
            "graph_count": len(self.contributions),
            "supported_graphs": sum(item.supported for item in self.contributions),
            "unsupported_graphs": sum(not item.supported for item in self.contributions),
            "collected_terms": len(self.polynomial.terms),
            "graphs": [item.report() for item in self.contributions],
        }


class PlaneCubicGraphContributionCompiler(SageObject):
    """Exact graph-to-polynomial compiler in the compact plane-cubic sector."""

    def __init__(self, rings=None, twisted_backend=None):
        self.rings = rings or PlaneCubicCoefficientRing(16)
        self.factors = PlaneCubicGraphFactors(self.rings)
        self.twisted_backend = twisted_backend or FullTwistedZeroVertexBackend(
            self.rings
        )

    def _ordinary_insertion(self, probe, marking):
        return TwistedInsertion.from_elliptic(probe.insertions[marking])

    def _class_for_insertion(self, insertion):
        coefficients = [0, 0, 0]
        coefficients[insertion.h_power] = insertion.scale
        return P2Class(self.rings, coefficients)

    def _incident_edges(self, graph, zero):
        return tuple(
            (edge_index, infinity, contact)
            for edge_index, (z, infinity, contact) in enumerate(graph.edges)
            if z == zero
        )

    def _stable_zero_options(self, probe, graph, zero):
        r"""Return ``(coefficient, {edge_index: infinity_H_power})`` options."""
        ordinary = tuple(
            self._ordinary_insertion(probe, marking)
            for marking in graph.markings_at_zero_vertex(zero)
        )
        edges = self._incident_edges(graph, zero)
        local_options = []
        for edge_index, _, contact in edges:
            remaining_dimension = PlaneCubicDimension.twisted_zero_virtual_dimension(
                graph.zero_genera[zero], graph.zero_degrees[zero],
                len(ordinary) + len(edges),
            ) - sum(item.codimension for item in ordinary)
            maximum_psi = max(ZZ(0), remaining_dimension)
            options_for_edge = []
            edge_class = self.factors.edge(contact)
            for psi_power in range(maximum_psi + 1):
                flag_class = self.factors.stable_flag_descendant(contact, psi_power)
                for diagonal_power in range(3):
                    dual_power = p2_dual_power(diagonal_power)
                    for flag_h_power, flag_coefficient in enumerate(
                            flag_class.coefficients):
                        zero_h_power = diagonal_power + flag_h_power
                        if zero_h_power > 2 or not flag_coefficient:
                            continue
                        for edge_h_power, edge_coefficient in enumerate(
                                edge_class.coefficients):
                            infinity_h_power = dual_power + edge_h_power
                            if infinity_h_power > 2 or not edge_coefficient:
                                continue
                            options_for_edge.append((
                                flag_coefficient * edge_coefficient,
                                TwistedInsertion(
                                    zero_h_power, psi_power, 1,
                                    label="edge_%s" % edge_index,
                                ),
                                edge_index, infinity_h_power,
                            ))
            local_options.append(tuple(options_for_edge))

        if not edges:
            local_products = (tuple(),)
        else:
            local_products = product(*local_options)
        answer = []
        for selected in local_products:
            coefficient = self.rings.full(1)
            zero_insertions = list(ordinary)
            infinity_insertions = {}
            for local_coefficient, zero_insertion, edge_index, infinity_power in selected:
                coefficient *= local_coefficient
                zero_insertions.append(zero_insertion)
                infinity_insertions[edge_index] = infinity_power
            request = TwistedZeroVertexRequest(
                graph.zero_genera[zero], graph.zero_degrees[zero], zero_insertions,
                label="CJR stable zero vertex",
            )
            value = self.twisted_backend.evaluate(request)
            if value:
                answer.append((coefficient * value, infinity_insertions))
        return tuple(answer)

    def _nonspecial_zero_options(self, graph, zero):
        edge_index, _, contact = self._incident_edges(graph, zero)[0]
        local_class = self.factors.nonspecial_edge_pair(contact)
        return tuple(
            (coefficient, {edge_index: h_power})
            for h_power, coefficient in enumerate(local_class.coefficients)
            if coefficient
        )

    def _marked_zero_options(self, probe, graph, zero):
        edge_index, _, contact = self._incident_edges(graph, zero)[0]
        marking = graph.markings_at_zero_vertex(zero)[0]
        ordinary = self._ordinary_insertion(probe, marking)
        ordinary_class = self._class_for_insertion(ordinary)
        local_class = self.factors.edge(contact) \
            * self.factors.zero_marked(contact, ordinary.psi_power) \
            * ordinary_class
        return tuple(
            (coefficient, {edge_index: h_power})
            for h_power, coefficient in enumerate(local_class.coefficients)
            if coefficient
        )

    def _nodal_zero_options(self, graph, zero):
        left, right = self._incident_edges(graph, zero)
        left_index, _, left_contact = left
        right_index, _, right_contact = right
        common_class = self.factors.edge(left_contact) \
            * self.factors.edge(right_contact) \
            * self.factors.zero_nodal(left_contact, right_contact)
        answer = []
        H = P2Class.H(self.rings)
        powers = (P2Class.one(self.rings), H, H * H)
        for left_basis in range(3):
            for right_basis in range(3):
                coefficient = (
                    common_class * powers[left_basis] * powers[right_basis]
                ).integrate()
                if coefficient:
                    answer.append((coefficient, {
                        left_index: p2_dual_power(left_basis),
                        right_index: p2_dual_power(right_basis),
                    }))
        return tuple(answer)

    def _zero_options(self, probe, graph, zero):
        vertex_type = graph.zero_vertex_type(zero)
        if vertex_type == "stable":
            return self._stable_zero_options(probe, graph, zero)
        if vertex_type == "nonspecial":
            return self._nonspecial_zero_options(graph, zero)
        if vertex_type == "marked":
            return self._marked_zero_options(probe, graph, zero)
        if vertex_type == "nodal":
            return self._nodal_zero_options(graph, zero)
        raise ValueError("invalid zero vertex cannot be compiled")

    def _effective_vertex(self, graph, infinity, edge_insertions, psi_power):
        paired = sorted(
            (-contact, edge_insertions[edge_index])
            for edge_index, (_, vertex, contact) in enumerate(graph.edges)
            if vertex == infinity
        )
        contacts = tuple(contact for contact, _ in paired)
        insertions = tuple(h_power for _, h_power in paired)
        vertex = EffectiveVertex(
            graph.infinity_genera[infinity],
            graph.infinity_degrees[infinity],
            contacts,
            psi_min=psi_power,
            insertions=insertions,
        )
        if not vertex.is_balanced():
            raise AssertionError("enumerated infinity vertex lost its balance")
        return vertex

    def _infinity_options(self, graph, infinity, edge_insertions):
        insertion_codimension = sum(
            edge_insertions[edge_index]
            for edge_index, (_, vertex, _) in enumerate(graph.edges)
            if vertex == infinity
        )
        valence = graph.infinity_edge_valences()[infinity]
        bound = PlaneCubicDimension.psi_min_power_bound(
            graph.infinity_genera[infinity], valence, insertion_codimension
        )
        if bound < 0:
            return tuple()
        return tuple(
            (
                self.factors.infinity_descendant_coefficient(psi_power),
                self._effective_vertex(
                    graph, infinity, edge_insertions, psi_power
                ),
            )
            for psi_power in range(bound + 1)
        )

    def compile_graph(self, probe, graph):
        if not isinstance(probe, ProbeSpec):
            raise TypeError("probe must be a ProbeSpec")
        if not graph.is_valid():
            raise ValueError("only valid localization graphs can be compiled")
        if (graph.total_genus(), len(graph.marking_vertices), graph.total_degree()) != (
                probe.genus, probe.marking_count, probe.ambient_degree):
            raise ValueError("graph discrete data do not match the probe")

        provenance = ["graph automorphism factor 1/%s" % graph.automorphism_order()]
        try:
            zero_option_lists = tuple(
                self._zero_options(probe, graph, zero)
                for zero in range(graph.zero_count)
            )
        except UnsupportedGeometryError as error:
            return GraphContribution(graph, unsupported=(error,), provenance=provenance)

        zero_products = product(*zero_option_lists) if zero_option_lists else (tuple(),)
        polynomial = EffectivePolynomial(self.rings.full_field)
        for selected_zero in zero_products:
            coefficient = self.rings.full(graph.localization_weight())
            edge_insertions = {}
            collision = False
            for local_coefficient, local_insertions in selected_zero:
                coefficient *= local_coefficient
                for edge_index, h_power in local_insertions.items():
                    if edge_index in edge_insertions:
                        collision = True
                    edge_insertions[edge_index] = h_power
            if collision or len(edge_insertions) != graph.edge_count:
                raise AssertionError("every localization edge needs one infinity insertion")

            infinity_option_lists = tuple(
                self._infinity_options(graph, infinity, edge_insertions)
                for infinity in range(graph.infinity_count)
            )
            infinity_products = product(*infinity_option_lists) \
                if infinity_option_lists else (tuple(),)
            for selected_infinity in infinity_products:
                total_coefficient = coefficient
                effective_vertices = []
                for infinity_coefficient, vertex in selected_infinity:
                    total_coefficient *= infinity_coefficient
                    effective_vertices.append(vertex)
                polynomial.add_term(total_coefficient, effective_vertices)

        return GraphContribution(
            graph, polynomial=polynomial, provenance=provenance
        )

    def compile_probe(self, probe):
        if not isinstance(probe, ProbeSpec):
            raise TypeError("probe must be a ProbeSpec")
        graphs = PlaneCubicGraphEnumerator(
            probe.genus, probe.marking_count, probe.ambient_degree
        ).graphs()
        return ProbeCompilation(
            probe,
            tuple(self.compile_graph(probe, graph) for graph in graphs),
            self.rings.full_field,
        )

    def laurent_coefficient(self, value, power=0, precision=None):
        series = self.rings.to_laurent(value, precision)
        return series[ZZ(power)]

    def extracted_polynomial(self, compilation, power=0, precision=None):
        if not isinstance(compilation, ProbeCompilation):
            raise TypeError("expected a ProbeCompilation")
        answer = {}
        for factors, coefficient in compilation.polynomial.terms.items():
            extracted = self.laurent_coefficient(coefficient, power, precision)
            if extracted:
                answer[factors] = extracted
        return answer
