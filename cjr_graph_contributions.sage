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
    r"""Exact graph-to-polynomial compiler in the compact plane-cubic sector.

    Both marking-cotangent conventions are supported.  For ``"log"``, CJR
    equation (8.21) evaluates the coarse log-domain cotangent line on an
    unstable marked zero tail.  For ``"stabilized"``, the pullback through
    ``st`` moves that descendant to the corresponding contact marking on the
    adjacent infinity vertex.  Markings on stable zero vertices are unchanged
    because their components survive stabilization.
    """

    def __init__(self, rings=None, twisted_backend=None,
                 base_weight_specialization="nonequivariant"):
        self.rings = rings or PlaneCubicCoefficientRing(16)
        self.factors = PlaneCubicGraphFactors(self.rings)
        self.twisted_backend = twisted_backend or FullTwistedZeroVertexBackend(
            self.rings
        )
        if base_weight_specialization == "nonequivariant":
            pass
        elif isinstance(base_weight_specialization, str):
            raise ValueError(
                "base weights must be 'nonequivariant', None, or a triple"
            )
        elif base_weight_specialization is not None:
            if len(base_weight_specialization) != 3 \
                    or len(set(base_weight_specialization)) != 3:
                raise ValueError("base weights must be three distinct values")
            base_weight_specialization = tuple(
                QQ(value) for value in base_weight_specialization
            )
        self.base_weight_specialization = base_weight_specialization
        self._specialized_coefficient_cache = {}
        self._graph_cache = {}

    def _specialize_coefficient(self, value):
        """Apply the auxiliary base-weight specialization once per value."""
        value = self.rings.full(value)
        if self.base_weight_specialization is None:
            return value
        specialized = self._specialized_coefficient_cache.get(value)
        if specialized is None:
            if self.base_weight_specialization == "nonequivariant":
                specialized = self.rings.nonequivariant_base_limit(value)
            else:
                specialized = self.rings.specialize_base_weights(
                    value, self.base_weight_specialization
                )
            self._specialized_coefficient_cache[value] = specialized
        return specialized

    def _graphs_for_probe(self, probe):
        """Reuse graph enumeration across probes with the same discrete data."""
        key = (
            ZZ(probe.genus), ZZ(probe.marking_count),
            ZZ(probe.ambient_degree),
        )
        if key not in self._graph_cache:
            self._graph_cache[key] = PlaneCubicGraphEnumerator(*key).graphs()
        return self._graph_cache[key]

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
        r"""Return ``(coefficient, {edge: (H_power, contact_psi)})`` options."""
        ordinary = tuple(
            self._ordinary_insertion(probe, marking)
            for marking in graph.markings_at_zero_vertex(zero)
        )
        edges = self._incident_edges(graph, zero)
        # The O(3) Euler class is equivariant.  Its powers of t can replace
        # part (or all) of its ordinary Chow codimension, so the flag-series
        # truncation is controlled by the virtual dimension of the ambient
        # stable-map space, not by the dimension after subtracting the
        # virtual Euler rank.  The old twisted-dimension bound discarded
        # genuine overdimensioned equivariant descendants.  Already on
        # Mbar_0,2(P2,1), the missing psi term changes the genus-one CJR
        # degree-one consistency equation by 1/24.
        plane_dimension = PlaneCubicDimension.plane_virtual_dimension(
            graph.zero_genera[zero], graph.zero_degrees[zero],
            len(ordinary) + len(edges),
        )
        ordinary_codimension = sum(item.codimension for item in ordinary)
        maximum_flag_codimension = plane_dimension - ordinary_codimension
        if maximum_flag_codimension < 0:
            return tuple()
        local_options = []
        for edge_index, _, contact in edges:
            options_for_edge = []
            edge_class = self.factors.edge(contact)
            for psi_power in range(maximum_flag_codimension + 1):
                flag_class = self.factors.stable_flag_descendant(contact, psi_power)
                for diagonal_power in range(3):
                    dual_power = p2_dual_power(diagonal_power)
                    for flag_h_power, flag_coefficient in enumerate(
                            flag_class.coefficients):
                        zero_h_power = diagonal_power + flag_h_power
                        if zero_h_power > 2 or not flag_coefficient:
                            continue
                        if psi_power + zero_h_power \
                                > maximum_flag_codimension:
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
                infinity_insertions[edge_index] = (infinity_power, ZZ.zero())
            if sum(item.codimension for item in zero_insertions) \
                    > plane_dimension:
                continue
            request = TwistedZeroVertexRequest(
                graph.zero_genera[zero], graph.zero_degrees[zero], zero_insertions,
                label="CJR stable zero vertex",
            )
            value = self.twisted_backend.evaluate(request)
            # This value has already summed all P2-fixed loci at the stable
            # zero vertex.  Evaluating the auxiliary base weights here keeps
            # all subsequent graph products in QQ(t), rather than first
            # constructing enormous expressions in three lambda variables.
            value = self._specialize_coefficient(value)
            if value:
                answer.append((coefficient * value, infinity_insertions))
        return tuple(answer)

    def _nonspecial_zero_options(self, graph, zero):
        edge_index, _, contact = self._incident_edges(graph, zero)[0]
        local_class = self.factors.nonspecial_edge_pair(contact)
        return tuple(
            (coefficient, {edge_index: (h_power, ZZ.zero())})
            for h_power, coefficient in enumerate(local_class.coefficients)
            if coefficient
        )

    def _marked_zero_options(self, probe, graph, zero):
        edge_index, _, contact = self._incident_edges(graph, zero)[0]
        marking = graph.markings_at_zero_vertex(zero)[0]
        ordinary = self._ordinary_insertion(probe, marking)
        ordinary_class = self._class_for_insertion(ordinary)
        if probe.psi_convention == "stabilized":
            # The fiber edge and the marked rational zero tail are contracted
            # by st.  The surviving stable-map marking is the infinity
            # half-edge, so st^*(bar_psi_i) restricts to psi at that contact.
            local_psi_power = ZZ.zero()
            contact_psi_power = ordinary.psi_power
        else:
            # CJR (8.21): psi on the coarse log-domain curve restricts to
            # H_infinity/c-t on the unstable marked tail.
            local_psi_power = ordinary.psi_power
            contact_psi_power = ZZ.zero()
        local_class = self.factors.edge(contact) \
            * self.factors.zero_marked(contact, local_psi_power) \
            * ordinary_class
        return tuple(
            (coefficient, {
                edge_index: (h_power, contact_psi_power)
            })
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
                        left_index: (p2_dual_power(left_basis), ZZ.zero()),
                        right_index: (p2_dual_power(right_basis), ZZ.zero()),
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
            (
                -contact,
                edge_insertions[edge_index][0],
                edge_insertions[edge_index][1],
            )
            for edge_index, (_, vertex, contact) in enumerate(graph.edges)
            if vertex == infinity
        )
        contacts = tuple(contact for contact, _, _ in paired)
        insertions = tuple(h_power for _, h_power, _ in paired)
        contact_psi = tuple(power for _, _, power in paired)
        vertex = EffectiveVertex(
            graph.infinity_genera[infinity],
            graph.infinity_degrees[infinity],
            contacts,
            psi_min=psi_power,
            insertions=insertions,
            contact_psi=contact_psi if any(contact_psi) else tuple(),
        )
        if not vertex.is_balanced():
            raise AssertionError("enumerated infinity vertex lost its balance")
        return vertex

    def _infinity_options(self, graph, infinity, edge_insertions):
        insertion_codimension = sum(
            sum(edge_insertions[edge_index])
            for edge_index, (_, vertex, _) in enumerate(graph.edges)
            if vertex == infinity
        )
        contacts = tuple(
            -contact
            for edge_index, (_, vertex, contact) in enumerate(graph.edges)
            if vertex == infinity
        )
        psi_power = PlaneCubicDimension.infinity_required_psi_min_power(
            graph.infinity_genera[infinity],
            graph.infinity_degrees[infinity],
            contacts,
            insertion_codimension,
        )
        if psi_power < 0:
            return tuple()
        vertex = self._effective_vertex(
            graph, infinity, edge_insertions, psi_power
        )
        if not vertex.is_dimension_zero():
            raise AssertionError(
                "enumerated infinity vertex violates virtual dimension"
            )
        return ((
            self.factors.infinity_descendant_coefficient(psi_power),
            vertex,
        ),)

    def _combined_zero_states(self, graph, zero_option_lists):
        r"""Aggregate equal edge states before constructing infinity data.

        Many stable-zero flag expansions produce the same collection of
        infinity insertions with different local coefficients.  Expanding the
        Cartesian product first used to rebuild the same EffectiveVertex and
        coerce the same universal descendant coefficient millions of times.
        Summation is linear, so combine those coefficients while the state is
        still the small hashable tuple ``(edge, (H-power, psi-power))``.
        """
        localization_weight = self.rings.full(graph.localization_weight())
        states = {tuple(): localization_weight}

        # Convolve one zero vertex at a time.  In particular, collapse equal
        # local edge states *before* taking products with the remaining zero
        # vertices.  Stable-zero expansions can contain many terms with the
        # same eventual (H, psi) data; constructing their full Cartesian
        # product before collecting like terms makes those multiplicities
        # grow exponentially.
        for options in zero_option_lists:
            local_states = {}
            for local_coefficient, local_insertions in options:
                local_signature = tuple(sorted(local_insertions.items()))
                local_states[local_signature] = local_states.get(
                    local_signature, self.rings.full_field.zero()
                ) + local_coefficient
                if not local_states[local_signature]:
                    del local_states[local_signature]

            next_states = {}
            for signature, coefficient in states.items():
                occupied_edges = frozenset(edge for edge, _ in signature)
                for local_signature, local_coefficient in local_states.items():
                    if any(edge in occupied_edges for edge, _ in local_signature):
                        raise AssertionError(
                            "an edge received insertions from two zero vertices"
                        )
                    combined_signature = tuple(sorted(
                        signature + local_signature
                    ))
                    next_states[combined_signature] = next_states.get(
                        combined_signature, self.rings.full_field.zero()
                    ) + coefficient * local_coefficient
                    if not next_states[combined_signature]:
                        del next_states[combined_signature]
            states = next_states

        if any(len(signature) != graph.edge_count for signature in states):
            raise AssertionError(
                "every localization edge needs one infinity insertion"
            )
        return tuple(
            (coefficient, dict(signature))
            for signature, coefficient in sorted(states.items())
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

        polynomial = EffectivePolynomial(self.rings.full_field)
        for coefficient, edge_insertions in self._combined_zero_states(
                graph, zero_option_lists):
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
        if not probe.is_stable_map_type():
            raise ValueError(
                "the underlying connected stable-map data (g,n,d)="
                "(%s,%s,%s) are unstable"
                % (probe.genus, probe.marking_count, probe.ambient_degree)
            )
        graphs = self._graphs_for_probe(probe)
        return ProbeCompilation(
            probe,
            tuple(self.compile_graph(probe, graph) for graph in graphs),
            self.rings.full_field,
        )

    def laurent_coefficient(self, value, power=0, precision=None):
        # ``precision`` remains in the public signature for compatibility.
        # Exact recurrence at t=infinity has no truncation boundary.
        if self.base_weight_specialization is not None:
            value = self._specialize_coefficient(value)
        return self.rings.laurent_coefficient_at_infinity(value, power)

    def extracted_polynomial(self, compilation, power=0, precision=None):
        if not isinstance(compilation, ProbeCompilation):
            raise TypeError("expected a ProbeCompilation")
        answer = {}
        for factors, coefficient in compilation.polynomial.terms.items():
            extracted = self.laurent_coefficient(coefficient, power, precision)
            if extracted:
                answer[factors] = extracted
        return answer
