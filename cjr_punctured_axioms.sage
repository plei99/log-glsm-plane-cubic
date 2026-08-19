r"""CJR unit-removal axioms for plane-cubic effective vertices.

Chen--Janda--Ruan, *Punctured logarithmic R-maps*, Theorem 8.1
(equations (8.3) and (8.4)) gives string- and divisor-type equations for
removing a unit sector.  For a hypersurface their input has ``r=d=1``; the
unit sector therefore has contact order ``-1``.  In the plane-cubic case

    infinity_X = P^2,       N = O(3),       psi_DF = c1(O(-infinity)) = -3H.

Consequently every correction term in Theorem 8.1 is again an effective
vertex with an ambient insertion in ``{1,H,H^2}``; terms containing ``H^3``
vanish.  This module translates the cycle formulas into exact scalar linear
relations understood by :class:`InfinityVertexOrchestrator`.

Corollary 8.3 is included as a redundant but useful dilaton consistency
relation.  The genus-zero vanishing and the two dimension-zero genus-one base
values follow from Proposition 1.21 and equations (9.17), (9.21).
"""

load("cjr_probe_factory.sage")


class PlaneCubicPuncturedAxioms(SageObject):
    r"""Generate CJR relations among dimension-zero ``EffectiveVertex`` data.

    Insertions are powers of the hyperplane class ``H``.  Contacts and
    insertions are kept paired and put into the same canonical order as the
    localization graph compiler.
    """

    UNIT_CONTACT = ZZ(-1)
    MAX_H_POWER = ZZ(2)
    PSI_DF_MULTIPLIER = ZZ(-3)

    def __init__(self, coefficient_field=QQ):
        self.coefficient_field = coefficient_field

    @staticmethod
    def _insertion_tuple(vertex):
        if not vertex.insertions:
            return (ZZ.zero(),) * len(vertex.contacts)
        if len(vertex.insertions) != len(vertex.contacts):
            raise ValueError(
                "an effective vertex must have one insertion per contact"
            )
        return tuple(ZZ(power) for power in vertex.insertions)

    @classmethod
    def _canonical_vertex(cls, template, contacts, insertions, psi_min):
        pairs = tuple(sorted(zip(contacts, insertions)))
        return EffectiveVertex(
            template.genus,
            template.ambient_degree,
            tuple(contact for contact, _ in pairs),
            psi_min=psi_min,
            insertions=tuple(power for _, power in pairs),
            label=template.label,
        )

    @classmethod
    def _remove_marking(cls, vertex, index, psi_min=None):
        insertions = cls._insertion_tuple(vertex)
        contacts = vertex.contacts[:index] + vertex.contacts[index + 1:]
        insertions = insertions[:index] + insertions[index + 1:]
        if psi_min is None:
            psi_min = vertex.psi_min
        return cls._canonical_vertex(
            vertex, contacts, insertions, ZZ(psi_min)
        )

    @classmethod
    def _raise_insertion(cls, vertex, index, amount):
        amount = ZZ(amount)
        insertions = list(cls._insertion_tuple(vertex))
        insertions[index] += amount
        if insertions[index] > cls.MAX_H_POWER:
            return None
        return cls._canonical_vertex(
            vertex, vertex.contacts, tuple(insertions), vertex.psi_min
        )

    @staticmethod
    def _stable_after_forgetting(vertex):
        base_markings = len(vertex.contacts) - 1
        return 2 * vertex.genus - 2 + base_markings > 0

    @classmethod
    def _unit_positions(cls, vertex, insertion_power):
        insertions = cls._insertion_tuple(vertex)
        return tuple(
            index for index, (contact, insertion) in enumerate(
                zip(vertex.contacts, insertions)
            )
            if contact == cls.UNIT_CONTACT
            and insertion == ZZ(insertion_power)
        )

    @staticmethod
    def _require_no_contact_descendants(vertex):
        if any(vertex.contact_psi):
            raise UnsupportedGeometryError(
                "the implemented punctured unit-removal equations do not "
                "include ordinary cotangent powers at retained contacts"
            )

    @staticmethod
    def _add_term(terms, vertex, coefficient):
        if vertex is None or not coefficient:
            return
        terms[vertex] = terms.get(vertex, 0) + coefficient
        if not terms[vertex]:
            del terms[vertex]

    def _compiled_relation(self, vertex, name, remove_index, terms):
        ordered = tuple(
            (self.coefficient_field(coefficient), (factor,))
            for factor, coefficient in sorted(
                terms.items(), key=lambda item: item[0].order_key()
            )
            if coefficient
        )
        probe = ProbeSpec(
            vertex.genus,
            # EffectiveVertex and ProbeSpec both store the ambient P^2
            # degree.  The factor of three converts an intrinsic elliptic
            # degree to ambient degree, but no such conversion is present in
            # this punctured relation.  Multiplying here used to corrupt the
            # relation metadata (and checkpoint/probe-degree bookkeeping) in
            # every positive infinity-degree sector.
            vertex.ambient_degree,
            (),
            label="CJR %s remove=%s vertex=%r" % (
                name, remove_index, vertex.signature()
            ),
            psi_convention="log",
        )
        return CompiledLocalizationRelation(
            probe,
            t_power=0,
            coefficient_field=self.coefficient_field,
            known_gw=0,
            terms=ordered,
            compilation=None,
            is_complete=True,
        )

    def string_relation(self, vertex, unit_index):
        r"""Return equation (8.3) after substituting ``psi_DF=-3H``.

        If ``k=vertex.psi_min``, the returned relation is

        ``V_{sigma+1,k} = sum_j |c_j| sum_{a=0}^{k-1}
        (-3)^a V_{sigma,k-1-a}(H_j^a)``.
        """
        self._require_no_contact_descendants(vertex)
        if unit_index not in self._unit_positions(vertex, 0):
            raise ValueError("string removal needs a (-1,1) unit marking")
        if not self._stable_after_forgetting(vertex):
            raise ValueError("the forgotten curve is outside the stable range")

        terms = {}
        self._add_term(terms, vertex, 1)
        k = ZZ(vertex.psi_min)
        for descendant in range(k):
            base = self._remove_marking(
                vertex, unit_index, psi_min=k - 1 - descendant
            )
            for marking, contact in enumerate(base.contacts):
                raised = self._raise_insertion(base, marking, descendant)
                coefficient = (
                    abs(contact) * self.PSI_DF_MULTIPLIER**descendant
                )
                self._add_term(terms, raised, -coefficient)
        return self._compiled_relation(
            vertex, "8.3 string", unit_index, terms
        )

    def divisor_relation(self, vertex, unit_index):
        r"""Return equation (8.4) for the divisor ``H`` at a unit marking."""
        self._require_no_contact_descendants(vertex)
        if unit_index not in self._unit_positions(vertex, 1):
            raise ValueError("divisor removal needs a (-1,H) unit marking")
        if not self._stable_after_forgetting(vertex):
            raise ValueError("the forgotten curve is outside the stable range")

        terms = {}
        self._add_term(terms, vertex, 1)
        k = ZZ(vertex.psi_min)
        base_same_psi = self._remove_marking(
            vertex, unit_index, psi_min=k
        )
        # int_beta H is the ambient degree stored on the infinity vertex.
        self._add_term(
            terms, base_same_psi, -ZZ(vertex.ambient_degree)
        )
        for descendant in range(k):
            base = self._remove_marking(
                vertex, unit_index, psi_min=k - 1 - descendant
            )
            for marking, contact in enumerate(base.contacts):
                # D_{I_mu} * psi_DF^a = H * (-3H)^a.
                raised = self._raise_insertion(
                    base, marking, descendant + 1
                )
                coefficient = (
                    abs(contact) * self.PSI_DF_MULTIPLIER**descendant
                )
                self._add_term(terms, raised, -coefficient)
        return self._compiled_relation(
            vertex, "8.4 divisor", unit_index, terms
        )

    def dilaton_relation(self, vertex, unit_index):
        r"""Return Corollary 8.3 in the plane-cubic basis.

        For ``k=vertex.psi_min>=1`` this is

        ``V_{sigma+1,k}(1) + 3 V_{sigma+1,k-1}(H)
          = (2g-2+n) V_{sigma,k-1}``.

        It follows from the string and divisor relations and the balancing
        equation, but retaining it catches convention or sign mistakes.
        """
        self._require_no_contact_descendants(vertex)
        if unit_index not in self._unit_positions(vertex, 0):
            raise ValueError("dilaton removal needs a (-1,1) unit marking")
        if vertex.psi_min < 1:
            raise ValueError("the coefficient dilaton relation needs psi_min>=1")
        if not self._stable_after_forgetting(vertex):
            raise ValueError("the forgotten curve is outside the stable range")

        terms = {}
        self._add_term(terms, vertex, 1)
        insertions = list(self._insertion_tuple(vertex))
        insertions[unit_index] = ZZ.one()
        divisor_vertex = self._canonical_vertex(
            vertex,
            vertex.contacts,
            tuple(insertions),
            vertex.psi_min - 1,
        )
        # Corollary 8.3 contains ``t-ev^* psi_DF``.  Keep this tied to
        # the convention above rather than duplicating its negative by hand.
        self._add_term(
            terms, divisor_vertex, -self.PSI_DF_MULTIPLIER
        )
        base = self._remove_marking(
            vertex, unit_index, psi_min=vertex.psi_min - 1
        )
        multiplier = 2 * vertex.genus - 2 + len(base.contacts)
        self._add_term(terms, base, -multiplier)
        return self._compiled_relation(
            vertex, "Corollary 8.3 dilaton", unit_index, terms
        )

    @staticmethod
    def _relation_signature(relation):
        return tuple(
            (coefficient, tuple(factor.signature() for factor in factors))
            for coefficient, factors in relation.terms
        )

    def relations_for(self, vertex):
        """Return every applicable, deduplicated CJR unit-removal relation."""
        if not isinstance(vertex, EffectiveVertex):
            raise TypeError("CJR axioms apply to EffectiveVertex objects")
        if any(vertex.contact_psi):
            # The implemented scalar forms of Theorem 8.1 involve psi_min
            # and evaluations, but no ordinary cotangent lines at the
            # punctures being retained.  Do not silently apply them to the
            # contact-descendant vertices introduced by stabilization.
            return tuple()
        if not vertex.is_balanced() or not vertex.is_dimension_zero():
            return tuple()
        if not self._stable_after_forgetting(vertex):
            return tuple()

        relations = []
        for index in self._unit_positions(vertex, 0):
            relations.append(self.string_relation(vertex, index))
            if vertex.psi_min:
                relations.append(self.dilaton_relation(vertex, index))
        for index in self._unit_positions(vertex, 1):
            relations.append(self.divisor_relation(vertex, index))

        seen = set()
        answer = []
        for relation in relations:
            signature = self._relation_signature(relation)
            if signature in seen:
                continue
            seen.add(signature)
            answer.append(relation)
        return tuple(answer)

    def cycle_unit_axiom_test_classes(self, vertex, unit_index):
        r"""Return the test classes for which CJR II (8.8) is numerical.

        Theorem 8.4 is an identity of cycles on ``Mbar_{g,n+1}``:

            p_* (psi_min^k [R_{s+1}]) = pi^* p_* (psi_min^k [R_s])
              + sum_j |c_j| delta_{j,n+1} . pi^* p_* (sum_k' ev_j^*(psi_DF^k')
                                                      psi_min^{k-1-k'} [R_s]).

        Pairing with a test class ``Theta`` and applying the projection formula
        replaces ``Theta`` by ``pi_* Theta`` on the first term and by
        ``Theta|_delta`` on the boundary terms.  The variables tracked by this
        repository are dimension-zero effective invariants, so those cycles are
        multiples of a point class and only ``deg Theta = 0`` survives.  The
        single admissible test class is therefore the fundamental class, and
        (8.8) then reduces to the already-implemented string equation (8.3).

        This method exists to make that collapse explicit and checkable rather
        than to generate rows: a caller looking for additional independent
        relations will get an empty tuple beyond the unit, which is the correct
        mathematical answer, not a missing feature.
        """
        self._require_no_contact_descendants(vertex)
        if unit_index not in self._unit_positions(vertex, 0):
            raise ValueError("the cycle unit axiom needs a (-1,1) unit marking")
        if not vertex.is_dimension_zero():
            raise UnsupportedGeometryError(
                "the cycle unit axiom becomes numerical only on dimension-zero "
                "effective invariants; higher cycles need tautological "
                "unknowns that this solver does not carry"
            )
        return (ZZ.one(),)

    def cycle_unit_axiom_relation(self, vertex, unit_index):
        """Return (8.8) paired with its unique admissible test class."""
        self.cycle_unit_axiom_test_classes(vertex, unit_index)
        return self.string_relation(vertex, unit_index)

    @staticmethod
    def genus_one_closed_form():
        r"""Derive the genus-one base values from CJR II (9.17) and (9.21).

        Section 9.4 of *Punctured logarithmic R-maps* computes the ``g=n=1``,
        ``beta=0`` case completely.  With ``H`` the Hodge line bundle on
        ``Mbar_{1,1}`` and ``h=c_1(H^vee)=-lambda_1``, the two inputs are

            (9.17)  psi_min = -c_1(H^vee |X| N),
            (9.21)  [R]^red = c_dim(H^vee |X| (Omega^vee + O - N)),

        the latter read off as the dimension-one part of

            c(H^vee |X| Omega^vee) * c(H^vee) / c(H^vee |X| N).

        For the plane cubic ``infinity_X = P^2``, ``N = O(3)`` and
        ``Omega^vee = T_{P^2}``.  Truncating by ``lambda_1^2 = 0`` and
        ``H^3 = 0`` gives ``[R]^red = 3H^2`` and ``psi_min = lambda_1 - 3H``.

        Returning these rather than hard-coding ``1/8`` and ``0`` makes the
        repository's ``psi_min`` sign, the ``3H^2`` normalization, and the
        ambient-insertion convention a checkable consequence of the paper.
        """
        ring = PolynomialRing(QQ, ["lam", "H"])
        lam, hyperplane = ring.gens()

        def truncate(value):
            # lambda_1^2 = 0 on Mbar_{1,1} and H^3 = 0 on P^2.
            value = ring(value)
            return ring(sum(
                (
                    coefficient * monomial
                    for coefficient, monomial in zip(
                        value.coefficients(), value.monomials()
                    )
                    if monomial.degree(lam) <= 1
                    and monomial.degree(hyperplane) <= 2
                ),
                ring.zero(),
            ))

        h = -lam
        chern_tangent = (
            1 + (2 * h + 3 * hyperplane)
            + (h**2 + 3 * h * hyperplane + 3 * hyperplane**2)
        )
        chern_hodge = 1 + h
        normal = h + 3 * hyperplane
        inverse_normal = 1 - normal + normal**2
        total = truncate(truncate(chern_tangent * chern_hodge) * inverse_normal)
        reduced_cycle = ring(sum(
            (
                coefficient * monomial
                for coefficient, monomial in zip(
                    total.coefficients(), total.monomials()
                )
                if monomial.degree() == 2
            ),
            ring.zero(),
        ))
        psi_min = -normal

        def integrate(value):
            # int_{Mbar_{1,1}} lambda_1 = 1/24 and int_{P^2} H^2 = 1.
            value = ring(value)
            answer = QQ.zero()
            for coefficient, monomial in zip(
                    value.coefficients(), value.monomials()):
                if monomial.degree(lam) == 1 and monomial.degree(hyperplane) == 2:
                    answer += QQ(coefficient) / 24
            return answer

        return {
            "reduced_cycle": reduced_cycle,
            "psi_min": psi_min,
            "psi_min_value": integrate(truncate(psi_min * reduced_cycle)),
            "divisor_value": integrate(truncate(hyperplane * reduced_cycle)),
        }

    def known_value(self, vertex):
        r"""Return a CJR closed value, or ``None`` when no formula applies.

        For ``E=O(3)`` the genus-one base cycle is equation (9.21) on
        ``Mbar_{1,1} x P^2`` and ``psi_min`` is equation (9.17).  Expanding
        the Chern class gives

        ``[R_{1,(-1)}]^red = 3 H^2``, hence
        ``integral H = 0`` and ``integral psi_min = 1/8``.
        """
        if not isinstance(vertex, EffectiveVertex):
            raise TypeError("CJR values apply to EffectiveVertex objects")
        if any(vertex.contact_psi):
            return None
        if not vertex.is_balanced() or not vertex.is_dimension_zero():
            return None
        if vertex.genus == 0:
            # Proposition 1.21(1), since O(3) is nef.
            return self.coefficient_field.zero()
        if vertex.genus != 1:
            return None
        if vertex.ambient_degree != 0 \
                or any(contact != -1 for contact in vertex.contacts):
            # Proposition 1.21(2) and Lemma 9.4, since O(3) is ample.
            return self.coefficient_field.zero()
        if len(vertex.contacts) != 1:
            return None
        insertions = self._insertion_tuple(vertex)
        if vertex.psi_min == 0 and insertions == (ZZ.one(),):
            return self.coefficient_field.zero()
        if vertex.psi_min == 1 and insertions == (ZZ.zero(),):
            return self.coefficient_field(QQ(1) / 8)
        return None
