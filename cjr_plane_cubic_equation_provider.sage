r"""
A concrete CJR equation provider for plane-cubic one-point stationary theory.

Bloch--Okounkov gives

    sum_{g>=0} <tau_(2g-2)(pt)>_g z^(2g)
      = exp(sum_{m>=1} 2 E_(2m)(q) z^(2m)/(2m)!).

For the plane cubic, split each primitive block into

    2 E_(2m)(q)/(2m)! = A_(2m)(q) + b_(2m),

where ``A_(2m)`` has strictly positive q-degree and is the resummed
O(3)-twisted zero-level block, while ``b_(2m)`` is degree-zero infinity
data.  The coefficient of z^(2g) is a finite sum indexed by integer
partitions of g.  This yields a genuine triangular equation provider for
``InfinityVertexDP`` in every genus in this one-point sector.

The inferred ``b_(2m)`` are primitive stationary cumulants.  They package
the sum of CJR contact-resolved infinity graphs contributing to this sector;
they are not automatically equal to a single punctured vertex integral.
"""

import argparse
import json
import os
import sys


load("log_glsm_infinity_dp.sage")
load("cjr_bipartite_graphs.sage")


def primitive_stationary_vertex(genus):
    """Return the aggregate degree-zero primitive unknown of genus ``genus``."""
    genus = ZZ(genus)
    if genus <= 0:
        raise ValueError("a primitive stationary vertex has positive genus")
    # The one-contact representative satisfies the degree-zero balance
    # c+1=-(2g-2), hence c=-(2g-1).  The label records that the value is the
    # aggregate primitive cumulant, not a claim about this profile alone.
    return EffectiveVertex(
        genus,
        0,
        (-(2 * genus - 1),),
        label="aggregate_one_point_primitive_cumulant",
    )


def _partition_multiplicity(partition):
    multiplicities = {}
    for part in partition:
        multiplicities[part] = multiplicities.get(part, 0) + 1
    return prod(factorial(value) for value in multiplicities.values())


def primitive_o3_twisted_block(genus, max_degree):
    r"""Return ``A_(2g)=2/(2g)! sum_d sigma_(2g-1)(d) q^d``."""
    genus = ZZ(genus)
    max_degree = ZZ(max_degree)
    if genus <= 0:
        raise ValueError("genus must be positive")
    if max_degree < 0:
        raise ValueError("max_degree must be nonnegative")
    Q = PowerSeriesRing(QQ, "q", default_prec=max_degree + 1)
    q = Q.gen()
    weight = 2 * genus
    scale = QQ(2) / factorial(weight)
    return Q(sum(
        scale
        * sum(QQ(k)^(weight - 1) for k in divisors(degree))
        * q^degree
        for degree in range(1, max_degree + 1)
    )).add_bigoh(max_degree + 1)


def basic_minus_two_vertex(genus, ambient_degree=0):
    r"""Return the all-contact-``-2`` plane-cubic vertex.

    Balance forces the number of punctures to be

        n_2 = 2g-2-3D.

    This is an *analogue* of the basic effective vertices of CJR III,
    Lemma 10.12, and not an instance of them.  CJR II, Definition 9.19, admits
    a basic effective invariant only for ``g>=2`` and only when

        (3 - dim X + rk E)*(g-1) - int_beta c_1(K_X (x) det E) = 0.

    For ``X=P^2`` and ``E=O(3)`` the bundle ``K_X (x) det E`` is trivial, so
    the integral vanishes for every ``beta`` and the left side is ``2(g-1)``,
    which is zero only at ``g=1`` -- excluded by that definition.  Assumption
    9.22, needed by Corollary 9.24 to determine effective invariants from basic
    ones, likewise requires ``2(g-1)<=0`` and fails for every ``g>=2``.  The
    controlling quantity is ``3-dim Z``, which vanishes for the Calabi--Yau
    threefolds the theory targets but equals ``2`` for an elliptic curve.

    The vertex returned here is therefore a bookkeeping device.  Its star
    coefficient ``1/n_2!`` is confirmed independently by the graph compiler,
    but it carries none of the determination properties of Definition 9.19.
    """
    genus = ZZ(genus)
    ambient_degree = ZZ(ambient_degree)
    number_of_contacts = 2 * genus - 2 - 3 * ambient_degree
    if number_of_contacts <= 0:
        raise ValueError("the basic -2 vertex requires 2g-2-3D > 0")
    vertex = EffectiveVertex(
        genus,
        ambient_degree,
        (-2,) * number_of_contacts,
        label="basic_all_contact_minus_two",
    )
    if not vertex.is_balanced():
        raise AssertionError("the constructed basic vertex must be balanced")
    return vertex


def unstable_leaf_edge_pair_factor(contact_order, moving_weight=None):
    r"""Product of a hypersurface edge and its nonspecial unstable 0-leaf.

    From CJR III, Section 9.3, with ``w=t-ev^*H_infinity``, the edge gives

        c^c/(c! w^c),

    while the unstable 0-vertex gives ``w^2/c``.  Their product is

        c^(c-1)/(c! w^(c-2)).

    It is exactly one for contact order two, explaining the cancellation in
    the dominant basic-effective star graph.
    """
    contact_order = ZZ(contact_order)
    if contact_order <= 0:
        raise ValueError("an edge contact order must be positive")
    if moving_weight is None:
        moving_weight = SR.var("w")
    return (
        contact_order^(contact_order - 1)
        / (factorial(contact_order) * moving_weight^(contact_order - 2))
    )


def basic_star_diagonal_coefficient(vertex):
    r"""Return the dominant star-graph coefficient for a basic -2 vertex.

    Every edge/unstable-leaf pair cancels to one.  Only permutations of the
    identical leaves remain, giving ``1/n_2!``.
    """
    if not isinstance(vertex, EffectiveVertex):
        raise TypeError("vertex must be an EffectiveVertex")
    if not vertex.is_balanced() or not vertex.contacts:
        raise ValueError("the star target must be a balanced punctured vertex")
    if any(contact != -2 for contact in vertex.contacts):
        raise NotImplementedError(
            "the simple star cancellation is currently proved only for -2 contacts"
        )
    if unstable_leaf_edge_pair_factor(2, QQ.one()) != 1:
        raise AssertionError("contact-two edge/leaf cancellation failed")
    return QQ(1) / factorial(len(vertex.contacts))


class PlaneCubicOnePointEquationProvider:
    r"""Integer-partition CJR provider for primitive stationary cumulants."""

    def __init__(self, max_genus, verification_degree=4):
        self.max_genus = ZZ(max_genus)
        self.verification_degree = ZZ(verification_degree)
        if self.max_genus <= 0:
            raise ValueError("max_genus must be positive")
        if self.verification_degree < 0:
            raise ValueError("verification_degree must be nonnegative")
        self.vertices = {
            genus: primitive_stationary_vertex(genus)
            for genus in range(1, self.max_genus + 1)
        }
        self._known = {}
        self._equations = {}

    def known_stationary_series(self, genus):
        genus = ZZ(genus)
        if genus not in self.vertices:
            raise ValueError("genus lies outside this provider's truncation")
        if genus not in self._known:
            self._known[genus] = connected_stationary_qseries(
                [2 * genus - 2], self.verification_degree
            )
        return self._known[genus]

    def graph_enumerator(self, genus, ambient_degree=0):
        r"""Return the finite CJR graph search for this one-point probe.

        The returned graphs retain contact orders and automorphisms.  The
        partition formula in ``equation_for_genus`` is a resummation of their
        contributions; assigning a contribution to every individual graph
        additionally requires the full twisted descendant and gluing rules.
        """
        genus = ZZ(genus)
        ambient_degree = ZZ(ambient_degree)
        if genus not in self.vertices:
            raise ValueError("genus lies outside this provider's truncation")
        return PlaneCubicGraphEnumerator(genus, 1, ambient_degree)

    def localization_graphs(self, genus, ambient_degree=0):
        """Return the one-point localization graph isomorphism classes."""
        return self.graph_enumerator(genus, ambient_degree).graphs()

    def combinatorial_infinity_dependencies(self, genus, ambient_degree=0):
        r"""List ``(g,D,contacts)`` occurring at infinity in the graph sum.

        These are the base contact profiles only.  Descendant powers and
        evaluation insertions arise when the localization contribution is
        expanded and are therefore deliberately not guessed here.
        """
        return tuple(sorted(set(
            data
            for graph in self.localization_graphs(genus, ambient_degree)
            for data in graph.all_infinity_vertex_data()
        )))

    def equation_for_genus(self, genus):
        r"""Return the degree-zero localization equation in genus ``genus``.

        Expanding ``exp(sum b_(2m) x^m)`` gives one monomial for every
        integer partition of ``genus``, with inverse multiplicity factorial.
        The singleton partition is the nonzero diagonal term.
        """
        genus = ZZ(genus)
        if genus not in self.vertices:
            raise ValueError("genus lies outside this provider's truncation")
        if genus in self._equations:
            return self._equations[genus]

        target = self.vertices[genus]
        polynomial, known_series = self.known_stationary_series(genus)
        terms = []
        for partition in Partitions(genus):
            factors = tuple(self.vertices[ZZ(part)] for part in partition)
            coefficient = QQ(1) / _partition_multiplicity(partition)
            terms.append((coefficient, factors))

        equation = LocalizationEquation(
            target,
            known_gw=known_series[0],
            twisted_zero_level=0,
            terms=terms,
            probe_label="degree-zero <tau_%s(pt)>_(g=%s) from %s"
                        % (2 * genus - 2, genus, polynomial),
        )
        self._equations[genus] = equation
        return equation

    def __call__(self, vertex):
        for genus, candidate in self.vertices.items():
            if vertex == candidate:
                return self.equation_for_genus(genus)
        raise KeyError("the vertex is not in this provider's truncation: %r" % vertex)

    def solver(self):
        return InfinityVertexDP(self)

    def solve(self):
        solver = self.solver()
        ordered = tuple(self.vertices[g] for g in range(1, self.max_genus + 1))
        return solver, solver.solve_up_to(ordered)

    def expected_bernoulli_value(self, genus):
        genus = ZZ(genus)
        return -bernoulli(2 * genus) / (
            2 * genus * factorial(2 * genus)
        )

    def reconstructed_stationary_series(self, genus, values):
        r"""Assemble all O(3)-twisted and inferred infinity graph terms."""
        genus = ZZ(genus)
        if genus not in self.vertices:
            raise ValueError("genus lies outside this provider's truncation")
        max_degree = self.verification_degree
        Q = PowerSeriesRing(QQ, "q", default_prec=max_degree + 1)
        primitives = [Q.zero()]
        for part in range(1, genus + 1):
            primitives.append(
                primitive_o3_twisted_block(part, max_degree)
                + values[self.vertices[part]]
            )

        # If F(x)=exp(P(x)), then g*F_g=sum_{m=1}^g m*P_m*F_(g-m).
        assembled = [Q.one()]
        for total_genus in range(1, genus + 1):
            assembled.append(sum(
                part * primitives[part] * assembled[total_genus - part]
                for part in range(1, total_genus + 1)
            ) / total_genus)
        return assembled[genus]

    def verify(self, values=None):
        if values is None:
            solver, values = self.solve()
        checks = {}
        for genus in range(1, self.max_genus + 1):
            vertex = self.vertices[genus]
            polynomial, known = self.known_stationary_series(genus)
            reconstructed = self.reconstructed_stationary_series(genus, values)
            checks[genus] = {
                "diagonal_nonzero": any(
                    factors == (vertex,) and coefficient != 0
                    for coefficient, factors in self.equation_for_genus(genus).terms
                ),
                "bernoulli_value": (
                    values[vertex] == self.expected_bernoulli_value(genus)
                ),
                "all_q_coefficients": reconstructed == known,
                "bo_polynomial": str(polynomial),
            }
        return checks

    def report(self):
        solver, values = self.solve()
        checks = self.verify(values)
        return {
            "sector": "one-point stationary plane-cubic theory",
            "max_genus": int(self.max_genus),
            "verification_degree": int(self.verification_degree),
            "vertices": [
                {
                    "genus": genus,
                    "representative_contacts": [
                        int(c) for c in self.vertices[genus].contacts
                    ],
                    "value": str(values[self.vertices[genus]]),
                    "equation_terms": len(
                        self.equation_for_genus(genus).terms
                    ),
                    "checks": checks[genus],
                }
                for genus in range(1, self.max_genus + 1)
            ],
        }


def print_plane_cubic_provider_report(max_genus=4, max_degree=6, as_json=False):
    provider = PlaneCubicOnePointEquationProvider(max_genus, max_degree)
    report = provider.report()
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print("Plane-cubic one-point CJR equation provider")
    print("genus | representative contact | infinity cumulant | graph terms")
    print("------+------------------------+-------------------+------------")
    for vertex in report["vertices"]:
        print(
            "%5s | %22s | %17s | %s"
            % (
                vertex["genus"],
                vertex["representative_contacts"],
                vertex["value"],
                vertex["equation_terms"],
            )
        )
    print("All q-series checks:", all(
        all(
            item["checks"][name]
            for name in ("diagonal_nonzero", "bernoulli_value", "all_q_coefficients")
        )
        for item in report["vertices"]
    ))


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-genus", type=int, default=4)
    parser.add_argument("--max-degree", type=int, default=6)
    parser.add_argument("--json", action="store_true")
    arguments = parser.parse_args()
    print_plane_cubic_provider_report(
        arguments.max_genus, arguments.max_degree, arguments.json
    )


_entrypoint = os.path.basename(sys.argv[0])
if _entrypoint in (
        "cjr_plane_cubic_equation_provider.sage",
        "cjr_plane_cubic_equation_provider.sage.py",
):
    _main()
