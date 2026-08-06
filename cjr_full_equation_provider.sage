r"""End-to-end CJR equation provider built from the enumerator and compiler."""

import argparse
import json
import os
import sys

load("cjr_graph_contributions.sage")
load("elliptic_probe_values.sage")
load("cjr_probe_factory.sage")


class PlaneCubicFullEquationProvider(SageObject):
    r"""Compile exact coefficient relations and feed them to the DP solver.

    Stable zero vertices are computed by full base-torus localization with
    admcycles integration.  A partial injected backend remains supported for
    diagnostics and externally registered exact values.
    """

    def __init__(self, rings=None, twisted_backend=None,
                 known_backend=None, laurent_precision=16):
        self.rings = rings or PlaneCubicCoefficientRing(laurent_precision)
        self.twisted_backend = twisted_backend or FullTwistedZeroVertexBackend(
            self.rings
        )
        self.compiler = PlaneCubicGraphContributionCompiler(
            self.rings, self.twisted_backend
        )
        self.known_backend = known_backend or EllipticProbeValueBackend()
        self._compilations = {}
        self._relations = {}
        self._assignments = {}

    def compilation(self, probe):
        if probe not in self._compilations:
            self._compilations[probe] = self.compiler.compile_probe(probe)
        return self._compilations[probe]

    def relation(self, probe, t_power=0, require_complete=True):
        key = probe, ZZ(t_power)
        if key in self._relations:
            relation = self._relations[key]
            if require_complete and not relation.is_complete:
                relation.compilation.require_complete()
            return relation

        compilation = self.compilation(probe)
        if require_complete:
            compilation.require_complete()
        extracted = self.compiler.extracted_polynomial(
            compilation, t_power, self.rings.laurent_precision
        )
        constant = extracted.pop(tuple(), self.rings.base_field.zero())
        known = self.known_backend.reduced_log_glsm_value(probe) \
            if ZZ(t_power) == 0 else QQ.zero()
        relation = CompiledLocalizationRelation(
            probe, t_power, self.rings.base_field,
            known_gw=known,
            twisted_zero_level=constant,
            terms=tuple((coefficient, factors)
                        for factors, coefficient in extracted.items()),
            compilation=compilation,
        )
        self._relations[key] = relation
        return relation

    def assign_probe(self, target, probe, t_power=0):
        if not isinstance(target, EffectiveVertex):
            raise TypeError("target must be an EffectiveVertex")
        self._assignments[target] = probe, ZZ(t_power)

    def __call__(self, target):
        if target not in self._assignments:
            raise KeyError("no probe has been assigned to %r" % target)
        probe, t_power = self._assignments[target]
        return self.relation(probe, t_power, require_complete=True).equation_for(target)

    def solver(self, initial_values=None):
        return InfinityVertexDP(
            self,
            initial_values=initial_values,
            coefficient_field=self.rings.base_field,
        )

    def candidate_relations(self, genus, ambient_degree, max_markings=4,
                            t_powers=(0,), require_complete=False,
                            include_unit_relatives=False):
        factory = ProbeFactory(self.rings.base_field)
        stationary = factory.stationary_candidates(
            genus, ambient_degree, max_markings=max_markings
        )
        probes = list(stationary)
        if include_unit_relatives:
            for probe in stationary:
                probes.extend(factory.unit_relatives(probe))
        relations = []
        for probe in probes:
            for t_power in t_powers:
                relations.append(self.relation(
                    probe, t_power, require_complete=require_complete
                ))
        return tuple(relations)

    def relation_report(self, probe, t_power=0, require_complete=False):
        relation = self.relation(
            probe, t_power, require_complete=require_complete
        )
        return {
            "probe": str(probe),
            "t_power": int(t_power),
            "complete": relation.is_complete,
            "known_gw": str(relation.known_gw),
            "twisted_zero_level": str(relation.twisted_zero_level),
            "terms": [
                {
                    "coefficient": str(coefficient),
                    "factors": [factor.to_record() for factor in factors],
                }
                for coefficient, factors in relation.terms
            ],
            "compilation": relation.compilation.report(),
        }


def _linear_vertex_coefficient(relation, target):
    """Collect the coefficient of the one-factor monomial ``target``."""
    return relation.coefficient_field(sum(
        coefficient
        for coefficient, factors in relation.terms
        if factors == (target,)
    ))


def genus_two_profile_split_rank_witness(provider=None):
    r"""Compile an exact rank witness separating genus-two contact sectors.

    The string and dilaton relatives of ``<tau_2(pt)>_2`` give independent
    rows on a contact-resolved ``(-3,-1)``/``(-2,-2)`` pair.  This is one
    block in the descending-Laurent reconstruction that ultimately separates
    the primitive ``(-3)`` and ``(-2,-2)`` contributions.

    This calculation compiles two 28-graph probes and can take several
    minutes.  It returns the provider as well as the exact relations so the
    caller can continue constructing the full dependency closure.
    """
    provider = provider or PlaneCubicFullEquationProvider(
        laurent_precision=8
    )
    factory = ProbeFactory(provider.rings.base_field)
    one_point = ProbeSpec.stationary(2, 0, (2,))
    string_probe, dilaton_probe = factory.unit_relatives(one_point)
    relations = (
        provider.relation(string_probe, t_power=0),
        provider.relation(dilaton_probe, t_power=0),
    )
    targets = (
        EffectiveVertex(
            2, 0, (-3, -1), psi_min=0, insertions=(0, 2)
        ),
        EffectiveVertex(
            2, 0, (-2, -2), psi_min=0, insertions=(1, 1)
        ),
    )
    rows = tuple(
        tuple(_linear_vertex_coefficient(relation, target)
              for target in targets)
        for relation in relations
    )
    matrix = Matrix(provider.rings.base_field, rows)
    return {
        "provider": provider,
        "probes": (string_probe, dilaton_probe),
        "relations": relations,
        "targets": targets,
        "matrix": matrix,
        "determinant": matrix.det(),
        "full_rank": matrix.rank() == len(targets),
    }


def genus_two_profile_split_report(provider=None):
    """Return the rank witness in a JSON-serializable form."""
    witness = genus_two_profile_split_rank_witness(provider)
    return {
        "targets": [target.to_record() for target in witness["targets"]],
        "probes": [str(probe) for probe in witness["probes"]],
        "rows": [[str(value) for value in row]
                 for row in witness["matrix"].rows()],
        "determinant": str(witness["determinant"]),
        "full_rank": bool(witness["full_rank"]),
    }


def genus_two_resummed_end_to_end(max_intrinsic_degree=8):
    r"""Complete aggregate genus-two reconstruction through a degree bound.

    This milestone uses the proven flat-coordinate O(3)-twisted blocks.  It
    verifies every coefficient of ``<pt psi^2>`` while keeping the primitive
    genus-two infinity value explicitly labeled as an aggregate.
    """
    max_intrinsic_degree = ZZ(max_intrinsic_degree)
    reconstruction = reconstruct_genus_two_infinity_vertices(
        max_intrinsic_degree
    )
    known = reconstruction["known_gw_series"]
    reconstructed = reconstruction["reconstructed_gw_series"]
    if known != reconstructed:
        raise ArithmeticError("genus-two aggregate reconstruction failed")
    return {
        "max_intrinsic_degree": max_intrinsic_degree,
        "known_series": known,
        "reconstructed_series": reconstructed,
        "b2": reconstruction["genus_one_basic_contact_minus_one"],
        "b4_aggregate": reconstruction["genus_two_primitive_profile_combination"],
        "coefficients": tuple(known[degree]
                              for degree in range(max_intrinsic_degree + 1)),
        "ambient_degrees": tuple(3 * degree
                                 for degree in range(max_intrinsic_degree + 1)),
    }


def genus_two_end_to_end_report(max_intrinsic_degree=8):
    result = genus_two_resummed_end_to_end(max_intrinsic_degree)
    return {
        "invariant": "<pt psi^2>_(2,1,d)",
        "degree_convention": "intrinsic degree; ambient P2 degree is 3d",
        "max_intrinsic_degree": int(result["max_intrinsic_degree"]),
        "b2": str(result["b2"]),
        "b4_aggregate": str(result["b4_aggregate"]),
        "values": [
            {
                "intrinsic_degree": degree,
                "ambient_degree": 3 * degree,
                "value": str(result["coefficients"][degree]),
            }
            for degree in range(result["max_intrinsic_degree"] + 1)
        ],
        "reconstruction_check": (
            result["known_series"] == result["reconstructed_series"]
        ),
        "contact_resolution": (
            "b4 is the aggregate primitive genus-two cumulant; the generic "
            "compiler keeps (-3) and (-2,-2) distinct"
        ),
    }


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=8,
                        help="maximum intrinsic elliptic-curve degree")
    parser.add_argument("--json", action="store_true")
    parser.add_argument(
        "--profile-split-rank", action="store_true",
        help="compile the slower string/dilaton contact-splitting rank witness",
    )
    arguments = parser.parse_args()
    if arguments.profile_split_rank:
        report = genus_two_profile_split_report()
        if arguments.json:
            print(json.dumps(report, indent=2, sort_keys=True))
            return
        print("Genus-two contact-profile rank witness")
        print("matrix:")
        for row in report["rows"]:
            print("  ", row)
        print("determinant =", report["determinant"])
        print("full rank =", report["full_rank"])
        return
    report = genus_two_end_to_end_report(arguments.max_degree)
    if arguments.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print("Plane-cubic log-GLSM genus-two end-to-end aggregate")
    print("intrinsic degree | ambient degree | <pt psi^2>_(2,1,d)")
    print("-----------------+----------------+--------------------")
    for item in report["values"]:
        print("%16s | %14s | %s" % (
            item["intrinsic_degree"], item["ambient_degree"], item["value"]
        ))
    print("b2 =", report["b2"])
    print("b4 (aggregate) =", report["b4_aggregate"])
    print("reconstruction check =", report["reconstruction_check"])


if __name__ == "__main__" and os.path.basename(sys.argv[0]) in (
        "cjr_full_equation_provider.sage",
        "cjr_full_equation_provider.sage.py"):
    _main()
