r"""End-to-end CJR equation provider built from the enumerator and compiler."""

import argparse
import json
import os
import sys

load("cjr_graph_contributions.sage")
load("elliptic_probe_values.sage")
load("cjr_probe_factory.sage")
load("cjr_infinity_chow.sage")
load("cjr_stabilization_boundary.sage")
load("o3_givental_teleman.sage")


class PlaneCubicFullEquationProvider(SageObject):
    r"""Compile exact coefficient relations and feed them to the DP solver.

    Stable zero vertices are computed by full base-torus localization with
    admcycles integration.  A partial injected backend remains supported for
    diagnostics and externally registered exact values.

    The graph compiler supports both CJR's log-domain marking psi classes and
    stabilization pullbacks of ordinary GW psi classes.  The latter transfer
    descendants on contracted marked zero tails to contact descendants on the
    adjacent infinity vertex, so ``EllipticProbeValueBackend`` supplies their
    exact compact side.  In the effective basis, the CJR-I virtual
    stabilization-boundary comparison instead supplies that same known
    number to a distinct log-domain probe, leaving only ``psi_min``
    descendants at infinity.
    """

    def __init__(self, rings=None, twisted_backend=None,
                 known_backend=None, chow_backend=None, laurent_precision=16,
                 stabilization_comparison=None,
                 zero_vertex_cache=None, autosave_zero_vertices=False,
                 hodge_cache=None, zero_vertex_strategy="localization",
                 base_weight_specialization="nonequivariant"):
        self.rings = rings or PlaneCubicCoefficientRing(laurent_precision)
        if base_weight_specialization == "nonequivariant":
            self.base_weight_specialization = "nonequivariant"
        elif isinstance(base_weight_specialization, str):
            raise ValueError(
                "base weights must be 'nonequivariant', None, or a triple"
            )
        elif base_weight_specialization is None:
            self.base_weight_specialization = None
        else:
            self.base_weight_specialization = tuple(
                QQ(value) for value in base_weight_specialization
            )
        self.zero_vertex_strategy = str(zero_vertex_strategy)
        if self.zero_vertex_strategy not in (
                "localization", "hybrid", "givental"):
            raise ValueError(
                "zero_vertex_strategy must be localization, hybrid, or givental"
            )
        if hodge_cache is None and zero_vertex_cache is not None:
            hodge_cache = str(zero_vertex_cache) + ".hodge.sqlite"
        if twisted_backend is not None:
            self.twisted_backend = twisted_backend
        else:
            localization_cache = zero_vertex_cache
            if localization_cache is not None \
                    and self.base_weight_specialization is not None:
                localization_cache = specialized_zero_vertex_cache_path(
                    localization_cache, self.base_weight_specialization
                )
            localization = FullTwistedZeroVertexBackend(
                self.rings, cache_path=localization_cache,
                autosave=autosave_zero_vertices,
                hodge_cache_path=hodge_cache,
                base_weight_specialization=self.base_weight_specialization,
            )
            if self.zero_vertex_strategy == "localization":
                self.twisted_backend = localization
            else:
                calibrated = O3CalibratedGiventalBackend(self.rings)
            if self.zero_vertex_strategy == "hybrid":
                self.twisted_backend = HybridTwistedZeroVertexBackend(
                    localization, calibrated,
                )
            elif self.zero_vertex_strategy == "givental":
                self.twisted_backend = calibrated
        self.compiler = PlaneCubicGraphContributionCompiler(
            self.rings, self.twisted_backend,
            base_weight_specialization=self.base_weight_specialization,
        )
        self.known_backend = known_backend or EllipticProbeValueBackend()
        self.stabilization_comparison = (
            stabilization_comparison
            or StabilizationBoundaryComparison(self.known_backend)
        )
        self.chow_backend = chow_backend or CJRInfinityChowBackend()
        self._compilations = {}
        self._relations = {}
        self._assignments = {}

    def compilation(self, probe):
        if probe not in self._compilations:
            self._compilations[probe] = self.compiler.compile_probe(probe)
        return self._compilations[probe]

    def relation(self, probe, t_power=0, require_complete=True):
        t_power = ZZ(t_power)
        chow_degree = self.chow_backend.require_scalar(probe, t_power)
        if t_power == 0 and not probe.is_dimension_zero():
            raise UnsupportedGeometryError(
                "a t^0 under-dimensioned probe is a positive-dimensional "
                "compact Chow class, not a numerical invariant"
            )
        key = probe, t_power
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
        if t_power != 0:
            known = QQ.zero()
        elif probe.has_descendants and probe.psi_convention == "log":
            # CJR I's virtual comparison includes the stabilization-boundary
            # correction.  Compile the graph with the log-domain psi class,
            # but evaluate its compact side using the corresponding
            # stabilized hypersurface descendant.
            known = self.stabilization_comparison.reduced_value(probe)
        else:
            known = self.known_backend.reduced_log_glsm_value(probe)
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

    def boundary_compared_relation(self, stabilized_probe, t_power=0,
                                   require_complete=True):
        r"""Compile a known stabilized descendant in the effective basis.

        The returned relation is deliberately keyed by the corresponding
        log-domain probe.  Its known side is the stabilized cubic invariant,
        transferred by CJR I (1.10), while its graph side uses CJR III
        (8.21) and therefore has no ordinary contact descendants.
        """
        log_probe = self.stabilization_comparison.log_probe(
            stabilized_probe
        )
        return self.relation(
            log_probe, t_power=t_power,
            require_complete=require_complete,
        )

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
                            include_unit_relatives=False,
                            max_unit_insertions=1,
                            include_chow_relations=False,
                            max_chow_unit_insertions=1,
                            chow_primary_only=False,
                            effective_basis_only=False):
        factory = ProbeFactory(self.rings.base_field)
        if effective_basis_only and (
                not include_chow_relations or not chow_primary_only
                or include_unit_relatives):
            raise ValueError(
                "effective-basis-only probes require primary Chow relations "
                "alongside log stationary rows, without unit relatives"
            )
        # CJR Theorem 10.1 uses log-domain descendants on the zero side.  The
        # stabilization-boundary comparison supplies their known compact
        # values, while (8.21) ensures that no ordinary contact descendant is
        # introduced at infinity.  The broader diagnostic path continues to
        # localize the stabilized class directly.
        probes = list(factory.stationary_candidates(
            genus, ambient_degree, max_markings=max_markings,
            psi_convention=("log" if effective_basis_only else "stabilized"),
        ))
        if include_unit_relatives:
            probes.extend(factory.mixed_unit_candidates(
                genus, ambient_degree,
                max_point_markings=max_markings,
                max_unit_insertions=max_unit_insertions,
            ))
        if include_chow_relations:
            probes.extend(factory.chow_candidates(
                genus, ambient_degree,
                max_markings=max_markings,
                max_unit_insertions=max_chow_unit_insertions,
                primary_only=chow_primary_only,
            ))
        relations = []
        for probe in probes:
            for t_power in self.chow_backend.scalar_powers(
                    probe, t_powers):
                relation = self.relation(
                    probe, t_power,
                    require_complete=require_complete,
                )
                if effective_basis_only and any(
                        factor.contact_psi
                        for _, factors in relation.terms
                        for factor in factors):
                    raise UnsupportedGeometryError(
                        "an effective-basis probe emitted an ordinary "
                        "contact cotangent descendant"
                    )
                relations.append(relation)
        return tuple(relations)

    def relation_report(self, probe, t_power=0, require_complete=False):
        relation = self.relation(
            probe, t_power, require_complete=require_complete
        )
        report = {
            "probe": str(probe),
            "t_power": int(t_power),
            "chow_degree": self.chow_backend.degree(
                probe, t_power
            ).to_record(),
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
        if probe.has_descendants and probe.psi_convention == "log":
            report["stabilization_boundary_comparison"] = (
                self.stabilization_comparison.comparison_record(probe)
            )
        return report


def genus_two_profile_split_rank_witness(provider=None):
    r"""The old witness must be rebuilt in the enlarged descendant basis.

    Stabilization pullback is now compiled correctly, but it produces ordinary
    cotangent powers at infinity contact markings.  The former two-column
    witness omitted those columns, so quoting its determinant would still be
    invalid.  Use ``candidate_relations`` and exact rank selection on the full
    emitted support instead.
    """
    raise UnsupportedGeometryError(
        "the former genus-two profile-split matrix omitted contact-descendant "
        "infinity vertices; rebuild the witness on the full emitted support"
    )


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
        try:
            report = genus_two_profile_split_report()
        except UnsupportedGeometryError as error:
            parser.error(str(error))
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
