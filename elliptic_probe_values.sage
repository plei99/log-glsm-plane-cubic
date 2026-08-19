r"""Known connected stable-map values for ``ProbeSpec`` objects.

Stationary invariants are read from ``bo_coefficient.sage`` through its
quasimodular-form wrapper.  Unit insertions with psi power zero or one are
reduced by the string and dilaton equations.  Unsupported nonstationary
probes are rejected explicitly.

These are descendants after stabilization.  They are not the descendants on
the coarse log-GLSM universal curve appearing in CJR III, equation (8.21),
when a marked zero tail is contracted by ``st``.
"""

load("log_glsm_conventions.sage")
load("log_glsm_infinity_vertices.sage")


class EllipticProbeValueBackend(SageObject):
    """Exact known-GW values with stabilized stable-map psi classes."""

    def __init__(self):
        self._series_cache = {}
        self._value_cache = {}

    def stationary_series(self, descendants, max_intrinsic_degree):
        descendants = tuple(ZZ(value) for value in descendants)
        max_intrinsic_degree = ZZ(max_intrinsic_degree)
        key = descendants, max_intrinsic_degree
        if key not in self._series_cache:
            self._series_cache[key] = connected_stationary_qseries(
                descendants, max_intrinsic_degree
            )
        return self._series_cache[key]

    def value(self, probe):
        if not isinstance(probe, ProbeSpec):
            raise TypeError("probe must be a ProbeSpec")
        if probe in self._value_cache:
            return self._value_cache[probe]
        if not probe.connected:
            raise UnsupportedGeometryError(
                "the current known-value wrapper accepts connected probes"
            )
        if not probe.is_stable_map_type():
            raise UnsupportedGeometryError(
                "the connected stable-map data (g,n,d)=(%s,%s,%s) are "
                "unstable"
                % (probe.genus, probe.marking_count, probe.ambient_degree)
            )
        if probe.has_descendants and probe.psi_convention != "stabilized":
            raise UnsupportedGeometryError(
                "Bloch--Okounkov computes stabilized stable-map descendants, "
                "not cotangent lines on the pre-stabilization log-GLSM curve"
            )
        if probe.ambient_degree % 3:
            # A curve contained in a plane cubic has ambient class 3r.  The
            # ambient log-GLSM localization problem is still defined in every
            # P2 degree, but equation (9.7) has an empty hypersurface side in
            # the other degrees.  Keeping these exact zero probes is useful:
            # they give substantially smaller homogeneous equations for
            # positive-degree infinity vertices.
            value = QQ.zero()
        elif not probe.is_dimension_zero():
            value = QQ.zero()
        elif all(item.kind == "point" for item in probe.insertions):
            degree = probe.intrinsic_degree
            if degree == 0 and probe.marking_count > 1:
                # A degree-zero map to the elliptic curve is constant, so all
                # evaluation maps factor through the same copy of E.  Two
                # stationary insertions contribute [pt]^2=0 in H*(E).  This
                # avoids constructing enormous Bloch--Okounkov expressions
                # for high-valence probes whose constant term is known to
                # vanish geometrically.
                value = QQ.zero()
            else:
                _, series = self.stationary_series(
                    probe.stationary_descendants(), degree
                )
                value = QQ(series[degree])
        else:
            value = self._reduce_unit_insertion(probe)
        self._value_cache[probe] = QQ(value)
        return self._value_cache[probe]

    def _reduce_unit_insertion(self, probe):
        for index, insertion in enumerate(probe.insertions):
            if insertion.kind != "unit":
                continue
            remaining = probe.insertions[:index] + probe.insertions[index + 1:]
            if insertion.psi_power == 0:
                # String equation.
                total = QQ.zero()
                for other_index, other in enumerate(remaining):
                    if other.psi_power == 0:
                        continue
                    reduced = EllipticInsertion(other.kind, other.psi_power - 1)
                    new_insertions = (
                        remaining[:other_index] + (reduced,)
                        + remaining[other_index + 1:]
                    )
                    total += self.value(ProbeSpec(
                        probe.genus, probe.ambient_degree, new_insertions,
                        connected=True, label="string reduction",
                        psi_convention=probe.psi_convention,
                    ))
                return total
            if insertion.psi_power == 1:
                # Dilaton equation, with n equal to the remaining markings.
                factor = 2 * probe.genus - 2 + len(remaining)
                return factor * self.value(ProbeSpec(
                    probe.genus, probe.ambient_degree, remaining,
                    connected=True, label="dilaton reduction",
                    psi_convention=probe.psi_convention,
                ))
            raise UnsupportedGeometryError(
                "unit descendants psi^%s require a nonstationary backend"
                % insertion.psi_power
            )
        raise UnsupportedGeometryError("probe is not reducible to stationary theory")

    def reduced_log_glsm_value(self, probe):
        r"""Return the reduced log-GLSM side with the requested psi pullback.

        Equation (9.7) identifies ``st_*[U]^red`` with the GW virtual class.
        Primary probes need no convention choice.  For descendants, the graph
        compiler now localizes ``st^*(bar_psi)``: on a contracted marked zero
        tail it records the descendant at the adjacent infinity contact.
        Bloch--Okounkov therefore supplies the correct stabilized known side.
        A genuine log-domain descendant still needs an independent backend.
        """
        if probe.has_descendants and probe.psi_convention != "stabilized":
            raise UnsupportedGeometryError(
                "Bloch--Okounkov does not determine a log-domain descendant"
            )
        return plane_cubic_reduced_sign(
            probe.genus, probe.ambient_degree
        ) * self.value(probe)

    @staticmethod
    def supports_log_glsm_probe(probe):
        """Whether this backend knows the compact side of this graph probe."""
        return not probe.has_descendants or probe.psi_convention == "stabilized"

    def reduced_stabilized_value(self, probe):
        """Apply the reduced-class sign to a stabilized GW invariant."""
        if probe.has_descendants and probe.psi_convention != "stabilized":
            raise UnsupportedGeometryError(
                "the requested descendant does not use stabilized psi"
            )
        return plane_cubic_reduced_sign(
            probe.genus, probe.ambient_degree
        ) * self.value(probe)


def known_elliptic_invariant(probe):
    return EllipticProbeValueBackend().value(probe)
