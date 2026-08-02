r"""Known connected elliptic-curve values for ``ProbeSpec`` objects.

Stationary invariants are read from ``bo_coefficient.sage`` through its
quasimodular-form wrapper.  Unit insertions with psi power zero or one are
reduced by the string and dilaton equations.  Unsupported nonstationary
probes are rejected explicitly.
"""

load("log_glsm_conventions.sage")
load("log_glsm_infinity_vertices.sage")


class EllipticProbeValueBackend(SageObject):
    """Exact known-GW right-hand side for connected compact probes."""

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
        if not probe.is_dimension_zero():
            value = QQ.zero()
        elif all(item.kind == "point" for item in probe.insertions):
            degree = probe.intrinsic_degree
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
                    ))
                return total
            if insertion.psi_power == 1:
                # Dilaton equation, with n equal to the remaining markings.
                factor = 2 * probe.genus - 2 + len(remaining)
                return factor * self.value(ProbeSpec(
                    probe.genus, probe.ambient_degree, remaining,
                    connected=True, label="dilaton reduction",
                ))
            raise UnsupportedGeometryError(
                "unit descendants psi^%s require a nonstationary backend"
                % insertion.psi_power
            )
        raise UnsupportedGeometryError("probe is not reducible to stationary theory")

    def reduced_log_glsm_value(self, probe):
        return plane_cubic_reduced_sign(
            probe.genus, probe.ambient_degree
        ) * self.value(probe)


def known_elliptic_invariant(probe):
    return EllipticProbeValueBackend().value(probe)

