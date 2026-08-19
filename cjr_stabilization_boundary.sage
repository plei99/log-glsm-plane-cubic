r"""Virtual stabilization-boundary comparison for compact cubic probes.

There are two genuinely different cotangent-line classes in the localization
calculation:

* ``bar_psi`` on the stabilized hypersurface stable-map space; and
* ``psi`` on the coarse universal log-GLSM curve.

They are not equal as classes: their difference is supported where
stabilization contracts a rational component.  It is therefore incorrect to
feed a Bloch--Okounkov descendant value directly to a graph formula using the
log-domain class merely by renaming the class.

For compact zero-sector markings, however, the virtual comparison used in
CJR I, equation (1.10) and Theorem 1.4, includes precisely these boundary
terms and gives equality of the resulting *numbers*.  In the plane-cubic
normalization it says

    integral_[U(P,D)]^red product ev_i^*(alpha_i) psi_i^a_i
      = (-1)^(1-g+3D)
        integral_[M(E,D/3)]^vir product ev_i^*(alpha_i) bar_psi_i^a_i.

This module packages that numerical comparison.  The returned graph probe
continues to use ``psi_convention='log'`` and is localized with CJR (8.21),
so its infinity factors involve evaluations and ``psi_min`` only.  The known
compact side is evaluated on a separate stabilized copy of the probe.  No
class-level identification of the two psi classes is made.
"""

load("elliptic_probe_values.sage")


class StabilizationBoundaryComparison(SageObject):
    r"""Transfer known stabilized descendants to CJR log-domain probes.

    The comparison is restricted to the compact plane-cubic sector encoded
    by :class:`ProbeSpec`: connected maps, zero-sector markings, and ambient
    insertions ``1`` or ``pt``.  These are exactly the hypotheses used by the
    current graph enumerator.
    """

    def __init__(self, stabilized_backend=None):
        self.stabilized_backend = (
            stabilized_backend or EllipticProbeValueBackend()
        )

    @staticmethod
    def stabilized_probe(log_probe):
        """Return the stabilized counterpart without identifying classes."""
        if not isinstance(log_probe, ProbeSpec):
            raise TypeError("the stabilization comparison needs a ProbeSpec")
        if log_probe.psi_convention != "log":
            raise ValueError("the input must use the log-domain psi convention")
        return ProbeSpec(
            log_probe.genus,
            log_probe.ambient_degree,
            log_probe.insertions,
            connected=log_probe.connected,
            label=log_probe.label,
            psi_convention="stabilized",
        )

    @staticmethod
    def log_probe(stabilized_probe):
        """Return the CJR graph probe paired with a stabilized known side."""
        if not isinstance(stabilized_probe, ProbeSpec):
            raise TypeError("the stabilization comparison needs a ProbeSpec")
        if stabilized_probe.psi_convention != "stabilized":
            raise ValueError("the input must use the stabilized psi convention")
        return ProbeSpec(
            stabilized_probe.genus,
            stabilized_probe.ambient_degree,
            stabilized_probe.insertions,
            connected=stabilized_probe.connected,
            label=stabilized_probe.label,
            psi_convention="log",
        )

    def reduced_value(self, log_probe):
        r"""Return the boundary-corrected reduced log-GLSM descendant.

        This method deliberately calls ``reduced_stabilized_value`` on a
        distinct probe.  Calling ``value(log_probe)`` would bypass the
        convention guard and incorrectly assert equality of line classes.
        """
        stabilized = self.stabilized_probe(log_probe)
        method = getattr(
            self.stabilized_backend, "reduced_stabilized_value", None
        )
        if method is None:
            raise UnsupportedGeometryError(
                "the known-value backend does not implement stabilized "
                "plane-cubic descendants"
            )
        return method(stabilized)

    def comparison_record(self, log_probe):
        """Return auditable metadata for reports and regression tests."""
        stabilized = self.stabilized_probe(log_probe)
        return {
            "kind": "CJR-I virtual stabilization-boundary comparison",
            "log_probe": log_probe.to_record(),
            "stabilized_probe": stabilized.to_record(),
            "reduced_value": str(self.reduced_value(log_probe)),
        }
