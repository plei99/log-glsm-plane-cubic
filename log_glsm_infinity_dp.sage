"""Compatibility entry point; implementation lives in log_glsm."""

from log_glsm._compat import export_legacy as _export_legacy
_export_legacy('log_glsm_infinity_dp', globals())
del _export_legacy
