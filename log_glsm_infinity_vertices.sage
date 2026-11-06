"""Compatibility entry point; implementation lives in log_glsm."""

from log_glsm._compat import export_legacy as _export_legacy
_export_legacy('log_glsm_infinity_vertices', globals())
del _export_legacy

if __name__ == "__main__":
    import os as _os
    import sys as _sys
    if _os.path.basename(_sys.argv[0]) in ('log_glsm_infinity_vertices.sage', 'log_glsm_infinity_vertices.sage.py'):
        from log_glsm.cli.log_glsm_infinity_vertices import _main as _run
        raise SystemExit(_run())
