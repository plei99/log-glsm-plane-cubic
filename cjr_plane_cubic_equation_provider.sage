"""Compatibility entry point; implementation lives in log_glsm."""

from log_glsm._compat import export_legacy as _export_legacy
_export_legacy('cjr_plane_cubic_equation_provider', globals())
del _export_legacy

if __name__ == "__main__":
    import os as _os
    import sys as _sys
    if _os.path.basename(_sys.argv[0]) in ('cjr_plane_cubic_equation_provider.sage', 'cjr_plane_cubic_equation_provider.sage.py'):
        from log_glsm.cli.cjr_plane_cubic_equation_provider import _main as _run
        raise SystemExit(_run())
