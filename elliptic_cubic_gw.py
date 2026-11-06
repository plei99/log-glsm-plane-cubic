"""Compatibility entry point; implementation lives in log_glsm."""

from log_glsm._compat import export_legacy as _export_legacy
_export_legacy('elliptic_cubic_gw', globals())
del _export_legacy

if __name__ == "__main__":
    import os as _os
    import sys as _sys
    if _os.path.basename(_sys.argv[0]) in ('elliptic_cubic_gw.py', 'elliptic_cubic_gw.py.py'):
        from log_glsm.cli.elliptic_cubic_gw import main as _run
        raise SystemExit(_run())
