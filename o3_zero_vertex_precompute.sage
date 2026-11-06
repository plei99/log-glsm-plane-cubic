"""Compatibility entry point; implementation lives in log_glsm."""

from log_glsm._compat import export_legacy as _export_legacy
_export_legacy('o3_zero_vertex_precompute', globals())
del _export_legacy

if __name__ == "__main__":
    import os as _os
    import sys as _sys
    if _os.path.basename(_sys.argv[0]) in ('o3_zero_vertex_precompute.sage', 'o3_zero_vertex_precompute.sage.py'):
        from log_glsm.cli.o3_zero_vertex_precompute import _main as _run
        raise SystemExit(_run())
