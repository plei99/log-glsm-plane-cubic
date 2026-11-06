"""Compatibility entry point; implementation lives in log_glsm."""

from log_glsm._compat import export_legacy as _export_legacy
_export_legacy('cjr_bipartite_graphs', globals())
del _export_legacy

if __name__ == "__main__":
    import os as _os
    import sys as _sys
    if _os.path.basename(_sys.argv[0]) in ('cjr_bipartite_graphs.sage', 'cjr_bipartite_graphs.sage.py'):
        from log_glsm.cli.cjr_bipartite_graphs import _command_line as _run
        raise SystemExit(_run(_sys.argv))
