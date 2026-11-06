"""Import this package's .sage sources through Sage's own preparser.

Only log_glsm modules are handled. Python's import cache owns module identity,
so loading a dependency cannot redefine another module's classes or caches.
"""

import importlib.abc
import importlib.util
from pathlib import Path
import sys


class _SageLoader(importlib.abc.Loader):
    def __init__(self, source):
        self.source = source

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        from sage.repl.preparse import preparse_file

        source = self.source.read_text(encoding="utf-8")
        # Inline Sage numeric constructors instead of hoisting constants before
        # the source's imports. This retains exact rational division and powers.
        code = compile(
            preparse_file(source, numeric_literals=False), str(self.source), "exec"
        )
        exec(code, module.__dict__)


class _SageFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if not fullname.startswith("log_glsm.") or path is None:
            return None
        name = fullname.rsplit(".", 1)[-1]
        for directory in path:
            source = Path(directory) / (name + ".sage")
            if source.is_file():
                return importlib.util.spec_from_file_location(
                    fullname, source, loader=_SageLoader(source)
                )
        return None


def install_importer():
    if not any(isinstance(finder, _SageFinder) for finder in sys.meta_path):
        sys.meta_path.append(_SageFinder())
