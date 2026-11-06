"""Compatibility for the historical root-level load() and script entry points.

New code imports named objects from their owning modules. This adapter alone
reconstructs the transitive namespace exposed by the old load() chains.
"""

from functools import lru_cache
from importlib import import_module
import json
from pathlib import Path


@lru_cache(maxsize=1)
def _manifest():
    return json.loads(Path(__file__).with_name("_legacy.json").read_text())


@lru_cache(maxsize=None)
def _bindings(name):
    record = _manifest()[name]
    result = {}
    for dependency in record["dependencies"]:
        result.update(_bindings(dependency))
    for symbol, origin in record["bindings"].items():
        kind, module, *tail = origin
        if kind == "module":
            value = getattr(import_module("log_glsm." + module), symbol)
        elif kind == "from":
            value = getattr(import_module(module), tail[0])
        else:
            value = import_module(module if tail[0] else module.split(".")[0])
        result[symbol] = value
    return result


def export_legacy(name, namespace):
    namespace.update(_bindings(name))
