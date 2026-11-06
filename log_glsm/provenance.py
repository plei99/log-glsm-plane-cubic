"""Hash the mathematical implementation for reports and reusable caches."""

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def source_hashes(extra=()):
    """Return repository-relative SHA-256 hashes, including caller-supplied files."""
    paths = {
        p for p in (ROOT / "log_glsm").rglob("*")
        if p.suffix in (".py", ".sage", ".json") and "__pycache__" not in p.parts
    }
    paths.update(ROOT / p for p in extra)
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(paths)
    }


def source_fingerprint():
    """A cache key that changes when any implementation source or its path changes."""
    records = source_hashes()
    payload = "\n".join(name + ":" + value for name, value in records.items())
    return hashlib.sha256(payload.encode()).hexdigest()
