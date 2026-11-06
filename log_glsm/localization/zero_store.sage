r"""localization / zero_store."""

from sage.all import *

import json
import os
import sqlite3
from log_glsm.core.conventions import (
    NONEQUIVARIANT_BASE_PATH,
)

__all__ = ['PersistentZeroVertexStore']


class PersistentZeroVertexStore(SageObject):
    """Concurrent SQLite store for normalized exact zero-vertex values."""

    FORMAT = "log-glsm-o3-zero-vertices-sqlite"
    VERSION = 1

    def __init__(self, path, include_twist, lift_strategy,
                 base_weight_specialization=None):
        self.path = os.path.abspath(path)
        directory = os.path.dirname(self.path)
        if directory and not os.path.isdir(directory):
            os.makedirs(directory)
        self.connection = sqlite3.connect(self.path, timeout=float(60))
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("PRAGMA busy_timeout=60000")
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS metadata ("
            "key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS zero_vertices ("
            "request_key TEXT PRIMARY KEY, request_json TEXT NOT NULL, "
            "value TEXT NOT NULL)"
        )
        expected = {
            "format": self.FORMAT,
            "version": str(int(self.VERSION)),
            "include_twist": str(bool(include_twist)),
            "lift_strategy": str(lift_strategy),
        }
        if base_weight_specialization is not None:
            if isinstance(base_weight_specialization, str):
                expected["base_weight_specialization"] = json.dumps(
                    base_weight_specialization
                )
                if base_weight_specialization == "nonequivariant":
                    expected["nonequivariant_path"] = json.dumps(
                        NONEQUIVARIANT_BASE_PATH
                    )
            else:
                expected["base_weight_specialization"] = json.dumps([
                    str(QQ(value)) for value in base_weight_specialization
                ])
        metadata = dict(self.connection.execute(
            "SELECT key, value FROM metadata"
        ).fetchall())
        if metadata:
            if any(metadata.get(key) != value
                   for key, value in expected.items()):
                raise ValueError("incompatible O(3) zero-vertex SQLite cache")
        else:
            self.connection.executemany(
                "INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)",
                tuple(expected.items()),
            )
            self.connection.commit()

    @staticmethod
    def request_key(request):
        return json.dumps(
            request.localization_record(),
            sort_keys=True, separators=(",", ":")
        )

    def get(self, request):
        row = self.connection.execute(
            "SELECT value FROM zero_vertices WHERE request_key = ?",
            (self.request_key(request),),
        ).fetchone()
        return None if row is None else row[0]

    def put(self, request, value):
        self.connection.execute(
            "INSERT OR IGNORE INTO zero_vertices"
            "(request_key, request_json, value) VALUES (?, ?, ?)",
            (
                self.request_key(request),
                json.dumps(request.to_record(), sort_keys=True),
                str(value),
            ),
        )
        self.connection.commit()

    def __len__(self):
        return int(self.connection.execute(
            "SELECT COUNT(*) FROM zero_vertices"
        ).fetchone()[0])

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
