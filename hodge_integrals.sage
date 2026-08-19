r"""Exact intersection-number backends needed by twisted fixed loci.

Pure psi intersections are implemented in every genus by the DVV recursion.
The lambda-g theorem and the top-Hodge-triple formula supply two closed
families, while Mumford relations simplify repeated Hodge factors through
genus three. ``AdmcyclesHodgeIntegralBackend`` extends this to arbitrary
remaining products of lambda and psi classes in the tautological ring.
"""

load("log_glsm_conventions.sage")

import os
import sys
import json
import sqlite3


def _load_vendored_admcycles():
    """Import admcycles, preferring the workspace-local vendored package."""
    candidates = (
        os.path.join(os.getcwd(), "vendor"),
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "vendor"),
    )
    for candidate in candidates:
        if os.path.isdir(os.path.join(candidate, "admcycles")) \
                and candidate not in sys.path:
            sys.path.insert(0, candidate)
    try:
        import admcycles
        return admcycles
    except ImportError as error:
        raise UnsupportedGeometryError(
            "admcycles is required for general Hodge monomials; install it "
            "or keep the vendored package in ./vendor"
        ) from error


def _odd_double_factorial(number):
    number = ZZ(number)
    if number <= -1:
        return ZZ.one()
    return prod(ZZ(k) for k in range(1, number + 1, 2))


_top_hodge_descendant_cache = {}


def _top_hodge_descendant_multiplier(genus, degrees):
    r"""Universal string/dilaton multiplier for a pulled-back top class.

    If ``alpha`` has top degree on ``Mbar_g``, this returns the multiplier of
    ``integral(alpha)`` in ``integral(alpha prod psi_i^degrees[i])``.  Only
    the dimension-matching case ``sum(degrees) == len(degrees)`` survives.
    """
    genus = ZZ(genus)
    degrees = tuple(sorted((ZZ(value) for value in degrees), reverse=True))
    key = genus, degrees
    if key in _top_hodge_descendant_cache:
        return _top_hodge_descendant_cache[key]
    if sum(degrees) != len(degrees):
        value = ZZ.zero()
    elif not degrees:
        value = ZZ.one()
    elif degrees[-1] == 0:
        rest = degrees[:-1]
        value = sum(
            _top_hodge_descendant_multiplier(
                genus, rest[:index] + (degree - 1,) + rest[index + 1:]
            )
            for index, degree in enumerate(rest) if degree > 0
        )
    else:
        # Positive degrees summing to their count are all one.  Repeatedly
        # apply the dilaton equation down to the unmarked integral.
        value = prod(
            2 * genus - 2 + index for index in range(len(degrees))
        )
    _top_hodge_descendant_cache[key] = ZZ(value)
    return _top_hodge_descendant_cache[key]


class HodgeIntegralRequest(SageObject):
    """A hashable unsupported/supported Hodge integral specification."""

    def __init__(self, genus, psi_powers, lambda_indices=()):
        self.genus = ZZ(genus)
        # The markings are interchangeable after integration.  Canonicalizing
        # here turns the many permutations produced by fixed-graph
        # localization into one admcycles request and one cache entry.
        self.psi_powers = tuple(sorted(
            (ZZ(value) for value in psi_powers), reverse=True
        ))
        # lambda_0 is the unit, so it should not create a distinct cache key.
        self.lambda_indices = tuple(sorted(
            ZZ(value) for value in lambda_indices if ZZ(value) != 0
        ))
        if self.genus < 0 or any(value < 0 for value in self.psi_powers):
            raise ValueError("genus and psi powers must be nonnegative")
        if any(value < 0 or value > self.genus for value in self.lambda_indices):
            raise ValueError("lambda indices must lie between zero and genus")

    def signature(self):
        return self.genus, self.psi_powers, self.lambda_indices

    def low_genus_normal_form(self):
        r"""Reduce Hodge monomials using Mumford's relation through genus 3.

        The identity ``lambda(t) lambda(-t) = 1`` gives

        - ``lambda_1^2 = 0`` in genus 1;
        - ``lambda_1^2 = 2 lambda_2`` and ``lambda_2^2 = 0`` in genus 2;
        - ``lambda_1^2 = 2 lambda_2``,
          ``lambda_2^2 = 2 lambda_1 lambda_3``, and
          ``lambda_3^2 = 0`` in genus 3.

        Returning the scalar separately lets the caller reuse the canonical
        request in both the memory and SQLite caches.  Higher genus requests
        are deliberately left unchanged rather than applying an incomplete
        rewriting system.
        """
        if self.genus > 3 or not self.lambda_indices:
            return QQ.one(), self

        indices = list(self.lambda_indices)
        factor = QQ.one()
        while True:
            if self.genus == 1:
                if indices.count(ZZ.one()) >= 2:
                    return QQ.zero(), self
                break

            if indices.count(ZZ.one()) >= 2:
                indices.remove(ZZ.one())
                indices.remove(ZZ.one())
                indices.append(ZZ(2))
                indices.sort()
                factor *= 2
                continue

            if self.genus == 2:
                if indices.count(ZZ(2)) >= 2:
                    return QQ.zero(), self
                break

            # genus 3
            if indices.count(ZZ(2)) >= 2:
                indices.remove(ZZ(2))
                indices.remove(ZZ(2))
                indices.extend((ZZ.one(), ZZ(3)))
                indices.sort()
                factor *= 2
                continue
            if indices.count(ZZ(3)) >= 2:
                return QQ.zero(), self
            break

        normalized = HodgeIntegralRequest(
            self.genus, self.psi_powers, tuple(indices)
        )
        return factor, normalized

    def to_record(self):
        return {
            "genus": int(self.genus),
            "psi_powers": [int(value) for value in self.psi_powers],
            "lambda_indices": [int(value) for value in self.lambda_indices],
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["genus"],
            record.get("psi_powers", ()),
            record.get("lambda_indices", ()),
        )

    def __hash__(self):
        return hash(self.signature())

    def __eq__(self, other):
        return isinstance(other, HodgeIntegralRequest) and self.signature() == other.signature()

    def _repr_(self):
        return "HodgeIntegral(g=%s,psi=%s,lambda=%s)" % self.signature()


class PersistentHodgeIntegralStore(SageObject):
    r"""Concurrent SQLite store for exact Hodge integral values.

    SQLite avoids rewriting a large JSON object after every admcycles call,
    and WAL mode lets independent precompute processes share completed
    values safely.
    """

    FORMAT = "log-glsm-hodge-integrals"
    VERSION = 1

    def __init__(self, path):
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
            "CREATE TABLE IF NOT EXISTS integrals ("
            "request_key TEXT PRIMARY KEY, request_json TEXT NOT NULL, "
            "value TEXT NOT NULL)"
        )
        metadata = dict(self.connection.execute(
            "SELECT key, value FROM metadata"
        ).fetchall())
        expected = {
            "format": self.FORMAT,
            "version": str(int(self.VERSION)),
        }
        if metadata:
            if any(metadata.get(key) != value
                   for key, value in expected.items()):
                raise ValueError("unsupported Hodge-integral SQLite cache")
        else:
            self.connection.executemany(
                "INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)",
                tuple(expected.items()),
            )
            self.connection.commit()

    @staticmethod
    def request_key(request):
        return json.dumps(
            request.to_record(), sort_keys=True, separators=(",", ":")
        )

    def get(self, request):
        row = self.connection.execute(
            "SELECT value FROM integrals WHERE request_key = ?",
            (self.request_key(request),),
        ).fetchone()
        return None if row is None else QQ(row[0])

    def put(self, request, value):
        record = json.dumps(request.to_record(), sort_keys=True)
        self.connection.execute(
            "INSERT OR IGNORE INTO integrals"
            "(request_key, request_json, value) VALUES (?, ?, ?)",
            (self.request_key(request), record, str(QQ(value))),
        )
        # A previous worker may have been interrupted while evaluating this
        # same request.  Completion by any worker makes every matching
        # diagnostic marker stale.
        self.connection.execute(
            "DELETE FROM metadata WHERE key LIKE 'in_progress:%' "
            "AND value = ?", (record,),
        )
        self.connection.commit()

    def begin_evaluation(self, request):
        r"""Persist the exact request currently entering a slow backend.

        A general admcycles multiplication is intentionally atomic and can
        take much longer than an orchestration time budget.  Recording the
        request before constructing the tautological class makes an
        interrupted run diagnosable and lets a later run specialize precisely
        the family that caused the delay.  The process-specific key is safe
        when several cache-producing workers share the same WAL database.
        """
        key = "in_progress:%s" % os.getpid()
        record = json.dumps(request.to_record(), sort_keys=True)
        self.connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            (key, record),
        )
        self.connection.commit()
        return key

    def finish_evaluation(self, key):
        """Clear a completed in-progress diagnostic record."""
        self.connection.execute(
            "DELETE FROM metadata WHERE key = ?", (str(key),)
        )
        self.connection.commit()

    def __len__(self):
        return int(self.connection.execute(
            "SELECT COUNT(*) FROM integrals"
        ).fetchone()[0])

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None


class PsiIntersectionBackend(SageObject):
    """Witten--Kontsevich intersections by exact DVV recursion."""

    def __init__(self):
        self._cache = {}

    def integral(self, genus, degrees):
        genus = ZZ(genus)
        degrees = tuple(sorted((ZZ(value) for value in degrees), reverse=True))
        key = genus, degrees
        if key in self._cache:
            return self._cache[key]
        value = self._integral_uncached(genus, degrees)
        self._cache[key] = QQ(value)
        return self._cache[key]

    def _integral_uncached(self, genus, degrees):
        marking_count = len(degrees)
        if genus < 0 or any(value < 0 for value in degrees):
            return QQ.zero()
        if 2 * genus - 2 + marking_count <= 0:
            return QQ.zero()
        if sum(degrees) != 3 * genus - 3 + marking_count:
            return QQ.zero()
        if genus == 0 and degrees == (0, 0, 0):
            return QQ.one()
        if genus == 1 and degrees == (1,):
            return QQ(1) / 24

        first = degrees[0]
        rest = degrees[1:]
        if first == 0:
            # String equation; this case is mostly the genus-zero base.
            return sum(
                self.integral(genus, rest[:j] + (rest[j] - 1,) + rest[j + 1:])
                for j in range(len(rest)) if rest[j] > 0
            )

        denominator = _odd_double_factorial(2 * first + 1)
        merge = QQ.zero()
        for j, degree in enumerate(rest):
            coefficient = QQ(
                _odd_double_factorial(2 * first + 2 * degree - 1)
            ) / _odd_double_factorial(2 * degree - 1)
            merged_degrees = rest[:j] + (first + degree - 1,) + rest[j + 1:]
            merge += coefficient * self.integral(genus, merged_degrees)

        boundary = QQ.zero()
        if first >= 2:
            for left_degree in range(first - 1):
                right_degree = first - 2 - left_degree
                coefficient = (
                    _odd_double_factorial(2 * left_degree + 1)
                    * _odd_double_factorial(2 * right_degree + 1)
                )
                nonseparating = self.integral(
                    genus - 1, (left_degree, right_degree) + rest
                )
                separating = QQ.zero()
                rest_count = len(rest)
                for mask in range(1 << rest_count):
                    left_rest = tuple(
                        rest[j] for j in range(rest_count) if mask & (1 << j)
                    )
                    right_rest = tuple(
                        rest[j] for j in range(rest_count) if not mask & (1 << j)
                    )
                    for left_genus in range(genus + 1):
                        separating += self.integral(
                            left_genus, (left_degree,) + left_rest
                        ) * self.integral(
                            genus - left_genus, (right_degree,) + right_rest
                        )
                boundary += coefficient * (nonseparating + separating)

        return (merge + boundary / 2) / denominator


class HodgeIntegralBackend(SageObject):
    """Closed psi, lambda-g, and top-Hodge-triple intersection backend."""

    def __init__(self, cache_path=None):
        self.psi = PsiIntersectionBackend()
        self._cache = {}
        self.store = (
            None if cache_path is None
            else PersistentHodgeIntegralStore(cache_path)
        )
        self._persistent_hits = ZZ.zero()
        self._persistent_writes = ZZ.zero()

    def _lookup(self, request):
        if request in self._cache:
            return self._cache[request]
        if self.store is not None:
            value = self.store.get(request)
            if value is not None:
                self._persistent_hits += 1
                self._cache[request] = QQ(value)
                return self._cache[request]
        return None

    def _remember(self, request, value):
        value = QQ(value)
        self._cache[request] = value
        if self.store is not None:
            self.store.put(request, value)
            self._persistent_writes += 1
        return value

    def cache_info(self):
        return {
            "memory_entries": len(self._cache),
            "persistent_entries": (
                int(0) if self.store is None else int(len(self.store))
            ),
            "persistent_hits": int(self._persistent_hits),
            "persistent_writes": int(self._persistent_writes),
            "path": None if self.store is None else self.store.path,
        }

    @staticmethod
    def lambda_g_constant(genus):
        genus = ZZ(genus)
        if genus < 0:
            raise ValueError("genus must be nonnegative")
        if genus == 0:
            return QQ.one()
        return QQ(
            (2 ** (2 * genus - 1) - 1) * abs(bernoulli(2 * genus))
        ) / (2 ** (2 * genus - 1) * factorial(2 * genus))

    @staticmethod
    def top_hodge_triple_constant(genus):
        r"""Return ``integral lambda_g lambda_(g-1) lambda_(g-2)``.

        This is the Faber--Pandharipande closed formula.  In genus two the
        ``lambda_0`` factor is the unit and is omitted from request keys.
        """
        genus = ZZ(genus)
        if genus < 2:
            raise ValueError("the top Hodge triple requires genus at least 2")
        return QQ.one() / (2 * factorial(2 * genus - 2)) \
            * QQ(abs(bernoulli(2 * genus - 2))) / (2 * genus - 2) \
            * QQ(abs(bernoulli(2 * genus))) / (2 * genus)

    @staticmethod
    def top_hodge_triple_indices(genus):
        genus = ZZ(genus)
        if genus < 2:
            return ()
        return tuple(ZZ(index) for index in range(
            max(1, int(genus) - 2), int(genus) + 1
        ))

    def integral(self, request_or_genus, psi_powers=None, lambda_indices=()):
        if isinstance(request_or_genus, HodgeIntegralRequest):
            request = request_or_genus
        else:
            request = HodgeIntegralRequest(
                request_or_genus, psi_powers or (), lambda_indices
            )
        cached = self._lookup(request)
        if cached is not None:
            return cached

        factor, normalized = request.low_genus_normal_form()
        if factor == 0:
            return self._remember(request, QQ.zero())
        if normalized != request:
            return self._remember(request, factor * self.integral(normalized))

        nonzero_lambda = tuple(index for index in request.lambda_indices if index)
        if not nonzero_lambda:
            value = self.psi.integral(request.genus, request.psi_powers)
        elif nonzero_lambda == (request.genus,):
            expected = 2 * request.genus - 3 + len(request.psi_powers)
            if sum(request.psi_powers) != expected:
                value = QQ.zero()
            else:
                multinomial = factorial(expected) / prod(
                    factorial(value) for value in request.psi_powers
                )
                value = QQ(multinomial) * self.lambda_g_constant(request.genus)
        elif nonzero_lambda == self.top_hodge_triple_indices(request.genus):
            value = QQ(_top_hodge_descendant_multiplier(
                request.genus, request.psi_powers
            )) * self.top_hodge_triple_constant(request.genus)
        else:
            raise UnsupportedGeometryError(
                "general Hodge monomial %s is not implemented; "
                "pure psi and one lambda_g factor are exact" % (request,)
            )
        return self._remember(request, value)


class AdmcyclesHodgeIntegralBackend(HodgeIntegralBackend):
    r"""All tautological Hodge/psi monomials via admcycles ``evaluate()``."""

    def __init__(self, cache_path=None):
        super().__init__(cache_path=cache_path)
        self.admcycles = _load_vendored_admcycles()
        self._lambda_product_cache = {}

    def cache_info(self):
        info = super().cache_info()
        info["lambda_product_entries"] = len(self._lambda_product_cache)
        return info

    def integral(self, request_or_genus, psi_powers=None, lambda_indices=()):
        if isinstance(request_or_genus, HodgeIntegralRequest):
            request = request_or_genus
        else:
            request = HodgeIntegralRequest(
                request_or_genus, psi_powers or (), lambda_indices
            )
        cached = self._lookup(request)
        if cached is not None:
            return cached

        factor, normalized = request.low_genus_normal_form()
        if factor == 0:
            return self._remember(request, QQ.zero())
        if normalized != request:
            return self._remember(request, factor * self.integral(normalized))

        # General Hodge monomials are pulled back by a forgetful morphism, so
        # the ordinary string and dilaton equations reduce markings before we
        # construct any tautological class.  This is decisive in genus three:
        # lambda_2 lambda_3 psi_1^4 on Mbar_3,3 reduces to
        # lambda_2 lambda_3 psi_1^2 on Mbar_3,1 instead of triggering a huge
        # boundary-stratum multiplication in admcycles.
        degrees = request.psi_powers
        forgetful_target_is_stable = bool(degrees) and (
            2 * request.genus - 2 + len(degrees) - 1 > 0
        )
        if forgetful_target_is_stable and degrees[-1] == 0:
            retained = degrees[:-1]
            value = QQ.zero()
            for marking, degree in enumerate(retained):
                if degree:
                    reduced = list(retained)
                    reduced[marking] -= 1
                    value += self.integral(HodgeIntegralRequest(
                        request.genus, tuple(reduced), request.lambda_indices
                    ))
            return self._remember(request, value)
        if forgetful_target_is_stable and degrees[-1] == 1:
            retained = degrees[:-1]
            multiplier = 2 * request.genus - 2 + len(retained)
            return self._remember(
                request,
                multiplier * self.integral(HodgeIntegralRequest(
                    request.genus, retained, request.lambda_indices
                )),
            )

        # Do not pay for tautological-ring construction when the exact DVV or
        # lambda-g formula already covers the request.  Genus-zero fixed loci
        # in particular contain only pure psi integrals.
        nonzero_lambda = tuple(
            index for index in request.lambda_indices if index
        )
        if not nonzero_lambda \
                or nonzero_lambda == (request.genus,) \
                or nonzero_lambda == self.top_hodge_triple_indices(
                    request.genus
                ):
            return super().integral(request)

        genus = request.genus
        markings = len(request.psi_powers)
        if 2 * genus - 2 + markings <= 0:
            value = QQ.zero()
        else:
            codimension = sum(request.psi_powers) + sum(request.lambda_indices)
            dimension = 3 * genus - 3 + markings
            if codimension != dimension:
                value = QQ.zero()
            else:
                in_progress_key = None
                if self.store is not None:
                    in_progress_key = self.store.begin_evaluation(request)
                lambda_key = (
                    ZZ(genus), ZZ(markings), tuple(request.lambda_indices)
                )
                if lambda_key not in self._lambda_product_cache:
                    lambda_product = self.admcycles.fundclass(genus, markings)
                    for index in request.lambda_indices:
                        if index:
                            lambda_product = (
                                lambda_product * self.admcycles.lambdaclass(
                                    int(index), genus, markings
                                )
                            )
                    self._lambda_product_cache[lambda_key] = lambda_product

                # Tautological classes are multiplied functionally here; the
                # cached lambda product itself is never modified.
                tautological_class = self._lambda_product_cache[lambda_key]
                for marking, power in enumerate(request.psi_powers, start=1):
                    if power:
                        tautological_class = (
                            tautological_class * self.admcycles.psiclass(
                                marking, genus, markings
                            ) ** int(power)
                        )
                value = QQ(tautological_class.evaluate())
                if in_progress_key is not None:
                    self.store.finish_evaluation(in_progress_key)
        return self._remember(request, value)
