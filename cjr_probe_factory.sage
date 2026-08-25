r"""Exact rank-aware probe selection for infinity-vertex blocks."""

load("log_glsm_conventions.sage")
load("log_glsm_infinity_dp.sage")


class CompiledLocalizationRelation(SageObject):
    r"""A target-independent extracted relation ``known=twisted+terms``."""

    def __init__(self, probe, t_power, coefficient_field, known_gw,
                 twisted_zero_level=0, terms=(), compilation=None,
                 is_complete=None):
        self.probe = probe
        self.t_power = ZZ(t_power)
        self.coefficient_field = coefficient_field
        self.known_gw = coefficient_field(known_gw)
        self.twisted_zero_level = coefficient_field(twisted_zero_level)
        self.terms = tuple(
            (coefficient_field(coefficient), tuple(factors))
            for coefficient, factors in terms
        )
        self.compilation = compilation
        self._is_complete_override = (
            None if is_complete is None else bool(is_complete)
        )

    @property
    def is_complete(self):
        if self._is_complete_override is not None:
            return self._is_complete_override
        return self.compilation is None or self.compilation.is_complete

    def equation_for(self, target):
        return LocalizationEquation(
            target,
            self.known_gw,
            twisted_zero_level=self.twisted_zero_level,
            terms=self.terms,
            probe_label="%s [t^%s]" % (self.probe, self.t_power),
            coefficient_field=self.coefficient_field,
        )

    def to_record(self):
        """Serialize an extracted relation without its heavy graph objects."""
        return {
            "probe": self.probe.to_record(),
            "t_power": int(self.t_power),
            "known_gw": str(self.known_gw),
            "twisted_zero_level": str(self.twisted_zero_level),
            "is_complete": self.is_complete,
            "terms": [
                {
                    "coefficient": str(coefficient),
                    "factors": [factor.to_record() for factor in factors],
                }
                for coefficient, factors in self.terms
            ],
        }

    @classmethod
    def from_record(cls, record, coefficient_field=QQ):
        return cls(
            ProbeSpec.from_record(record["probe"]),
            record["t_power"],
            coefficient_field,
            known_gw=coefficient_field(record["known_gw"]),
            twisted_zero_level=coefficient_field(
                record.get("twisted_zero_level", "0")
            ),
            terms=tuple(
                (
                    coefficient_field(term["coefficient"]),
                    tuple(EffectiveVertex.from_record(factor)
                          for factor in term.get("factors", ())),
                )
                for term in record.get("terms", ())
            ),
            compilation=None,
            is_complete=record.get("is_complete", True),
        )


class ProbeBlockSelection(SageObject):
    """Independent relation rows selected for a same-rank target block."""

    def __init__(self, targets, relations, rows, coefficient_field,
                 rejected=(), rows_independent=False,
                 rejections_certified=True):
        self.targets = tuple(targets)
        self.relations = tuple(relations)
        self.rows = tuple(tuple(coefficient_field(value) for value in row)
                          for row in rows)
        self.coefficient_field = coefficient_field
        self.rejected = tuple(rejected)
        self.matrix = Matrix(coefficient_field, self.rows) if self.rows else Matrix(
            coefficient_field, 0, len(self.targets)
        )
        self._rank_cache = None
        # ``select_relations`` only ever appends rows its incremental basis
        # proved independent, so the exact fraction-field rank recomputation
        # can be skipped for selections built there.
        if rows_independent:
            self._rank_cache = ZZ(len(self.rows))
        # False when a fast selection rejected rows on specialization
        # evidence alone; such a selection must not be quoted as a certified
        # rank deficiency.  Full-rank selections are certified either way,
        # because acceptance always carries an exact certificate.
        self.rejections_certified = bool(rejections_certified)

    @property
    def rank(self):
        if self._rank_cache is None:
            self._rank_cache = self.matrix.rank()
        return self._rank_cache

    @property
    def is_full_rank(self):
        return self.rank == len(self.targets)

    def kernel_basis(self):
        return tuple(tuple(vector) for vector in self.matrix.right_kernel().basis())


class ExactIncrementalRowBasis(SageObject):
    r"""Incremental echelon basis without rebuilding a matrix per row."""

    def __init__(self, coefficient_field, column_count):
        self.coefficient_field = coefficient_field
        self.column_count = ZZ(column_count)
        self._pivots = []
        self._rows = []

    @property
    def rank(self):
        return ZZ(len(self._rows))

    def _reduce(self, row):
        """Return the reduction of ``row`` and its pivot, without inserting."""
        if len(row) != self.column_count:
            raise ValueError("incremental row has the wrong column count")
        vector = [self.coefficient_field(value) for value in row]
        for pivot, basis_row in zip(self._pivots, self._rows):
            coefficient = vector[pivot]
            if coefficient:
                for column in range(pivot, self.column_count):
                    vector[column] -= coefficient * basis_row[column]
        pivot = next(
            (column for column, value in enumerate(vector) if value), None
        )
        return vector, pivot

    def residue_is_nonzero(self, row):
        """Test independence against the basis without mutating it."""
        return self._reduce(row)[1] is not None

    def add(self, row):
        vector, pivot = self._reduce(row)
        if pivot is None:
            return False
        leading = vector[pivot]
        for column in range(pivot, self.column_count):
            vector[column] /= leading
        self._pivots.append(ZZ(pivot))
        self._rows.append(tuple(vector))
        return True


class CertifiedIncrementalRowBasis(SageObject):
    r"""Row-independence filter that certifies accepts by specialization.

    Exact echelon elimination over ``QQ(lambda0,lambda1,lambda2)`` dominated
    large orchestration runs: coefficients grow into tens of thousands of
    digits and integer multiplication enters the FFT range.  Specializing the
    auxiliary weights to rational numbers makes elimination cheap, but naive
    filtering is unsound in both directions, because a dependence whose
    coefficients have poles at the point need not specialize.

    The sound protocol maintains, per specialization point, the invariant
    that every *accepted* row specializes there and the specialized rows are
    linearly independent.  Under that invariant, dependence coefficients are
    ratios of minors that are regular and nonvanishing at the point, so exact
    dependence would survive specialization.  Hence:

    * a row whose specialization is independent of the specialized accepted
      rows is **certainly** exactly independent -- accepted cheaply;
    * a row whose specialization reduces to zero at every maintained point is
      ambiguous.  With ``certify_rejections=True`` (default) it is decided by
      exact elimination, reproducing the fully exact behavior.  With
      ``certify_rejections=False`` it is rejected on the specialization
      evidence alone and the basis records ``rejections_certified=False``;
      accepted rows remain certified either way.

    Whenever an exact decision accepts a row that no point certified, every
    maintained point has lost its invariant and is rebuilt at fresh points.
    If the coefficient field does not support rational specialization, or
    rebuilding fails repeatedly, the filter degrades to pure exact
    elimination.
    """

    def __init__(self, coefficient_field, column_count,
                 certify_rejections=True, specialization_points=None,
                 view_count=2, max_rebuild_attempts=8):
        self.coefficient_field = coefficient_field
        self.column_count = ZZ(column_count)
        self.certify_rejections = bool(certify_rejections)
        self.rejections_certified = True
        self._accepted = []
        self._exact = None
        self._exact_synced = 0
        self._max_rebuild_attempts = int(max_rebuild_attempts)
        self._specializer = self._build_specializer(coefficient_field)
        if specialization_points is not None:
            self._points = iter(specialization_points)
        else:
            self._points = self._default_points()
        self._views = []
        if self._specializer is not None:
            for _ in range(int(view_count)):
                view = self._fresh_view()
                if view is None:
                    self._views = []
                    self._specializer = None
                    break
                self._views.append(view)

    @property
    def rank(self):
        return ZZ(len(self._accepted))

    @staticmethod
    def _build_specializer(field):
        """Return ``(specialize(value, point), arity)`` or ``None``."""
        try:
            gens = field.gens()
            arity = len(gens)
            if arity == 0:
                return None

            def specialize(value, point):
                numerator = value.numerator()
                denominator = value.denominator()
                if arity == 1:
                    top = numerator(point[0])
                    bottom = denominator(point[0])
                else:
                    top = numerator(*point)
                    bottom = denominator(*point)
                bottom = QQ(bottom)
                if not bottom:
                    raise ZeroDivisionError
                return QQ(top) / bottom

            # Trial run: the field must actually specialize to QQ, and the
            # generators must be genuine variables rather than constants.
            trial = tuple(QQ(index + 2) for index in range(arity))
            for generator in gens:
                image = specialize(field(generator), trial)
                if image != trial[gens.index(generator)]:
                    return None
            specialize(field.one(), trial)
            return specialize, arity
        except (AttributeError, TypeError, ValueError, ZeroDivisionError):
            return None

    def _default_points(self):
        import random as _random
        generator = _random.Random(int(20260825))
        _, arity = self._specializer
        while True:
            yield tuple(QQ(generator.randint(int(2), int(10**6)))
                        for _ in range(arity))

    def _specialize_row(self, row, point):
        specialize, _ = self._specializer
        try:
            return tuple(specialize(value, point) for value in row)
        except ZeroDivisionError:
            return None

    def _fresh_view(self):
        """Build a point whose specialized accepted rows are independent."""
        for _ in range(self._max_rebuild_attempts):
            try:
                point = next(self._points)
            except StopIteration:
                return None
            basis = ExactIncrementalRowBasis(QQ, self.column_count)
            for row in self._accepted:
                spec = self._specialize_row(row, point)
                if spec is None or not basis.add(spec):
                    basis = None
                    break
            if basis is not None:
                return [point, basis]
        return None

    def _rebuild_views(self):
        rebuilt = []
        for _ in range(len(self._views) or 1):
            view = self._fresh_view()
            if view is None:
                break
            rebuilt.append(view)
        if not rebuilt:
            self._views = []
            self._specializer = None
        else:
            self._views = rebuilt

    def _sync_exact(self):
        if self._exact is None:
            self._exact = ExactIncrementalRowBasis(
                self.coefficient_field, self.column_count
            )
            self._exact_synced = 0
        while self._exact_synced < len(self._accepted):
            accepted = self._exact.add(self._accepted[self._exact_synced])
            if not accepted:
                raise AssertionError(
                    "certified rows must stay exactly independent"
                )
            self._exact_synced += 1

    def _record_accept(self, row):
        # A materialized exact basis is synchronized lazily by
        # ``_sync_exact`` through the ``_exact_synced`` counter.
        self._accepted.append(row)
        broken = False
        for view in self._views:
            point, basis = view
            spec = self._specialize_row(row, point)
            if spec is None or not basis.add(spec):
                broken = True
                break
        if broken:
            self._rebuild_views()

    def _decide_exactly(self, row):
        self._sync_exact()
        accepted = self._exact.add(row)
        if accepted:
            self._accepted.append(row)
            self._exact_synced += 1
            # No maintained point certified this row, so every invariant is
            # broken; rebuild them all at fresh points.
            if self._views:
                self._rebuild_views()
        return accepted

    def add(self, row):
        if len(row) != self.column_count:
            raise ValueError("incremental row has the wrong column count")
        exact_row = tuple(self.coefficient_field(value) for value in row)
        for view in self._views:
            point, basis = view
            spec = self._specialize_row(exact_row, point)
            if spec is None:
                continue
            if basis.residue_is_nonzero(spec):
                self._record_accept(exact_row)
                return True
        if not self._views:
            return self._decide_exactly(exact_row)
        if not self.certify_rejections:
            self.rejections_certified = False
            return False
        return self._decide_exactly(exact_row)


class ProbeFactory(SageObject):


    """Generate stationary candidates and choose exact independent rows."""

    def __init__(self, coefficient_field=QQ):
        self.coefficient_field = coefficient_field

    def stationary_candidates(self, genus, ambient_degree, max_markings=4,
                              psi_convention="stabilized"):
        genus = ZZ(genus)
        ambient_degree = ZZ(ambient_degree)
        max_markings = ZZ(max_markings)
        if genus < 1 or max_markings < 1:
            raise ValueError("stationary candidate generation needs g,n positive")
        if psi_convention not in ProbeSpec._PSI_CONVENTIONS:
            raise ValueError("unknown stationary psi convention")
        total_psi = 2 * genus - 2
        candidates = []
        for markings in range(1, max_markings + 1):
            profiles = set(
                tuple(sorted(profile, reverse=True))
                for profile in IntegerVectors(total_psi, markings)
            )
            for profile in sorted(profiles, reverse=True):
                candidate = ProbeSpec.stationary(
                    genus, ambient_degree, profile,
                    label="stationary candidate",
                )
                candidates.append(
                    candidate.with_psi_convention(psi_convention)
                )
        return tuple(candidates)

    @staticmethod
    def _symmetric_psi_profiles(total, count):
        """Return partitions of ``total`` among indistinguishable markings."""
        total = ZZ(total)
        count = ZZ(count)
        if count < 0 or total < 0:
            return tuple()
        if count == 0:
            return (tuple(),) if total == 0 else tuple()
        return tuple(sorted(set(
            tuple(sorted(profile, reverse=True))
            for profile in IntegerVectors(total, count)
        ), reverse=True))

    def chow_candidates(self, genus, ambient_degree, max_markings=1,
                        max_unit_insertions=1, primary_only=False):
        r"""Generate positive-dimensional probes for graded Chow extraction.

        These probes are not themselves numerical invariants.  A caller must
        use only Laurent powers ``k <= -probe.dimension_defect()``; the
        :class:`CJRInfinityChowBackend` enforces that condition.  Markings of
        the same insertion kind are indistinguishable, so their psi profiles
        are stored in decreasing order.
        """
        genus = ZZ(genus)
        ambient_degree = ZZ(ambient_degree)
        max_markings = ZZ(max_markings)
        max_unit_insertions = ZZ(max_unit_insertions)
        if genus < 1:
            raise ValueError("Chow probes need positive genus")
        if max_markings < 0 or max_unit_insertions < 0:
            raise ValueError("Chow-probe bounds must be nonnegative")

        candidates = []
        unmarked = ProbeSpec(
            genus, ambient_degree, (), connected=True,
            label="CJR Chow probe", psi_convention="stabilized",
        )
        if unmarked.dimension_defect() > 0:
            candidates.append(unmarked)

        for marking_count in range(1, max_markings + 1):
            for unit_count in range(
                    min(marking_count, max_unit_insertions) + 1):
                point_count = marking_count - unit_count
                top_psi = 2 * genus - 2 + unit_count
                psi_totals = (ZZ.zero(),) if primary_only else range(top_psi)
                for total_psi in psi_totals:
                    for unit_total in range(total_psi + 1):
                        point_total = total_psi - unit_total
                        for unit_profile in self._symmetric_psi_profiles(
                                unit_total, unit_count):
                            for point_profile in self._symmetric_psi_profiles(
                                    point_total, point_count):
                                insertions = tuple(
                                    EllipticInsertion("unit", power)
                                    for power in unit_profile
                                ) + tuple(
                                    EllipticInsertion("point", power)
                                    for power in point_profile
                                )
                                probe = ProbeSpec(
                                    genus, ambient_degree, insertions,
                                    connected=True, label="CJR Chow probe",
                                    psi_convention="stabilized",
                                )
                                if probe.dimension_defect() > 0:
                                    candidates.append(probe)
        return tuple(candidates)

    def unit_relatives(self, probe):
        r"""Return dimension-zero string and dilaton relatives of a probe.

        For a stationary probe, adding ``tau_1(1)`` gives its dilaton
        relative.  Adding ``tau_0(1)`` and raising one point-descendant by
        one gives a string relative.  With stabilized psi their known
        elliptic-curve values are computed by
        :class:`EllipticProbeValueBackend`.  CJR graph compilation uses the
        log-domain convention instead; only convention-compatible sides may
        be combined in a localization equation.
        """
        if not isinstance(probe, ProbeSpec):
            raise TypeError("unit relatives require a ProbeSpec")
        if any(insertion.kind != "point" for insertion in probe.insertions):
            raise ValueError("unit relatives are generated from stationary probes")

        relatives = []
        for index, insertion in enumerate(probe.insertions):
            raised = EllipticInsertion("point", insertion.psi_power + 1)
            insertions = (
                (EllipticInsertion("unit", 0),)
                + probe.insertions[:index]
                + (raised,)
                + probe.insertions[index + 1:]
            )
            relatives.append(ProbeSpec(
                probe.genus, probe.ambient_degree, insertions,
                connected=probe.connected, label="string relative",
                psi_convention=probe.psi_convention,
            ))
        relatives.append(ProbeSpec(
            probe.genus, probe.ambient_degree,
            (EllipticInsertion("unit", 1),) + probe.insertions,
            connected=probe.connected, label="dilaton relative",
            psi_convention=probe.psi_convention,
        ))
        return tuple(relatives)

    def mixed_unit_candidates(self, genus, ambient_degree,
                              max_point_markings=4,
                              max_unit_insertions=1,
                              psi_convention="stabilized"):
        r"""Generate canonical string/dilaton-reducible mixed probes.

        There are ``p`` point insertions and ``u`` unit insertions, with every
        unit carrying psi power zero or one.  If ``a`` of the units carry
        psi power one, the dimension-zero condition is

        ``sum(point psi powers) = 2*g - 2 + u - a``.

        Only one representative of each permutation orbit is needed because
        the compact invariant and the complete localization sum are symmetric
        in the markings.  At unit depth one this is the deduplicated union of
        :meth:`unit_relatives` over all stationary candidates; larger depths
        provide genuinely new string/dilaton equations.
        """
        genus = ZZ(genus)
        ambient_degree = ZZ(ambient_degree)
        max_point_markings = ZZ(max_point_markings)
        max_unit_insertions = ZZ(max_unit_insertions)
        if genus < 1 or max_point_markings < 1:
            raise ValueError("mixed candidates need positive genus and point count")
        if max_unit_insertions < 1:
            raise ValueError("max_unit_insertions must be positive")

        candidates = []
        for point_count in range(1, max_point_markings + 1):
            for unit_count in range(1, max_unit_insertions + 1):
                for dilaton_count in range(unit_count + 1):
                    total_point_psi = (
                        2 * genus - 2 + unit_count - dilaton_count
                    )
                    profiles = set(
                        tuple(sorted(profile, reverse=True))
                        for profile in IntegerVectors(
                            total_point_psi, point_count
                        )
                    )
                    units = (
                        (EllipticInsertion("unit", 0),)
                        * (unit_count - dilaton_count)
                        + (EllipticInsertion("unit", 1),) * dilaton_count
                    )
                    if unit_count == 1:
                        label = (
                            "dilaton relative" if dilaton_count
                            else "string relative"
                        )
                    else:
                        label = "mixed unit depth %s" % unit_count
                    for profile in sorted(profiles, reverse=True):
                        points = tuple(
                            EllipticInsertion("point", power)
                            for power in profile
                        )
                        candidates.append(ProbeSpec(
                            genus, ambient_degree, units + points,
                            connected=True, label=label,
                            psi_convention=psi_convention,
                        ))
        return tuple(candidates)

    def diagonal_row(self, relation, targets, lower_values=None):
        targets = tuple(targets)
        target_set = set(targets)
        lower_values = dict(lower_values or {})
        row = [relation.coefficient_field.zero()] * len(targets)
        nonlinear = False
        unresolved = set()
        for coefficient, factors in relation.terms:
            block_factors = [factor for factor in factors if factor in target_set]
            if len(block_factors) > 1:
                nonlinear = True
                continue
            multiplier = coefficient
            skipped = False
            for factor in factors:
                if factor in target_set and not skipped:
                    skipped = True
                    continue
                if factor not in lower_values:
                    unresolved.add(factor)
                    multiplier = None
                    break
                multiplier *= relation.coefficient_field(lower_values[factor])
            if multiplier is not None and block_factors:
                row[targets.index(block_factors[0])] += multiplier
        return tuple(row), nonlinear, tuple(sorted(
            unresolved, key=lambda item: item.order_key()
        ))

    ROW_FILTERS = ("exact", "certified", "fast")

    def select_relations(self, targets, relations, lower_values=None,
                         row_filter="exact",
                         specialization_points=None):
        r"""Choose exact independent rows for a same-rank target block.

        ``row_filter`` selects the independence test:

        * ``"exact"`` (default): incremental echelon elimination over the
          coefficient field.  Measured fastest on the genus-three frontier
          component (343 targets, 576 rows: 1.5s, against 5.7s certified and
          4.3s fast), because most candidate rows there carry small
          coefficients and specializing every entry costs more than
          eliminating it.
        * ``"certified"``: specialization certificates for acceptance with
          exact fallback for ambiguous rows.  Mathematically equivalent to
          ``"exact"``; intended for row families whose entries are large
          enough that fraction-field elimination enters the FFT range.
        * ``"fast"``: like ``"certified"``, but ambiguous rows are rejected
          on specialization evidence alone.  Accepted rows remain exactly
          certified, so a full-rank selection is sound; a deficient one
          carries ``rejections_certified=False`` and must not be quoted as a
          certified kernel.

        ``specialization_points`` injects deterministic points for tests.
        """
        targets = tuple(targets)
        if not targets:
            raise ValueError("a probe block needs at least one target")
        if row_filter not in self.ROW_FILTERS:
            raise ValueError("unknown row filter %r" % (row_filter,))
        selected_relations = []
        selected_rows = []
        rejected = []
        if row_filter == "exact":
            row_basis = ExactIncrementalRowBasis(
                self.coefficient_field, len(targets)
            )
        else:
            row_basis = CertifiedIncrementalRowBasis(
                self.coefficient_field, len(targets),
                certify_rejections=(row_filter == "certified"),
                specialization_points=specialization_points,
            )
        for relation in relations:
            if relation.coefficient_field is not self.coefficient_field:
                raise TypeError("all candidate relations must use the factory field")
            if not relation.is_complete:
                rejected.append((relation, "incomplete graph compilation"))
                continue
            row, nonlinear, unresolved = self.diagonal_row(
                relation, targets, lower_values
            )
            if nonlinear:
                rejected.append((relation, "nonlinear in target block"))
                continue
            if unresolved:
                rejected.append((relation, "unresolved lower factors %s" % (unresolved,)))
                continue
            if row_basis.add(row):
                selected_relations.append(relation)
                selected_rows.append(row)
            else:
                rejected.append((relation, "row is dependent"))
            if row_basis.rank == len(targets):
                break
        return ProbeBlockSelection(
            targets, selected_relations, selected_rows,
            self.coefficient_field, rejected=rejected,
            rows_independent=True,
            rejections_certified=getattr(
                row_basis, "rejections_certified", True
            ),
        )
