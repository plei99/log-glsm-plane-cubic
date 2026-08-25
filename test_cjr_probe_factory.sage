"""Tests for certified specialization-filtered probe-row selection."""

load("cjr_probe_factory.sage")


def _field():
    return PolynomialRing(
        QQ, names=("lambda0", "lambda1", "lambda2")
    ).fraction_field()


def test_certificate_matches_exact_on_generic_rows():
    """Accept/reject decisions must equal fully exact elimination."""
    field = _field()
    l0, l1, l2 = field.gens()
    rows = (
        (1, l0, l1 * l2),
        (l0, l0**2, l0 * l1 * l2),          # multiple of the first
        (0, 1, l2),
        (1, l0 + 1, l1 * l2 + l2),          # row1 + row3
        (l1 / (l0 + 7), 0, 1),
        (0, 0, 0),
        (2, 2 * l0 + 1, 2 * l1 * l2 + l2),  # 2*row1 + row3
    )
    certified = CertifiedIncrementalRowBasis(field, 3)
    exact = ExactIncrementalRowBasis(field, 3)
    for row in rows:
        assert certified.add(row) == exact.add(row)
    assert certified.rank == exact.rank == 3
    assert certified.rejections_certified


def test_pole_coefficient_dependence_is_still_rejected():
    r"""The classical trap: dependence whose coefficients have poles.

    ``r = (r1 - r2) / (2*(l0-5))`` is exactly dependent on ``r1, r2``, but
    the dependence coefficients blow up on ``l0 = 5``.  At a point with
    ``l0 = 5`` the specializations of ``r1`` and ``r2`` coincide, so a naive
    one-point filter would see ``spec(r) = (0,0,1)`` as independent and
    wrongly accept an exactly dependent row.  The invariant machinery must
    rebuild the poisoned view and reject exactly.
    """
    field = _field()
    l0, _, _ = field.gens()
    r1 = (field(1), field(0), l0 - 5)
    r2 = (field(1), field(0), 5 - l0)
    dependent = (field(0), field(0), field(1))

    points = iter((
        (QQ(5), QQ(7), QQ(11)),     # poisoned: l0 = 5
        (QQ(13), QQ(17), QQ(19)),   # sound second view
        (QQ(23), QQ(29), QQ(31)),   # rebuild target for the poisoned view
        (QQ(37), QQ(41), QQ(43)),
    ))
    basis = CertifiedIncrementalRowBasis(
        field, 3, specialization_points=points
    )
    assert basis.add(r1)
    assert basis.add(r2)
    assert basis.rank == 2
    assert not basis.add(dependent), \
        "an exactly dependent row must never be accepted"
    assert basis.rank == 2
    assert basis.rejections_certified

    # The same trap in fully exact arithmetic, as the reference answer.
    exact = ExactIncrementalRowBasis(field, 3)
    assert exact.add(r1) and exact.add(r2) and not exact.add(dependent)


def test_point_exhaustion_degrades_to_exact_mode():
    """Running out of specialization points must fall back to exact rows."""
    field = _field()
    l0, _, _ = field.gens()
    points = iter(((QQ(5), QQ(7), QQ(11)),))  # a single poisoned point
    basis = CertifiedIncrementalRowBasis(
        field, 3, specialization_points=points, view_count=1
    )
    assert basis.add((field(1), field(0), l0 - 5))
    # Ambiguous at the only point; exact accept breaks the invariant and the
    # rebuild exhausts the injected points, forcing pure exact mode.
    assert basis.add((field(1), field(0), 5 - l0))
    assert not basis._views
    assert not basis.add((field(0), field(0), field(1)))
    assert basis.rank == 2


def test_uncertified_rejections_are_flagged():
    """Fast mode may under-select but must say so, and never over-select."""
    field = _field()
    l0, _, _ = field.gens()
    # Both views poisoned on l0 = 5: the exactly independent second row is
    # spec-dependent everywhere the filter looks.
    points = iter((
        (QQ(5), QQ(7), QQ(11)),
        (QQ(5), QQ(13), QQ(17)),
    ))
    basis = CertifiedIncrementalRowBasis(
        field, 3, certify_rejections=False, specialization_points=points
    )
    assert basis.add((field(1), field(0), l0 - 5))
    assert not basis.add((field(1), field(0), 5 - l0))
    assert not basis.rejections_certified
    assert basis.rank == 1


def test_rational_field_disables_certificates():
    """Fields without genuine variables must use exact elimination."""
    basis = CertifiedIncrementalRowBasis(QQ, 2)
    assert not basis._views
    assert basis.add((1, 0))
    assert basis.add((0, 1))
    assert not basis.add((2, 3))
    assert basis.rank == 2
    assert basis.rejections_certified


def test_select_relations_certification_metadata():
    """Selections built by the factory carry the certification flags."""
    field = _field()
    l0, _, _ = field.gens()
    x = EffectiveVertex(1, 0, (-1,), label="x")
    y = EffectiveVertex(1, 0, (-1, -1), label="y")
    probe_1 = ProbeSpec.stationary(1, 0, (0,), label="p1")
    probe_2 = ProbeSpec.stationary(1, 0, (1,), label="p2")
    relation_1 = CompiledLocalizationRelation(
        probe_1, 0, field, 1, terms=((l0, (x,)), (1, (y,))),
    )
    relation_2 = CompiledLocalizationRelation(
        probe_2, 0, field, 2, terms=((1, (x,)), (l0, (y,))),
    )
    factory = ProbeFactory(field)
    for row_filter in ProbeFactory.ROW_FILTERS:
        selection = factory.select_relations(
            (x, y), (relation_1, relation_2), row_filter=row_filter
        )
        assert selection.is_full_rank
        assert selection.rank == 2
        assert selection.rejections_certified
        assert selection.relations == (relation_1, relation_2)

    # Deficient block: a single relation for two targets.
    deficient = factory.select_relations((x, y), (relation_1,))
    assert not deficient.is_full_rank
    assert deficient.rank == 1
    assert deficient.rejections_certified
    kernel = deficient.kernel_basis()
    assert len(kernel) == 1


test_certificate_matches_exact_on_generic_rows()
test_pole_coefficient_dependence_is_still_rejected()
test_point_exhaustion_degrades_to_exact_mode()
test_uncertified_rejections_are_flagged()
test_rational_field_disables_certificates()
test_select_relations_certification_metadata()
print("all certified probe-factory tests passed")
