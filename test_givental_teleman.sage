"""Tests for the generic R-action and the calibrated O(3) connection."""

load("o3_fixed_locus_graphs.sage")
load("o3_givental_teleman.sage")


def run_tests():
    # Identity R recovers a rank-one TFT exactly.
    identity_r = TruncatedRMatrix(QQ, (
        identity_matrix(QQ, 1), zero_matrix(QQ, 1, 1),
    ))
    tft = SemisimpleTFTData(QQ, (QQ(2),), identity_r)
    identity = GiventalTelemanCorrelator(tft)
    assert identity.evaluate(0, (((QQ.one(),), 0),) * 3) == 2
    assert identity.evaluate(1, (((QQ(3),), 1),)) == QQ(1) / 8

    rings = PlaneCubicCoefficientRing()
    localization = P2FixedLocusEvaluator(rings, include_twist=True)
    givental = O3DegreeZeroGiventalBackend(rings)

    requests = (
        TwistedZeroVertexRequest(0, 0, (
            TwistedInsertion(1), TwistedInsertion(1), TwistedInsertion(0),
        )),
        TwistedZeroVertexRequest(1, 0, (TwistedInsertion(0),)),
        TwistedZeroVertexRequest(1, 0, (TwistedInsertion(1, 1),)),
        TwistedZeroVertexRequest(2, 0, ()),
    )
    for request in requests:
        assert givental.evaluate(request) == localization.evaluate(request)

    data = O3DegreeZeroGiventalData(rings, max_power=5)
    assert all(not defect for defect in data.r_matrix.symplectic_defects(
        data.tft.metric
    ))

    quantum = O3TwistedQuantumRing(rings, max_degree=3)
    momenta = quantum.canonical_momenta()
    assert len(momenta) == 3
    assert all(quantum._truncate(quantum.relation(root)) == 0
               for root in momenta)
    assert tuple(metric[0] for metric in quantum.canonical_metrics()) \
        == data.tft.metric
    multiplication = quantum.quantum_multiplication_by_h()
    relation_coefficients = quantum.relation.list()
    relation_matrix = sum(
        relation_coefficients[power] * multiplication ** power
        for power in range(len(relation_coefficients))
    )
    assert not relation_matrix

    hybrid = HybridTwistedZeroVertexBackend(
        FullTwistedZeroVertexBackend(rings), givental=givental, verify=True
    )
    assert hybrid.evaluate(requests[1]) == localization.evaluate(requests[1])

    positive = TwistedZeroVertexRequest(0, 1, (
        TwistedInsertion(2), TwistedInsertion(2),
    ))
    assert hybrid.evaluate(positive) \
        == FullTwistedZeroVertexBackend(rings).evaluate(positive)

    # The positive-q calibration uses standard-lift genus-zero three-point
    # invariants, the q=0 QRR integration constants, and the S-matrix
    # ancestor--descendant transformation.  Primary and descendant checks
    # exercise different parts of this construction.
    standard_localization = FullTwistedZeroVertexBackend(
        rings, lift_strategy="standard"
    )
    calibrated = O3CalibratedGiventalBackend(
        rings, lift_strategy="standard"
    )
    positive_requests = (
        TwistedZeroVertexRequest(1, 1, (TwistedInsertion(1),)),
        TwistedZeroVertexRequest(1, 1, (TwistedInsertion(1, 1),)),
        TwistedZeroVertexRequest(1, 2, (TwistedInsertion(1, 1),)),
    )
    for request in positive_requests:
        assert calibrated.evaluate(request) \
            == standard_localization.evaluate(request)

    calibrated_data, _ = calibrated._evaluator(positive_requests[-1])
    assert all(not defect for defect in
               calibrated_data.r_matrix.symplectic_defects(
                   calibrated_data.canonical_metric
               ))
    assert calibrated_data.quantum_product.transpose() \
        * calibrated_data._series_matrix(calibrated_data.flat_metric) \
        == calibrated_data._series_matrix(calibrated_data.flat_metric) \
        * calibrated_data.quantum_product

    verified_hybrid = HybridTwistedZeroVertexBackend(
        standard_localization, givental=calibrated, verify=True
    )
    assert verified_hybrid.evaluate(positive_requests[1]) \
        == standard_localization.evaluate(positive_requests[1])


run_tests()
print("all Givental-Teleman reconstruction tests passed")
