r"""Givental--Teleman backend for the equivariant degree-zero O(3) theory.

At ambient degree zero the auxiliary base torus makes the state algebra a
direct sum of three fields.  At the fixed point ``p_i`` its metric is

``eta_i = (t-3 lambda_i) / prod_(j != i)(lambda_i-lambda_j)``.

The normalized Hodge class is a product of two tangent Euler twists and the
inverse fibre twist. Mumford GRR therefore gives a diagonal R-matrix.  This
module reconstructs the resulting CohFT by the general stable-graph engine in
``givental_teleman.sage``.  In positive degree it calibrates the R-matrix from
the small quantum connection and applies the S-matrix ancestor--descendant
transformation before extracting an individual Novikov coefficient.
"""

import itertools

if "TwistedZeroVertexBackend" not in globals():
    load("o3_twisted_plane_vertices.sage")
if "FullTwistedZeroVertexBackend" not in globals():
    load("o3_fixed_locus_graphs.sage")
if "GiventalTelemanCorrelator" not in globals():
    load("givental_teleman.sage")


class O3TwistedQuantumRing(SageObject):
    r"""Truncated small quantum Frobenius algebra of the dual Euler twist.

    The dual-twisted Picard--Fuchs principal symbol is

    ``prod_i(p-lambda_i) = q (t-3p)^3``.

    Its three roots are constructed by Hensel lifting from ``p=lambda_i``.
    This avoids algebraic root expressions and keeps every coefficient in the
    exact auxiliary-equivariant fraction field.
    """

    def __init__(self, rings=None, max_degree=3):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.field = self.rings.full_field
        self.lambdas = self.rings.base_weights()
        self.max_degree = ZZ(max_degree)
        if self.max_degree < 0:
            raise ValueError("the quantum-degree bound must be nonnegative")
        self.series_ring = PowerSeriesRing(
            self.field, "q", default_prec=int(self.max_degree + 1)
        )
        self.q = self.series_ring.gen()
        self.polynomial_ring = PolynomialRing(self.series_ring, "p")
        self.p = self.polynomial_ring.gen()
        self.relation = self.polynomial_ring(
            prod(self.p - weight for weight in self.lambdas)
            - self.q * (self.rings.t - 3 * self.p) ** 3
        )
        self._momenta = None
        self._metrics = None

    def _truncate(self, value):
        return self.series_ring(value).add_bigoh(int(self.max_degree + 1))

    def canonical_momenta(self):
        if self._momenta is None:
            derivative = self.relation.derivative(self.p)
            roots = []
            for weight in self.lambdas:
                root = self._truncate(weight)
                # Newton--Hensel iteration doubles q-adic precision; a few
                # extra iterations are harmless for the small truncations used
                # in higher-genus reconstruction.
                for _ in range(int(ceil(log(max(2, self.max_degree + 1), 2))) + 1):
                    numerator = self._truncate(self.relation(root))
                    denominator = self._truncate(derivative(root))
                    root = self._truncate(root - numerator / denominator)
                roots.append(root)
            self._momenta = tuple(roots)
        return self._momenta

    def canonical_metrics(self):
        r"""Diagonal residue metric ``(t-3p_i)/partial_p relation(p_i)``."""
        if self._metrics is None:
            derivative = self.relation.derivative(self.p)
            self._metrics = tuple(
                self._truncate(
                    (self.rings.t - 3 * root) / derivative(root)
                )
                for root in self.canonical_momenta()
            )
        return self._metrics

    def quantum_multiplication_by_h(self):
        r"""Return multiplication by ``H`` in the basis ``1,H,H^2``."""
        coefficients = self.relation.list()
        coefficients += [self.series_ring.zero()] * (4 - len(coefficients))
        leading = coefficients[3]
        if not leading.is_unit():
            raise ArithmeticError("the quantum relation is not monic up to a unit")
        reduction = tuple(-coefficients[index] / leading for index in range(3))
        return Matrix(self.series_ring, (
            (0, 0, reduction[0]),
            (1, 0, reduction[1]),
            (0, 1, reduction[2]),
        ))

    def planned_insertions(self, insertions, lift_strategy="sparse"):
        return O3DegreeZeroGiventalData(
            self.rings, max_power=1, lift_strategy=lift_strategy
        ).planned_insertions(insertions)

    def insertion_vector(self, insertion):
        values = []
        for root in self.canonical_momenta():
            if insertion.h_power == 0:
                value = self.series_ring.one()
            elif insertion.lift_kind == "standard":
                value = root ** insertion.h_power
            elif insertion.lift_kind == "vanish":
                value = root - self.lambdas[int(insertion.lift_index)]
            else:
                index = int(insertion.lift_index)
                value = prod(
                    root - self.lambdas[other]
                    for other in range(3) if other != index
                )
            values.append(self._truncate(insertion.scale * value))
        return tuple(values)

    def tft_data(self, r_matrix, label="positive-degree O(3) calibration"):
        """Combine these Frobenius data with a calibrated q-series R-matrix."""
        return SemisimpleTFTData(
            self.series_ring, self.canonical_metrics(), r_matrix,
            label=label,
        )


class O3DegreeZeroGiventalData(SageObject):
    """Exact semisimple TFT and diagonal R-matrix at ``q=0``."""

    def __init__(self, rings=None, max_power=8, include_twist=True,
                 lift_strategy="sparse"):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.field = self.rings.full_field
        self.lambdas = self.rings.base_weights()
        self.max_power = ZZ(max_power)
        self.include_twist = bool(include_twist)
        self.lift_strategy = str(lift_strategy)
        if self.max_power < 1:
            raise ValueError("R-matrix precision must be positive")
        if self.lift_strategy not in ("sparse", "standard"):
            raise ValueError("lift strategy must be sparse or standard")

        metric = []
        diagonal_coefficients = []
        for colour in range(3):
            tangent_weights = tuple(
                self.lambdas[colour] - self.lambdas[other]
                for other in range(3) if other != colour
            )
            tangent_euler = prod(tangent_weights)
            fibre_weight = self.rings.t - 3 * self.lambdas[colour]
            metric.append(
                fibre_weight / tangent_euler
                if self.include_twist else tangent_euler ** (-1)
            )

            logarithm = {}
            for degree in range(1, self.max_power + 1):
                bernoulli_value = bernoulli(degree + 1)
                if not bernoulli_value:
                    continue
                characteristic = -factorial(degree - 1) * sum(
                    weight ** (-degree) for weight in tangent_weights
                )
                if self.include_twist:
                    characteristic -= (
                        factorial(degree - 1) * fibre_weight ** (-degree)
                    )
                logarithm[ZZ(degree)] = (
                    characteristic * bernoulli_value
                    / factorial(degree + 1)
                )
            diagonal_coefficients.append(
                exponential_series_coefficients(
                    self.field, logarithm, self.max_power
                )
            )

        matrices = []
        for degree in range(self.max_power + 1):
            matrices.append(diagonal_matrix(
                self.field,
                tuple(diagonal_coefficients[colour][degree]
                      for colour in range(3)),
            ))
        self.r_matrix = TruncatedRMatrix(self.field, tuple(matrices))
        self.tft = SemisimpleTFTData(
            self.field, tuple(metric), self.r_matrix,
            label="degree-zero equivariant O(3)-twisted P2",
        )

    def planned_insertions(self, insertions):
        insertions = tuple(insertions)
        if self.lift_strategy == "standard":
            return insertions
        counters = {1: ZZ.zero(), 2: ZZ.zero()}
        planned = []
        for insertion in insertions:
            if insertion.h_power in (1, 2) \
                    and insertion.lift_kind == "standard":
                index = counters[int(insertion.h_power)] % 3
                counters[int(insertion.h_power)] += 1
                planned.append(insertion.with_sparse_lift(index))
            else:
                planned.append(insertion)
        return tuple(planned)

    def insertion_vector(self, insertion):
        return tuple(
            self.field(insertion.restriction(colour, self.lambdas))
            for colour in range(3)
        )


class O3DegreeZeroGiventalBackend(TwistedZeroVertexBackend):
    """Evaluate stable degree-zero requests by the diagonal R-action."""

    def __init__(self, rings=None, include_twist=True,
                 lift_strategy="sparse"):
        super().__init__(rings=rings)
        self.include_twist = bool(include_twist)
        self.lift_strategy = str(lift_strategy)
        self._evaluators = {}
        self._cache = {}

    def _evaluator(self, request):
        dimension = 3 * request.genus - 3 + request.valence
        precision = max(ZZ.one(), ZZ(dimension + 1))
        if precision not in self._evaluators:
            data = O3DegreeZeroGiventalData(
                self.rings, precision,
                include_twist=self.include_twist,
                lift_strategy=self.lift_strategy,
            )
            self._evaluators[precision] = (
                data, GiventalTelemanCorrelator(data.tft)
            )
        return self._evaluators[precision]

    def is_cached(self, request):
        scale, normalized = request.normalized_scales()
        return not scale or normalized in self._cache

    def supports(self, request):
        return request.degree == 0

    def provenance(self, request):
        if request.degree == 0:
            return "Givental-Teleman diagonal R-action at q=0"
        return "unsupported: positive-degree R-matrix not calibrated"

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        if request.degree != 0:
            raise UnsupportedGeometryError(
                "the positive-degree O(3) R-matrix is not calibrated for %s"
                % request
            )
        if 2 * request.genus - 2 + request.valence <= 0:
            raise ValueError("the degree-zero stable-map type is unstable")
        scale, normalized = request.normalized_scales()
        if not scale:
            return self.rings.full_field.zero()
        if normalized not in self._cache:
            data, evaluator = self._evaluator(normalized)
            planned = data.planned_insertions(normalized.insertions)
            inputs = tuple(
                (data.insertion_vector(insertion), insertion.psi_power)
                for insertion in planned
            )
            self._cache[normalized] = evaluator.evaluate(
                normalized.genus, inputs
            )
        return self.rings.full(scale * self._cache[normalized])


class O3CalibratedGiventalData(SageObject):
    r"""Positive-degree calibration from the small quantum connection.

    The raw hypergeometric principal symbol does not by itself give flat
    coordinates: its residue pairing varies with ``q``.  We instead compute
    multiplication by ``H`` from exact genus-zero three-point invariants in
    the flat basis ``1,H,H^2``.  If ``C`` is that multiplication matrix,
    ``E`` is the matrix of (unnormalized) idempotents, ``P`` is the diagonal
    matrix of eigenvalues, and ``A=E^(-1) D E`` for ``D=q d/dq``, flatness of
    the canonical fundamental solution gives

    ``[P,R_(k+1)] = D R_k + A R_k``.

    Its diagonal integration constants are precisely the ambiguity left by
    the quantum differential equation.  They are fixed here by the exact
    degree-zero quantum-Riemann--Roch matrix.
    """

    def __init__(self, rings=None, max_degree=1, max_power=3,
                 genus_zero_backend=None, lift_strategy="sparse"):
        self.rings = rings or PlaneCubicCoefficientRing()
        self.field = self.rings.full_field
        self.lambdas = self.rings.base_weights()
        self.max_degree = ZZ(max_degree)
        self.max_power = ZZ(max_power)
        self.lift_strategy = str(lift_strategy)
        if self.max_degree < 0:
            raise ValueError("the quantum-degree bound must be nonnegative")
        if self.max_power < 1:
            raise ValueError("R-matrix precision must be positive")
        if self.lift_strategy not in ("sparse", "standard"):
            raise ValueError("lift strategy must be sparse or standard")

        self.series_ring = PowerSeriesRing(
            self.field, "q", default_prec=int(self.max_degree + 1)
        )
        self.q = self.series_ring.gen()
        # The calibration must use the actual flat cohomology basis.  Sparse
        # lifts are an evaluation optimization and would change the
        # equivariant representatives of the three-point tensor.
        self.genus_zero_backend = genus_zero_backend or \
            FullTwistedZeroVertexBackend(
                self.rings, lift_strategy="standard"
            )

        degree_zero = O3DegreeZeroGiventalData(
            self.rings,
            max_power=self.max_power + self.max_degree,
            include_twist=True,
            lift_strategy=self.lift_strategy,
        )
        self._degree_zero = degree_zero
        vandermonde = Matrix(self.field, (
            tuple(weight ** power for power in range(3))
            for weight in self.lambdas
        ))
        self.flat_metric = (
            vandermonde.transpose()
            * diagonal_matrix(self.field, degree_zero.tft.metric)
            * vandermonde
        )
        self.quantum_product = self._quantum_product_by_h()
        self.canonical_momenta = self._canonical_momenta()
        self.idempotents = self._idempotent_matrix()
        self.canonical_metric = self._canonical_metric()
        self.connection_matrix = self.idempotents.inverse() \
            * self._derivative_matrix(self.idempotents)
        self.metric_connection = diagonal_matrix(
            self.series_ring,
            tuple(
                self._truncate(
                    self._derivative(value) / (2 * value)
                )
                for value in self.canonical_metric
            ),
        )
        self.r_matrix = self._calibrated_r_matrix()
        self._s_matrices = [identity_matrix(self.series_ring, 3)]
        self.tft = SemisimpleTFTData(
            self.series_ring, self.canonical_metric, self.r_matrix,
            label="calibrated positive-degree equivariant O(3)-twisted P2",
        )

    def _truncate(self, value):
        return self.series_ring(value).add_bigoh(int(self.max_degree + 1))

    def _series_matrix(self, matrix):
        return Matrix(
            self.series_ring, matrix.nrows(), matrix.ncols(),
            lambda row, column: self._truncate(matrix[row, column]),
        )

    def _coefficient_matrix(self, matrix, degree):
        degree = ZZ(degree)
        return Matrix(
            self.field, matrix.nrows(), matrix.ncols(),
            lambda row, column: self.field(matrix[row, column][degree]),
        )

    def _derivative(self, value):
        return self._truncate(sum(
            degree * value[degree] * self.q ** degree
            for degree in range(1, self.max_degree + 1)
        ))

    def _derivative_matrix(self, matrix):
        return Matrix(
            self.series_ring, matrix.nrows(), matrix.ncols(),
            lambda row, column: self._derivative(matrix[row, column]),
        )

    def _three_point_series(self, left_power, right_power):
        answer = self.series_ring.zero()
        for degree in range(self.max_degree + 1):
            request = TwistedZeroVertexRequest(0, degree, (
                TwistedInsertion(1),
                TwistedInsertion(left_power),
                TwistedInsertion(right_power),
            ))
            answer += self.series_ring(
                self.genus_zero_backend.evaluate(request)
            ) * self.q ** degree
        return self._truncate(answer)

    def _quantum_product_by_h(self):
        tensor = zero_matrix(self.series_ring, 3, 3)
        for left in range(3):
            for right in range(left, 3):
                value = self._three_point_series(left, right)
                tensor[left, right] = value
                tensor[right, left] = value
        product = self._series_matrix(self.flat_metric.inverse()) * tensor

        # At the large-radius point this must be ordinary equivariant
        # multiplication by H.  This catches lift/sign mistakes before they
        # contaminate every coefficient of the R-matrix.
        elementary_two = sum(
            self.lambdas[left] * self.lambdas[right]
            for left in range(3) for right in range(left + 1, 3)
        )
        classical = Matrix(self.field, (
            (0, 0, prod(self.lambdas)),
            (1, 0, -elementary_two),
            (0, 1, sum(self.lambdas)),
        ))
        if self._coefficient_matrix(product, 0) != classical:
            raise ArithmeticError(
                "the genus-zero calibration has the wrong large-radius product"
            )
        if product.transpose() * self._series_matrix(self.flat_metric) \
                != self._series_matrix(self.flat_metric) * product:
            raise ArithmeticError(
                "twisted quantum multiplication is not metric self-adjoint"
            )
        return product

    def _canonical_momenta(self):
        polynomial_ring = PolynomialRing(self.series_ring, "p_cal")
        variable = polynomial_ring.gen()
        characteristic = self.quantum_product.charpoly(variable)
        derivative = characteristic.derivative(variable)
        roots = []
        iterations = int(ceil(log(max(2, self.max_degree + 1), 2))) + 1
        for weight in self.lambdas:
            root = self._truncate(weight)
            for _ in range(iterations):
                root = self._truncate(
                    root - self._truncate(characteristic(root))
                    / self._truncate(derivative(root))
                )
            if self._truncate(characteristic(root)):
                raise ArithmeticError("Hensel lifting of a canonical root failed")
            roots.append(root)
        return tuple(roots)

    def _idempotent_matrix(self):
        unit = vector(self.series_ring, (1, 0, 0))
        first = self.quantum_product * unit
        second = self.quantum_product * first
        cyclic = Matrix(self.series_ring, (
            tuple(unit), tuple(first), tuple(second)
        )).transpose()
        evaluation = Matrix(self.series_ring, (
            tuple(root ** power for power in range(3))
            for root in self.canonical_momenta
        ))
        idempotents = cyclic * evaluation.inverse()
        diagonal = diagonal_matrix(
            self.series_ring, self.canonical_momenta
        )
        if self.quantum_product * idempotents != idempotents * diagonal:
            raise ArithmeticError("the canonical idempotents do not diagonalize H")
        return idempotents

    def _canonical_metric(self):
        metric = (
            self.idempotents.transpose()
            * self._series_matrix(self.flat_metric)
            * self.idempotents
        )
        for row in range(3):
            for column in range(3):
                if row != column and metric[row, column]:
                    raise ArithmeticError("canonical idempotents are not orthogonal")
        return tuple(metric[index, index] for index in range(3))

    @staticmethod
    def _commutator(left, right):
        return left * right - right * left

    def _calibrated_r_matrix(self):
        # A coefficient q^d of R_k can depend on the large-radius value of
        # R_(k+d), so compute a triangular precision window and discard its
        # auxiliary upper-right edge afterwards.
        internal_power = int(self.max_power + self.max_degree)
        initial = self._degree_zero.r_matrix
        zero = zero_matrix(self.field, 3, 3)
        coefficients = [
            [Matrix(self.field, zero)
             for _ in range(self.max_degree + 1)]
            for _ in range(internal_power + 1)
        ]
        for power in range(internal_power + 1):
            coefficients[power][0] = Matrix(
                self.field, initial.coefficient(power)
            )

        connection = tuple(
            self._coefficient_matrix(self.connection_matrix, degree)
            for degree in range(self.max_degree + 1)
        )
        metric_connection = tuple(
            self._coefficient_matrix(self.metric_connection, degree)
            for degree in range(self.max_degree + 1)
        )
        momenta = tuple(diagonal_matrix(
            self.field,
            tuple(root[degree] for root in self.canonical_momenta),
        ) for degree in range(self.max_degree + 1))
        if connection[0]:
            raise ArithmeticError("the canonical connection must vanish at q=0")

        for degree in range(1, self.max_degree + 1):
            largest_power = internal_power - degree
            for power in range(largest_power + 1):
                current = coefficients[power][degree]
                base = sum(
                    connection[shift] * coefficients[power][degree - shift]
                    - coefficients[power][degree - shift]
                    * metric_connection[shift]
                    for shift in range(1, degree + 1)
                ) - sum(
                    self._commutator(
                        momenta[shift],
                        coefficients[power + 1][degree - shift],
                    )
                    for shift in range(1, degree + 1)
                )
                for index in range(3):
                    current[index, index] = -base[index, index] / degree
                residual = degree * current + base
                if any(residual[index, index] for index in range(3)):
                    raise ArithmeticError(
                        "the diagonal R-matrix flatness recursion failed"
                    )
                if power < largest_power:
                    following = coefficients[power + 1][degree]
                    for row in range(3):
                        for column in range(3):
                            if row == column:
                                continue
                            following[row, column] = (
                                residual[row, column]
                                / (self.lambdas[row] - self.lambdas[column])
                            )

        matrices = []
        for power in range(self.max_power + 1):
            matrices.append(Matrix(
                self.series_ring, 3, 3,
                lambda row, column: self._truncate(sum(
                    coefficients[power][degree][row, column]
                    * self.q ** degree
                    for degree in range(self.max_degree + 1)
                )),
            ))
        calibrated = TruncatedRMatrix(
            self.series_ring, tuple(matrices)
        )
        if any(calibrated.symplectic_defects(self.canonical_metric)):
            raise ArithmeticError(
                "the calibrated R-matrix violates the symplectic condition"
            )
        return calibrated

    def planned_insertions(self, insertions):
        return O3DegreeZeroGiventalData(
            self.rings, max_power=1, lift_strategy=self.lift_strategy
        ).planned_insertions(insertions)

    def _flat_insertion_vector(self, insertion):
        if insertion.h_power == 0:
            flat = vector(self.series_ring, (1, 0, 0))
        elif insertion.lift_kind == "standard":
            flat = vector(
                self.series_ring,
                tuple(1 if power == insertion.h_power else 0
                      for power in range(3)),
            )
        elif insertion.lift_kind == "vanish":
            weight = self.lambdas[int(insertion.lift_index)]
            flat = vector(self.series_ring, (-weight, 1, 0))
        else:
            index = int(insertion.lift_index)
            others = tuple(
                self.lambdas[other] for other in range(3)
                if other != index
            )
            flat = vector(
                self.series_ring,
                (prod(others), -sum(others), 1),
            )
        return vector(
            self.series_ring,
            tuple(self._truncate(insertion.scale * value) for value in flat),
        )

    def insertion_vector(self, insertion):
        coordinates = self.idempotents.inverse() \
            * self._flat_insertion_vector(insertion)
        return tuple(self._truncate(value) for value in coordinates)

    def s_matrices(self, maximum):
        r"""Return the reduced fundamental solution through ``z^-maximum``.

        Removing the classical factor ``q^(H/z)`` gives

        ``D S_(m+1) = C S_m - S_m C(0)``, ``S_0=1``.

        The zero constant term for ``m>0`` is the large-radius calibration.
        """
        maximum = ZZ(maximum)
        if maximum < 0:
            raise ValueError("S-matrix precision must be nonnegative")
        classical = self._series_matrix(
            self._coefficient_matrix(self.quantum_product, 0)
        )
        while len(self._s_matrices) <= maximum:
            previous = self._s_matrices[-1]
            derivative = self.quantum_product * previous \
                - previous * classical
            integrated = Matrix(
                self.series_ring, 3, 3,
                lambda row, column: self._truncate(sum(
                    derivative[row, column][degree]
                    / degree * self.q ** degree
                    for degree in range(1, self.max_degree + 1)
                )),
            )
            self._s_matrices.append(integrated)
        return tuple(self._s_matrices[:int(maximum + 1)])

    def ancestor_terms(self, insertion):
        r"""Expand one stable-map descendant into ancestor insertions.

        ``tau_k(gamma)`` becomes
        ``sum_(m=0)^k bar_tau_(k-m)(S_m gamma)``.
        """
        flat = self._flat_insertion_vector(insertion)
        inverse_idempotents = self.idempotents.inverse()
        return tuple(
            (
                tuple(self._truncate(value) for value in (
                    inverse_idempotents * matrix * flat
                )),
                ZZ(insertion.psi_power - power),
            )
            for power, matrix in enumerate(
                self.s_matrices(insertion.psi_power)
            )
        )


class O3CalibratedGiventalBackend(TwistedZeroVertexBackend):
    """Evaluate every stable degree through a calibrated small R-matrix."""

    def __init__(self, rings=None, genus_zero_backend=None,
                 lift_strategy="sparse"):
        super().__init__(rings=rings)
        self.genus_zero_backend = genus_zero_backend or \
            FullTwistedZeroVertexBackend(
                self.rings, lift_strategy="standard"
            )
        self.lift_strategy = str(lift_strategy)
        self._evaluators = {}
        self._cache = {}

    def _evaluator(self, request):
        dimension = 3 * request.genus - 3 + request.valence
        precision = max(ZZ.one(), ZZ(dimension + 1))
        key = ZZ(request.degree), precision
        covering = tuple(
            candidate for candidate in self._evaluators
            if candidate[0] >= key[0] and candidate[1] >= key[1]
        )
        if covering:
            return self._evaluators[min(
                covering, key=lambda candidate: (
                    candidate[0] + candidate[1], candidate
                )
            )]
        if key not in self._evaluators:
            data = O3CalibratedGiventalData(
                self.rings,
                max_degree=request.degree,
                max_power=precision,
                genus_zero_backend=self.genus_zero_backend,
                lift_strategy=self.lift_strategy,
            )
            self._evaluators[key] = (
                data, GiventalTelemanCorrelator(data.tft)
            )
        return self._evaluators[key]

    def is_cached(self, request):
        scale, normalized = request.normalized_scales()
        return not scale or normalized in self._cache

    def supports(self, request):
        return 2 * request.genus - 2 + request.valence > 0

    def provenance(self, request):
        return (
            "Givental-Teleman reconstruction from the genus-zero twisted "
            "quantum connection and q=0 quantum-Riemann-Roch calibration"
        )

    def evaluate(self, request):
        if not isinstance(request, TwistedZeroVertexRequest):
            raise TypeError("request must be a TwistedZeroVertexRequest")
        if 2 * request.genus - 2 + request.valence <= 0:
            raise ValueError("Givental reconstruction needs a stable (g,n)")
        scale, normalized = request.normalized_scales()
        if not scale:
            return self.rings.full_field.zero()
        if normalized not in self._cache:
            data, evaluator = self._evaluator(normalized)
            planned = data.planned_insertions(normalized.insertions)
            term_families = tuple(
                data.ancestor_terms(insertion) for insertion in planned
            )
            series = data.series_ring.zero()
            for inputs in itertools.product(*term_families):
                series += evaluator.evaluate(normalized.genus, inputs)
            self._cache[normalized] = self.rings.full_field(
                series[normalized.degree]
            )
        return self.rings.full(scale * self._cache[normalized])


class HybridTwistedZeroVertexBackend(TwistedZeroVertexBackend):
    r"""Route supported requests to Givental and everything else to fallback.

    ``verify=True`` evaluates supported requests by both methods and raises on
    disagreement.  It is intended for calibration tests, not production runs.
    """

    def __init__(self, fallback, givental=None, verify=False):
        # Sage ``load`` may create a fresh Python class object for the same
        # source module. Use the backend protocol rather than brittle class
        # identity so this wrapper composes with already-loaded localization.
        if not hasattr(fallback, "evaluate") or not hasattr(fallback, "rings"):
            raise TypeError("fallback must implement the zero-vertex protocol")
        self.fallback = fallback
        self.givental = givental or O3DegreeZeroGiventalBackend(
            rings=fallback.rings
        )
        self.verify = bool(verify)
        super().__init__(
            rings=fallback.rings,
            hodge_backend=getattr(fallback, "hodge", None),
        )

    def _supported(self, request):
        if hasattr(self.givental, "supports"):
            return bool(self.givental.supports(request))
        return request.degree == 0

    def is_cached(self, request):
        backend = self.givental if self._supported(request) else self.fallback
        return hasattr(backend, "is_cached") and backend.is_cached(request)

    def provenance(self, request):
        backend = self.givental if self._supported(request) else self.fallback
        return backend.provenance(request)

    def evaluate(self, request):
        if self._supported(request):
            value = self.givental.evaluate(request)
            if self.verify:
                reference = self.fallback.evaluate(request)
                if value != reference:
                    raise ArithmeticError(
                        "Givental/localization mismatch for %s: %s != %s"
                        % (request, value, reference)
                    )
            return value
        return self.fallback.evaluate(request)
