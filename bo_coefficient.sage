"""
Coefficient extraction for the Bloch--Okounkov n-point determinant.

For the notation and Eisenstein-series normalization of Appendix A in
arXiv:2205.12777, use

    disconnected_coefficient([d1,...,dn])
    connected_coefficient([d1,...,dn])

The second function applies the set-partition Möbius inversion of (A.1).
The lower-level ``bo_coefficient`` uses powers of z directly and normalized
Eisenstein series (constant term 1); it is retained for compatibility.
"""

def _normalized_eisenstein_polynomial(weight, E4, E6):
    """The normalized E_weight, written as a polynomial in E4 and E6."""
    if weight == 4:
        return E4
    if weight == 6:
        return E6
    if weight == 2:
        raise ValueError("E2 is supplied separately")

    # The monomials E4^u E6^v of weight ``weight`` are a basis of M_weight.
    basis = [(u, v) for v in range(weight // 6 + 1)
             for u in range(weight // 4 + 1) if 4*u + 6*v == weight]
    d = len(basis)
    if not d:
        raise ValueError("there is no holomorphic modular form of weight %s" % weight)

    Q = PowerSeriesRing(QQ, 'q', default_prec=d)
    q = Q.gen()

    def eisenstein_qexp(k):
        return Q(1) + sum(-QQ(2*k)/bernoulli(k) *
                          sum(QQ(e)^(k-1) for e in divisors(m)) * q^m
                          for m in range(1, d))

    e4q, e6q, ekq = eisenstein_qexp(4), eisenstein_qexp(6), eisenstein_qexp(weight)
    A = Matrix(QQ, d, d,
               lambda r, c: (e4q^basis[c][0] * e6q^basis[c][1])[r])
    b = vector(QQ, [ekq[r] for r in range(d)])
    coeffs = A.solve_right(b)
    return sum(coeffs[i] * E4^u * E6^v for i, (u, v) in enumerate(basis))


def bo_coefficient(a):
    r"""
    Return [z_1^a_1 ... z_n^a_n] (q)_infinity * RHS of the formula.

    Here theta is normalized by

        theta(z) = z exp(sum_{m>=1} B_(2m) E_(2m) z^(2m)
                            / (2m (2m)!)).

    Consequently the return value lies in QQ[E2,E4,E6].
    """
    a = tuple(ZZ(x) for x in a)
    if not a:
        raise ValueError("a must be nonempty")
    if any(x < 0 for x in a):
        raise ValueError("all entries of a must be non-negative")

    n = len(a)
    total_degree = sum(a)
    # In F(t*z), the desired homogeneous coefficient occurs at t^total_degree.
    # The product of the n theta denominators contributes t^n.
    M = total_degree + n

    R = PolynomialRing(QQ, names=('E2', 'E4', 'E6'))
    E2, E4, E6 = R.gens()
    S = PolynomialRing(R, names=tuple('z%s' % (i+1) for i in range(n)))
    z = S.gens()
    T = PowerSeriesRing(S, 't', default_prec=M+1)
    t = T.gen()

    # theta must be known M+n degrees before differentiation (r <= n).
    theta_precision = M + n
    U = PowerSeriesRing(R, 'x', default_prec=theta_precision+1)
    x = U.gen()
    eisenstein = {2: E2}
    for m in range(2, theta_precision // 2 + 1):
        eisenstein[2*m] = _normalized_eisenstein_polynomial(2*m, E4, E6)
    Hx = R.zero()
    for m in range(1, theta_precision // 2 + 1):
        weight = 2*m
        Ek = eisenstein[weight]
        Hx += bernoulli(2*m) * Ek * x^(2*m) / (2*m * factorial(2*m))
    theta = x * Hx.exp()
    theta_coeff = theta.list()

    def theta_derivative_over_factorial(r, L):
        # theta^(r)(tL)/r! = sum_d binomial(d,r) theta_d (tL)^(d-r).
        if r < 0:
            return T.zero()
        return sum(T(theta_coeff[d] * binomial(d, r) * (t*L)^(d-r))
                   for d in range(r, len(theta_coeff)))

    def H_at(L):
        return sum(T(bernoulli(2*m) *
                     eisenstein[2*m] *
                     (t*L)^(2*m) / (2*m * factorial(2*m)))
                   for m in range(1, M // 2 + 1))

    K = S.fraction_field()
    homogeneous_part = K.zero()
    for p in Permutations(range(n)):
        partial = [S.zero()]
        for i in p:
            partial.append(partial[-1] + z[i])

        entries = []
        for i in range(n):
            entries.append([theta_derivative_over_factorial(j-i+1, partial[n-j-1])
                            for j in range(n)])
        numerator = matrix(T, entries).det()
        inverse_exponent = (-sum(H_at(partial[j]) for j in range(1, n+1))).exp()
        coeff_tM = (numerator * inverse_exponent)[M]
        homogeneous_part += K(coeff_tM) / K(prod(partial[1:]))

    # All non-coordinate poles cancel in the permutation sum.  The remaining
    # pole is exactly z1...zn; multiplying by it gives a polynomial.
    regular_numerator = homogeneous_part * K(prod(z))
    if regular_numerator.denominator() != 1:
        raise ArithmeticError("unexpected uncancelled pole; increase precision")
    # We multiplied by z_1...z_n to remove the coordinate poles above.
    target = prod(z[i]^(a[i] + 1) for i in range(n))
    return regular_numerator.numerator().monomial_coefficient(target)


def _descendant_degrees(d):
    """Validate and freeze a nonempty tuple of descendant degrees."""
    d = tuple(ZZ(x) for x in d)
    if not d:
        raise ValueError("d must be nonempty")
    if any(x < 0 for x in d):
        raise ValueError("all descendant degrees must be non-negative")
    return d


def _in_paper_eisenstein_convention(f):
    r"""
    Convert normalized Eisenstein series to the convention of 2205.12777.

    The paper uses

        E_k = zeta(1-k)/2 + sum_{m>=1} sigma_{k-1}(m) q^m,

    whereas ``bo_coefficient`` internally uses constant-term-one series.
    Thus E2_norm=-24*E2, E4_norm=240*E4, E6_norm=-504*E6.
    """
    P = PolynomialRing(QQ, names=('E2', 'E4', 'E6'))
    E2, E4, E6 = P.gens()
    source = f.parent()
    conversion = source.hom([-24*E2, 240*E4, -504*E6], P)
    return conversion(f)


def disconnected_coefficient(d):
    r"""
    Return the vacuum-stripped disconnected stationary series D_d(q).

    More precisely, for d=(d_1,...,d_n), this returns

        [z_1^(d_1+1) ... z_n^(d_n+1)] (q)_infinity F_E(z_1,...,z_n),

    expressed using the E2,E4,E6 convention of arXiv:2205.12777.
    By (A.1), D_d is the set-partition sum of products of connected
    series C_(d restricted to a block).
    """
    d = _descendant_degrees(d)
    return _in_paper_eisenstein_convention(
        bo_coefficient(tuple(x + 1 for x in d))
    )


def connected_coefficient(d):
    r"""
    Return the connected stationary generating series C_d(q).

    Equation (A.1) says

        D_d = sum_{pi partition of {1,...,n}} prod_{B in pi} C_(d|B).

    Möbius inversion on the lattice of set partitions therefore gives

        C_d = sum_pi (-1)^(|pi|-1) (|pi|-1)!
                         prod_{B in pi} D_(d|B).

    The result is in QQ[E2,E4,E6] of weight sum_i(d_i+2), with the
    Eisenstein-series normalization used in arXiv:2205.12777.
    """
    d = _descendant_degrees(d)
    n = len(d)
    P = PolynomialRing(QQ, names=('E2', 'E4', 'E6'))

    # F_E is symmetric, so sorting a block gives a useful cache key.
    disconnected_cache = {}

    def D(block_degrees):
        key = tuple(sorted(block_degrees))
        if key not in disconnected_cache:
            disconnected_cache[key] = disconnected_coefficient(key)
        return disconnected_cache[key]

    answer = P.zero()
    for partition in SetPartitions(n):
        blocks = [tuple(d[i-1] for i in block) for block in partition]
        number_of_blocks = len(blocks)
        mobius = (-1)^(number_of_blocks - 1) * factorial(number_of_blocks - 1)
        answer += mobius * prod(D(block) for block in blocks)
    return answer
