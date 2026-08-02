import sys

from itertools import product

from sage.rings.rational_field import QQ  # pylint: disable=import-error
from sage.arith.functions import lcm  # pylint: disable=import-error
from sage.arith.misc import gcd  # pylint: disable=import-error
from sage.symbolic.ring import SR  # pylint: disable=import-error
from sage.matrix.constructor import Matrix  # pylint: disable=import-error
from sage.modules.free_module_element import free_module_element  # pylint: disable=import-error
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors  # pylint: disable=import-error
from sage.misc.cachefunc import cached_method  # pylint: disable=import-error

from admcycles.diffstrata.generalisedstratum import GeneralisedStratum, Stratum
from admcycles.diffstrata.sig import Signature

import admcycles.admcycles


def test_calL(sig):
    """
    Compare calL and the 'error term' of cnb.

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import test_calL
        sage: test_calL((1,1,1,1,-6))
        sage: test_calL((4,))
    """
    X = GeneralisedStratum([Signature(sig)])
    for i, B in enumerate(X.bics):
        ll = lcm(B.LG.prongs.values())
        assert X.calL(((i,), 0), 0) + ll * X.cnb(((i,), 0), ((i,), 0)) + X.gen_pullback_taut(X.xi_at_level(0, ((i,), 0)),
                                                                                             ((i,), 0), ((i,), 0)) + (-1) * X.gen_pullback_taut(X.xi_at_level(1, ((i,), 0)), ((i,), 0), ((i,), 0)) == X.ZERO


def meromorphic_tests():
    """
    EXAMPLES::

        sage: from admcycles.diffstrata import *
        sage: X=GeneralisedStratum([Signature((1,1,1,1,-6))])
        sage: (X.xi^2).evaluate(quiet=True)
        25

        sage: X=GeneralisedStratum([Signature((-2,-2,-2,-2,6))])
        sage: (X.xi^X.dim()).evaluate(quiet=True)
        30

        Testing normal bundles:

        sage: X=GeneralisedStratum([Signature((2,2,-2))])
        sage: td0 = X.taut_from_graph((0,))
        sage: td1 = X.taut_from_graph((1,))
        sage: td4 = X.taut_from_graph((4,))
        sage: td5 = X.taut_from_graph((5,))
        sage: td8 = X.taut_from_graph((8,))
        sage: assert X.ZERO == td0^2*td1 - td0*td1*td0  # not 'safe'  # doctest:+SKIP
        sage: assert X.ZERO == td0^2*td5 - td0*td5*td0  # not 'safe'  # doctest:+SKIP
        sage: assert X.ZERO == td0^2*td8 - td0*td8*td0  # not 'safe'  # doctest:+SKIP
        sage: assert (td8^3*td4).evaluate(quiet=True) == (td8^2*td4*td8).evaluate(quiet=True) == (td8*td4*td8^2).evaluate(quiet=True)

    """
    pass


def gengenustwotests():
    """
    EXAMPLES::

        sage: from admcycles.diffstrata import *
        sage: X=GeneralisedStratum([Signature((2,))])
        sage: assert (X.xi^X.dim()).evaluate(quiet=True) == -1/640
        sage: ct=[i for i, B in enumerate(X.bics) if len(B.LG.edges) == 1]
        sage: ct_index = ct[0]  # compact type graph
        sage: ct_taut = X.taut_from_graph((ct_index,))
        sage: assert (X.c2_E*ct_taut).evaluate(quiet=True) == 1/48
        sage: banana_index = 1 - ct_index  # only 2 BICs!
        sage: banana_taut = X.taut_from_graph((banana_index,))
        sage: assert (X.c2_E*banana_taut).evaluate(quiet=True) == 1/24
        sage: assert (X.c1_E*banana_taut**2).evaluate(quiet=True) == -1/16
        sage: assert (X.c1_E*banana_taut*ct_taut).evaluate(quiet=True) == 1/24
        sage: assert (X.c1_E*ct_taut**2).evaluate(quiet=True) == -1/48
        sage: assert (ct_taut**3).evaluate(quiet=True) == 1/96
        sage: assert (banana_taut**3).evaluate(quiet=True) == 1/48
        sage: assert (ct_taut*banana_taut**2).evaluate(quiet=True) == 0
        sage: assert (ct_taut*ct_taut*banana_taut).evaluate(quiet=True) == -1/48

        sage: X=GeneralisedStratum([Signature((1,1))])
        sage: at1 = X.taut_from_graph((1,))
        sage: at2 = X.taut_from_graph((2,))
        sage: X.ZERO == at1^2*at2 + (-1)*at1*at2*at1
        True
        sage: X.ZERO == (-1)*(at1+at2)^3*at2 + at1^3*at2 + 3*at1^2*at2*at2  + 3*at1*at2^2*at2 + at2^3*at2
        True
        sage: X.ZERO == (-1)*(at1+at2)^4 + at1^4 + 4*at1^3*at2 + 6*at1^2*at2^2 + 4*at1^1*at2^3 + at2^4
        True
        sage: (X.xi^X.dim()).evaluate(quiet=True)
        0
        sage: psi1 = AdditiveGenerator(X,((),0),{1:1}).as_taut()
        sage: (X.xi^3*psi1).evaluate(quiet=True)
        -1/360
    """
    pass


def genusthreetests():
    """
    Testing Normal bundles in min stratum (4):

    EXAMPLES::

        sage: from admcycles.diffstrata import *
        sage: X=GeneralisedStratum([Signature((4,))])

        sage: v_banana = [i for i, B in enumerate(X.bics) if B.LG.genera == [1,1,0]]
        sage: v_banana_taut=X.taut_from_graph((v_banana[0],))
        sage: g1_banana = [i for i, B in enumerate(X.bics) if B.LG.genera == [1,1]]
        sage: g1_banana_taut=X.taut_from_graph((g1_banana[0],))

        sage: assert g1_banana_taut**2*v_banana_taut + (-1)*g1_banana_taut*v_banana_taut*g1_banana_taut == X.ZERO
        sage: assert v_banana_taut**2*g1_banana_taut + (-1)*v_banana_taut*g1_banana_taut*v_banana_taut == X.ZERO

        sage: (X.xi_with_leg(quiet=True)^X.dim()).evaluate(quiet=True)  # long time  # optional - local
        305/580608
    """
    pass


def genusfourtests():
    """
    Tests in H(6): (think of smart NB test...)

    EXAMPLES::

        sage: from admcycles.diffstrata import *

    """
    pass


class BananaSuite:
    """
    A frontend for the Stratum (k, 1, -k-1).

    This class models the situation of sec 10.2 of CMZ20
    with a method D for accessing the divisors in the
    notation of loc. cit. (either D(i), i=2,3,4 or D(i,a)
    for i=1,5).
    """

    def __init__(self, k):
        """
        Initialise the genus 1 stratum (k, 1, -k-1).

        Args:
            k (int): order of zero
        """
        self._k = k
        self._X = GeneralisedStratum([Signature((k, 1, -k - 1))])

    def D(self, i, a=1):
        """
        The divisor using the notation of Sec. 10.2.

        Args:
            i (int): index 1,2,3,4,5
            a (int, optional): prong for i=1 or 5. Defaults to 1.

        Returns:
            ELGTautClass: Tautological class of the divisor D_{i,a}.
        """
        if i == 1:
            for b, B in enumerate(self._X.bics):
                v = B.LG.verticesonlevel(0)[0]
                if (B.LG.genus(v) == 0 and
                    len(B.LG.ordersonvertex(v)) == 3 and
                    len(B.LG.edges) == 2 and
                        a in B.LG.prongs.values()):
                    assert -self._k - 1 in B.LG.ordersonvertex(v),\
                        "%r" % B
                    return self._X.taut_from_graph((b,))
        elif i == 2:
            for b, B in enumerate(self._X.bics):
                v = B.LG.verticesonlevel(0)[0]
                if (B.LG.genus(v) == 1 and len(B.LG.ordersonvertex(v)) == 2):
                    assert -self._k - 1 in B.LG.ordersonvertex(v),\
                        "%r" % B
                    return self._X.taut_from_graph((b,))
        elif i == 3:
            for b, B in enumerate(self._X.bics):
                v = B.LG.verticesonlevel(0)[0]
                if (B.LG.genus(v) == 1 and len(B.LG.ordersonvertex(v)) == 1):
                    return self._X.taut_from_graph((b,))
        elif i == 4:
            for b, B in enumerate(self._X.bics):
                v = B.LG.verticesonlevel(-1)[0]
                if (B.LG.genus(v) == 1 and len(B.LG.ordersonvertex(v)) == 2):
                    return self._X.taut_from_graph((b,))
        elif i == 5:
            for b, B in enumerate(self._X.bics):
                v = B.LG.verticesonlevel(0)[0]
                if (B.LG.genus(v) == 0 and
                    len(B.LG.ordersonvertex(v)) == 4 and
                    len(B.LG.edges) == 2 and
                        a in B.LG.prongs.values()):
                    assert -self._k - 1 in B.LG.ordersonvertex(v),\
                        "%r" % B
                    return self._X.taut_from_graph((b,))
        else:
            return None

    def check(self, quiet=False):
        """
        Check Prop 10.1

        Args:
            quiet (bool, optional): No output. Defaults to False.

        Returns:
            bool: Should always return True.
        """
        check_list = []

        def delta(k, a):
            if a == k / 2:
                return QQ(1) / QQ(2)
            else:
                return 1
        for a in range(1, self._k + 1):
            si = (self.D(1, a)**2).evaluate(quiet=True)
            rhs = -QQ(delta(self._k + 1, a) * self._k *
                      gcd(a, self._k + 1 - a)) / QQ(lcm(a, self._k + 1 - a))
            check_list.append(si == rhs)
            if not quiet:
                print("D(1,%r)^2 = %r, RHS = %r" % (a, si, rhs))
        for a in range(1, self._k):
            si = (self.D(5, a)**2).evaluate(quiet=True)
            rhs = -QQ(delta(self._k, a) * (self._k + 1) *
                      gcd(a, self._k - a)) / QQ(lcm(a, self._k - a))
            check_list.append(si == rhs)
            if not quiet:
                print("D(5,%r)^2 = %r, RHS = %r" % (a, si, rhs))
        return all(check_list)

    def banana_tests(self):
        """
        EXAMPLES::

            sage: from admcycles.diffstrata.tests import BananaSuite
            sage: B=BananaSuite(2)
            sage: B.check(quiet=True)
            True
            sage: (B._X.xi**2).evaluate(quiet=True) == QQ(2**4 - 1)/QQ(24)
            True
            sage: assert((B._X.xi_with_leg(quiet=True)*B.D(5,1)).evaluate(quiet=True) == \
                          B._X.xi_at_level(0,B.D(5,1).psi_list[0][1].enh_profile,quiet=True).evaluate(quiet=True) == -1)
            sage: assert(B._X.xi_at_level(0,B.D(5,1).psi_list[0][1].enh_profile,leg=3,quiet=True).evaluate(quiet=True) == -1)

            sage: B=BananaSuite(5)
            sage: B.check(quiet=True)
            True
            sage: B=BananaSuite(5)
            sage: (B._X.xi**2).evaluate(quiet=True) == QQ(5**4 - 1)/QQ(24)
            True

            sage: B=BananaSuite(6)
            sage: assert((B._X.xi_with_leg(quiet=True)*B.D(5,2)).evaluate(quiet=True) == \
                          B._X.xi_at_level(0,B.D(5,2).psi_list[0][1].enh_profile,leg=1,quiet=True).evaluate(quiet=True) == \
                          B._X.xi_at_level(0,B.D(5,2).psi_list[0][1].enh_profile,leg=2,quiet=True).evaluate(quiet=True) == \
                          B._X.xi_at_level(0,B.D(5,2).psi_list[0][1].enh_profile,leg=3,quiet=True).evaluate(quiet=True) == \
                          B._X.xi_at_level(0,B.D(5,2).psi_list[0][1].enh_profile,leg=4,quiet=True).evaluate(quiet=True) == -12)

            sage: B=BananaSuite(10)
            sage: B.check(quiet=True)
            True
            sage: (B._X.xi**2).evaluate(quiet=True) == QQ(10**4 - 1)/QQ(24)
            True
        """
        pass


def rGRCtests():
    """
    Test the surjectivity of the _to_bic maps.

    EXAMPLES::

        sage: from admcycles.diffstrata import *
        sage: X=GeneralisedStratum([Signature((2,2,-2))])
        sage: assert all(set(X.DG.top_to_bic(i).keys()) == set(range(len(X.bics[i].level(0).bics))) for i in range(len(X.bics)))
        sage: assert all(set(X.DG.bot_to_bic(i).keys()) == set(range(len(X.bics[i].level(1).bics))) for i in range(len(X.bics)))
    """
    pass


def middle_level_degeneration(sig):
    """
    Check if gluing in middle bics gives (at least as a set) all length three profiles.

    Maybe think of some more sophisticated test...

    Args:
        sig (tuple): Signature tuple.

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import middle_level_degeneration
        sage: middle_level_degeneration((1,1))
        sage: middle_level_degeneration((2,2,-2))
        sage: middle_level_degeneration((4,))
        sage: middle_level_degeneration((2,2,2,-6))
    """
    X = GeneralisedStratum([Signature(sig)])
    three_level_graphs = X.enhanced_profiles_of_length(2)
    four_level_profiles_set = set(X.lookup_list[3])
    seen = set()
    for ep in three_level_graphs:
        p, _i = ep
        for b in X.DG.middle_to_bic(ep).values():
            seen.add((p[0], b, p[1]))
    assert seen == four_level_profiles_set


def leg_test(sig, quiet=False):
    """
    Tests on each dimension 1 graphs of the stratum with signature sig:
    We test on each one-dimensional level if the evaluation of the xi glued in at
    that level is the same (for every leg!) as the product of the graph with xi on
    the whole stratum.

    Args:
        sig (tuple): Signature of a stratum.
        quiet (bool, optional): No output. Defaults to False.

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import leg_test
        sage: leg_test((6,1,-7),quiet=True)
        sage: leg_test((3,-1,-2),quiet=True)
        sage: leg_test((1,1),quiet=True)
        sage: leg_test((2,2,-2),quiet=True)  # long time
    """
    X = GeneralisedStratum([Signature(sig)])
    d = X.dim() - 1  # codim
    for p in X.lookup_list[d]:
        components = X.lookup(p)
        for i, B in enumerate(components):
            enh_profile = (p, i)
            global_xi_prod = (X.xi_with_leg(quiet=True) *
                              X.taut_from_graph(p, i)).evaluate(quiet=True)
            top_dim = X.lookup_graph(*enh_profile).level(0).dim()
            if not quiet:
                print("Graph %r: xi evaluated: %r (dim of Level 0: %r)" %
                      (enh_profile, global_xi_prod, top_dim))
            if top_dim == 0:
                assert global_xi_prod == 0
            for l in range(d):
                L = B.level(l)
                if L.dim() != 1:
                    continue
                first = None
                for leg in L.leg_dict:
                    level_xi_prod = X.xi_at_level(
                        l, (p, i), leg=leg, quiet=True).evaluate(quiet=True)
                    if not first:
                        first = level_xi_prod
                    if not quiet:
                        print("level: %r, leg: %r, xi ev: %r" %
                              (l, leg, level_xi_prod))
                    if quiet:
                        assert first == level_xi_prod
                        if l == 0:
                            assert global_xi_prod == level_xi_prod


def stratumclasstests():
    """
    Tests of stratum class calculations.

    EXAMPLES::

        sage: from admcycles.diffstrata import *

        sage: X=GeneralisedStratum([Signature((23,5,-13,-17))])
        sage: assert X.res_stratum_class([(0,2)]).evaluate(quiet=True) == 5
    """
    pass


def commutativity_check(sig):
    """
    Run a (large) commutativity check on Stratum(sig)
    to check the normal bundle.

    More precisely, we check all top-degree products
    of BICs in this stratum, multiplying them from
    right-to-left and from left-to-right and checking
    that the evaluations agree.

    Args:
        sig (tuple): signature tuple

    Raises:
        RuntimeError: raised if a commutator doesn't
            evaluate to 0.
    """
    X = GeneralisedStratum([Signature(sig)])
    n = X.dim()
    for T in product(range(len(X.bics)), repeat=n):
        print("Starting IPs")
        print(T)
        PR = X.taut_from_graph((T[0],))
        RP = X.taut_from_graph((T[-1],))
        for i in range(1, n):
            PR *= X.taut_from_graph((T[i],))
            RP *= X.taut_from_graph((T[-1 - i],))
        print(T[0], T[1])
        RP = RP.evaluate(quiet=True)
        PR = PR.evaluate(quiet=True)
        if PR - RP != 0:
            print(T, " gives ", PR - RP)
            raise RuntimeError


def c2_banana_tests(k):
    """
    Check c2 evaluation for the bananas.

    Args:
        k (int): zero of the banana

    Returns:
        RationalNumber: difference of evaluate of c2 and c2
            formula in terms of c1 and ch2 (should be 0!!!)

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import c2_banana_tests
        sage: for k in range(1,10): assert c2_banana_tests(k) == 0
    """
    X = GeneralisedStratum([Signature((-k - 1, k, 1))])
    assert (X.c2_E).evaluate(quiet=True) == QQ(k * (k + 1)) / QQ(6)
    return (X.c2_E).evaluate(quiet=True) - QQ(1) / QQ(2) * \
        (X.c1_E**2 + (-2) * X.ch2_E).evaluate(quiet=True)


def c2_test(sig):
    """
    Compare c2_E to (ch_1^2 - 2ch_2)/2.

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import c2_test
        sage: assert c2_test((1,1,1,1,-6)) == 2
        sage: assert c2_test((-2,-2,-2,-2,6)) == 2
    """
    X = GeneralisedStratum([Signature(sig)])
    c2ev = (X.c2_E).evaluate(quiet=True)
    diff = c2ev - QQ(1) / QQ(2) * (X.c1_E**2 + (-2) * X.ch2_E).evaluate(quiet=True)
    assert diff == 0
    return c2ev


def c1_test(k):
    """
    Check the Euler characteristic of (k,-k) (modular curves).

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import c1_test
        sage: for k in range(2,20): assert c1_test(k)
    """
    X = GeneralisedStratum([Signature((k, -k))])
    return (X.c1_E).evaluate(quiet=True) == QQ(k * k - 1) / QQ(12)


def ch1_pow_test(sig_tuple, deg=3):
    """
    Compare ch1_pow(deg) to c1^deg in Stratum(sig_tuple).

    We multiply with classes to obtain top-degree and then
    evaluate.

    This should produce a series of 0s.

    Args:
        sig_tuple (tuple): signature tuple
        deg (int, optional): degree. Defaults to 3.

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import ch1_pow_test
        sage: ch1_pow_test((1,1), 2)
        Calculating difference...
        Products of BICs:
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        Products with graphs of codim 2:
        0 0 0 0
    """
    X = GeneralisedStratum([Signature(sig_tuple)])
    codim = X.dim() - deg
    print('Calculating difference...')
    c = X.chern_poly(upto=1)
    diff = c[1]**deg - X.ch1_pow(deg)
    if codim == 0:
        print(diff.evaluate(quiet=True))
        return
    print('Products of BICs:')
    for pr in product(range(len(X.bics)), repeat=codim):
        expr = diff
        for b in pr:
            expr *= X.taut_from_graph((b,))
        ev = expr.evaluate(quiet=True)
        print(ev, end=' ')
        sys.stdout.flush()
    print()
    if codim > 1:
        print('Products with graphs of codim %r:' % codim)
        for ep in X.enhanced_profiles_of_length(codim):
            print((diff * X.taut_from_graph(*ep)).evaluate(quiet=True), end=' ')
            sys.stdout.flush()
    print()


def chern_poly_test(sig):
    """
    Compare chern_poly and chern_class.

    We multiply with classes to obtain top-degree and then
    evaluate.

    Args:
        sig (tuple): signature tuple

    Returns:
        bool: True

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import chern_poly_test
        sage: chern_poly_test((2,))  # doctest:+ELLIPSIS
        Calculating Chern Polynomial via character...
        Calculating difference...
        Comparing top classes: 0
        Comparing Products of 1 BICs and psi classes:
        BIC 0 0 BIC 1 0 Psi 1 0
        Comparing Products of 2 BICs and psi classes:
        BIC 0 BIC 0 0 BIC 0 BIC 1 0 BIC 0 Psi 1 0 BIC 1 BIC 0 0 BIC 1 BIC 1 0 BIC 1 Psi 1 0 Psi 1 BIC 0 0 Psi 1 BIC 1 0 Psi 1 Psi 1 0
        Products with graphs of codim 2:
        Profile (..., 0) 0
        All tests passed: True
        True
    """
    X = GeneralisedStratum([Signature(sig)])
    print('Calculating Chern Polynomial via character...')
    c = X.chern_poly()
    print('Calculating difference...')
    diff = [c[i] - X.chern_class(i) for i in reversed(range(1, X.dim() + 1))]
    evs = []
    for codim in range(X.dim()):
        assert diff[codim].is_equidimensional()
        if diff[codim].psi_list:
            assert diff[codim].psi_list[0][1].degree == X.dim() - codim
        if codim == 0:
            print("Comparing top classes:", end=' ')
            sys.stdout.flush()
            ev = diff[codim].evaluate(quiet=True)
            print(ev)
            evs.append(ev)
            continue
        print('Comparing Products of %r BICs and psi classes:' % codim)
        for pr in product(range(len(X.bics) + len(sig)), repeat=codim):
            expr = diff[codim]
            for b in pr:
                if b < len(X.bics):
                    print('BIC %r' % b, end=' ')
                    expr *= X.taut_from_graph((b,))
                else:
                    psi_num = b - len(X.bics) + 1
                    print('Psi %r' % psi_num, end=' ')
                    expr *= X.psi(psi_num)
            ev = expr.evaluate(quiet=True)
            print(ev, end=' ')
            evs.append(ev)
            sys.stdout.flush()
        print()
        if codim > 1:
            print('Products with graphs of codim %r:' % codim)
            for ep in X.enhanced_profiles_of_length(codim):
                print('Profile %r' % (ep,), end=' ')
                ev = (diff[codim] * X.taut_from_graph(*ep)).evaluate(quiet=True)
                print(ev, end=' ')
                sys.stdout.flush()
                evs.append(ev)
            print()
    passed = all(ev == 0 for ev in evs)
    print('All tests passed:', passed)
    return passed


def chern_poly_tests(sig_list):
    """
    Apply chern_poly_test to a list of signatures.

    Args:
        sig_list (iterable): list of signature tuples.

    EXAMPLES::

        sage: from admcycles.diffstrata.tests import chern_poly_tests
        sage: chern_poly_tests([(0,),(2,-2)])
        Entering Stratum (0,)
        Calculating Chern Polynomial via character...
        Calculating difference...
        Comparing top classes: 0
        All tests passed: True
        Done.
        Entering Stratum (2, -2)
        Calculating Chern Polynomial via character...
        Calculating difference...
        Comparing top classes: 0
        All tests passed: True
        Done.
        All strata tests passed: True
    """
    check_vec = []
    for sig in sig_list:
        print('Entering Stratum', sig)
        check_vec.append(chern_poly_test(sig))
        print('Done.')
    print('All strata tests passed:', all(check_vec))


class C3_coefficient_hunter:
    """
    A class that illustrates how to use symbolic variables to
    test coefficients in explicit formulas of c_k.
    """

    def __init__(self, sig=None, fct=None):
        self.NUM_VARS = 33
        self.DEG = 3
        if fct is None:
            self.fct = 'c3_E'
        else:
            self.fct = fct
        self.t = t = SR.var('t', self.NUM_VARS)
        self.var_list = list(t)
        self.M = None
        self.constants = []
        if not (sig is None):
            self.add_stratum(sig)

    def add_stratum(self, sig):
        X = GeneralisedStratum([Signature(sig)])
        print('Calculating difference...')
        expr = getattr(X, self.fct)() - X.chern_poly(upto=self.DEG)[self.DEG]
        codim = X.dim() - self.DEG
        if codim == 0:
            self._add_eq(expr.evaluate(quiet=True))
            return
        print('Products of BICs and Psis:')
        for pr in product(range(len(X.bics) + len(sig)), repeat=codim):
            diff = expr
            for b in pr:
                if b < len(X.bics):
                    print('BIC %r' % b, end=' ')
                    diff *= X.taut_from_graph((b,))
                else:
                    psi_num = b - len(X.bics) + 1
                    print('Psi %r' % psi_num, end=' ')
                    diff *= X.psi(psi_num)
            ev = diff.evaluate(quiet=True)
            self._add_eq(ev)
        if codim > 1:
            print('Products with graphs of codim %r:' % codim)
            for ep in X.enhanced_profiles_of_length(codim):
                print('Profile %r' % (ep,), end=' ')
                self._add_eq(
                    (expr * X.taut_from_graph(*ep)).evaluate(quiet=True))

    def _add_eq(self, expr):
        print("Adding equation: %r" % expr)
        if expr == 0:
            return
        eqn = [expr.coefficient(v) for v in self.var_list]
        cst = expr.substitute({v: 0 for v in expr.free_variables()})
        if self.M is None:
            self.M = Matrix(QQ, [eqn])
        else:
            self.M = self.M.stack(Matrix(QQ, [eqn]))
        self.constants.append(cst)

    def solve(self):
        self.solution = self.M.solve_right(
            free_module_element(QQ, self.constants))

    def __str__(self):
        s = ["Coefficient Matrix:\n"]
        s.append(str(self.M))
        s.append('\nRank: %r' % self.M.rank())
        s.append("\nConstants:\n")
        s.append(str(self.constants))
        self.solve()
        s.append('\nSolution:\n')
        s.append(str(self.solution))
        return ''.join(s)


class IntersectionMatrix:
    """
    The intersection matrix of a stratum.
    """

    def __init__(self, sig):
        """
        Initialise stratum and cache.

        Args:
            sig (tuple): signature tuple
        """
        self.X = Stratum(sig)
        self.info_vecs = {}
        self.codim_one_summary()

    @cached_method
    def codim_xis(self, k):
        """
        Classes of codimension k, that graphs with <= k levels
        with xi powers on levels to reach deg k.

        Args:
            k (int): degree

        Returns:
            list: list of ELGTautClasses of deg k
        """
        print('Calculating classes of codim %r...' % k, end=' ')
        sys.stdout.flush()
        classes = []
        info_vec = []
        for l in range(k + 1):
            xi_deg = k - l
            for ep in self.X.enhanced_profiles_of_length(l):
                AG = self.X.additive_generator(ep)
                # distribute xi powers:
                if xi_deg == 0:
                    # no xis to distribute
                    classes.append(AG.as_taut())
                    info_vec.append((ep, {}))
                    continue
                level_dims = [AG.level_dim(l) for l in range(AG.codim + 1)]
                # we number the positive-dimensional levels:
                pos_dim_level_dict = {}
                i = 0
                for l, dim in enumerate(level_dims):
                    if dim > 0:
                        pos_dim_level_dict[i] = l
                        i += 1
                # not the most efficient way, but good enough for now:
                for exponents in WeightedIntegerVectors(
                        xi_deg, [1] * len(pos_dim_level_dict)):
                    if any(exponents[i] > level_dims[pos_dim_level_dict[i]]
                           for i in range(len(pos_dim_level_dict))):
                        continue
                    prod = AG.as_taut()
                    for i, e in enumerate(exponents):
                        curr_xi = self.X.xi_at_level_pow(
                            pos_dim_level_dict[i], ep, e)
                        prod = self.X.intersection(prod, curr_xi, ep)
                    classes.append(prod)
                    info_vec.append(
                        (ep, {l: exponents[i] for i, l in pos_dim_level_dict.items()}))
        print('%r found' % len(classes))
        self.info_vecs[k] = info_vec
        return classes

    @cached_method
    def int_matrix(self, k=1):
        """
        Matrix of evaluations of top-degree classes (products of codim_xis(k)
        and codim_xis(dim-k)).

        Args:
            k (int, optional): degree. Defaults to 1.

        Returns:
            list: list of lists of rational numbers.
        """
        x_classes = self.codim_xis(k)
        y_classes = self.codim_xis(self.X.dim() - k)
        M = [[self.X.ZERO for _ in range(len(y_classes))]
             for _ in range(len(x_classes))]
        print('Calculating intersection matrix for codim %r...' % k, end=' ')
        sys.stdout.flush()
        for i, x in enumerate(x_classes):
            for j, y in enumerate(y_classes):
                prod = x * y
                assert prod.is_equidimensional()
                assert prod == self.X.ZERO or prod.psi_list[0][1].degree == self.X.dim(
                )
                M[i][j] = (prod).evaluate()
        print('Done!')
        return M

    def codim_one_summary(self):
        """
        Display codim 1 matrix summary.
        """
        print(self.X)
        M = self.int_matrix(1)
        rk = Matrix(M).rank()
        print('Codim 1: Rank of %r x %r matrix: %r' % (len(M), len(M[0]), rk))

    def summary(self):
        """
        Summary of matrices in all codimensions.
        """
        print(self.X)
        for k in range(self.X.dim() + 1):
            M = self.int_matrix(k)
            rk = Matrix(M).rank()
            print('Codim %r: Rank of %r x %r matrix: %r' %
                  (k, len(M), len(M[0]), rk))

    def print_matrix(self, k=1):
        """
        Human readable output of int_matrix(k) values.

        Args:
            k (int, optional): degree. Defaults to 1.
        """
        for row in range(len(self.int_matrix(k))):
            self.print_row(row + 1, k)

    def info(self, i, k=1):
        """
        A string representation of the i-th codim k class.

        Args:
            i (int): index of codim_xis(k)
            k (int, optional): degree. Defaults to 1.

        Returns:
            String: string description of which levels carry xis.
        """
        if k not in self.info_vecs:
            self.codim_xis(k)
        ep, d = self.info_vecs[k][i]
        s = str(ep)
        for l, e in d.items():
            s += ' level %r: xi^%r,' % (l, e)
        return s

    def entry(self, row, col, k=1):
        """
        Human-readable representation of the entry row, col of int_matrix(k).

        Note that this prints the classes being multiplied, not the values.

        Args:
            row (int): row index (math notation, i.e. starting at 1)
            col (int): col index (math notation, i.e. starting at 1)
            k (int, optional): degree. Defaults to 1.

        Returns:
            String: info of the two factors at this entry
        """
        # math notation, i.e. starting at 1!!
        return self.info(row - 1, k) + ' * ' + self.info(col - 1, self.X.dim() - k)

    def print_row(self, row, k=1):
        """
        Human-readable output of the row of int_matrix(k)

        Args:
            row (int): row (math notation, i.e. starting at 1)
            k (int, optional): degree. Defaults to 1.
        """
        # use math notation, i.e. starting at 1!!
        M = self.int_matrix(k)
        print("Row %r" % row)
        for col, value in enumerate(M[row - 1]):
            print('Col %r: %s value: %r' %
                  (col + 1, self.entry(row, col + 1, k), value))

    def print_col(self, col, k=1):
        """
        Human-readable output of the col of int_matrix(k)

        Args:
            col (int): column (math notation, i.e. starting at 1)
            k (int, optional): degree. Defaults to 1.
        """
        # use math notation, i.e. starting at 1!!
        M = self.int_matrix(k)
        print("Col %r" % col)
        for row in range(len(M)):
            print('Row %r: %s value: %r' %
                  (row + 1, self.entry(row + 1, col, k), M[row][col - 1]))

#########################################


def push_pull(sig, T_name=None):
    print('Testing non-horizontal push-pull for stratum %r' % (sig,), end=' ')
    sys.stdout.flush()
    X = Stratum(sig)
    if T_name == 'xi':
        print('using xi:')
        T = X.xi
    else:
        T = X.ONE
        print('using [1]:')
    g = X._g[0]
    n = X._n
    for G in admcycles.admcycles.list_strata(g, n, 1):
        if len(G.genera()) == 1:
            # horizontal
            continue
        print("%r: " % G, end='')
        sys.stdout.flush()
        pullback = X.boundary_pullback(G)
        if pullback == X.ZERO:
            print('ZERO')
        else:
            LHS = (T * pullback).to_prodtautclass().pushforward()
            RHS = admcycles.admcycles.tautclass(
                [admcycles.admcycles.decstratum(G)]) * T.to_prodtautclass().pushforward()
            print((LHS - RHS).is_zero())
