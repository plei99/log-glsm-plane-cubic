class kSignature:
    """
    A signature of a stratum of k-differentials.

    Attributes:
        sig (tuple): signature tuple
        g (int): genus
        n (int): total number of points
        p (int): number of poles
        z (int): number of zeroes
        poles (tuple): tuple of pole orders
        zeroes (tuple): tuple of zero orders
        pole_ind (tuple): tuple of indices of poles
        zero_ind (tuple): tuple of indices of zeroes

    EXAMPLES::

        sage: from admcycles.diffstrata.sig import kSignature
        sage: sig=kSignature((2,1,-1,0), k=1)
        sage: sig.g
        2
        sage: sig.n
        4
        sage: sig.poles
        (-1,)
        sage: sig.zeroes
        (2, 1)
        sage: sig.pole_ind
        (2,)
        sage: sig.zero_ind
        (0, 1)
        sage: sig.p
        1
        sage: sig.z
        2
    """

    def __init__(self, sig, k):
        """
        Initialise signature

        Args:
            sig (tuple): signature tuple of integers adding up to k*(2g-2)
            k (int): order of the differential.
        """
        self.sig = tuple(sig)
        sum_sig = sum(sig)
        if sum_sig % (2 * k) != 0:
            raise ValueError("Error! Illegal signature: Genus not an integer")
        self.g = int(sum_sig / (2 * k)) + 1
        self.n = len(sig)
        self.k = k
        self.poles = tuple(p for p in sig if p < 0)
        self.zeroes = tuple(z for z in sig if z > 0)
        self.marked_points = tuple(k for k in sig if k == 0)
        self.p = len(self.poles)
        self.z = len(self.zeroes)
        self.pole_ind = tuple(i for i, p in enumerate(sig) if p < 0)
        self.zero_ind = tuple(i for i, z in enumerate(sig) if z > 0)

    def __repr__(self):
        if self.k == 1:
            return "Signature(%r)" % (self.sig,)
        else:
            return "Signature(%r, k=%r)" % (self.sig, self.k)

    def __hash__(self):
        return hash((self.sig, self.k))

    def __eq__(self, other):
        try:
            return self.sig == other.sig
        except AttributeError:
            return False


class Signature(kSignature):
    """
    A signature of an abelian stratum.

    Attributes:
        sig (tuple): signature tuple
        g (int): genus
        n (int): total number of points
        p (int): number of poles
        z (int): number of zeroes
        poles (tuple): tuple of pole orders
        zeroes (tuple): tuple of zero orders
        pole_ind (tuple): tuple of indices of poles
        zero_ind (tuple): tuple of indices of zeroes

    EXAMPLES::

        sage: from admcycles.diffstrata.sig import Signature
        sage: sig=Signature((2,1,-1,0))
        sage: sig.g
        2
        sage: sig.n
        4
        sage: sig.poles
        (-1,)
        sage: sig.zeroes
        (2, 1)
        sage: sig.pole_ind
        (2,)
        sage: sig.zero_ind
        (0, 1)
        sage: sig.p
        1
        sage: sig.z
        2
    """

    def __init__(self, sig):
        super().__init__(sig, k=1)


class QuadraticSignature(kSignature):
    """
    A signature of an quadratic stratum.

    Attributes:
        sig (tuple): signature tuple
        difftype (string): either "p" for primitive or "gs" of global square
        g (int): genus
        n (int): total number of points
        p (int): number of poles
        z (int): number of zeroes
        poles (tuple): tuple of pole orders
        zeroes (tuple): tuple of zero orders
        pole_ind (tuple): tuple of indices of poles
        zero_ind (tuple): tuple of indices of zeroes

    EXAMPLES::

        sage: from admcycles.diffstrata.sig import QuadraticSignature
        sage: sig=QuadraticSignature((4,1,-1,0), "p")
        sage: sig.difftype
        'p'
        sage: sig.g
        2
        sage: sig.n
        4
        sage: sig.poles
        (-1,)
        sage: sig.zeroes
        (4, 1)
        sage: sig.pole_ind
        (2,)
        sage: sig.zero_ind
        (0, 1)
        sage: sig.p
        1
        sage: sig.z
        2
    """

    def __init__(self, sig, difftype):
        super().__init__(sig, k=2)
        assert difftype in ["p", "gs"]
        if difftype == "gs" and any(a % 2 != 0 for a in sig):
            raise ValueError(
                "Quadratic differentials with odd-order points can't be global squares.")
        self.difftype = difftype

    def __repr__(self):
        return "QuadraticSignature(%r, %r)" % (self.sig, self.difftype)

    def __hash__(self):
        return hash((self.sig, self.k, self.difftype))
