from collections import Counter

# pylint does not know sage
from sage.structure.sage_object import SageObject  # pylint: disable=import-error
from sage.misc.cachefunc import cached_method  # pylint: disable=import-error

import admcycles.diffstrata.additivegenerator
import admcycles.admcycles


class ELGTautClass (SageObject):
    """
    A Tautological class of a stratum X, i.e. a formal sum of of psi classes on
    EmbeddedLevelGraphs.

    This is encoded by a list of summands.

    Each summand corresponds to an AdditiveGenerator with coefficient.

    Thus an ELGTautClass is a list with entries tuples (coefficient, AdditiveGenerator).

    These can be added, multiplied and reduced (simplified).

    INPUT :

      * X : GeneralisedStratum that we are on
      * psi_list : list of tuples (coefficient, AdditiveGenerator) as
            described above.
      * reduce=True : call self.reduce() on initialisation
    """

    def __init__(self, X, psi_list, reduce=True):
        self._psi_list = psi_list
        self._X = X
        if reduce:
            self.reduce()

    @classmethod
    def from_hash_list(cls, X, hash_list):
        # This does not reduce!
        return cls(X, [(c, X.additive_generator_from_hash(h))
                       for c, h in hash_list], reduce=False)

    @property
    def psi_list(self):
        return self._psi_list

    def __repr__(self):
        return "ELGTautClass(X=%r,psi_list=%r)"\
            % (self._X, self._psi_list)

    def __str__(self):
        str = "Tautological class on %s\n" % self._X
        for coeff, psi in self._psi_list:
            str += "%s * %s + \n" % (coeff, psi)
        return str

    def __eq__(self, other):
        if isinstance(
                other,
                admcycles.diffstrata.additivegenerator.AdditiveGenerator):
            return self == other.as_taut()
        try:
            return self._psi_list == other._psi_list
        except AttributeError:
            return False

    def __add__(self, other):
        # for sum, we need to know how to add '0':
        if other == 0:
            return self
        try:
            if not self._X == other._X:
                return NotImplemented
            new_psis = self._psi_list + other._psi_list
            return ELGTautClass(self._X, new_psis)
        except AttributeError:
            return NotImplemented

    def __iadd__(self, other):
        return self.__add__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return (-1) * self

    def __sub__(self, other):
        return self + (-1) * other

    def __mul__(self, other):
        if 0 == other:
            return 0
        elif self._X.ONE == other:
            return self
        # convert AdditiveGenerators to Tautclasses:
        if isinstance(
                other,
                admcycles.diffstrata.additivegenerator.AdditiveGenerator):
            return self * other.as_taut()
        try:
            # check if other is a tautological class
            other._psi_list
        except AttributeError:
            # attempt scalar multiplication:
            new_psis = [(coeff * other, psi) for coeff, psi in self._psi_list]
            return ELGTautClass(self._X, new_psis, reduce=False)
        if not self._X == other._X:
            return NotImplemented
        else:
            return self._X.intersection(self, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, exponent):
        if exponent == 0:
            return self._X.ONE
        # TODO: quick check for going over top degree?
        prod = self
        for _ in range(1, exponent):
            prod = self * prod
        return prod

    @cached_method
    def is_equidimensional(self):
        """
        Determine if all summands of self have the same degree.

        Note that the empty empty tautological class (ZERO) gives True.

        Returns:
            bool: True if all AdditiveGenerators in self.psi_list are of same degree,
                False otherwise.
        """
        if self.psi_list:
            first_deg = self.psi_list[0][1].degree
            return all(AG.degree == first_deg for _c, AG in self.psi_list)
        return True

    def reduce(self):
        """
        Reduce self.psi_list by combining summands with the same AdditiveGenerator
        and removing those with coefficient 0 or that die for dimension reasons.
        """
        # we use the hash of the AdditiveGenerators to group:
        hash_dict = Counter()
        for c, AG in self._psi_list:
            hash_dict[AG] += c
        self._psi_list = [(c, AG) for AG, c in hash_dict.items()
                          if c != 0 and AG.dim_check()]

    # To evaluate, we go through the AdditiveGenerators and
    # take the (weighted) sum of the AdditiveGenerators.
    def evaluate(self, quiet=True, warnings_only=False,
                 admcycles_output=False):
        """
        Evaluation of self, i.e. cap with fundamental class of X.

        This is the sum of the evaluation of the AdditiveGenerators in psi_list
        (weighted with their coefficients).

        Each AdditiveGenerator is (essentially) the product of its levels,
        each level is (essentially) evaluated by admcycles.

        Args:
            quiet (bool, optional): No output. Defaults to True.
            warnings_only (bool, optional): Only warnings output. Defaults to False.
            admcycles_output (bool, optional): admcycles debugging output. Defaults to False.

        Returns:
            QQ: integral of self on X

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: X=GeneralisedStratum([Signature((0,0))])
            sage: assert (X.xi^2).evaluate() == 0

            sage: X=GeneralisedStratum([Signature((1,1,1,1,-6))])
            sage: assert set([(X.cnb(((i,),0),((i,),0))).evaluate() for i in range(len(X.bics))]) == {-2, -1}
        """
        if warnings_only:
            quiet = True
        DS_list = []
        for c, AG in self.psi_list:
            value = AG.evaluate(
                quiet=quiet,
                warnings_only=warnings_only,
                admcycles_output=admcycles_output)
            DS_list.append(c * value)
        if not quiet:
            print("----------------------------------------------------")
            print("In summary: We sum")
            for i, summand in enumerate(DS_list):
                print("Contribution %r from AdditiveGenerator" % summand)
                print(self.psi_list[i][1])
                print("(With coefficient %r)" % self.psi_list[i][0])
            print("To obtain a total of %r" % sum(DS_list))
            print("----------------------------------------------------")
        return sum(DS_list)

    def extract(self, i):
        """
        Return the i-th component of self.

        Args:
            i (int): index of self._psi_list

        Returns:
            ELGTautClass: coefficient * AdditiveGenerator at position i of self.
        """
        return ELGTautClass(self._X, [self._psi_list[i]], reduce=False)

    @cached_method
    def degree(self, d):
        """
        The degree d part of self.

        Args:
            d (int): degree

        Returns:
            ELGTautClass: degree d part of self
        """
        new_psis = []
        for c, AG in self.psi_list:
            if AG.degree == d:
                new_psis.append((c, AG))
        return ELGTautClass(self._X, new_psis, reduce=False)

    @cached_method
    def list_by_degree(self):
        """
        A list of length X.dim with the degree d part as item d

        Returns:
            list: list of ELGTautClasses with entry i of degree i.
        """
        deg_psi_list = [[] for _ in range(self._X.dim() + 1)]
        for c, AG in self.psi_list:
            deg_psi_list[AG.degree].append((c, AG))
        return [ELGTautClass(self._X, piece, reduce=False)
                for piece in deg_psi_list]

    def is_pure_psi(self):
        """
        Check if self is ZERO or a psi-product on the stratum.

        Returns:
            boolean: True if self has at most one summand and that is of the form
                AdditiveGenerator(((), 0), psis).

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: X=Stratum((2,))
            sage: X.ZERO.is_pure_psi()
            True
            sage: X.ONE.is_pure_psi()
            True
            sage: X.psi(1).is_pure_psi()
            True
            sage: X.xi.is_pure_psi()
            False
        """
        if not self.psi_list:
            return True
        return len(
            self.psi_list) == 1 and self.psi_list[0][1].enh_profile == ((), 0)

    def to_prodtautclass(self):
        """
        Transforms self into an admcycles prodtautclass on the stable graph of the smooth
        graph of self._X.

        Note that this is essentially the pushforward to M_g,n, i.e. we resolve residues
        and multiply with the correct Strataclasses along the way.

        Returns:
            prodtautclass: admcycles prodtautclass corresponding to self pushed forward
                to the stable graph with one vertex.

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: X=Stratum((2,))
            sage: X.ONE.to_prodtautclass()
            Outer graph : [2] [[1]] []
            Vertex 0 :
            Graph :      [2] [[1]] []
            Polynomial : -7/24*(kappa_1)_0 + 79/24*psi_1
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1] [[1, 2, 3]] [(2, 3)]
            Polynomial : -1/48
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1, 1] [[2], [1, 3]] [(2, 3)]
            Polynomial : -19/24
            sage: (X.xi^X.dim()).evaluate() == (X.xi^X.dim()).to_prodtautclass().pushforward().evaluate()
            True
        """
        G = self._X.smooth_LG
        stgraph = G.LG.stgraph
        total = admcycles.admcycles.prodtautclass(stgraph, terms=[])
        for c, AG in self.psi_list:
            ptc = AG.to_prodtautclass()
            # sort vertices by connected component:
            vertex_map = {}
            # note that every vertex of G has at least one leg (that is a
            # marked point of _X):
            for v, _ in enumerate(G.LG.genera):
                mp_on_stratum = G.dmp[G.LG.legs[v][0]]
                # find this marked point on AG:
                l_AG = AG._G.dmp_inv[mp_on_stratum]
                # get the corresponding vertex
                LG = AG._G.LG
                v_AG = LG.vertex(l_AG)
                # we use the underlying graph:
                UG_v = LG.UG_vertex(v_AG)
                for w, g, kind in LG.UG_without_infinity().connected_component_containing_vertex(
                        UG_v):
                    if kind != 'LG':
                        continue
                    vertex_map[w] = v
            # map legs of AG._G to smooth_LG
            # CAREFUL: This goes in the OTHER direction!
            leg_map = {G.dmp_inv[mp]: ldeg for ldeg, mp in AG._G.dmp.items()}
            pf = ptc.partial_pushforward(stgraph, vertex_map, leg_map)
            total += c * pf
        return total
