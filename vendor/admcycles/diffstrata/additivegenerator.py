
# pylint does not know sage
from sage.structure.sage_object import SageObject  # pylint: disable=import-error
from sage.misc.cachefunc import cached_method  # pylint: disable=import-error
from sage.rings.rational_field import QQ  # pylint: disable=import-error

import admcycles.admcycles
import admcycles.stratarecursion

import admcycles.diffstrata.elgtautclass

from admcycles.diffstrata.auxiliary import hash_AG


class AdditiveGenerator (SageObject):
    """
    Product of Psi classes on an EmbeddedLevelGraph (of a stratum X).

    The information of a product of psi-class on an EmbeddedLevelGraph, i.e. a
    leg_dict and an enhanced_profile, where leg_dict is a dictionary on the legs
    leg -> exponent of the LevelGraph associated to the enhanced profile, i.e.
    (profile,index) or None if we refer to the class of the graph.

    We (implicitly) work inside some stratum X, where the enhanced profile
    makes sense.

    This class should be considered constant (hashable)!
    """

    def __init__(self, X, enh_profile, leg_dict=None):
        """
        AdditiveGenerator for psi polynomial given by leg_dict on graph
        corresponding to enh_profile in X.

        Args:
            X (GeneralisedStratum): enveloping stratum

            enh_profile (tuple): enhanced profile (in X)

            leg_dict (dict, optional): dictionary leg of enh_profile -> exponent
            encoding a psi monomial. Defaults to None.
        """
        self._X = X
        self._hash = hash_AG(leg_dict, enh_profile)
        self._enh_profile = (tuple(enh_profile[0]), enh_profile[1])
        self._leg_dict = leg_dict
        self._G = self._X.lookup_graph(*enh_profile)
        # dictionary leg -> level
        # Careful! These are leg numbers on the whole graph, not on
        # the graphs inside the LevelStrata!!
        self._level_dict = {}
        if leg_dict is not None:
            for l in leg_dict:
                self._level_dict[l] = self._G.LG.level_number(
                    self._G.LG.levelofleg(l))
        self._inv_level_dict = {}
        for leg in self._level_dict:
            try:
                self._inv_level_dict[self._level_dict[leg]].append(leg)
            except KeyError:
                self._inv_level_dict[self._level_dict[leg]] = [leg]

    @classmethod
    def from_hash(cls, X, hash):
        """
        AdditiveGenerator from a hash generated with hash_AG.

        Args:
            X (GeneralisedStratum): Enveloping stratum.
            hash (tuple): hash from hash_AG

        Returns:
            AdditiveGenerator: AG from hash.
        """
        if hash[0] is None:
            leg_dict = None
        else:
            leg_dict = dict(hash[0])
        return cls(X, (hash[1], hash[2]), leg_dict)

    def __hash__(self):
        return hash(self._hash)

    def __eq__(self, other):
        try:
            return self._hash == other._hash
        except AttributeError:
            return NotImplemented

    def __repr__(self):
        return "AdditiveGenerator(X=%r,enh_profile=%r,leg_dict=%r)"\
            % (self._X, self._enh_profile, self._leg_dict)
        # Better, but destroys tests:
        # return "AdditiveGenerator(enh_profile=%r,leg_dict=%r)"\
        #      % (self._enh_profile, self._leg_dict)

    def __str__(self):
        str = ""
        if self._leg_dict is not None:
            for l in self._leg_dict:
                str += "Psi class %r with exponent %r on level %r * "\
                    % (l, self._leg_dict[l], self._level_dict[l])
        str += "Graph %r" % (self._enh_profile,)
        return str

    def __mul__(self, other):
        """
        Multiply to psi products on the same graph (add dictionaries).

        Args:
            other (AdditiveGenerator): Product of psi classes on same graph.

        Returns:
            AdditiveGenerator: Product of psi classes on same graph.

        EXAMPLES::


            Also works without legs.

        """
        # Check that other is an AdditiveGenerator for the same graph:
        try:
            if self._X != other._X or self._enh_profile != other._enh_profile:
                return NotImplemented
            other_leg_dict = other._leg_dict
        except AttributeError:
            return NotImplemented
        # "unite" the leg_dicts:
        if self._leg_dict is None:
            self_leg_dict = {}
        else:
            self_leg_dict = self._leg_dict
        if other_leg_dict is None:
            other_leg_dict = {}
        new_leg_dict = {l: self_leg_dict.get(l, 0) + other_leg_dict.get(l, 0)
                        for l in set(self_leg_dict) | set(other_leg_dict)}
        return self._X.additive_generator(self._enh_profile, new_leg_dict)

    def __rmul__(self, other):
        self.__mul__(other)

    def __pow__(self, n):
        return self.pow(n)

    @property
    def enh_profile(self):
        return self._enh_profile

    @property
    def psi_degree(self):
        """
        Sum of powers of psi classes of self.
        """
        if self._leg_dict is None:
            return 0
        else:
            return sum(self._leg_dict.values())

    @cached_method
    def dim_check(self):
        """
        Check if, on any level, the psi degree is higher than the dimension.

        Returns:
            bool: False if the class is 0 for dim reasons, True otherwise.
        """
        # remove if degree > dim(X)
        if self.degree > self._X.dim():
            return False
        if self.codim == 0:
            # Avoid crazy infinite recursion for smooth graph :-)
            return True
        # for each level, check if psi product on level exceeds level dimension
        for level_number in range(self.codim + 1):
            assert self.level_dim(level_number) >= 0
            if self.degree_on_level(
                    level_number) > self.level_dim(level_number):
                return False
        return True

    @property
    def codim(self):
        """
        The codimension of the graph (number of levels)

        Returns:
            int: length of the profile
        """
        return len(self._enh_profile[0])

    @property
    def degree(self):
        """
        Degree of class, i.e. codimension of graph + psi-degree

        Returns:
            int: codim + psi_degree
        """
        # degree = codim of graph + powers of psi classes
        return self.codim + self.psi_degree

    @property
    def leg_dict(self):
        return self._leg_dict

    @property
    def level_dict(self):
        """
        The dictionary mapping leg -> level
        """
        return self._level_dict

    @property
    def inv_level_dict(self):
        """
        The dictionary mapping level -> list of legs on level.

        Returns:
            dict: level -> list of legs.
        """
        return self._inv_level_dict

    @cached_method
    def degree_on_level(self, level):
        """
        Total degree of psi classes on level.

        Args:
            level (int): (relative) level number (i.e. 0...codim)

        Raises:
            RuntimeError: Raised for level number out of range.

        Returns:
            int: sum of exponents of psis appearing on this level.
        """
        if level not in range(self.codim + 1):
            raise RuntimeError(
                "Illegal level number: %r on %r" % (level, self))
        try:
            return sum(self._leg_dict[leg]
                       for leg in self._inv_level_dict[level])
        except KeyError:
            # no psis on this level
            return 0

    def level(self, level_number):
        """
        Level of underlying graph.

        Args:
            level_number (int): (relative) level number (0...codim)

        Returns:
            LevelStratum: Stratum at level level_number of self._G.
        """
        return self._G.level(level_number)

    @cached_method
    def level_dim(self, level_number):
        """
        Dimension of level level_number.

        Args:
            level_number (int): (relative) level number (i.e. 0...codim)

        Returns:
            int: dimension of GeneralisedLevelStratum
        """
        level = self._G.level(level_number)
        return level.dim()

    @property
    def stack_factor(self):
        """
        The stack factor, that is the product of the prongs of the underlying graph
        divided by the product of the ells of the BICs and the automorphisms.

        Returns:
            QQ: stack factor
        """
        try:
            return self._stack_factor
        except AttributeError:
            # to get g_Gamma, we have to take the product of prongs/lcm for
            # each bic:
            prod = 1
            for k in self._G.LG.prongs.values():
                prod *= k

            p, _ = self.enh_profile

            bic_contr = 1
            for i in p:
                bic_contr *= self._X.bics[i].ell

            stack_factor = QQ(prod) / QQ(bic_contr *
                                         len(self._G.automorphisms))

            self._stack_factor = stack_factor
            return self._stack_factor

    @cached_method
    def as_taut(self):
        """
        Helper method, returns [(1,self)] as default input to ELGTautClass.
        """
        return admcycles.diffstrata.elgtautclass.ELGTautClass(self._X, [
                                                              (1, self)])

    @cached_method
    def is_in_ambient(self, ambient_enh_profile):
        """
        Check if ambient_enh_profile is an ambient graph, i.e. self is a degeneration
        of ambient_enh_profile.

        INPUT:

        ambient_enh_profile: tuple
        An enhanced profile.

        OUTPUT:

        True if there exists a leg map, False otherwise.
        """
        return self._X.is_degeneration(self._enh_profile, ambient_enh_profile)

    @cached_method
    def pow(self, n, amb=None):
        """
        Recursively calculate the n-th power of self (in amb), caching all results.

        Args:
            n (int): exponent
            amb (tuple, optional): enhanced profile. Defaults to None.

        Returns:
            ELGTautClass: self^n in CH(amb)
        """
        if amb is None:
            ONE = self._X.ONE
            amb = ((), 0)
        else:
            ONE = self._X.taut_from_graph(*amb)
        if n == 0:
            return ONE
        return self._X.intersection(self.as_taut(), self.pow(n - 1, amb), amb)

    @cached_method
    def exp(self, c, amb=None, stop=None):
        """
        exp(c * self) in CH(amb), calculated via exp_list.

        Args:
            c (QQ): coefficient
            amb (tuple, optional): enhanced profile. Defaults to None.
            stop (int, optional): cut-off. Defaults to None.

        Returns:
            ELGTautClass: the tautological class associated to the
                graded list exp_list.
        """
        # graded pieces are already reduced:
        new_taut_list = []
        for T in self.exp_list(c, amb, stop):
            new_taut_list.extend(T.psi_list)
        return admcycles.diffstrata.elgtautclass.ELGTautClass(
            self._X, new_taut_list, reduce=False)

    @cached_method
    def exp_list(self, c, amb=None, stop=None):
        """
        Calculate exp(c * self) in CH(amb).

        We calculate exp as a sum of powers (using self.pow, i.e. cached)
        and check at each step if the power vanishes (if yes, we obviously stop).

        The result is returned as a list consisting of the graded pieces.

        Optionally, one may specify the cut-off degree using stop (by
        default this is dim + 1).

        Args:
            c (QQ): coefficient
            amb (tuple, optional): enhanced profile. Defaults to None.
            stop (int, optional): cut-off. Defaults to None.

        Returns:
            list: list of ELGTautClasses
        """
        c = QQ(c)
        if amb is None:
            ONE = self._X.ONE
            amb = ((), 0)
        else:
            ONE = self._X.taut_from_graph(*amb)
        e = [ONE]
        f = ONE
        coeff = QQ(1)
        k = QQ(0)
        if stop is None:
            stop = self._X.dim() + 1
        while k < stop and f != self._X.ZERO:
            k += 1
            coeff *= c / k
            f = self.pow(k, amb)
            e.append(coeff * f)
        return e

    def pull_back(self, deg_enh_profile):
        """
        Pull back self to the graph associated to deg_enh_profile.

        Note that this returns an ELGTautClass as there could be several maps.

        More precisely, we return the sum over the pulled back classes divided
        by the number of undegeneration maps.

        Args:
            deg_enh_profile (tuple): enhanced profile of graph to pull back to.

        Raises:
            RuntimeError: raised if deg_enh_profile is not a degeneration of the
                underlying graph of self.

        Returns:
            ELGTautClass: sum of pullbacks of self to deg_enh_profile for each
                undegeneration map divided by the number of such maps.

        """
        if self._leg_dict is None:
            # trivial pullback
            return admcycles.diffstrata.elgtautclass.ELGTautClass(
                self._X, [(1, self._X.additive_generator(deg_enh_profile))])
        else:
            leg_maps = self._X.explicit_leg_maps(
                self._enh_profile, deg_enh_profile)
            if leg_maps is None:
                raise RuntimeError("Pullback failed: %r is not a degeneration of %r")\
                    % (deg_enh_profile, self._enh_profile)
            psi_list = []
            aut_factor = QQ(1) / QQ(len(leg_maps))
            for leg_map in leg_maps:
                new_leg_dict = {leg_map[l]: e for l,
                                e in self._leg_dict.items()}
                psi_list.append(
                    (aut_factor, self._X.additive_generator(
                        deg_enh_profile, new_leg_dict)))
            return admcycles.diffstrata.elgtautclass.ELGTautClass(
                self._X, psi_list)

    def psis_on_level(self, l):
        """
        The psi classes on level l of self.

        Args:
            l (int): level, i.e. 0,...,codim

        Returns:
            dict: psi dictionary on self.level(l).smooth_LG
        """
        L = self.level(l)
        # The psi classes on this level should be expressed in terms of the legs
        # of the smooth_LG of L:
        EG = L.smooth_LG
        try:
            # Careful: the legs of the smooth_LG are numbered 1,...,n
            # The psi classes are still numbered inside the whole graph
            # The conversion runs through the embedding of the LevelStratum
            # and back through the embedding of smooth_LG (dmp_inv)
            psis = {EG.dmp_inv[L.leg_dict[leg]]: self.leg_dict[leg]
                    for leg in self.inv_level_dict[l]}
        except KeyError:
            # no psis on this level
            psis = {}
        return psis

    def evaluate(self, quiet=False, warnings_only=False,
                 admcycles_output=False):
        """
        Evaluate self (cap with the fundamental class of self._X).

        Note that this gives 0 if self is not a top-degree class.

        Evaluation works by taking the product of the evaluation of each level
        (i.e. evaluating, for each level, the psi monomial on this level) and
        multiplying this with the stack factor.

        The psi monomials on the levels are evaluated using admcycles (after
        removing residue conditions).

        Raises a RuntimeError if there are inconsistencies with the psi degrees
        on the levels.

        INPUT:

        quiet: boolean (optional)
        If set to true, then get no output. Defaults to False.

        warnings_only: boolean (optional)
        If set to true, then output warnings. Defaults to False.

        admcycles_output: boolean (optional)
        If set to true, prints debugging info (used when evaluating levels). Defaults to False.

        OUTPUT:

        the integral of self on X as a rational number.
        """
        if self.degree < self._X.dim():
            if not quiet or warnings_only:
                print("Warning: %r is not of top degree: %r (instead of %r)" %
                      (self, self.degree, self._X.dim()))
            return 0
        level_list = []
        for l in range(self.codim + 1):
            if self.degree_on_level(l) < self.level_dim(l):
                raise RuntimeError(
                    "%r is of top degree, but not on level %r" % (self, l))
            L = self.level(l)
            value = L.evaluate(
                psis=self.psis_on_level(l),
                quiet=quiet,
                warnings_only=warnings_only,
                admcycles_output=admcycles_output)
            if value == 0:
                return 0
            level_list.append(value)
        # product over levels:
        prod = 1
        for p in level_list:
            prod *= p
        if not quiet:
            print("----------------------------------------------------")
            print("Contribution of Additive generator:")
            print(self)
            print("Product of level-wise integrals: %r" % prod)
            print("Stack factor: %r" % self.stack_factor)
            print("Total: %r" % (prod * self.stack_factor))
            print("----------------------------------------------------")
        return self.stack_factor * prod

    def to_prodtautclass(self, relabel=False):
        """
        Transform self into an admcycles prodtautclass on the underlying stgraph of self.

        Note that this gives the pushforward to M_g,n in the sense that we multiply with
        Strataclass and remove all residue conditions.

        Returns:
            prodtautclass: the prodtautclass of self, multiplied with the Strataclasses of
                the levels and all residue conditions removed.

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: X=Stratum((2,))
            sage: X.additive_generator(((),0)).to_prodtautclass()
            Outer graph : [2] [[1]] []
            Vertex 0 :
            Graph :      [2] [[1]] []
            Polynomial : -7/24*(kappa_1)_0 + 79/24*psi_1
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1] [[1, 3, 4]] [(3, 4)]
            Polynomial : -1/48
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1, 1] [[3], [1, 4]] [(3, 4)]
            Polynomial : -19/24
            sage: from admcycles.stratarecursion import Strataclass
            sage: X=GeneralisedStratum([Signature((4,-2,-2))], res_cond=[[(0,1)], [(0,2)]])
            sage: (X.additive_generator(((),0)).to_prodtautclass().pushforward() - Strataclass(1, 1, [4,-2,-2], res_cond=[2])).is_zero()
            True


        TESTS::

            sage: from admcycles import diffstrata, psiclass
            sage: X = diffstrata.generalisedstratum.GeneralisedStratum(sig_list = [diffstrata.sig.Signature(tuple([8,-3,-2,-3]))], res_cond = [[(0,1)],[(0,2)]])
            sage: X.psi(1).evaluate()
            9
            sage: v = X.ONE.to_prodtautclass().pushforward()
            sage: (v*psiclass(1,1,4)).evaluate()
            9

        We noticed the problem that the markings of the prodtautclass of an AdditiveGenerator are usually not the standard one,
        and hence its pushforward to M_g,n is not comparable to that of other AdditiveGenerator. Thus we solve this
        by adding the keyword "relabel"::

            sage: X=diffstrata.generalisedstratum.Stratum((4,-2))
            sage: X.additive_generator(((1,),0)).to_prodtautclass(relabel=True)
            Outer graph : [1, 0] [[3, 4], [1, 2, 5, 6]] [(3, 5), (4, 6)]
            Vertex 0 :
            Graph :      [1] [[1, 2]] []
            Polynomial : 1/2
            Vertex 1 :
            Graph :      [0] [[1, 2, 3, 4]] []
            Polynomial : psi_4
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1] [[1, 2]] []
            Polynomial : 1/2
            Vertex 1 :
            Graph :      [0, 0] [[2, 3, 9], [1, 4, 12]] [(9, 12)]
            Polynomial : 3
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1] [[1, 2]] []
            Polynomial : 1/2
            Vertex 1 :
            Graph :      [0, 0] [[3, 4, 9], [1, 2, 12]] [(9, 12)]
            Polynomial : -3


        We fixed the problem that when the ELGTautClass arised from removing residue conditions of a level stratum has only one term, the
        original algorithm would break down. We can test it::

            sage: X=diffstrata.generalisedstratum.GeneralisedStratum([diffstrata.sig.Signature((0,-1,-1)),diffstrata.sig.Signature((0,-1,-1))], res_cond=[[(0,1),(1,1)]])
            sage: X.additive_generator(((),0)).to_prodtautclass()
            Outer graph : [0, 0] [[1, 2, 3], [4, 5, 6]] []
            Vertex 0 :
            Graph :      [0] [[1, 2, 3]] []
            Polynomial : 1
            Vertex 1 :
            Graph :      [0] [[1, 2, 3]] []
            Polynomial : 1

        """
        LG = self._G.LG
        stgraph = LG.stgraph
        if any(self.level(l).zeroStratumClass()
               for l in range(self.codim + 1)):
            return admcycles.admcycles.prodtautclass(stgraph, terms=[])  # ZERO

        if len(LG.genera) == 1 and self._X.res_cond == []:  # just some stratum class without res_cond
            adm_psis = admcycles.admcycles.decstratum(stgraph, psi=self.leg_dict)
            adm_psis_taut = admcycles.admcycles.tautclass([adm_psis])
            sig = self._X._sig_list[0].sig
            g = self._X._sig_list[0].g
            stratum_class = admcycles.stratarecursion.Strataclass(g, 1, sig)
            result = admcycles.admcycles.prodtautclass(stgraph, protaut=[adm_psis_taut * stratum_class])
            return result

        # Now, we are dealing with nontrivial graphs
        alpha = []  # prodtautclasses on levels
        vertices = []  # embedding of level into stgraph
        for l in range(self.codim + 1):
            psis = self.psis_on_level(l)  # the psi classes on level l
            T = self.level(l).remove_res_cond(psis)

            # turn to the to_prodtautclass of ELGTautClass which will recursively be solved
            alpha.append(T.to_prodtautclass())
            vertices.append(LG.verticesonlevel(LG.internal_level_number(l)))
        prod = self.stack_factor * admcycles.admcycles.prodtautclass(stgraph)

        for l, ptc in enumerate(alpha):
            prod = prod.factor_pullback(vertices[l], ptc)  # returns product (!)

        if relabel:
            # we standardise the marking of the EmbeddedLevelGraph and leg_dict for psi
            standard_legdict = self._G.standard_markings()
            newEmbLG = self._G.relabel(standard_legdict, tidyup=False)
            stgraph = newEmbLG.LG.stgraph

            # get the prodtaut under the standardised stable graph
            prod = admcycles.admcycles.prodtautclass(stgraph, prod.terms)

        return prod
