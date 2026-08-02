import itertools

# pylint does not know sage
from sage.structure.sage_object import SageObject  # pylint: disable=import-error
from sage.matrix.constructor import matrix  # pylint: disable=import-error
from sage.misc.flatten import flatten  # pylint: disable=import-error
from sage.misc.cachefunc import cached_method  # pylint: disable=import-error
from sage.rings.rational_field import QQ  # pylint: disable=import-error
from sage.arith.functions import lcm  # pylint: disable=import-error


class EmbeddedLevelGraph(SageObject):
    """
    LevelGraph inside a generalised stratum.

    Note that the points of the enveloping GeneralisedStratum are of the form
    (i,j) where i is the component and j the index of sig of that component,
    while the points of the level graph are numbers 1,...,n.

    Thus, dmp is a dictionary mapping integers to tuples of integers.

    Attributes:

    * LG (LevelGraph): underlying LevelGraph
    * X (GeneralisedStratum): enveloping stratum
    * dmp (dict): (bijective!) dictionary marked points of LG -> points of stratum
    * dmp_inv (dict): inverse of dmp
    * dlevels (dict): (bijective!) dictionary levels of LG -> new level numbering
    * dlevels_inv (dict): inverse of dlevels
    * top (GeneralisedStratum): (if self is a BIC) top component
    * bot (GeneralisedStratum): (if self is a BIC) bottom component
    * clutch_dict (dict): (if self is a BIC) dictionary mapping points of top
      stratum to points of bottom stratum where there is an edge in self.
    * emb_top (dict): (if self is a BIC) dictionary mapping points of stratum top
      to the corresponding points of the enveloping stratum.
    * emb_bot (dict): (if self is a BIC) dictionary mapping points of stratum bot
      to the corresponding points of the enveloping stratum.
    * automorphisms (list of list of dicts): automorphisms
    * codim (int): codimension of LevelGraph in stratum
      number_of_levels (int): Number of levels of self.

    Note that attempting to access any of the attributes top, bot, clutch_dict,
    emb_top or emb_bot will raise a ValueError if self is not a BIC.
    """
    def __init__(self, X, LG, dmp, dlevels):
        """
        Initialises EmbeddedLevelGraph.

        Args:
            LG (LevelGraph): underlying LevelGraph
            X (GeneralisedStratum): enveloping stratum
            dmp (dictionary): (bijective!) dictionary marked points of LG -> points of stratum
            dlevels (dictionary): (bijective!) dictionary levels of LG -> new level numbering
        """
        self.LG = LG
        self.X = X
        self.dmp = dmp
        self.dmp_inv = {value: key for key, value in dmp.items()}
        self.add_vertices_at_infinity()
        self.dlevels = dlevels
        self.dlevels_inv = {value: key for key, value in dlevels.items()}
        self._top = None
        self._bot = None
        self._clutch_dict = None
        self._emb_top = None
        self._emb_bot = None
        self._automorphisms = None
        self._level = {}
        self._ell = None
        self.codim = self.LG.codim()
        self.number_of_levels = len(set(self.dlevels.keys()))

    def __repr__(self):
        return "EmbeddedLevelGraph(LG=%r,dmp=%r,dlevels=%r)" % (
            self.LG, self.dmp, self.dlevels)

    def __str__(self):
        return (
            "Embedded Level Graph consisting of %s with point dictionary %s and level dictionary %s" %
            (self.LG, self.dmp, self.dlevels))

    def explain(self):
        """
        A more user-friendly display of __str__ :-)
        """
        def _list_print(L):
            if len(L) > 1:
                s = ['s ']
                for x in L[:-2]:
                    s.append('%r, ' % x)
                s.append('%r ' % L[-2])
                s.append('and %r.' % L[-1])
                return ''.join(s)
            else:
                return ' %r.' % L[0]

        def _num(i):
            if i == 1:
                return 'one edge'
            else:
                return '%r edges' % i
        print("LevelGraph embedded into stratum %s with:" % self.X)
        LG = self.LG
        for l in range(LG.numberoflevels()):
            internal_l = LG.internal_level_number(l)
            print("On level %r:" % l)
            for v in LG.verticesonlevel(internal_l):
                print("* A vertex (number %r) of genus %r" % (v, LG.genus(v)))
        levels_of_mps = list(
            set(LG.level_number(LG.levelofleg(leg)) for leg in self.dmp))
        print("The marked points are on level%s" %
              _list_print(sorted(levels_of_mps)))
        print("More precisely, we have:")
        for leg in self.dmp:
            print(
                "* Marked point %r of order %r on vertex %r on level %r" %
                (self.dmp[leg],
                 LG.orderatleg(leg),
                 LG.vertex(leg),
                    LG.level_number(
                    LG.levelofleg(leg))))
        print("Finally, we have %s. More precisely:" % _num(len(LG.edges)))
        edge_dict = {e: (LG.vertex(e[0]), LG.vertex(e[1])) for e in LG.edges}
        edge_dict_inv = {}
        for k, v in edge_dict.items():
            if v in edge_dict_inv:
                edge_dict_inv[v].append(k)
            else:
                edge_dict_inv[v] = [k]
        for e in edge_dict_inv:
            print("* %s between vertex %r (on level %r) and vertex %r (on level %r) with prong%s" %
                  (_num(len(edge_dict_inv[e])),
                   e[0], LG.level_number(LG.levelofvertex(e[0])),
                   e[1], LG.level_number(LG.levelofvertex(e[1])),
                   # _write_prongs()
                   _list_print([LG.prong(ee) for ee in edge_dict_inv[e]])))

    def __eq__(self, other):
        if not isinstance(other, EmbeddedLevelGraph):
            return False
        return self.LG == other.LG and self.dmp == other.dmp and self.dlevels == other.dlevels

    @cached_method
    def is_bic(self):
        return self.LG.is_bic()

    @property
    def ell(self):
        """
        If self is a BIC: the lcm of the prongs.

        Raises:
            RuntimeError: raised if self is not a BIC.

        Returns:
            int: lcm of the prongs.
        """
        if self._ell is None:
            if not self.is_bic():
                raise RuntimeError("ell only defined for BICs!")
            self._ell = lcm(self.LG.prongs.values())
        return self._ell

    # Dawei's positivity coefficients
    @property
    def b(self):
        if self.X._h0 > 1:
            raise ValueError('Cannot compute b on disconnected stratum.')
        g = self.X._sig_list[0].g
        val = 0
        # super inefficient, but probably good enough for now:
        # take the underlying StableGraph and, for each edge,
        # contract all other edges and check which graph we end up with:
        stgraph = self.LG.stgraph
        for e in self.LG.edges:
            ee = self.LG.edges[:]
            ee.remove(e)
            curr_graph = stgraph.copy()
            for contr in ee:
                curr_graph.contract_edge(contr)
            assert curr_graph.edges() == [e]
            if len(curr_graph.genera()) == 2:
                # compact type:
                i = min(curr_graph.genera())
                val += QQ(6 * i * (g - i)) / QQ((g + 3) * self.LG.prong(e))
            else:
                # irreducible
                assert len(curr_graph.genera()) == 1
                val += QQ(g + 1) / QQ((g + 3) * self.LG.prong(e))
        return self.ell * val

    @property
    def c(self):
        return self.ell * (self.bot.kappa_EKZ -
                           self.X.kappa_EKZ * QQ(self.bot.N) / QQ(self.X.N))
    ######

    @property
    def top(self):
        if self._top is None:
            self.split()
        return self._top

    @property
    def bot(self):
        if self._bot is None:
            self.split()
        return self._bot

    @property
    def clutch_dict(self):
        if self._clutch_dict is None:
            self.split()
        return self._clutch_dict

    @property
    def emb_bot(self):
        if self._emb_bot is None:
            self.split()
        return self._emb_bot

    @property
    def emb_top(self):
        if self._emb_top is None:
            self.split()
        return self._emb_top

    def add_vertices_at_infinity(self):
        """
        We add the vertices at infinity to the underlying_graph of self.LG.

        These are given by the residue conditions.

        More precisely: Recall that the underlying_graph of self.LG has vertices
        and edges of self.LG stored in the form UG_vertex = (vertex number, genus, 'LG')
        and edges of the underlying graph are of the form: (UG_vertex, UG_vertex, edge name)
        We now add vertices 'at level infinity' by adding, for each res_cond of self.X

        * a UG_vertex called (i, 0, 'res') (where i is the index of the condition in res_cond
          we are currently considering) and edges so that each residue
          condition corresponds to an edge from the corresponding pole to some
          residue at 'level infinity'. We store these in the form:
        * (res_vertex, UG_vertex, resiedgej)
          Here UG_vertex is the vertex of self.LG, in the form (vertex number, genus, 'LG'),
          that res_vertex is attached to and j is the leg of that vertex (as a leg of self.LG!)
          corresponding to the pole that resi should be attached to.
        """
        # remove any that might already be there:
        existing_residues = [v for v in self.LG.underlying_graph.vertices(sort=True)
                             if v[2] == 'res']
        for v in existing_residues:
            self.LG.underlying_graph.delete_vertex(v)
        # add a vertex for every residue condition:
        # TODO: remove duplicates?
        edges = []
        for i, rc in enumerate(self.X.res_cond):
            v_name = (i, 0, 'res')
            for p in rc:
                e_name = 'res%redge%r' % (i, self.dmp_inv[p])
                v_on_graph = self.LG.vertex(self.dmp_inv[p])
                edges.append((self.LG.UG_vertex(v_on_graph), v_name, e_name))
        self.LG.underlying_graph.add_edges(edges)

    @property
    @cached_method
    def residue_matrix_from_RT(self):
        """
        The matrix associated to the residue conditions imposed by the residue theorem
        on each vertex of self.

        Returns:
            SAGE Matrix: matrix of residue conditions given by RT
        """
        poles_by_vertex = {}
        for p in self.X._polelist:
            vertex = self.LG.vertex(self.dmp_inv[p])
            try:
                poles_by_vertex[vertex].append(p)
            except KeyError:
                poles_by_vertex[vertex] = [p]
        rows = []
        for v in poles_by_vertex:
            rows.append([int(p in poles_by_vertex[v])
                         for p in self.X._polelist])
        return matrix(QQ, rows)

    @property
    @cached_method
    def full_residue_matrix(self):
        """
        Residue matrix with GRC conditions and RT conditions (for each vertex).

        OUTPUT:

        A matrix with number of poles columns and a row for each condition.
        """
        M = self.X.residue_matrix()
        if M:
            M = M.stack(self.residue_matrix_from_RT)
        else:
            M = self.residue_matrix_from_RT
        return M

    def residue_zero(self, pole):
        """
        Check if the residue at pole is forced zero by residue conditions.

        NOTE: We DO include the RT on the vertices in this check!

        Args:
            pole (tuple): pole (as a point (i,j) of self.X)

        Returns:
            bool: True if forced zero, False otherwise.
        """
        # add the equation corresponding to the residue at pole to the residue matrix
        # and see if the rank changes:
        i = self.X._polelist.index(pole)
        res_vec = [[int(i == j) for j in range(len(self.X._polelist))]]
        RM = self.full_residue_matrix
        # RM = self.X.residue_matrix()
        if RM:
            stacked = RM.stack(matrix(res_vec))
            return stacked.rank() == self.full_residue_matrix.rank()
            # return stacked.rank() == self.X.residue_matrix().rank()
        else:
            return False

    def level(self, l):
        """
        The generalised stratum on level l.

        Note that this is cached, i.e. on first call, it is stored in the dictionary
        _level.

        INPUT:

        l: integer
        The relative level number (0,...,codim)

        OUTPUT:

        The LevelStratum that is

        * a list of Signatures (one for each vertex on the level)
        * a list of residue conditions, i.e. a list [res_1,...,res_r]
          where each res_k is a list of tuples [(i_1,j_1),...,(i_n,j_n)]
          where each tuple (i,j) refers to the point j (i.e. index) on the
          component i and such that the residues at these points add up
          to 0.
        * a dictionary of legs, i.e. n -> (i,j) where n is the original
          number of the point (on the LevelGraph self) and i is the
          number of the component, j the index of the point in the signature tuple.

        Note that LevelStratum is a GeneralisedStratum together with
        a leg dictionary. Here, we provide an additional attribute:

        * leg_orbits, a nested list giving the orbits of the points on the level
          under the automorphism group of self.
        """
        try:
            LS = self._level[l]
        except KeyError:
            # for the residue conditions: We add the residue conditions from
            # the enveloping stratum:
            # We do this by passing those poles with residue forced
            # zero as those to be ignored in the residue calculations
            # performed by the LevelGraph:
            # We have to translate them to points on self:
            # Note that self.LG knows the "level at infinity"
            excluded_poles = tuple(self.dmp_inv[p] for p in flatten(
                self.X.res_cond, max_level=1))
            LS = self.LG.stratum_from_level(l, excluded_poles=excluded_poles)
            # add automorphism info
            LS.leg_orbits = []
            seen = set()
            for leg in LS._leg_dict:
                if leg in seen:
                    continue
                curr_orbit = [LS._leg_dict[leg]]
                for _v_map, l_map in self.automorphisms:
                    curr_orbit.append(LS._leg_dict[l_map[leg]])
                    seen.update([l_map[leg]])
                LS.leg_orbits.append(list(set(curr_orbit))
                                     )  # remove duplicates
            self._level[l] = LS
        return LS

    def legs_are_isomorphic(self, leg, other_leg):
        """
        Check if leg and other_leg are in the same Aut-orbit of self.

        Args:
            leg (int): leg on self.LG
            other_leg (int): leg on self.LG

        Raises:
            RuntimeError: If leg is not in any aut-orbit of the level it should be on.

        Returns:
            bool: True if they are in the same orbit of self.level(levelofleg),
                    False, otherwise.

        EXAMPLES::

            Note the asymmetric banana graph.

            The symmetric one has isomorphisms.

            Legs are isomorphic to themselves.

            It's symmetric.

        """
        level = self.LG.level_number(self.LG.levelofleg(leg))
        other_level = self.LG.level_number(self.LG.levelofleg(other_leg))
        if level != other_level:
            return False
        assert level == other_level
        emb_leg = self.level(level)._leg_dict[leg]
        emb_other_leg = self.level(level)._leg_dict[other_leg]
        for orbit in self.level(level).leg_orbits:
            if emb_leg in orbit:
                return emb_other_leg in orbit
        else:
            raise RuntimeError(
                "leg %r not in any orbit %r of %r" %
                (leg, self.level(level).leg_orbits, self.level(level)))

    @cached_method
    def edge_orbit(self, edge):
        """
        The edge orbit of edge in self.

        raises a ValueError if edge is not an edge of self.LG

        INPUT:

        edge: tuple
        An edge of ``self.LG``, i.e. tuple (start leg, end leg), where start
        leg should not be on a lower level than end leg.

        OUTPUT:

        The set of edges in automorphism orbit of ``edge``.
        """
        if edge not in self.LG.edges:
            raise ValueError("%r is not an edge of %r!" % (edge, self))
        s = set([edge])
        for v_map, l_map in self.automorphisms:
            new_edge = (l_map[edge[0]], l_map[edge[1]])
            s.add(new_edge)
        return s

    def len_edge_orbit(self, edge):
        """
        Length of the edge orbit of edge in self.

        Args:
            edge (tuple): edge of self.LG, i.e. tuple (start leg, end leg), where
                start leg should not be on a lower level than end leg.

        Raises:
            ValueError: if edge is not an edge of self.LG

        Returns:
            int: length of the aut-orbit of edge.

        EXAMPLES::


            Prongs influence the orbit length.

        """
        return len(self.edge_orbit(edge))

    def automorphisms_stabilising_legs(self, leg_tuple):
        stabs = []
        for v_map, l_map in self.automorphisms:
            for l in leg_tuple:
                if l_map[l] != l:
                    break
            else:  # no break
                stabs.append(l_map)
        return stabs

    def delta(self, i):
        """
        Squish all levels except for i.

        Note that delta(1) contracts everything except top-level and that the
        argument is interpreted via internal_level_number (i.e. a relative level number).

        Moreover, dlevels is set to map to 0 and -1(!).

        Args:
            i (int): Level not to be squished.

        Returns:
            EmbeddedLevelGraph: Embedded BIC (result of applying delta to the
                underlying LevelGraph)
        """
        newLG = self.LG.delta(i, quiet=True)
        newdmp = self.dmp.copy()
        # level_number is (positive!) relative level number.
        newdlevels = {l: -newLG.level_number(l) for l in newLG.levels}
        return EmbeddedLevelGraph(self.X, newLG, newdmp, newdlevels)

    def squish_vertical(self, level):
        """
        Squish level crossing below level 'level'.

        Note that in contrast to the levelgraph method, we work with relative
        level numbers here!

        Args:
            level (int): relative (!) level number.

        Returns:
            EmbeddedLevelGraph: Result of squishing.

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: X=GeneralisedStratum([Signature((4,))])
            sage: p = X.enhanced_profiles_of_length(4)[0][0]
            sage: g = X.lookup_graph(p)

            lookup_graph uses the sorted profile (note that these do not have to be reduced!):

            sage: assert any(g.squish_vertical(0).is_isomorphic(G) for G in X.lookup(p[1:]))
            sage: assert any(g.squish_vertical(1).is_isomorphic(G) for G in X.lookup(p[:1]+p[2:]))
            sage: assert any(g.squish_vertical(2).is_isomorphic(G) for G in X.lookup(p[:2]+p[3:]))
            sage: assert any(g.squish_vertical(3).is_isomorphic(G) for G in X.lookup(p[:3]))

            Squishing outside the range of levels does nothing:

            sage: assert g.squish_vertical(4) == g

            Recursive squishing removes larger parts of the profile:

            sage: assert any(g.squish_vertical(3).squish_vertical(2).is_isomorphic(G) for G in X.lookup(p[:2]))
        """
        newLG = self.LG.squish_vertical(
            self.LG.internal_level_number(level), quiet=True)
        newdmp = self.dmp.copy()
        # level_number is (positive!) relative level number.
        newdlevels = {l: -newLG.level_number(l) for l in newLG.levels}
        return EmbeddedLevelGraph(self.X, newLG, newdmp, newdlevels)

    def split(self):
        """
        Splits embedded BIC self into top and bottom component.

        Raises a ValueError if self is not a BIC.

        OUTPUT:

        A dictionary consising of

        * the GeneralisedStratum self.X
        * the top component LevelStratum
        * the bottom component LevelStratum
        * the clutching dictionary mapping ex-half-edges on
          top to their partners on bottom (both as points in the
          respective strata!)
        * a dictionary embedding top into the stratum of self
        * a dictionary embedding bot into the stratum of self

        Note that clutch_dict, emb_top and emb_bot are dictionaries between
        points of strata, i.e. after applying dmp to the points!
        """
        if not self.is_bic():
            raise ValueError(
                "Error: %s is not a BIC! Cannot be split into Top and Bottom component!" %
                self)
        self._top = self.level(0)
        self._bot = self.level(1)
        # To construct emb_top and emb_bot, we have to combine self.dmp with the
        # the leg_dicts of top and bot.
        # More precisely: emb_top is the composition of the inverse of the leg_dict
        # of top, i.e. top.stratum_number, and self.dmp
        # (giving a map from the points of top to the points of the enveloping
        # stratum of self) and the same for bot.
        # We implement this by iterating over the marked points of self on top level,
        # which are exactly the keys of self.dmp that are on top level.
        # Note that we make extra sure that we didn't mess up the level numbering by
        # using the relative level numbering (where the top level is guaranteed to be 0
        # and the bottom level is 1 (positive!)).
        self._emb_top = {self._top.stratum_number(l): self.dmp[l]
                         for l in iter(self.dmp)
                         if self.LG.level_number(self.LG.levelofleg(l)) == 0}
        self._emb_bot = {self._bot.stratum_number(l): self.dmp[l]
                         for l in iter(self.dmp)
                         if self.LG.level_number(self.LG.levelofleg(l)) == 1}
        # Because this is a BIC, all edges of self are cut in this process
        # and this is exactly the dictionary we must remember
        # WARNING: Here we assume that e[0] is on top level and e[1] is on bottom
        #   This is assured by tidy_up, e.g. after initialisation!
        # Note that all these dictionaries map points of GeneralisedStrata to each
        # other so we must take the corresponding stratum_number!
        self._clutch_dict = {
            self._top.stratum_number(
                e[0]): self._bot.stratum_number(
                e[1]) for e in self.LG.edges}
        return {'X': self.X, 'top': self._top, 'bottom': self._bot,
                'clutch_dict': self._clutch_dict,
                'emb_dict_top': self._emb_top, 'emb_dict_bot': self._emb_bot, }

    def is_legal(self):
        """
        Check the R-GRC for self.

        Returns:
            bool: result of R-GRC.
        """
        # Check if any levels are empty:
        # Note that this can only happen if self.X has simple poles (as we never
        # have horizontal edges)
        if list(self.X.simple_poles()):
            if any(self.level(l).is_empty()
                   for l in range(self.number_of_levels)):
                return False
        # poles are excluded if they are contained in _any_ residue condition of the stratum.
        # In particular, they are _not_ excluded if they are only restrained by
        # the RT on some component!
        poles_in_rc_stratum = flatten(self.X.res_cond, max_level=1)
        poles_in_rc_graph = tuple(self.dmp_inv[p] for p in poles_in_rc_stratum)
        return self.LG.is_legal(excluded_poles=poles_in_rc_graph, quiet=True)

    def standard_markings(self):
        r"""
        Construct a dictionary for relabelling the markings. A standard labelling will label the legs
        of markings first and then the half edges. If the generalised stratum has only one component,
        the standard label of a marking will be exactly the position of that marking in the signature.

        EXAMPLES::

            sage: from admcycles.diffstrata.generalisedstratum import Stratum
            sage: X=Stratum((1,1))
            sage: X.bics[0]
            EmbeddedLevelGraph(LG=LevelGraph([1, 0],[[1, 2], [3, 4, 5, 6]],[(1, 5), (2, 6)],{1: 0, 2: 0, 3: 1, 4: 1, 5: -2, 6: -2},[0, -1],True),dmp={3: (0, 0), 4: (0, 1)},dlevels={0: 0, -1: -1})
            sage: X.bics[0].standard_markings()
            {1: 3, 2: 4, 3: 1, 4: 2, 5: 5, 6: 6}

        """
        n_list = [0 for i in range(self.X._h0)]  # list of number of markings on each component of the stratum
        for t in self.dmp_inv:
            n_list[t[0]] += 1

        # list such that the j-th entry is the sum of numbers of
        # markings on the components of smaller indices
        n_sums = [0] + [sum(n_list[i] for i in range(j))
                        for j in range(1, self.X._h0)]
        new_leg_dict = {}  # the mapping dict for relabelling the legs
        h = 1
        for i in range(1, len(self.LG.poleorders) + 1):
            if i in self.dmp:
                new_leg_dict[i] = (self.dmp[i][1] +
                                   n_sums[self.dmp[i][0]] + 1)
            else:
                new_leg_dict[i] = len(self.dmp) + h
                h = h + 1
        return new_leg_dict

    def relabel(self, legdict, tidyup=True):
        r"""
        Relabel the EmbeddedLevelGraph by a given dictionary.

        INPUT:

            - legdict (dict): A dictionary indicating the map from old markings
              to the new ones

        EXAMPLES::
            sage: from admcycles.diffstrata.generalisedstratum import Stratum
            sage: X = Stratum((1,1))
            sage: dict1={1: 3, 2: 4, 3: 1, 4: 2, 5: 5, 6: 6}
            sage: X.bics[0].relabel(dict1)
            EmbeddedLevelGraph(LG=LevelGraph([1, 0],[[3, 4], [1, 2, 5, 6]],[(3, 5), (4, 6)],{1: 1, 2: 1, 3: 0, 4: 0, 5: -2, 6: -2},[0, -1],True),dmp={1: (0, 0), 2: (0, 1)},dlevels={0: 0, -1: -1})

        """
        newLG = self.LG.relabel(legdict, tidyup)
        new_dmp = {legdict[i]: j for i, j in self.dmp.items()}
        if tidyup:
            new_dmp = {a: b for a, b in sorted(new_dmp.items())}
        newEmbLG = EmbeddedLevelGraph(self.X, newLG, new_dmp, self.dlevels)

        return newEmbLG

    def is_isomorphic(self, other):
        """
        Check if self and other are isomorphic (as EmbeddedLevelGraphs).

        Args:
            other (EmbeddedLevelGraph): Graph to check isomorphism.

        Returns:
            bool: True if there exists at least one isomorphism.
        """
        # TODO: Maybe include a way to check against unembedded LGs
        # TODO: Check embedding!
        if not isinstance(other, EmbeddedLevelGraph):
            return False
        try:
            next(self.isomorphisms(other))
            return True
        except StopIteration:
            return False

    @property
    def automorphisms(self):
        """
        The automorphisms of self (as automorphisms of the underlying LevelGraph,
        respecting the embedding, see doc of isomorphisms).

        OUTPUT:

        A list of pairs ``(map_of_vertices, map_of_legs)``. Both ``maps_of_vertices``
        and ``map_of_legs`` are dictionaries.
        """
        if not self._automorphisms:
            self._automorphisms = list(self.isomorphisms(self))
        return self._automorphisms

    def isomorphisms(self, other):
        """
        Generator yielding the "next" isomorphism of self and other.

        Note that while this gives an "isomorphism" from self.LG to other.LG, this
        is not necessarily an isomorphism of the LevelGraphs (the numbered points may
        be permuted if this is "fixed" by the embedding).

        INPUT:

        other: EmbeddedLevelGraph
        The (potentially) isomorphic EmbeddedLevelGraph.

        OUTPUT:

        An iterator going through isomorphisms. Each isomorphism is encoded by a
        pair of dictionaries ``(vertex_map, leg_map)`` The dictionaries
        ``vertex_map`` (respectively ``leg_map``) is to the mapping of the
        vertices (resp. legs) of ``self.LG`` to the vertices (resp. legs) of
        ``other.LG``.
        """
        # Isomorphisms of EmbeddedLevelGraphs:
        # An isomorphism of EmbeddedLevelGraph is a set of compatible level isomorphisms.
        # We iterate through the isomorphisms on each level and yield whenever we find
        # compatible levelisomorphisms for all levels.
        # Note that we use dlevels for this, as these should be compatible.
        # There are (at least) two ways in which this can be optimised:
        # * don't go through the entire product before checking edge compatibility!
        # * choose a smart ordering of levels (e.g. number of vertices)
        isom_vertices = {}
        isom_legs = {}
        for level_isos in itertools.product(
                *[self._level_iso(other, l) for l in self.dlevels.values()]):
            for level_iso_v, level_iso_l in level_isos:
                isom_vertices.update(level_iso_v)
                isom_legs.update(level_iso_l)
            # check edge compatibility
            for e in self.LG.edges:
                if (isom_legs[e[0]], isom_legs[e[1]]) not in other.LG.edges:
                    break
            else:  # no break
                yield isom_vertices.copy(), isom_legs.copy()

    def _level_iso(self, other, l):
        """
        Generator yielding the "next" isomorphism of level l of self and other.

        Here, l is a value of dlevels (this should be compatible).

        Note that we require the graph to have no horizontal edges, i.e. the level
        contains no edges!

        TODO: * Maybe add future horizontal support?
              * Maybe use relative level number instead? (this seems to give weird behaviour
                right now...)

        Args:
            other (EmbeddedLevelGraph): The embedded graph we are checking for isomorphism.
            l (int): Level number (embedded into the stratum, i.e. value of dlevels).

        Yields:
            tuple: the next isomorphism of levels:
                dict: vertices of self.LG -> vertices of other.LG
                dict: legs of self.LG -> legs of other.LG
        """
        # Isomorphisms of levels:
        # An isomorphism of levels consist of
        # * a map: vertices to vertices
        # * a map: legs to legs
        # respecting:
        # * the genus,
        # * the number of legs on every vertex,
        # * the order at every leg,
        # * the marked points of the stratum (via dmp).
        ####
        # First, we extract the data for level l from self and other.
        # Note that we do not use stratum_from_level to avoid all the overhead.
        # TODO: All this should be cached!!
        l_self = self.LG.internal_level_number(l)
        l_other = other.LG.internal_level_number(l)
        # we need to be careful to distinguish the indices in the list of genera
        # of the LevelGraph from the actual genera.
        vv_self_idx = self.LG.verticesonlevel(l_self)  # list of indices
        vv_other_idx = other.LG.verticesonlevel(l_other)  # list of indices
        if len(vv_self_idx) != len(vv_other_idx):
            return
        vv_self = [self.LG.genus(i) for i in vv_self_idx]  # list of genera
        vv_other = [other.LG.genus(i) for i in vv_other_idx]  # list of genera
        # extract the legs: (nested lists)
        legs_self = [self.LG.legsatvertex(v) for v in vv_self_idx]
        legs_other = [other.LG.legsatvertex(v) for v in vv_other_idx]
        # build dictionary: leg -> index in vv
        leg_dict_self = {l: i for i, legs in enumerate(
            legs_self) for l in legs}
        leg_dict_other = {l: i for i, legs in enumerate(
            legs_other) for l in legs}
        if len(leg_dict_self) != len(leg_dict_other):
            return
        # for quick checks, we create sorted lists of the orders at each vertex
        order_sorted_self = [sorted([self.LG.orderatleg(l) for l in legs])
                             for legs in legs_self]
        order_sorted_other = [sorted([other.LG.orderatleg(l) for l in legs])
                              for legs in legs_other]
        # We create the two maps as dictionaries:
        # index of vv_self -> index of vv_other
        isom_vert = {}
        # leg number (on self.LG) -> leg number (on other.LG)
        isom_legs = {}
        # We also want to keep track of whom we've already mapped:
        # source is a dictionary: genus -> list of indices of vv_self
        source = {}
        for i, g in enumerate(vv_self):
            try:
                source[g].append(i)
            except KeyError:
                source[g] = [i]
        # target is a dictionary: genus -> list of indices of vv_other
        target = {}
        for i, g in enumerate(vv_other):
            try:
                target[g].append(i)
            except KeyError:
                target[g] = [i]
        # for the legs we build copies of the nested lists to manipulate
        legs_source = [legs[:] for legs in legs_self]
        legs_target = [legs[:] for legs in legs_other]
        # Next, we exclude some deal-breakers:
        # * The same genera must appear.
        if sorted(vv_self) != sorted(vv_other):
            return
        # * The same embedded points have to be on this level (they have to be
        # mapped to each other!)
        # In particular, this gives a part of the leg map (and thus also of the
        # vertex map).
        for p_self, p in self.dmp.items(
        ):  # p is the "shared" point of the stratum
            p_other = other.dmp_inv[p]
            # If neither point is on this level, we continue:
            if not (p_self in leg_dict_self or p_other in leg_dict_other):
                continue
            # The vertex of p_self must map to that of p_other.
            # Three things can fail here:
            # * only one of the two points is on this level.
            if ((p_self in leg_dict_self and p_other not in leg_dict_other) or (
                    p_self not in leg_dict_self and p_other in leg_dict_other)):
                return
            v_self = leg_dict_self[p_self]
            v_other = leg_dict_other[p_other]
            # * the points are on incompatible vertices (genus or numbers/orders of legs!)
            if (vv_self[v_self] != vv_other[v_other] or
                len(legs_self[v_self]) != len(legs_other[v_other]) or
                    order_sorted_self[v_self] != order_sorted_other[v_other]):
                return
            # * two points are on the same vertex in one case, but on different vertices
            #   in the other. I.e. v_self is already being mapped somewhere other than v_other
            #   or v_other is already being mapped to (by someone else)
            try:
                if isom_vert[v_self] != v_other:
                    return
            except KeyError:  # v_self not being mapped yet, i.e. still in source
                g = vv_other[v_other]
                if v_other in target[g]:  # make sure v_other is still a possible target
                    isom_vert[v_self] = v_other
                    source[g].remove(v_self)
                    target[g].remove(v_other)
                else:
                    return
            # now we can safely map the legs:
            isom_legs[p_self] = p_other
            # and remove them from source and target (so they won't be
            # reassigned later)
            legs_source[v_self].remove(p_self)
            legs_target[v_other].remove(p_other)
        # Next, we construct maps of the remaining vertices.
        # For this, we use a small recursive function:
        curr_v_map = {}
        legal_v_maps = []

        def vertex_maps(sl, tl):
            if not sl:
                # all entries of tl should be None at this point:
                assert all(tv is None for tv in tl)
                legal_v_maps.append(curr_v_map.copy())
                return
            curr_v = sl.pop()
            curr_legs = len(legs_self[curr_v])
            # try to map curr_v to tl:
            for i, tv in enumerate(tl):
                # we temporarily set "hit" targets to None so we don't have to worry
                # about indexing...
                if tv is None:
                    continue
                # check if legs _can_ be compatible:
                if (curr_legs != len(legs_other[tv]) or
                        order_sorted_self[curr_v] != order_sorted_other[tv]):
                    continue
                tl[i] = None
                curr_v_map[curr_v] = tv
                vertex_maps(sl, tl)
                # undo
                del curr_v_map[curr_v]
                tl[i] = tv
            # undo
            sl.append(curr_v)
        # the function for the legs is almost the same, just the condition is
        # different:
        curr_l_map = {}
        legal_l_maps = []

        def leg_maps(sl, tl):
            if not sl:
                # all entries of tl should be None at this point:
                assert all(tleg is None for tleg in tl)
                legal_l_maps.append(curr_l_map.copy())
                return
            curr_l = sl.pop()
            # try to map curr_l to tl:
            for i, tleg in enumerate(tl):
                # we temporarily set "hit" targets to None so we don't have to worry
                # about indexing...
                if tleg is None:
                    continue
                # check if orders are compatible:
                if self.LG.orderatleg(curr_l) == other.LG.orderatleg(tleg):
                    tl[i] = None
                    curr_l_map[curr_l] = tleg
                    leg_maps(sl, tl)
                    # undo
                    del curr_l_map[curr_l]
                    tl[i] = tleg
            # undo
            sl.append(curr_l)
        # Now we build the list of all vertex isomorphisms going through the
        # vertices by genus
        v_isom_list = []
        for g in source:
            legal_v_maps = []  # will get filled by vertex_maps
            vertex_maps(source[g], target[g])
            v_isom_list.append(legal_v_maps[:])  # copy!
        # v_isom_list is now a nested list of maps for each genus.
        # the product consists of tuples, one map for every genus.
        for v_maps in itertools.product(*v_isom_list):
            for v_map in v_maps:
                # this overwrites exactly the vertices in source.
                isom_vert.update(v_map)
            # Finally, the returned vertex map should use the indexing of the
            # LevelGraph, not of the level:
            return_isom_vert = {vv_self_idx[k]: vv_other_idx[v]
                                for k, v in isom_vert.items()}
            # Now we build all leg maps, again as a product of all maps at vertices.
            # Note: This also included the previously assigned vertices (with
            # marked points...)
            l_isom_list = []
            for v in isom_vert:
                # Construct leg maps:
                # We work with legs_source and legs_target, i.e. the list
                # of legs with the marked points removed.
                legal_l_maps = []
                leg_maps(legs_source[v], legs_target[isom_vert[v]])
                l_isom_list.append(legal_l_maps[:])  # copy!
            for l_maps in itertools.product(*l_isom_list):
                for l_map in l_maps:
                    isom_legs.update(l_map)
                yield return_isom_vert.copy(), isom_legs.copy()
