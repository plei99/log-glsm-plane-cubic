# pylint does not know sage
from sage.structure.sage_object import SageObject  # pylint: disable=import-error
from sage.misc.flatten import flatten  # pylint: disable=import-error

import admcycles.diffstrata.generalisedstratum
import admcycles.diffstrata.levelgraph
import admcycles.diffstrata.embeddedlevelgraph

#######################################################################
#######################################################################
# Point Dictionaries
#######################################################################
# * Points on generalised strata are of the form (i,j) where i is the
# index of the component and j is the index in the signature of
# the component.
# * An EmbeddedLevelGraph is a LevelGraph together with a dictionary
# mapping all non-half-edges (i.e. marked points) to the the points
# of the associated stratum.
# Note that the numbering of points starts with 1 while the indices
# of the signatures start with 0. I.e. in the case of a gen. stratum
# with only one component, the dictionary is i -> (0,i-1).
# * For clutching, we need dictionaries that map points of strata to
# points of other strata. This allows us to embed clutched LevelGraphs
# by composition.
# * When splitting a BIC into top and bottom, we consequently need to
# remember the dictionary that maps marked points of the top and
# bottom  gen. stratum into the points of stratum the BIC was
# embedded into.
##
# Thus, a splitting of a BIC inside stratum X must result in two generalised
# strata, top and bottom, as well as three dictionaries:
# * emb_top: points of top -> points of X
# * emb_bot: points of bottom -> points of X
# * clutching: points of top -> points of bottom
# All of these dictionaries are injective, the keys of emp_top and clutching
# give a partition of the points of top, the keys of emb_bot and the values
# of clutching give a partition of the points of bottom and the values of
# emb_top and emb_bot give a partition of the points of X.
#######################################################################
#######################################################################


def clutch(
        X,
        top,
        bottom,
        clutch_dict,
        emb_dict_top,
        emb_dict_bot,
        middle=None,
        clutch_dict_lower={},
        clutch_dict_long={},
        emb_dict_mid={}):
    """
    Clutch EmbeddedLevelGraphs top and bottom to get a new EmbeddedLevelGraph.

    If middle is specified, we do a "double clutching" top - middle - bottom to
    give a graph in X.

    .. NOTE::

        * this cannot be split into two steps, as the intermediate graph is
          not a graph in X!
        * We include the possibility to clutch points in the top component *past*
          the middle to the bottom component (long edges!)

    This is the inverse procedure to EmbeddedLevelGraph.split() and
    GeneralisedStratum.doublesplit().

    Note that if top or bottom is None, we return None.

    INPUT:

    X (GeneralisedStratum): Stratum that the clutched graph will live in

    top (EmbeddedLevelGraph): Top component to be clutched

    bottom (EmbeddedLevelGraph): Bottom component to be clutched

    middle (EmbeddedLevelGraph, optional, defaults to None): Middle component to
    be clutched.

    clutch_dict (dict): The dictionary giving the clutching information,
    (i.e. the edges that will be inserted).
    In the case of a "simple clutch", i.e. top to bottom with no middle,
    a dictionary:

    * point on top stratum -> point on bottom stratum
      In the case of a middle graph, i.e. a "double clutching":
    * point on top stratum -> point on middle stratum.
      In this case, we also need clutch_dict_lower for the other clutching
      and clutch_dict_long for the long edges.

    clutch_dict_lower (dict, optional): In the case of a "double clutching", the
    info for the "lower" clutching, i.e. a dictionary
    point on middle stratum -> point on bottom stratum

    clutch_dict_long (dict, optional): In the case of a "double clutching", the
    info for the "long edges", i.e. a dictionary
    point on top stratum -> point on bottom stratum

    emb_dict_top (dict): dictionary giving the embedding of the resulting
    EmbeddedLevelGraph, i.e. points on top stratum -> point on
    enveloping stratum.

    emb_dict_bot (dict): dictionary giving the embedding of the resulting
    EmbeddedLevelGraph, i.e. points on bottom stratum -> point on
    enveloping stratum.

    emb_dict_mid (dict, optional): dictionary giving the embedding of the
    resulting EmbeddedLevelGraph, i.e. points on middle stratum -> point on
    enveloping stratum.

    OUTPUT:

    EmbeddedLevelGraph: Clutched LevelGraph together with embedding.

    Returns:
        tuple: EmbeddedLevelGraph: as above

    EXAMPLES::

        sage: from admcycles.diffstrata import *
        sage: X=GeneralisedStratum([Signature((1,1))])
        sage: assert clutch(**X.bics[1].split()).is_isomorphic(X.bics[1])
        sage: assert all(clutch(**B.split()).is_isomorphic(B) for B in X.bics)

        The same works for 3-level graphs and doublesplit:

        sage: X=GeneralisedStratum([Signature((1,1))])
        sage: assert all(X.lookup_graph(*ep).is_isomorphic(clutch(**X.doublesplit(ep))) for ep in X.enhanced_profiles_of_length(2))
        sage: X=GeneralisedStratum([Signature((2,2,-2))])
        sage: assert all(X.lookup_graph(*ep).is_isomorphic(clutch(**X.doublesplit(ep))) for ep in X.enhanced_profiles_of_length(2))  # long time (4 seconds)
    """
    if top is None or bottom is None:
        return None
    # Build new EmbeddedLevelGraph:
    # This is somewhat similar to levelstratum's unite_embedded_graphs
    # but here the result is embedded into a _different_ stratum!
    # In case top, middle or bottom are LevelStrata, we replace them by their smooth
    # component (this is mainly to allow clutch to eat the output of split...)
    if isinstance(top, admcycles.diffstrata.generalisedstratum.LevelStratum):
        top = top.smooth_LG
    if isinstance(
            bottom, admcycles.diffstrata.generalisedstratum.LevelStratum):
        bottom = bottom.smooth_LG
    if isinstance(
            middle, admcycles.diffstrata.generalisedstratum.LevelStratum):
        # in particular, middle is not None!
        middle = middle.smooth_LG
    # genera simply get combined
    if middle is None:
        newgenera = top.LG.genera + bottom.LG.genera
    else:
        newgenera = top.LG.genera + middle.LG.genera + bottom.LG.genera
    # We renumber the levels consecutively starting at 0.
    top_levels = top.LG.numberoflevels()
    top_level_dict = {top.LG.internal_level_number(
        l): -l for l in range(top_levels)}
    newlevels = [top_level_dict[l] for l in top.LG.levels]
    if middle is None:
        mid_levels = 0
        mid_level_dict = {}
    else:
        # the levels of the middle component have to be shifted:
        mid_levels = middle.LG.numberoflevels()
        mid_level_dict = {middle.LG.internal_level_number(l): -l - top_levels
                          for l in range(mid_levels)}
        newlevels += [mid_level_dict[l] for l in middle.LG.levels]
    bot_levels = bottom.LG.numberoflevels()
    # The levels of the bottom component have to be shifted by both top and
    # middle:
    bot_level_dict = {
        bottom.LG.internal_level_number(l): -
        l -
        top_levels -
        mid_levels for l in range(bot_levels)}
    newlevels += [bot_level_dict[l] for l in bottom.LG.levels]
    # TODO: What is a sensible level dictionary?
    # At the moment, we choose the identity.
    newdlevels = {-l: -l for l in range(top_levels + mid_levels + bot_levels)}
    # The legs have to be renumbered in analogy to unite_embedded_graphs:
    newlegs = []
    newedges = []
    newpoleorders = {}
    newdmp = {}  # the new embedding is given by the emb_dicts
    # For adding the new edges, we have to record the new leg numbers and this will
    # happen at different stages (on each component).
    # The idea is, for each set of edges, to create a list of the form
    # [(starting leg 1, starting leg 2,....), (ending leg 1, ending leg 2,...)]
    # and then zip these lists to create the edges.
    # The points are currently stored (as points on the respective stratum) in the
    # various clutching dictionaries.
    # We want to move them to lists of the form described above.
    if middle is None:
        # If this is a "simple clutch", we only have to add edges between top
        # and bottom.
        clutch_edges = []  # the new edges to be added
    else:
        # Otherwise, we have to distinguish between:
        # * edges from top to middle:
        clutch_edges_tm = []
        # * edges from middle to bottom:
        clutch_edges_mb = []
        # * and edges from top to bottom (going past middle):
        clutch_edges_tb = []
    leg_number_shift = 0  # the shift for renumbering legs
    legs = 0  # leg counter
    # We now go through the components, starting with the top:
    #
    # Note that we also have to choose the appropriate embedding dictionary.
    # We will also have to record the clutching information, as the points are all
    # going to be renamed, so we have to remember who will be clutched.
    #
    if middle is None:
        # in the "simple clutch" case, it makes sense to also pass
        # the clutch points to each component; these are stored in clutch_dict
        comp_data = ((top, emb_dict_top, clutch_dict.keys()),
                     (bottom, emb_dict_bot, clutch_dict.values()))
    else:
        # We have to sort out the clutching in the loop anyway, so we just pass
        # a dummy.
        comp_data = ((top, emb_dict_top, {}),
                     (middle, emb_dict_mid, {}),
                     (bottom, emb_dict_bot, {}))
    # emb_g: graph to be clutched
    # emb_dict: dictionary
    #   marked points of the enveloping stratum of emb_g -> points of X
    for emb_g, emb_dict, clutch_points in comp_data:
        leg_dict = {}  # old number -> new number
        for i, l in enumerate(flatten(emb_g.LG.legs)):
            # the point numbers on bottom component are shifted by the top
            # number
            newlegnumber = leg_number_shift + i + 1
            leg_dict[l] = newlegnumber
            # while we're at it, we add the pole orders:
            newpoleorders[newlegnumber] = emb_g.LG.poleorders[l]
            legs += 1
            # For the dictionary of marked points (for the embedding), we
            # must distinguish if this is a marked point or a half-edge.
            # Marked points are simply the ones for which we have a key
            # in emb_dict :-)
            # The newdmp must then map this point to the image of its dmp in
            # the new stratum, i.e. the image under emb_dict:
            # Note that emb_dict eats stratum points, i.e. the image of the
            # leg under dmp.
            try:
                newdmp[newlegnumber] = emb_dict[emb_g.dmp[l]]
            except KeyError:
                pass
        leg_number_shift += legs
        # build (nested) list of legs:
        newlegs += [[leg_dict[l] for l in comp] for comp in emb_g.LG.legs]
        # finally, the edges are renumbered accordingly:
        newedges += [(leg_dict[e[0]], leg_dict[e[1]]) for e in emb_g.LG.edges]
        # We now record the clutching information:
        # Note that the points to be clutched are points on the Stratum!
        #
        # Because everything was renamed with leg_dict, we need to *now* remember
        # the edges to clutch (see discussion above).
        #
        # For this, we now have to distinguish if we're working with
        # two or three components:
        # * In the two component case, the points to be clutched are
        #   the keys (top) and values (bottom) of clutch_dict.
        #   Said differently: clutch_dict encodes the edges that will be added.
        if middle is None:
            # We save "half" of the clutch points in each step, afterwards we just
            # have to combine these lists to obtain all new edges.
            clutch_edges.append(
                tuple(leg_dict[emb_g.dmp_inv[p]] for p in clutch_points))
        # * In the three component case,
        #       * the *keys* of clutch_dict are points of the *top* component
        #       * the *values* clutch_dict are those points on the *middle* component that
        #           are clutched to the *top* component
        #       * the *keys* of clutch_dict_lower are points of the *middle* component
        #       * the *values* of clutch_dict_lower are points of the *bottom* component that
        #           are clutched to the *middle* component
        #       * the *keys* of clutch_dict_long are points of the *top* component
        #       * the *values* of clutch_dict_long are points of the *bottom* component that
        #           are clutched to the *top* component (long edges going past the middle)
        #   Said differently:
        #       * clutch_dict encodes edges top - middle
        #       * clutch_dict_lower encodes edges middle - bottom
        #       * clutch_dict_long encodes edges top - bottom (long edges going past middle)
        else:
            # We save the tuple of legs on the current component:
            if emb_g is top:
                clutch_edges_tm.append(tuple(
                    leg_dict[emb_g.dmp_inv[p]] for p in clutch_dict.keys()
                ))
                clutch_edges_tb.append(tuple(
                    leg_dict[emb_g.dmp_inv[p]] for p in clutch_dict_long.keys()
                ))
            elif emb_g is middle:
                clutch_edges_tm.append(tuple(
                    leg_dict[emb_g.dmp_inv[p]] for p in clutch_dict.values()
                ))
                clutch_edges_mb.append(tuple(
                    leg_dict[emb_g.dmp_inv[p]] for p in clutch_dict_lower.keys()
                ))
            else:
                assert emb_g is bottom
                clutch_edges_tb.append(tuple(
                    leg_dict[emb_g.dmp_inv[p]] for p in clutch_dict_long.values()
                ))
                clutch_edges_mb.append(tuple(
                    leg_dict[emb_g.dmp_inv[p]] for p in clutch_dict_lower.values()
                ))
    # Now the actual clutching happens, i.e. adding in the extra edges between
    # the points of clutch_dict:
    if middle is None:
        assert len(clutch_edges) == 2
        newedges += list(zip(*clutch_edges))
    else:
        assert len(clutch_edges_tm) == 2
        assert len(clutch_edges_tb) == 2
        assert len(clutch_edges_mb) == 2
        newedges += list(zip(*clutch_edges_tm)) + \
            list(zip(*clutch_edges_tb)) + \
            list(zip(*clutch_edges_mb))

    LG = admcycles.diffstrata.levelgraph.LevelGraph(
        newgenera, newlegs, newedges, newpoleorders, newlevels
    )

    ELG = admcycles.diffstrata.embeddedlevelgraph.EmbeddedLevelGraph(
        X, LG, newdmp, newdlevels)

    if not ELG.is_legal():
        # This should not happen if R-GRC is checked!
        # print("Illegal clutch attempted: %r, %r" % (top,bottom))
        raise RuntimeError("Illegal clutch, this should not happen anymore!")

    return ELG


class GenDegenerationGraph(SageObject):
    """
    The degeneration graph of the Generalised Stratum X.

    A degeneration graph of a (generalised) stratum is built recursively.
    For that we need to know:

    * The BICs inside X
    * For each BIC its top and bottom component (generalised strata)
    * For each BIC the clutching map of these components
    * For each top and bottom component of each BIC an embedding of their
      respective BICs into BIC(X)

    This allows us to translate products of BICs on X into clutchings of
    products of BICs on top and bottom components and gives a recursive
    procedure to compute the LevelGraph(s) for any product of BICs in X.

    INPUT:

    X (GeneralisedStratum): The GeneralisedStratum we give the degeneration
    graph of.
    """
    def __init__(self, X):
        self.X = X
        self.n = len(self.X.bics)  # possibly generates all bics of X (!)
        self._top_to_bic = [None for _ in range(
            self.n)]  # initialise on demand
        self._bot_to_bic = [None for _ in range(
            self.n)]  # initialise on demand
        self._middle_to_bic = {}  # initialise on demand
        self._top_to_bic_inv = [None for _ in range(
            self.n)]  # initialise on demand
        self._bot_to_bic_inv = [None for _ in range(
            self.n)]  # initialise on demand
        self._middle_to_bic_inv = {}  # initialise on demand
        # Note that while top_to_bic and bot_to_bic aren't injective, their
        # images are disjoint, i.e. we can remove the guys we found for top
        # from the ones we're searching through for bot (and vice versa).
        # This happens for every bic, though.
        self._unused_bics_top = [
            list(range(self.n)) for _ in range(self.n)]
        self._unused_bics_bot = [
            list(range(self.n)) for _ in range(self.n)]
        # To build all LevelGraphs inside a (generalised) stratum, we store all
        # tuples of indices of bics that give a (non-empty!) 3-level graph:
        # Note that the multiplicities tell us something about automorphisms
        # so we keep track of them and the set in parallel (use
        # _append_triple_prod() !!)
        self._triple_prods = []
        self._triple_prod_set = set()
        self._triple_prod_set_list = None

    def __repr__(self):
        return "GenDegenerationGraph(X=%r)" % self.X

    def _append_triple_prod(self, t):
        # append the tuple t to _triple_prods (both the list and the set...)
        self._triple_prods.append(t)
        self._triple_prod_set.add(t)

    @property
    def triple_prods(self):
        """
        A list of profiles of the 3-level graphs in X.

        Note that this is generated on first call.

        Returns:
            list: list of tuples of length two of indices of X.bics.
        """
        if self._triple_prod_set_list is None:
            # the easiest way to initialise is to initialise all bic maps:
            for i in range(self.n):
                self.top_to_bic(i)
                self.bot_to_bic(i)
            self._triple_prod_set_list = list(self._triple_prod_set)
        return self._triple_prod_set_list

    def top_to_bic(self, i):
        """
        A dictionary indices of BICs of X.bic[i].top -> indices of X.bics.

        Note that this is created once and then reused.

        Args:
            i (int): Index of X.bics

        Raises:
            RuntimeError: Raised if the squish of the clutching of a degeneration of
                the top component with X.bic[i].bot is not found in X.bics.

        Returns:
            dict: Dictionary mapping BICs of X.bic[i].top to indices of X.bics.
        """
        if self._top_to_bic[i] is None:
            B = self.X.bics[i]
            top_to_bic = {}  # index of B.top.bics -> index of X.bics
            for j, Bt in enumerate(B.top.bics):
                # Note that Bt is embedded (as a Graph) into B.top
                # B.top is embedded (as a stratum) into B (via B.emb_top)
                G = clutch(self.X, Bt, B.bot, B.clutch_dict,
                           B.emb_top, B.emb_bot)
                if G is None:
                    # illegal clutch:
                    raise RuntimeError(
                        "Illegal clutch, this should not happen anymore!")
                # Because the degeneration is clutched on top, squishing the
                # bottom level gives a new divisor (i.e. delta(1))
                squished_bic = G.delta(1)
                for im_j in self._unused_bics_top[i]:
                    if im_j == i:
                        continue  # B^2 is empty
                    if squished_bic.is_isomorphic(self.X.bics[im_j]):
                        top_to_bic[j] = im_j
                        # Record profile of the 3-level graph:
                        self._append_triple_prod((im_j, i))
                        if self._bot_to_bic[i] is None:
                            # otherwise bot_to_bic has already been constructed and
                            # this is futile...
                            try:
                                self._unused_bics_bot[i].remove(im_j)
                            except ValueError:
                                # this is a case where the map is non-injective
                                pass
                        break
                else:  # no break
                    raise RuntimeError(
                        "delta(1) of %s, %s, not found in list of BICs of %s" %
                        (G, squished_bic, self.X))
            self._top_to_bic[i] = top_to_bic.copy()
        return self._top_to_bic[i]

    def bot_to_bic(self, i):
        """
        A dictionary indices of BICs of X.bic[i].bot -> indices of X.bics.

        Note that this is created once and then reused.

        Args:
            i (int): Index of X.bics

        Raises:
            RuntimeError: Raised if the squish of the clutching of a degeneration of
                the bottom component with X.bic[i].top is not found in X.bics.

        Returns:
            dict: Dictionary mapping BICs of X.bic[i].bot to indices of X.bics.
        """
        if self._bot_to_bic[i] is None:
            B = self.X.bics[i]
            bot_to_bic = {}  # index of B.top.bics -> index of X.bics
            for j, Bb in enumerate(B.bot.bics):
                # Note that Bb is embedded (as a Graph) into B.bot
                # B.bot is embedded (as a stratum) into B (via B.emb_bot)
                G = clutch(self.X, B.top, Bb, B.clutch_dict,
                           B.emb_top, B.emb_bot)
                if G is None:
                    # illegal clutch:
                    raise RuntimeError(
                        "Illegal clutch, this should not happen anymore!")
                # Because the degeneration is clutched on bottom, squishing the
                # top level gives a new divisor (i.e. delta(2))
                squished_bic = G.delta(2)
                for im_j in self._unused_bics_bot[i]:
                    if im_j == i:
                        continue  # B^2 is empty
                    if squished_bic.is_isomorphic(self.X.bics[im_j]):
                        bot_to_bic[j] = im_j
                        # Record profile of the 3-level graph:
                        self._append_triple_prod((i, im_j))
                        if self._top_to_bic[i] is None:
                            # otherwise top_to_bic has already been constructed and
                            # this is futile...
                            try:
                                self._unused_bics_top[i].remove(im_j)
                            except ValueError:
                                # this is a case where the map is non-injective
                                pass
                        break
                else:  # no break
                    raise RuntimeError(
                        "delta(2) of %s, %s, not found in list of BICs of %s" %
                        (G, squished_bic, self.X))
            self._bot_to_bic[i] = bot_to_bic.copy()
        return self._bot_to_bic[i]

    def middle_to_bic(self, enh_profile):
        """
        A dictionary for each 3-level graph mapping (indices of) BICs of the middle level
        to (indices of) BICs of the enveloping stratum.

        Raises a RuntimeError if

        * either enh_profile is not a 3-level graph
        * or if a degeneration of the middle level is isomorphic
          (after squishing) to more than one BIC of X.
        * or if no isomorphic BIC of X is found.

        INPUT:

        enh_profile (tuple): Enhanced profile of a 3-level graph.

        OUTPUT:

        dictionary: indices of bics of the middle level -> indices of X.bics
        """
        if enh_profile not in self._middle_to_bic:
            p, i = enh_profile
            if len(p) != 2:
                raise RuntimeError(
                    "Only 3-level graphs have a well-defined middle! %r" %
                    (enh_profile,))
            # The dictionary we want to create:
            mid_to_bic = {}
            # The easiest thing is to split the graph and reuse all the clutching
            # info this yields, replacing the middle level by one of its BICs:
            clutching_info = self.X.doublesplit(enh_profile)
            L = clutching_info['middle']
            # candidates are those BICs that appear as bottom degenerations
            # of the top level *and* as top degenerations of the bottom level:
            candidates = set()
            top_deg = set(self.top_to_bic(p[1]).values())
            for b in self.bot_to_bic(p[0]).values():
                if b in top_deg:
                    candidates.add(b)
            # We now clutch in degenerations of the middle level:
            for i, B in enumerate(L.bics):
                clutching_info['middle'] = B
                H = clutch(**clutching_info)
                # We now apply to delta to find the new BIC.
                # Note that delta is shifted with regard to the profile numbering.
                # It also allows us to recover p:
                assert H.delta(1).is_isomorphic(self.X.bics[p[0]])
                assert H.delta(3).is_isomorphic(self.X.bics[p[1]])
                # delta(2) is the one we are interested in:
                XB = H.delta(2)
                for b in candidates:
                    if XB.is_isomorphic(self.X.bics[b]):
                        # check if more than one is found
                        if i in mid_to_bic:
                            raise RuntimeError(
                                "BIC %r of %r seems to be isomorphic to several BICs of %r!" %
                                (i, L, self.X))
                        mid_to_bic[i] = b
                if i not in mid_to_bic:
                    raise RuntimeError(
                        "BIC %r of %r not found in BICs of %r!" %
                        (i, L, self.X))
            self._middle_to_bic[enh_profile] = mid_to_bic.copy()
        return self._middle_to_bic[enh_profile]

    def top_to_bic_inv(self, i):
        """
        Inverse of top_to_bic.

        More Precisely: A dictionary assigning indices of X.bics a list of indices of
        X.top.bics. Note that this can be more than one (if the intersection is not
        irreducible).

        Note also that not all indices of X.bics are keys (but if not, they should be
        keys of bot_to_bic_inv)!

        Args:
            i (int): index of X.bics

        Returns:
            dict: index of X.bics -> list of indices of X.top.bics.
        """
        if self._top_to_bic_inv[i] is None:
            self._top_to_bic_inv[i] = {}
            d = self.top_to_bic(i)  # possibly built here (!)
            for j, im_j in d.items():
                try:
                    self._top_to_bic_inv[i][im_j].append(j)
                except KeyError:
                    self._top_to_bic_inv[i][im_j] = [j]
        return self._top_to_bic_inv[i]

    def bot_to_bic_inv(self, i):
        """
        Inverse of :meth:`bot_to_bic`.

        More Precisely: A dictionary assigning indices of X.bics a list of indices of
        X.bot.bics. Note that this can be more than one (if the intersection is not
        irreducible).

        Note also that not all indices of X.bics are keys (but if not, they should be
        keys of top_to_bic_inv)!

        Args:
            i (int): index of X.bics

        Returns:
            dict: index of X.bics -> list of indices of X.bot.bics.
        """
        if self._bot_to_bic_inv[i] is None:
            self._bot_to_bic_inv[i] = {}
            d = self.bot_to_bic(i)  # possibly built here (!)
            for j, im_j in d.items():
                try:
                    self._bot_to_bic_inv[i][im_j].append(j)
                except KeyError:
                    self._bot_to_bic_inv[i][im_j] = [j]
        return self._bot_to_bic_inv[i]

    def middle_to_bic_inv(self, enh_profile):
        """
        Inverse of middle_to_bic.

        More Precisely: A dictionary assigning indices of X.bics a list of indices of
        {middle level of enh_profile}.bics.
        Note that this can be more than one (if the intersection is not
        irreducible).

        INPUT:

        enh_profile (tuple): enhanced profile of a 3-level graph.

        OUTPUT:

        index of X.bics -> list of indices of {middle level of enh_profile}.bics.
        """
        if enh_profile not in self._middle_to_bic_inv:
            self._middle_to_bic_inv[enh_profile] = {}
            d = self.middle_to_bic(enh_profile)  # possibly built here (!)
            for j, im_j in d.items():
                try:
                    self._middle_to_bic_inv[enh_profile][im_j].append(j)
                except KeyError:
                    self._middle_to_bic_inv[enh_profile][im_j] = [j]
        return self._middle_to_bic_inv[enh_profile]
