from copy import deepcopy

# pylint does not know sage
from sage.modules.free_module import FreeModule  # pylint: disable=import-error
from sage.rings.integer_ring import ZZ  # pylint: disable=import-error
from sage.arith.functions import lcm  # pylint: disable=import-error
from sage.misc.flatten import flatten  # pylint: disable=import-error
from sage.misc.cachefunc import cached_method  # pylint: disable=import-error

import admcycles.diffstrata.generalisedstratum
from admcycles.diffstrata.sig import Signature
from admcycles.diffstrata.klevelgraph import KLevelGraph


def smooth_LG(sig):
    """
    The smooth (i.e. one vertex) LevelGraph in the stratum (sig).

    INPUT:

    sig (Signature): signature of the stratum.

    OUTPUT:

    LevelGraph: The smooth graph in the stratum.
    """
    if sig.k != 1:
        raise ValueError("The signature is not for an abelian differential")
    return LevelGraph.fromPOlist(
        [sig.g], [list(range(1, sig.n + 1))], [], sig.sig, [0])


class LevelGraph(KLevelGraph):
    r"""
    Create a (stable) level graph.

    ..NOTE::

        This is a low-level class and should NEVER be invoced directly!
        Preferably, EmbeddedLevelGraphs should be used and these should be
        generated automatically by Stratum (or GeneralisedStratum).

    .. NOTE::

        We don't inherit from stgraph/StableGraph anymore, as LevelGraphs
        should be static objects!

    Extends admcycles stgraph to represent a level graph as a stgraph,
    i.e. a list of genera, a list of legs and a list of edges, plus a list of
    poleorders for each leg and a list of levels for each vertex.

    Note that the class will warn if the data is not admissible, i.e. if the
    graph is not stable or the pole orders at separating nodes across levels do
    not add up to -2 or -1 on the same level (unless this is suppressed with
    quiet=True).

    INPUT:

    genera : list
    List of genera of the vertices of length m.

    legs : list
    List of length m, where ith entry is list of legs attached to vertex i.
    By convention, legs are unique positive integers.

    edges : list
    List of edges of the graph. Each edge is a 2-tuple of legs.

    poleorders : dictionary
    Dictionary of the form leg number : poleorder

    levels : list
    List of length m, where the ith entry corresponds to the level of
    vertex i.
    By convention, top level is 0 and levels go down by 1, in particular
    they are NEGATIVE numbers!

    quiet : optional boolean (default = False)
    Suppresses most error messages.

    ALTERNATIVELY, the pole orders can be supplied as a list by calling
    LevelGraph.fromPOlist:

    poleorders : list
    List of order of the zero (+) or pole (-) at each leg. The ith element
    of the list corresponds to the order at the leg with the marking i+1
    (because lists start at 0 and the legs are positive integers).

    EXAMPLES:

    Creating a level graph with three components on different levels of genus 1,
    3 and 0. The bottom level has a horizontal node.::

        sage: from admcycles.diffstrata import *
        sage: G = LevelGraph.fromPOlist([1,3,0],[[1,2],[3,4,5],[6,7,8,9]],[(2,3),(5,6),(7,8)],[2,-2,0,6,-2,0,-1,-1,0],[-2,-1,0]); G  #  doctest: +SKIP
        LevelGraph([1, 3, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9]],[(3, 2), (6, 5), (7, 8)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0},[-2, -1, 0],True)

    or alternatively::

        sage: LevelGraph([1, 3, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9]],[(3, 2), (6, 5), (7, 8)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0},[-2, -1, 0],quiet=True)
        LevelGraph([1, 3, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9]],[(3, 2), (6, 5), (7, 8)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0},[-2, -1, 0],True)

    We get a warning if the graph has non-stable components: (not any more ;-))::

        sage: G = LevelGraph.fromPOlist([1,3,0,0],[[1,2],[3,4,5],[6,7,8,9],[10]],[(2,3),(5,6),(7,8),(9,10)],[2,-2,0,6,-2,0,-1,-1,0,-2],[-3,-2,0,-1]); G  # doctest: +SKIP
        Warning: Component 3 is not stable: g = 0 but only 1 leg(s)!
        Warning: Graph not stable!
        LevelGraph([1, 3, 0, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9], [10]],[(3, 2), (6, 5), (7, 8), (9, 10)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0, 10: -2},[-3, -2, 0, -1],True)
    """
    def __init__(self, genera, legs, edges, poleorders, levels, quiet=False):
        super().__init__(
            genera, legs, edges, poleorders, levels, 1, quiet)
        # the prong orders are stored in a dictionary edge : prong order
        self.prongs = self._gen_prongs()

    @classmethod
    def fromPOlist(cls, genera, legs, edges,
                   poleordersaslist, levels, quiet=False):
        """
        This gives a LevelGraph where the poleorders are given as a list, not
        directly as a dictionary.
        """
        # translate poleorder list to dictionary:
        sortedlegs = sorted(flatten(legs))
        if len(sortedlegs) != len(poleordersaslist):
            raise ValueError("Numbers of legs and pole orders do not agree!" +
                             " Legs: " + repr(len(sortedlegs)) +
                             " Pole orders: " + repr(len(poleordersaslist)))
        poleorderdict = {sortedlegs[i]: poleordersaslist[i]
                         for i in range(len(sortedlegs))}
        return cls(genera, legs, edges, poleorderdict, levels, quiet)

    @classmethod
    def fromKLevelGraph(cls, kLG):
        assert kLG.sig.k == 1
        return cls(kLG.genera, kLG.legs, kLG.edges,
                   kLG.poleorders, kLG.levels, quiet=True)

    def __repr__(self):
        return "LevelGraph(" + self.input_as_string() + ")"

    def __hash__(self):
        return hash(repr(self))

    def __str__(self):
        return ("LevelGraph " + repr(self.genera) + ' ' + repr(self.legs) + ' '
                + repr(self.edges) + ' ' + repr(self.poleorders) + ' '
                + repr(self.levels))

    def input_as_string(self):
        """
        return a string that can be given as argument to __init__ to give self
        """
        return (repr(self.genera) + ',' + repr(self.legs) + ','
                + repr(self.edges) + ',' + repr(self.poleorders) + ','
                + repr(self.levels) + ',True')

    ##########################################################
    # Interface:
    # --- provide a series of methods that translate various data
    # (maybe this can be optimised by smart caching later, so only
    # these should be used...)
    # up to now, we have:
    # * prongs_list :
    # returns : list of tuples (edge,prong)
    # * prong :
    #            INPUT : edge
    # returns : prong order at edge
    # * _pointtype :
    #            INPUT : order (integer)
    # returns : pole/marked point/zero (string)
    # * is_meromorphic :
    # returns : boolean (check for pole among marked points)
    # * highestorderpole :
    #            INPUT : vertex
    # returns : leg name (at v) or -1 if none found
    #########################################################
    @cached_method
    def prongs_list(self):
        return list(self.prongs.items())

    @cached_method
    def prong(self, e):
        """
        The prong order is the pole order of the higher-level pole +1
        """
        if self.levelofleg(e[0]) > self.levelofleg(e[1]):
            return self.orderatleg(e[0]) + 1
        else:
            return self.orderatleg(e[1]) + 1

    @cached_method
    def _pointtype(self, order):
        if order < 0:
            return "pole"
        elif order == 0:
            return "marked point"
        else:
            return "zero"

    @cached_method
    def is_meromorphic(self):
        """
        Returns True iff at least one of the MARKED POINTS is a pole.
        """
        return any(self.orderatleg(l) < 0 for l in self.list_markings())

    @cached_method
    def highestorderpole(self, v):
        """
        Returns the leg with the highest order free pole at v, -1 if no poles found.
        """
        minorder = 0
        leg = -1
        for l in self.list_markings(v):
            if self.orderatleg(l) < minorder:
                minorder = self.orderatleg(l)
                leg = l
        return leg
    #########################################################
    # end interface
    #########################################################

    #################################################################
    # Cleanup functions
    # USE ONLY IN __init__!!! KLevelGraphs should be static!!!
    #################################################################
    def _gen_prongs(self):
        """
        Generate the dictionary edge : prong order.
        The prong order is the pole order of the higher-level pole +1
        """
        return {e: self.prong(e) for e in self.edges}

    #################################################################
    # Extract a subgraph
    #################################################################

    def extract(self, vertices, edges):
        """
        Extract the subgraph of self (as a LevelGraph) consisting of vertices
        (as list of indices) and edges (as list of edges).

        Returns the levelgraph consisting of vertices, edges and all legs on
        vertices (of self) with their original poleorders and the original
        level structure.
        """
        return self.fromKLevelGraph(
            super().extract(vertices, edges))

    def levelGraph_from_subgraph(self, G):
        """
        Returns the LevelGraph associated to a subgraph of underlying_graph
        (with the level structure induced by self)
        """
        return self.fromKLevelGraph(
            super().levelGraph_from_subgraph(G))

    def stratum_from_level(self, l, excluded_poles=None):
        """
        Return the LevelStratum at (relative) level l.

        INPUT:

        l (int): relative level number (0,...,codim)

        excluded_poles (tuple, defaults to None): a list of poles (legs of the graph,
        that are marked points of the stratum!) to be ignored for r-GRC.

        OUTPUT:

        LevelStratum: the LevelStratum, i.e.

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
        a leg dictionary.
        """
        if self.is_horizontal():
            print("Error: Cannot extract levels of a graph with horizontal edges!")
            return None
        internal_l = self.internal_level_number(l)
        # first, we extract the obvious information from self at this level:
        vv = self.verticesonlevel(internal_l)
        legs_on_level = []
        sigs = []
        res_cond = []
        leg_dict = {}
        for i, v in enumerate(vv):
            legs_on_level += [deepcopy(self.legsatvertex(v))]
            leg_dict.update([(n, (i, j))
                             for j, n in enumerate(legs_on_level[-1])])
            sig_v = Signature(tuple(self.orderatleg(leg)
                                    for leg in legs_on_level[-1]))
            sigs += [sig_v]
        # Now we hunt for residue conditions.
        # We work with the underlying graph to include the "level at infinity":
        # We "cheat" here, to catch the edges going up from level l: These are the
        # edges of the graph above level l-1 that have a vertex on level l:
        # (see below for parsing of these edges)
        ee = []
        for e in self.UG_above_level(internal_l - 1).edges(sort=True):
            for UG_vertex in e[:2]:
                if UG_vertex[2] == 'res':
                    continue
                if self.levelofvertex(UG_vertex[0]) == internal_l:
                    ee.append(e)
                    break
        G = self.UG_above_level(internal_l)
        conn_comp = G.connected_components_subgraphs()
        # Residue conditions are obtained by sorting the poles on this level
        # by connected component of G they are attached to (note that this
        # includes the "level at infinity"!)
        for c in conn_comp:
            # if there is a free pole in c, we ignore this component
            if self._redeemed_by_merom_in_comp(
                    c, excluded_poles=excluded_poles):
                continue
            res_cond_c = []
            # Otherwise we record all poles to attached to this component in one
            # residue condition:
            for e in ee:
                # We work with edges of the underlying graph here, so we have to
                # be a bit more careful. Recall that edges of the underlying graph are of
                # the form:
                # * (UG_vertex, UG_vertex, edge name) if it's an edge of LG or
                # * (UG_vertex, UG_vertex, resiedgej) if it's an edge to level infinity
                # Here UG_vertex is either of the form (vertex number, genus, 'LG') or
                # (i, 0, 'res') where
                # i is the number of the vertex at infinity and resiedgej is the edge
                # connecting the residue i to the leg j (of self).
                if e[0] in c or e[1] in c:
                    # res_cond contains tuples (i,j) where i is the index of the
                    # component and j the index of the pole in sig_i.
                    #
                    # We need find the vertex on this level (i.e. not in c)
                    if e[0] in c:
                        UG_vertex = e[1]
                    else:
                        UG_vertex = e[0]
                    # to find the pole on v, we have to distinguish if this is an
                    # edge of self or if it connects to a residue at infinity:
                    edge = e[2]
                    if 'res' in edge:
                        # the leg is the number after 'edge' in the edge
                        # string:
                        _, leg = edge.split('edge')
                        leg = int(leg)
                    else:
                        # edge consists of (leg_top,leg_bot) and leg_bot is on
                        # level i:
                        leg = edge[1]
                    # We retrieve the stratum point from leg_dict:
                    res_cond_c += [leg_dict[leg]]
            if res_cond_c:
                res_cond += [res_cond_c]
        return admcycles.diffstrata.generalisedstratum.LevelStratum(
            sigs, res_cond, leg_dict)

    #################################################################
    # Check Residue condition (via inconvenient vertices/edges)
    #################################################################

    @cached_method
    def is_inconvenient_vertex(self, v):
        """
        Check if vertex is inconvenient, i.e.
        * g = 0
        * no simple poles
        * there exists a pole order m_i such that m_i > sum(m_j)-p-1
        Return boolean
        """
        # inconvenient vertices are of genus 0 and have no simple poles.
        if self.genus(v) > 0 or len(self.simplepolesatvertex(v)) > 0:
            return False
        # check inequality
        ll = self.legsatvertex(v)
        poles = [l for l in ll if self.orderatleg(l) < 0]
        polesum = sum([-self.orderatleg(l) for l in poles])
        p = len(poles)
        for l in ll:
            if self.orderatleg(l) > polesum - p - 1:
                return True
        return False

    def _redeemed_by_merom_in_comp(self, G, excluded_poles=None):
        """
        Check if there is a pole in the subgraph G (intended to be a connected
        component above a vertex).

        excluded_poles are ignored (for r-GRC).

        Returns boolean (True = exists a pole)
        """
        if excluded_poles is None:
            excluded_poles = []
        for w in G.vertices(sort=False):
            if w[2] == 'res':  # vertex at infinity
                continue
            for l in self.list_markings(w[0]):
                if (self.orderatleg(l) < 0) and (l not in excluded_poles):
                    return True
        return False

    @cached_method
    def is_illegal_vertex(self, v, excluded_poles=None):
        """
        Check if vertex is inconvenient and is not redeemed, i.e.
        * v is inconvenient
        * there are no two separate edges to the same connected component above (i.e. loop above)
        * if meromorphic: v is not connected to higher level marked poles

        We can also pass a tuple (hashing!) of poles that are to be excluded because of
        residue conditions.

        Return boolean
        """
        if not self.is_inconvenient_vertex(v):
            return False
        l = self.levelofvertex(v)
        # in the underlying_graph, vertices are of the form (index in genera,
        # genus, 'LG'):
        v_graph = self.UG_vertex(v)
        # edges not going down from v:
        ee = [e for e in self.edgesatvertex(v) if self.levelofleg(e[0]) >= l
              and self.levelofleg(e[1]) >= l]
        # in the holomorphic case, legal inconvenient vertices need at least
        # two edges not going down
        if len(ee) < 1 and not self.is_meromorphic():
            return True
        # check if v has two edges into one connected component that don't lie below v,
        # i.e. if v can be redeemed:
        abovegraph = self.UG_above_level(l - 1)
        cc = self.underlying_graph.subgraph([w for w in abovegraph.connected_component_containing_vertex(
            v_graph) if w[2] == 'res' or self.levelofvertex(w[0]) >= l])
        cc.delete_vertex(v_graph)
        conn_comp = cc.connected_components_subgraphs()
        if excluded_poles is None:
            excluded_poles = []
        freepoles = len([l for l in self.list_markings(v) if (
            self.orderatleg(l) < 0) and (l not in excluded_poles)])
        for G in conn_comp:
            # edges from v to G:
            # (Careful: cc does not have any edges with v anymore, so we must use abovegraph!)
            eeG = [e for e in abovegraph.edges(sort=True)
                   if (e[0] == v_graph and e[1] in G.vertices(sort=False))
                   or (e[1] == v_graph and e[0] in G.vertices(sort=False))]
            if len(eeG) > 1:
                # redeemed :-)
                return False
            # in the meromorphic case, we also check if a "free" pole exists in
            # a connected component above v
            if self.is_meromorphic() and self._redeemed_by_merom_in_comp(G,
                                                                         excluded_poles=excluded_poles):
                freepoles += 1
        if freepoles >= 2:
            # redeemed :-)
            return False
        return True

    @cached_method
    def is_illegal_edge(self, e, excluded_poles=None):
        """
        Check if edge is illegal, i.e. if
        * e is horizontal (not loop) and
        * there no simple loop over e
        * there is no "free" pole over each end point

        excluded_poles may be a tuple (hashing!) of poles to be excluded for r-GRC.

        Return boolean
        """
        if not self.is_horizontal(e) or (e[0] == e[1]):
            return False
        # check if there is a simple loop above e, i.e. if the subgraph
        # above e is still connected after e is removed (e is not a bridge)
        l = self.levelofleg(e[0])
        # note that verticesabove checks for >l
        abovegraph = self.UG_above_level(l - 1)
        # edges are encoded by (self.vertex(e[0]),self.vertex(e[1]),e) in
        # underlying_graph!
        cut = abovegraph.is_cut_edge(
            (self.UG_vertex(
                self.vertex(
                    e[0])), self.UG_vertex(
                self.vertex(
                    e[1])), e))
        # True = e is a cut-edge, i.e. connected components increase when removing e, i.e. no loop above,
        #       i.e. illegal unless also meromorphic
        # False = graph has a loop, i.e. not illegal
        if not cut:
            return False
        else:
            if self.is_meromorphic():
                # if there are "free" poles above the endpoints, we are still
                # fine
                abovegraph.delete_edge(
                    (self.UG_vertex(
                        self.vertex(
                            e[0])), self.UG_vertex(
                        self.vertex(
                            e[1])), e))
                polesfound = 0
                for l in e:
                    v = self.vertex(l)
                    G = abovegraph.connected_component_containing_vertex(
                        self.UG_vertex(v))
                    if self._redeemed_by_merom_in_comp(
                            abovegraph.subgraph(G), excluded_poles=excluded_poles):
                        polesfound += 1
                if polesfound == 2:
                    # we are fine
                    return False
        # no redemption :/
        return True

    @cached_method
    def is_legal(self, quiet=False, excluded_poles=None):
        """
        Check if self has illegal vertices or edges.

        excluded_poles may be a tuple (hashing!) of poles to be excluded for r-GRC.

        Return boolean
        """
        for v in range(len(self.genera)):
            if self.is_illegal_vertex(v, excluded_poles=excluded_poles):
                if not quiet:
                    print("Vertex", v, "is illegal!")
                return False
        for e in self.edges:
            if self.is_illegal_edge(e, excluded_poles=excluded_poles):
                if not quiet:
                    print("Edge", e, "is illegal!")
                return False
        return True

    #################################################################
    # Squishing ---- i.e. contracting horizontal edges / levels
    #################################################################

    def squish_horizontal(self, e):
        """
        Squish the horizontal edge e and return the new graph.

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: G = LevelGraph.fromPOlist([1,1], [[1,2],[3,4]], [(2,3)], [1,-1,-1,1],[0,0])
            sage: G.squish_horizontal((2,3))
            LevelGraph([2],[[1, 4]],[],{1: 1, 4: 1},[0],True)
        """
        return LevelGraph(*self._squish_horizontal(e), quiet=True)

    def squish_vertical_slow(self, l, adddata, quiet=False):
        """
        Squish the level l (and the next lower level!) and return the new graph.
        If addata=True is specified, we additionally return a boolean that tells
        us if a level was squished and the legs that no longer exist.

        More precisely, adddata=True returns a tuple (G,boolean,legs) where
        * G is the (possibly non-contracted) graph and
        * boolean tells us if a level was squished
        * legs is the (possibly empty) list of legs that don't exist anymore

        Implementation:
        At the moment, we remember all the edges (ee) that connect level l to
        level l-1. We then create a new LevelGraph with the same data as self,
        the only difference being that all vertices of level l-1 have now been
        assigned level l. We then remove all (now horizontal!) edges in ee.

        In particular, all points and edges retain their numbering (only the level might have changed).

        WARNING: Level l-1 does not exist anymore afterwards!!!

        Downside: At the moment we get a warning for each edge in ee, as these
        don't have legal pole orders for being horizontal (so checkedgeorders
        complains each time the constructor is called :/).
        """
        return self.fromKLevelGraph(
            super().squish_vertical_slow(l, adddata, quiet))

    def squish_vertical(self, l, quiet=True):
        """
        Squish the level l (and the next lower level!) and return the new graph.

        WARNING: Level l-1 does not exist anymore afterwards!!!

        Args:
            l (int): (internal) level number
            quiet (bool, optional): No output. Defaults to True.

        Returns:
            LevelGraph: new levelgraph with one level less.

        EXAMPLES::

            sage: from admcycles.diffstrata import *
            sage: G = LevelGraph.fromPOlist([1,2], [[1],[2,3,4]], [(1,2)], [0,-2,1,3],[0,-1])
            sage: G.squish_vertical(0)
            LevelGraph([3],[[3, 4]],[],{3: 1, 4: 3},[0],True)
        """
        return LevelGraph(*self._squish_vertical(l, quiet=quiet), quiet=True)

    def delta(self, k, quiet=False):
        """
        Squish all levels except for the k-th.

        Note that delta(1) contracts everything except top-level and that
        the argument is interpreted via internal_level_number

        WARNING: Currently giving an out of range level (e.g.
            0 or >= maxlevel) squishes the entire graph!!!

        Return the corresponding divisor.
        """
        return self.fromKLevelGraph(super().delta(k, quiet))

    #################################################################
    # end Squishing
    #################################################################

    def relabel(self, legdict, tidyup=True):
        new_legs = []
        for vlegs in self.legs:
            list1 = [legdict[i] for i in vlegs]
            if tidyup:
                list1.sort()
            new_legs.append(list1)
        new_edges = [(legdict[e[0]], legdict[e[1]]) for e in self.edges]
        new_poleorders = {legdict[i]: j for i, j in self.poleorders.items()}
        if tidyup:
            new_edges.sort()
            new_poleorders = {a: b for a, b in sorted(new_poleorders.items())}
        newLG = LevelGraph(self.genera, new_legs, new_edges,
                           new_poleorders, self.levels)
        return newLG

    #################################################################
    # Twist groups
    #################################################################

    def twist_group(self, quiet=True, with_degree=False):
        """
        This should not be needed! The stack factor should be computed
        using AdditiveGenerator.stack_factor!!!

        Calculate the index of the simple twist group inside the twist group.

        with_degree=True: return a tuple (e,g) where
            * e = index of simple twist in twist
            * g = index of twist group in level rotation group

        """
        # horizontal edges => trivial twist group :-)
        if self.is_horizontal():
            return 1
        N = len(self.edges)
        M = FreeModule(ZZ, N)
        # kernel of projections to Z/prong Z
        K = M.submodule([M.basis()[i] * self.prong(e)
                         for i, e in enumerate(self.edges)])
        if not quiet:
            print("kernel:", K)
        # submodule generated by edges crossing level i
        # i.e. image of R^Gamma (level rotation group) in ZZ^edges
        t_basis = [sum([M.basis()[i] for i, e in enumerate(self.edges)
                        if self.crosseslevel(e, self.internal_level_number(l))])
                   for l in range(1, self.numberoflevels())]
        E = M.submodule(t_basis)
        if not quiet:
            print("t-module:", E)
        # simple twist group: lcm of delta(l) edge prongs*t_i
        deltas = [self.delta(l, True) for l in range(1, self.numberoflevels())]
        if not quiet:
            print("deltas:", deltas)
        ll = [lcm([d.prong(e) for e in d.edges]) for d in deltas]
        if not quiet:
            print("lcms:", ll)
        st_basis = [sum([ll[l - 1] * M.basis()[i] for i, e in enumerate(self.edges)
                         if self.crosseslevel(e, self.internal_level_number(l))])
                    for l in range(1, self.numberoflevels())]
        S = M.submodule(st_basis)
        if not quiet:
            print("simple twist group:")
            print(S)
            print("Intersection:")
            print(K.intersection(E))
        # K.intersection(E) = Twist group (in ZZ^edges)
        tw_gr = K.intersection(E)
        if with_degree:
            return (S.index_in(tw_gr), tw_gr.index_in(E))
        else:
            return S.index_in(tw_gr)
