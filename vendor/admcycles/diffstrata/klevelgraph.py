from copy import deepcopy

# pylint does not know sage
from sage.structure.sage_object import SageObject  # pylint: disable=import-error
from sage.rings.integer_ring import ZZ  # pylint: disable=import-error
from sage.misc.flatten import flatten  # pylint: disable=import-error
from sage.functions.generalized import sign  # pylint: disable=import-error
from sage.misc.cachefunc import cached_method  # pylint: disable=import-error
from sage.graphs.graph import Graph  # pylint: disable=import-error
from sage.plot.text import text  # pylint: disable=import-error

from admcycles.admcycles import GraphIsom
from admcycles.stable_graph import StableGraph as stgraph
from admcycles.diffstrata.sig import kSignature


class KLevelGraph(SageObject):
    r"""
    Create a (stable) level graph for a k-differential.

    .. NOTE::

        This is a low-level class and should NEVER be invoked directly!
        Preferably, ``EmbeddedLevelGraphs`` should be used and these should be
        generated automatically by ``Stratum`` (or ``GeneralisedStratum``).

    .. NOTE::

        We do not inherit from stgraph/StableGraph anymore, as KLevelGraphs
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

    k : int
    The order of the differential.

    quiet : optional boolean (default = False)
    Suppresses most error messages.

    ALTERNATIVELY, the pole orders can be supplied as a list by calling
    KLevelGraph.fromPOlist:
    poleorders : list
    List of order of the zero (+) or pole (-) at each leg. The ith element
    of the list corresponds to the order at the leg with the marking i+1
    (because lists start at 0 and the legs are positive integers).

    EXAMPLES:

    Creating a level graph with three components on different levels of genus 1,
    3 and 0. The bottom level has a horizontal node.::

        sage: from admcycles.diffstrata.klevelgraph import KLevelGraph
        sage: G = KLevelGraph.fromPOlist([1,3,0],[[1,2],[3,4,5],[6,7,8,9]],[(2,3),(5,6),(7,8)],[2,-2,0,6,-2,0,-1,-1,0],[-2,-1,0],1); G  #  doctest: +SKIP
        KLevelGraph([1, 3, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9]],[(3, 2), (6, 5), (7, 8)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0},[-2, -1, 0],1,True)

    or alternatively::

        sage: KLevelGraph([1, 3, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9]],[(3, 2), (6, 5), (7, 8)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0},[-2, -1, 0],1,quiet=True)
        KLevelGraph([1, 3, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9]],[(3, 2), (6, 5), (7, 8)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0},[-2, -1, 0],1,True)

    Creating a level graph with two components on different levels of genus 1 and 1 and a quadratic differential.::

        sage: KLevelGraph.fromPOlist([1,1],[[1],[2,3]],[(1,2)],[0,-4,4],[0,-1],2)
        KLevelGraph([1, 1],[[1], [2, 3]],[(1, 2)],{1: 0, 2: -4, 3: 4},[0, -1],2,True)

    We get a warning if the graph has non-stable components: (not any more ;-))::

        sage: G = KLevelGraph.fromPOlist([1,3,0,0],[[1,2],[3,4,5],[6,7,8,9],[10]],[(2,3),(5,6),(7,8),(9,10)],[2,-2,0,6,-2,0,-1,-1,0,-2],[-3,-2,0,-1],1); G  # doctest: +SKIP
        Warning: Component 3 is not stable: g = 0 but only 1 leg(s)!
        Warning: Graph not stable!
        KLevelGraph([1, 3, 0, 0],[[1, 2], [3, 4, 5], [6, 7, 8, 9], [10]],[(3, 2), (6, 5), (7, 8), (9, 10)],{1: 2, 2: -2, 3: 0, 4: 6, 5: -2, 6: 0, 7: -1, 8: -1, 9: 0, 10: -2},[-3, -2, 0, -1],1,True)
    """
    def __init__(self, genera, legs, edges,
                 poleorders, levels, k, quiet=False):
        checks = False
        if len(genera) != len(legs):
            raise ValueError('genera and legs must have the same length')
        self.genera = genera
        self.legs = legs
        self.edges = edges
        self.poleorders = poleorders
        self.levels = levels
        # the signature consists of the marked points that are not half-edges
        sig_list = tuple(poleorders[l] for l in self.list_markings())
        self.sig = kSignature(sig_list, k)
        # associated stgraph...
        self._stgraph = None
        self._init_more()  # initiate some things and make sure all is in order
        if checks:
            self.checkadmissible(quiet)
        self._has_long_edge = None
        self._is_long = {}
        self._is_bic = None

    @classmethod
    def fromPOlist(cls, genera, legs, edges,
                   poleordersaslist, levels, k, quiet=False):
        """
        This gives a KLevelGraph where the poleorders are given as a list, not
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
        return cls(genera, legs, edges, poleorderdict, levels, k, quiet)

    def __repr__(self):
        return "KLevelGraph(" + self.input_as_string() + ")"

    def __hash__(self):
        return hash(repr(self))

    def __str__(self):
        return (
            "KLevelGraph " +
            repr(
                self.genera) +
            ' ' +
            repr(
                self.legs) +
            ' ' +
            repr(
                self.edges) +
            ' ' +
            repr(
                self.poleorders) +
            ' ' +
            repr(
                self.levels) +
            ' ' +
            repr(
                self.k))

    def input_as_string(self):
        """
        return a string that can be given as argument to __init__ to give self
        """
        return (repr(self.genera) + ',' + repr(self.legs) + ',' +
                repr(self.edges) + ',' + repr(self.poleorders) + ',' +
                repr(self.levels) + ',' + repr(self.sig.k) + ',True')

    def __eq__(self, other):
        if not isinstance(other, KLevelGraph):
            return False
        return (
            self.genera == other.genera) and (
            self.levels == other.levels) and (
            self.legs == other.legs) and (
                set(
                    self.edges) == set(
                        other.edges)) and (
                            self.poleorders == other.poleorders) and (
                                self.sig.k == other.sig.k)

    # reimplementing the admcycles stuff we use:
    @cached_method
    def g(self):
        genus = sum(self.genera) + len(self.edges) - \
            len(self.genera) + ZZ.one()
        assert genus == self.sig.g, "Signature %r does not match %r!" % (
            self.sig, self)
        return genus

    @property
    def stgraph(self):
        if self._stgraph is None:
            self._stgraph = stgraph([int(g) for g in self.genera],
                                    [[int(l) for l in leg_list]
                                     for leg_list in self.legs],
                                    self.edges)
        return self._stgraph

    @cached_method
    def list_markings(self, v=None):
        r"""
        Return the list of markings (non-edge legs) of self at vertex v.
        """
        if v is None:
            return tuple([j for v in range(len(self.genera))
                          for j in self.list_markings(v)])
        s = set(self.legs[v])
        for e in self.edges:
            s -= set(e)
        return tuple(s)

    ##########################################################
    # Interface:
    #   --- provide a series of methods that translate various data
    #   (maybe this can be optimised by smart caching later, so only
    #    these should be used...)
    # up to now, we have:
    #    * vertex : (overloaded from admcycles)
    #            INPUT : leg
    #            returns : vertex of leg (index of genera)
    #    * vertex_list :
    #            returns : list of vertices (index of genera)
    #    * edges_list :
    #            returns : list of edges
    #    * levelofvertex :
    #            INPUT : index of genera
    #            returns : level number (negative)
    #    * levelofleg :
    #            INPUT : leg label (given as positive integer)
    #            returns : level number (negative)
    #    * orderatleg :
    #            INPUT : leg label (given as positive integer)
    #            returns : poleorder (integer)
    #    * ordersonvertex :
    #            INPUT : index of genera
    #            returns : list of poleorders (integers)
    #    * verticesonlevel :
    #            INPUT : level number (negative)
    #            returns : list of indices of genera
    #    * edgesatvertex :
    #            INPUT : index of genera
    #            returns : list of edges
    #    * legsatvertex :
    #            INPUT : index of genera
    #            returns : list of legs
    #    * simplepolesatvertex :
    #            INPUT : index of genera
    #            returns : list of legs with pole order -1
    #    * genus :
    #            INPUT : index of genera
    #            returns : genus
    #    * codim :
    #            returns : codimension
    #    * is_bic :
    #            returns : boolean
    #    * edgesatlevel :
    #            INPUT : level number (negative)
    #            returns : list of edges with at least one node at that level (or horizontal)
    #    * horizontaledgesatlevel :
    #            INPUT : level number (negative)
    #            returns : list of horizontal edges
    #    * nextlowerlevel :
    #            INPUT : level number
    #            returns : level number of next lower level
    #    * lowest_level :
    #            returns : level number of lowest level
    #    * is_horizontal :
    #            INPUT : edge (or none for whole graph)
    #            returns : boolean success
    #    * has_loop :
    #            INPUT : vertex
    #            returns : boolean success
    #    * edgesgoingupfromlevel :
    #            INPUT : level
    #            returns : list of edges with e[1] on level
    #    * verticesabove :
    #            INPUT : level
    #            returns : list of vertices with level > l
    #    * edgesabove :
    #            INPUT : level
    #            returns : list of edges with level of each vertex > l
    #    * crosseslevel :
    #            INPUT : edge, level
    #            returns : boolean
    #    * edgesgoingpastlevel :
    #            INPUT : level
    #            returns : list of edges with start level > l and
    #                      end level < l
    #########################################################
    @cached_method
    def vertex(self, leg):
        """
        The vertex (as index of genera) where leg is located.

        Args:
            leg (int): leg

        Returns:
            int: index of genera
        """
        return self.legdict[leg]

    @cached_method
    def vertex_list(self):
        return list(range(len(self.genera)))

    @cached_method
    def edges_list(self):
        return self.edges

    @cached_method
    def levelofvertex(self, v):
        return self.levels[v]

    @cached_method
    def levelofleg(self, leg):
        return self.levelofvertex(self.vertex(leg))

    @cached_method
    def orderatleg(self, leg):
        return self.poleorders[leg]

    @cached_method
    def ordersonvertex(self, v):
        return [self.orderatleg(leg) for leg in self.legsatvertex(v)]

    @cached_method
    def verticesonlevel(self, level):
        return [v for v in range(len(self.genera))
                if self.levelofvertex(v) == level]

    @cached_method
    def edgesatvertex(self, v):
        return [e for e in self.edges if self.vertex(e[0]) == v or
                self.vertex(e[1]) == v]

    @cached_method
    def legsatvertex(self, v):
        return self.legs[v]

    @cached_method
    def is_halfedge(self, leg):
        """
        Check if leg has an edge attached to it.
        """
        return any(e[0] == leg or e[1] == leg for e in self.edges)

    @cached_method
    def simplepolesatvertex(self, v):
        return [l for l in self.legsatvertex(v) if self.orderatleg(l) == -1]

    @cached_method
    def genus(self, v=None):
        """
        Return the genus of vertex v.

        If v is None, return genus of the complete KLevelGraph.
        """
        if v is None:
            return self.g()  # from admcycles
        else:
            return self.genera[v]

    @cached_method
    def numberoflevels(self):
        return len(set(self.levels))

    def is_bic(self):
        if self._is_bic is None:
            self._is_bic = self.numberoflevels() == 2 and not self.is_horizontal()
        return self._is_bic

    @cached_method
    def codim(self):
        """
        Return the codim = No. of levels of self + horizontal edges.
        """
        return (self.numberoflevels() - 1 +
                sum(1 for e in self.edges if self.is_horizontal(e)))

    @cached_method
    def edgesatlevel(self, level):
        """
        Return list of edges with at least one node at level.
        """
        return [e for e in self.edges if self.levelofleg(e[0]) == level or
                self.levelofleg(e[1]) == level]

    @cached_method
    def horizontaledgesatlevel(self, level):
        return [e for e in self.edgesatlevel(level) if self.is_horizontal(e)]

    @cached_method
    def nextlowerlevel(self, l):
        """
        Return the next lower level number or  False if no legal level
        (e.g. lowest level) is given as input.

        Point of discussion: Is it better to tidy up the level numbers
        to be directly ascending or not?

        Pro tidy: * this becomes obsolete ;)
                  * easier to check isomorphisms?
        Con tidy: * we "see" where we came from
                  * easier to undo/glue back in after cutting
        """
        try:
            llindex = self.sortedlevels.index(l) - 1
        except ValueError:
            return False
        if llindex == -1:
            return False
        else:
            return self.sortedlevels[llindex]

    @cached_method
    def internal_level_number(self, i):
        """
        Return the internal i-th level, e.g.

        self.levels = [0,-2,-5,-3]

        then

        internal_level_number(0) is 0
        internal_level_number(1) is -2
        internal_level_number(2) is -3
        internal_level_number(3) is -4

        Returns None if the level does not exist.
        """
        reference_levels = list(reversed(sorted(list(set(self.levels)))))
        i = abs(i)
        if i >= len(reference_levels):
            return None
        else:
            return reference_levels[i]

    @cached_method
    def level_number(self, l):
        """
        Return the relative level number (positive!) of l, e.g.

        self.levels = [0,-2,-5,-3]

        then

        level_number(0) is 0
        level_number(-2) is 1
        level_number(-3) is 2
        level_number(-5) is 3

        Returns None if the level does not exist.
        """
        reference_levels = list(reversed(sorted(list(set(self.levels)))))
        try:
            return reference_levels.index(l)
        except ValueError:
            return None

    @cached_method
    def is_long(self, e):
        """
        Check if edge e is long, i.e. passes through more than one level.
        """
        try:
            return self._is_long[e]
        except KeyError:
            il = abs(self.level_number(self.levelofleg(
                e[0])) - self.level_number(self.levelofleg(e[1]))) > 1
            self._is_long[e] = il
            return il

    @cached_method
    def has_long_edge(self):
        if self._has_long_edge is None:
            for e in self.edges:
                if self.is_long(e):
                    self._has_long_edge = True
                    break
            else:
                self._has_long_edge = False
        return self._has_long_edge

    @cached_method
    def lowest_level(self):
        return self.sortedlevels[0]

    @cached_method
    def is_horizontal(self, e=None):
        """
        Check if edge e is a horizontal edge or if self has a horizontal edge.
        """
        if e is None:
            return any(self.is_horizontal(e) for e in self.edges)
        if e not in self.edges:
            print("Warning: " + repr(e) + " is not an edge of this graph!")
            return False
        return self.levelofleg(e[0]) == self.levelofleg(e[1])

    @cached_method
    def has_loop(self, v):
        """
        Check if vertex v has a loop.
        """
        return any(self.vertex(e[0]) == self.vertex(e[1])
                   for e in self.edgesatvertex(v))

    @cached_method
    def edgesgoingupfromlevel(self, l):
        """
        Return the edges going up from level l.

        This uses that e[0] is not on a lower level than e[1]!
        """
        return [e for e in self.edges if self.levelofleg(e[1]) == l and
                not self.is_horizontal(e)]

    @cached_method
    def verticesabove(self, l):
        """
        Return list of all vertices above level l.

        If l is None, return all vertices.
        """
        if l is None:
            return list(range(len(self.genera)))
        return [v for v in range(len(self.genera))
                if self.levelofvertex(v) > l]

    @cached_method
    def edgesabove(self, l):
        """
        Return a list of all edges above level l, i.e. start and end vertex
        have level > l.
        """
        return [e for e in self.edges
                if self.levelofleg(e[0]) > l and self.levelofleg(e[1]) > l]

    @cached_method
    def crosseslevel(self, e, l):
        """
        Check if e crosses level l (i.e. starts > l and ends <= l)
        """
        return self.levelofleg(e[0]) > l and self.levelofleg(e[1]) <= l

    @cached_method
    def edgesgoingpastlevel(self, l):
        """
        Return list of edges that go from above level l to below level l.
        """
        return [e for e in self.edges
                if self.levelofleg(e[0]) > l > self.levelofleg(e[1])]
    #########################################################
    # end interface
    #########################################################

    #########################################################
    # Checks
    #########################################################
    @cached_method
    def checkadmissible(self, quiet=False):
        """
        Run all kinds of checks on the level graph. Currently:
         * Check if orders on each komponent add up to k*(2g-2)
         * Check if graph is stable
         * Check if orders at edges sum to -k*2
         * Check if orders respect level crossing
        """
        admissible = True
        if not self.checkorders(quiet):
            if not quiet:
                print("Warning: Illegal orders!")
            admissible = False
        if not self.is_stable(quiet):
            if not quiet:
                print("Warning: Graph not stable!")
            admissible = False
        if not self.checkedgeorders(quiet):
            if not quiet:
                print("Warning: Illegal orders at a node!")
            admissible = False
        return admissible

    @cached_method
    def checkorders(self, quiet=False):
        """
        Check if the orders add up to k*(2g-2) on each component.
        """
        for v, g in enumerate(self.genera):
            if sum(self.ordersonvertex(v)) != self.sig.k * (2 * g - 2):
                if not quiet:
                    print("Warning: Component " +
                          repr(v) +
                          " orders add up to " +
                          repr(sum(self.ordersonvertex(v))) +
                          " but component is of genus " +
                          repr(g))
                return False
        return True

    @cached_method
    def is_stable(self, quiet=False):
        """
        Check if graph is stable.
        """
        # total graph:
        e = 2 * self.genus() - 2 + len(self.leglist)
        if e < 0:
            if not quiet:
                print("Warning: Total graph not stable! 2g-2+n = " + repr(e))
            return False

        # components:
        for v in range(len(self.genera)):
            if 3 * self.genus(v) - 3 + len(self.legsatvertex(v)) < 0:
                if not quiet:
                    print("Warning: Component " +
                          repr(v) +
                          " is not stable: g = " +
                          repr(self.genus(v)) +
                          " but only " +
                          repr(len(self.legsatvertex(v))) +
                          " leg(s)!")
                return False
        return True

    @cached_method
    def checkedgeorders(self, quiet=False):
        """
        Check that the orders at nodes (i.e. edges) sum to -k*2 and that lower
        level has lower order.
        """
        for e in self.edges:
            orders = self.orderatleg(e[0]) + self.orderatleg(e[1])
            if orders != -self.sig.k * 2:
                if not quiet:
                    print(
                        "Warning: Orders at edge " +
                        repr(e) +
                        " add up to " +
                        repr(orders) +
                        " instead of -k*2!")
                return False
            # iff the pole order at e[0] is > the poleorder at e[1], e[0]
            # should be on a higher level than e[1]
            if (sign(self.orderatleg(e[0]) - self.orderatleg(e[1])) !=
                    sign(self.levelofleg(e[0]) - self.levelofleg(e[1]))):
                if not quiet:
                    print("Warning: Orders at edge " + repr(e) +
                          " do not respect level crossing!")
                    print("Poleorders are",
                          (self.orderatleg(e[0]), self.orderatleg(e[1])),
                          "but levels are",
                          (self.levelofleg(e[0]), self.levelofleg(e[1]))
                          )
                return False
        return True

    #################################################################
    # End checks
    #################################################################

    #################################################################
    # Cleanup functions
    # USE ONLY IN __init__!!! KLevelGraphs should be static!!!
    #################################################################
    def _init_more(self):
        """
        This should be used _only_ for initialisation!

        Make sure everything is in order, in particular:
        * sortedlevels is fine
        * leglist is fine
        * legdict is fine
        * maxleg is fine
        * maxlevel is fine
        * underlying graph is fine
        """
        self.sortedlevels = sorted(self.levels)
        # legs is a nested list:
        self.leglist = flatten(self.legs)
        # legs as dictionary
        self.legdict = {l: v for v in range(len(self.genera))
                        for l in self.legs[v]}
        self.maxleg = max([max(j + [0]) for j in self.legs])
        # maxlevel is named a bit misleading, as it is the LOWEST level
        # (the max of the absolute values...)
        self.maxlevel = max([abs(l) for l in self.levels])
        # construct the "underlying graph" (as a sage graph)
        # the edges are labeled by their "real" name,
        # the vertices are of the form (i,g_i,'LG')
        # NOTE: In the case of an EmbeddedLevelGraph, vertices "at infinity" might
        # be added to this!
        self.underlying_graph = Graph([[self.UG_vertex(i) for i in range(len(self.genera))],
                                       [(self.UG_vertex(self.vertex(e[0])), self.UG_vertex(self.vertex(e[1])), e)
                                        for e in self.edges]],
                                      format='vertices_and_edges', loops=True, multiedges=True)

    def UG_vertex(self, v):
        """
        Vertex of the underlying graph corresponding to v.

        Args:
            v (int): vertex number (index of self.genera)

        Returns:
            tuple: tuple encoding the vertex of the Underlying Graph.
        """
        return (v, self.genera[v], 'LG')

    #################################################################
    # Extract a subgraph
    #################################################################
    def extract(self, vertices, edges):
        """
        Extract the subgraph of self (as a KLevelGraph) consisting of vertices
        (as list of indices) and edges (as list of edges).

        Returns the levelgraph consisting of vertices, edges and all legs on
        vertices (of self) with their original poleorders and the original
        level structure.
        """
        newvertices = deepcopy([self.genera[i] for i in vertices])
        newlegs = deepcopy([self.legs[i] for i in vertices])
        newedges = deepcopy(edges)
        newpoleorders = deepcopy({l: self.orderatleg(l)
                                  for l in flatten(newlegs)})
        newlevels = deepcopy([self.levels[i] for i in vertices])
        return KLevelGraph(newvertices, newlegs, newedges,
                           newpoleorders, newlevels, self.sig.k)

    def levelGraph_from_subgraph(self, G):
        """
        Returns the KLevelGraph associated to a subgraph of underlying_graph
        (with the level structure induced by self)
        """
        vertex_list = [v[0] for v in G.vertices(sort=True) if v[2] == 'LG']
        return self.extract(vertex_list, G.edge_labels())

    def UG_vertices_above_level(self, l):
        """
        The vertices above (internal) level l (including possible vertices at infinity).

        Used for checking residue conditions.

        Args:
            l (int): internal level number

        Returns:
            list: list of vertices of self.underlying_graph
        """
        # note that verticesabove checks for >l
        vertices_above = [self.UG_vertex(i) for i in self.verticesabove(l)]
        vertices_above += [v for v in self.underlying_graph.vertices(sort=True)
                           if v[2] == 'res']
        return vertices_above

    def UG_above_level(self, l):
        """
        Subgraph of Underlying Graph (including vertices at infinity) (strictly) above
        (relative!) level l.

        Args:
            l (int): relative level number

        Returns:
            SAGE Graph: subgraph of self.underlying_graph with all vertices strictly
                above level l.
        """
        internal_l = self.internal_level_number(l)
        vertices = self.UG_vertices_above_level(internal_l)
        abovegraph = self.underlying_graph.subgraph(vertices)
        return abovegraph

    def UG_without_infinity(self):
        """
        Subgraph of Underlying Graph without the vertices at infinity.

        Returns:
            SAGE graph
        """
        vertices = [self.UG_vertex(i) for i, _ in enumerate(self.genera)]
        return self.underlying_graph.subgraph(vertices)

    #################################################################
    # Squishing ---- i.e. contracting horizontal edges / levels
    #################################################################

    def _squish_horizontal(self, e):
        """
        Squish the horizontal edge e and returns the raw data for the new graph.
        """
        # make sure v <= w
        [v, w] = sorted([self.vertex(e[0]), self.vertex(e[1])])

        if not self.is_horizontal(e):
            print("Warning: " + repr(e) +
                  " is not a horizontal edge -- Not contracting.")
            return (self.genera, self.legs, self.edges,
                    self.poleorders, self.levels)

        newgenera = deepcopy(self.genera)
        newlegs = deepcopy(self.legs)
        newedges = deepcopy(self.edges)
        newpoleorders = deepcopy(self.poleorders)
        newlevels = deepcopy(self.levels)
        if v == w:
            # --horizontal loop getting contracted---
            # add genus to node:
            newgenera[v] += 1
        else:
            # --horizontal edge getting contracted---
            # transfer genus and legs from w to v:
            newgenera[v] += newgenera[w]
            newlegs[v] += newlegs[w]
            # remove all traces of w
            newgenera.pop(w)
            newlegs.pop(w)
            newlevels.pop(w)
        # remove edge:
        newedges.remove(e)
        # remove legs
        newlegs[v].remove(e[0])
        newlegs[v].remove(e[1])
        # remove poleorders
        del newpoleorders[e[1]]
        del newpoleorders[e[0]]

        return [newgenera, newlegs, newedges, newpoleorders, newlevels]

    def squish_horizontal(self, e):
        """
        Squish the horizontal edge e and return the new graph.

        EXAMPLES::

            sage: from admcycles.diffstrata.klevelgraph import KLevelGraph
            sage: G = KLevelGraph.fromPOlist([1,1], [[1,2],[3,4]], [(2,3)], [1,-1,-1,1],[0,0],1)
            sage: G.squish_horizontal((2,3))
            KLevelGraph([2],[[1, 4]],[],{1: 1, 4: 1},[0],1,True)
        """
        # Python2 only allows keyword arguments after *() arguments, so we need
        # to repack
        dat = self._squish_horizontal(e) + [self.sig.k]
        return KLevelGraph(*dat, quiet=True)

    def squish_vertical_slow(self, l, adddata=False, quiet=False):
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
        level l-1. We then create a new KLevelGraph with the same data as self,
        the only difference being that all vertices of level l-1 have now been
        assigned level l. We then remove all (now horizontal!) edges in ee.

        In particular, all points and edges retain their numbering (only the level might have changed).

        WARNING: Level l-1 does not exist anymore afterwards!!!

        Downside: At the moment we get a warning for each edge in ee, as these
        don't have legal pole orders for being horizontal (so checkedgeorders
        complains each time the constructor is called :/).
        """
        if l is None:
            return self
        ll = self.nextlowerlevel(l)
        if ll is False:
            if not quiet:
                print("Warning: Illegal level to contract: " + repr(l))
            if adddata:
                return (self, False, None)
            else:
                return self

        if not quiet:
            print("Squishing levels", l, "and", ll)

        vv = self.verticesonlevel(ll)
        # edges that go to next lower level, i.e. will be contracted
        ee = [e for e in self.edgesatlevel(l)
              if self.levelofleg(e[1]) == ll and not self.is_horizontal(e)]

        if not quiet:
            print("Contracting edges", ee)

        newgenera = deepcopy(self.genera)
        newlegs = deepcopy(self.legs)
        newedges = deepcopy(self.edges)
        newpoleorders = deepcopy(self.poleorders)
        newlevels = deepcopy(self.levels)

        # raise level
        for v in vv:
            newlevels[v] = l
        returngraph = KLevelGraph(
            newgenera,
            newlegs,
            newedges,
            newpoleorders,
            newlevels,
            self.sig.k,
            True)
        # contract edges that are now horizontal:
        for e in ee:
            returngraph = returngraph.squish_horizontal(e)

        if adddata:
            return (returngraph, True, flatten(ee))
        else:
            return returngraph

    def _squish_vertical(self, l, quiet=True):
        """
        Squish the level l (and the next lower level!) and return the raw data for the new graph.

        Implementation:
        At the moment, we remember all the edges (ee) that connect level l to
        level l-1. We then create a new KLevelGraph with the same data as self,
        the only difference being that all vertices of level l-1 have now been
        assigned level l. We then remove all edges in ee.

        (In contrast to the old squish_vertical_slow, this is now done in
        one step and not recursively so we avoid a lot of the copying of graphs.)

        In particular, all points and edges retain their numbering (only the level might have changed).

        WARNING: Level l-1 does not exist anymore afterwards!!!

        Args:
            l (int): (internal) level number
            quiet (bool, optional): No output. Defaults to True.

        Returns:
            tuple: raw data of the new graph
        """
        if l is None:
            return self
        ll = self.nextlowerlevel(l)
        levelset = set([l, ll])
        if ll is False:
            if not quiet:
                print("Warning: Illegal level to contract: " + repr(l))
            else:
                return (self.genera, self.legs, self.edges,
                        self.poleorders, self.levels)
        vertices = self.verticesonlevel(l) + self.verticesonlevel(ll)
        # edges that go to next lower level, i.e. will be contracted
        edges = [e for e in self.edges if set(
            [self.levelofleg(e[0]), self.levelofleg(e[1])]) == levelset]
        # we use two dictionaries for quick lookup while we reorder the info:
        genus_of_vertex = {v: self.genus(v) for v in vertices}
        vertex_of_leg = {
            leg: v for v in vertices for leg in self.legsatvertex(v)}
        legs_of_vertex = {v: self.legsatvertex(v)[:] for v in vertices}
        deleted_legs = []
        while edges:
            start, end = edges.pop()
            v = vertex_of_leg[start]
            w = vertex_of_leg[end]
            del vertex_of_leg[start]
            del vertex_of_leg[end]
            legs_of_vertex[v].remove(start)
            legs_of_vertex[w].remove(end)
            deleted_legs.extend([start, end])
            if v == w:
                # if e (now) connects a vertex to itself, increase its genus
                genus_of_vertex[v] += 1
            else:
                # if e connects two different vertices, combine their data
                # (moving everything from w to v)
                genus_of_vertex[v] += genus_of_vertex[w]
                del genus_of_vertex[w]
                for leg in legs_of_vertex[w]:
                    vertex_of_leg[leg] = v
                    legs_of_vertex[v].append(leg)
                del legs_of_vertex[w]
        # Build new graph with this data:
        # We use sets for more efficient lookup:
        vertices = set(vertices)
        deleted_legs = set(deleted_legs)

        newgenera = []
        newlegs = []
        newlevels = []
        # we copy the untouched vertices and insert the new ones:
        for v in range(len(self.genera)):
            if v in vertices:
                if v in genus_of_vertex:
                    newgenera.append(genus_of_vertex[v])
                    newlegs.append(legs_of_vertex[v][:])
                    newlevels.append(l)
                # otherwise it has been deleted
            else:
                newgenera.append(self.genus(v))
                newlegs.append(self.legsatvertex(v)[:])
                newlevels.append(self.levelofvertex(v))

        # for the edges, we remove those with a deleted edge:
        newedges = []
        for start, end in self.edges:
            if start in deleted_legs:
                assert end in deleted_legs
                continue
            assert end not in deleted_legs
            newedges.append((start, end))
        # for the new poleorders, we only have to remove the deleted legs:
        newpoleorders = {
            leg: p for leg,
            p in self.poleorders.items() if leg not in deleted_legs}

        return [newgenera, newlegs, newedges, newpoleorders, newlevels]

    def squish_vertical(self, l, quiet=True):
        """
        Squish the level l (and the next lower level!) and return the new graph.

        In particular, all points and edges retain their numbering (only the level might have changed).

        WARNING: Level l-1 does not exist anymore afterwards!!!

        Args:
            l (int): (internal) level number
            quiet (bool, optional): No output. Defaults to True.

        Returns:
            KLevelGraph: new levelgraph with one level less.

        EXAMPLES::

            sage: from admcycles.diffstrata.klevelgraph import KLevelGraph
            sage: G = KLevelGraph.fromPOlist([1,2], [[1],[2,3,4]], [(1,2)], [0,-2,1,3],[0,-1],1)
            sage: G.squish_vertical(0)
            KLevelGraph([3],[[3, 4]],[],{3: 1, 4: 3},[0],1,True)
        """
        # Python2 only allows keyword arguments after *() arguments, so we need
        # to repack
        dat = self._squish_vertical(l, quiet=quiet) + [self.sig.k]
        return KLevelGraph(*dat, quiet=True)

    def delta(self, k, quiet=False):
        """
        Squish all levels except for the k-th.

        Note that delta(1) contracts everything except top-level and that
        the argument is interpreted via internal_level_number

        WARNING: Currently giving an out of range level (e.g.
            0 or >= maxlevel) squishes the entire graph!!!

        Return the corresponding divisor.
        """
        G = self
        if self.is_horizontal():
            print("Error: Cannot delta a graph with horizontal edges!")
            return deepcopy(G)
        # recursive squishing forces us to go backwards, as the lower squished
        # level is always removed.
        for i in reversed(range(self.numberoflevels() - 1)):
            # unfortunately delta and squish use slightly incompatible level
            # numbering, as the "illegal" squish here is k-1
            if abs(k) - 1 == i:
                continue
            if not quiet:
                print("Squishing level", i)
            # note that it is important that internal_level_number of
            # self and not G is called in the argument, as the relative
            # numbering of levels in G changes, but the internal names don't
            G = G.squish_vertical(self.internal_level_number(i), quiet=quiet)
        return G

    #################################################################
    # end Squishing
    #################################################################

    #################################################################
    # Isomorphisms
    # DEPRECATED!!!! USE EmbeddedLevelGraphs instead!
    ####
    # This only remains for the BIC generation comparison tests.
    #################################################################
    def _init_invariant(self):
        """
        DEPRECATED!! DO NOT USE!!

        Compute a bunch of invariants (as a quick check for not isommorphic).
        Up to now we have:
        * sorted list of genera (by level)
        * number of edges
        * number of levels
        * names and pole orders of marked points
        * sorted list of pole orders
        """
        self._invariant = (
            len(self.edges),
            self.codim(),
            [(l, self.orderatleg(l))
             for l in sorted(self.list_markings())],
            sorted([self.orderatleg(l) for l in self.leglist]))

    def invariant(self):
        """
        DEPRECATED!!! DO NOT USE!!!
        """
        self._init_invariant()
        return self._invariant

    def _genisom_preserves_levels(self, other, dV):
        """
        Check if a "vertex isomorphism" preserves the (relative) level structure
        """
        return all(self.sortedlevels.index(self.levelofvertex(sv)) ==
                   other.sortedlevels.index(other.levelofvertex(ov))
                   for sv, ov in dV.items())

    def _legisom_preserves_poleorders(self, other, dL):
        """
        Check if a "leg isomorphism" preserves pole orders.
        """
        return all(self.orderatleg(sl) == other.orderatleg(ol)
                   for sl, ol in dL.items())

    def _leveliso(self, other, dV):
        """
        Give a dictionary translating the levels of self to other in accordance
        with the dictionary of vertices dV.
        """
        return {l: other.levelofvertex(dV[self.verticesonlevel(l)[0]])
                for l in self.levels}

    def is_isomorphic(self, other, adddata=False):
        """
        DEPRECATED!!! Instead, use EmbeddedLevelGraph!!!

        Check if self is isomorphic to other.
        Return boolean + list of isomorphisms if adddata=True.

        Note that our isomorphisms are stgraph isomorphisms plus a dictionary
        translating the level structure of self into that of other.

        At the moment we use admcycles.GraphIsom and check which of these
        preserve the KLevelGraph structure (using the above invariants, though!)
        Probably it would be much faster to reimplement this for KLevelGraphs
        as it seems quite redundant as is...
        """
        isoms = GraphIsom(self.stgraph, other.stgraph)
        # check if the stable graph isomorphisms preserve level structure
        # GraphIsoms consist of a tuple (dV,dL), where dV is a dictionary
        # translating the vertices and dL is one for the legs.
        levelisoms = [[dV, dL, self._leveliso(other, dV)] for [dV, dL] in isoms
                      if self._genisom_preserves_levels(other, dV) and
                      self._legisom_preserves_poleorders(other, dL)]
        if adddata:
            return (bool(levelisoms), levelisoms)
        else:
            return (bool(levelisoms))

    #################################################################
    # Plotting
    #################################################################
    def plot_obj(self):
        """
        Return a sage graphics object generated from a sage graph.

        TODO: Some things that are still ugly:
        * No legs!!
        * currently the vertices have labels (index,genus)
        * loops are vertical not horizontal
        """
        # Create the underlying graph with edge label prong order.
        G = Graph([self.genera,
                   [(self.vertex(e[0]), self.vertex(e[1]), self.prong(e))
                    for e in self.edges]],
                  format='vertices_and_edges', loops=True, multiedges=True)
        # unfortunately, sage insists of displaying the vertex name (index)
        # and not only the label (genus), so we display both...
        G.relabel(list(enumerate(self.genera)))
        h = {l: [(v, self.genus(v)) for v in self.verticesonlevel(l)]
             for l in self.levels}
        # at the moment we just display a list of legs and orders at the bottom
        legtext = "Marked points: "
        for l in sorted(set(self.levels)):
            points_at_level = [leg for leg in self.list_markings()
                               if self.levelofleg(leg) == l]
            if not points_at_level:
                continue
            legtext += "\nAt level %s: " % l
            for leg in points_at_level:
                legtext += "\n" + self._pointtype(leg)
                legtext += " of order " + repr(self.orderatleg(leg))
        print(legtext)
        return G.plot(edge_labels=True, vertex_size=1000,
                      heights=h, layout='ranked') + text(legtext, (0, -0.1),
                                                         vertical_alignment='top', axis_coords=True, axes=False)
