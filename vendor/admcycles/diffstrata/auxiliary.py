from collections import defaultdict

# pylint does not know sage
from sage.misc.flatten import flatten  # pylint: disable=import-error

import admcycles.diffstrata.levelgraph
import admcycles.diffstrata.embeddedlevelgraph

#################################################################
#################################################################
#   Auxiliary functions:
#################################################################
#################################################################


def unite_embedded_graphs(gen_LGs):
    return admcycles.diffstrata.embeddedlevelgraph.EmbeddedLevelGraph(
        *unite_embedded_k_graphs(
            gen_LGs,
            admcycles.diffstrata.levelgraph.LevelGraph))


def unite_embedded_k_graphs(gen_LGs, clsLG):
    """
    Create a (disconnected) EmbeddedLevelGraph from a tuple of tuples that generate EmbeddedLevelGraphs.

    (The name is slightly misleading, but of course it does not make sense to actually unite two complete
    EmbeddedLevelGraphs, as the checks would (and do!) throw errors otherwise! Therefore, this essentially
    takes the data of a LevelGraph embedded into each connected componenent of a GeneralisedStratum and
    returns an EmbeddedLevelGraph on the product.)

    This should be used on (products) of BICs in generalised strata.

    INPUT:

    gen_LGs: tuple
    A tuple of tuples that generate EmbeddedLevelGraphs.
    More precisely, each tuple is of the form:

    * X (GeneralisedStratum): Enveloping stratum (should be the same for all tuples!)
    * LG (LevelGraph): Underlying LevelGraph
    * dmp (dict): (partial) dictionary of marked points
    * dlevels (dict): (partial) dictionary of levels

    clsELG: class of the graph to return, either EmbeddedLevelGraph or EmbeddedQuadraticLevelGraph

    clsLG: class of the underlying levelgraph, either LevelGraph or QuadraticLevelGraph

    OUTPUT:

    The (disconnected) LevelGraph obtained from the input with the legs
    renumbered (continuously, starting with 1), and the levels numbered
    according to the embedding.
    """
    newgenera = []
    newlevels = []
    newlegs = []
    newpoleorders = {}
    newedges = []
    newdmp = {}
    newdlevels = {}
    max_leg_number = 0
    oldX = gen_LGs[0][0]  # for check that all belong to the same stratum:
    for emb_g in gen_LGs:
        # Unpack tuple:
        X, LG, dmp, dlevels = emb_g
        if X != oldX:
            raise RuntimeError(
                "Can't unite graphs on different Strata! %r" % gen_LGs)
        # the genera are just appended
        newgenera += LG.genera
        # same for the levels, but while we're at it, we might just as well
        # replace them by their embedding (then newdlevels will be trivial)
        # and these should be consistens for all graphs in the tuple.
        # Thus, newdlevels will be the identity.
        newlevels += [dlevels[l] for l in LG.levels]
        # the legs will have to be renumbered
        leg_dict = {}  # old number -> new number
        legs = 0
        for i, l in enumerate(flatten(LG.legs)):
            newlegnumber = max_leg_number + i + 1
            leg_dict[l] = newlegnumber
            # while we're at it, we add the pole orders:
            newpoleorders[newlegnumber] = LG.poleorders[l]
            # For the dictionary of marked points (for the embedding), we
            # must distinguish if this is a marked point or a half-edge.
            # Marked points are simply the ones for which we have a key
            # in dmp :-)
            try:
                newdmp[newlegnumber] = dmp[l]
            except KeyError:
                pass
            legs += 1
        max_leg_number += legs
        # append (nested) list of legs:
        newlegs += [[leg_dict[l] for l in comp] for comp in LG.legs]
        # finally, the edges are renumbered accordingly:
        newedges += [(leg_dict[e[0]], leg_dict[e[1]]) for e in LG.edges]
    # the levels are already numbered according to the embedding dict
    newdlevels = {l: l for l in newlevels}
    newLG = clsLG(
        newgenera, newlegs, newedges, newpoleorders, newlevels
    )
    return (X, newLG, newdmp, newdlevels)


def sort_with_dict(l):
    """
    Sort a list and provide a dictionary relating old and new indices.

    If x had index i in l, then x has index sorted_dict[i] in the sorted l.

    Args:
        l (list): List to be sorted.

    Returns:
        tuple: A tuple consisting of:
            list: The sorted list l.
            dict: A dictionary old index -> new index.
    """
    sorted_list = []
    sorted_dict = {}
    for i, (j, v) in enumerate(sorted(enumerate(l), key=lambda w: w[1])):
        sorted_list.append(v)
        sorted_dict[j] = i
    return sorted_list, sorted_dict


def get_squished_level(deg_ep, ep):
    """
    Get the (relative) level number of the level squished in ep.

    This is the index of the corresponding BIC in the profile.

    Args:
        deg_ep (tuple): enhanced profile
        ep (tuple): enhanced profile

    Raises:
        RuntimeError: raised if deg_ep is not a degeneration of ep

    Returns:
        int: relative level number
    """
    deg_p = deg_ep[0]
    p = set(ep[0])
    for i, b in enumerate(deg_p):
        if b not in p:
            break
    else:
        raise RuntimeError("%r is not a degeneration of %r!" % (deg_ep, p))
    return i

#################################################################
#################################################################
#    Auxiliary functions for caching:
#################################################################
#################################################################


def hash_AG(leg_dict, enh_profile):
    """
    The hash of an AdditiveGenerator, built from the psis and the enhanced profile.

    The hash-tuple is (leg-tuple,profile,index), where profile is
    changed to a tuple and leg-tuple is a nested tuple consisting of
    tuples (leg,exponent) (or None).

    Args:
        leg_dict (dict): dictioary for psi powers (leg -> exponent)
        enh_profile (tuple): enhanced profile

    Returns:
        tuple: nested tuple
    """
    if leg_dict is None:
        leg_hash = ()
    else:
        leg_hash = tuple(sorted(leg_dict.items()))
    return (leg_hash, tuple(enh_profile[0]), enh_profile[1])


def adm_key(sig, psis):
    """
    The hash of a psi monomial on a connected stratum without residue conditions.

    This is used for caching the values computed using admcycles (using
    GeneralisedStratum.adm_evaluate)

    The signature is sorted, the psis are renumbered accordingly and also
    sorted (with the aim of computing as few duplicates as possible).

    Args:
        sig (tuple): signature tuple
        psis (dict): psi dictionary

    Returns:
        tuple: nested tuple
    """
    sorted_psis = {}
    sorted_sig = []
    psi_by_order = defaultdict(list)
    # sort signature and relabel psis accordingly:
    # NOTE: Psis are labelled 'mathematically', i.e. 1,...,len(sig)
    for new_i, (old_i, order) in enumerate(
            sorted(enumerate(sig), key=lambda k: k[1])):
        psi_new_i = new_i + 1
        psi_old_i = old_i + 1
        sorted_sig.append(order)
        if psi_old_i in psis:
            assert not (psi_new_i in sorted_psis)
            psi_exp = psis[psi_old_i]
            sorted_psis[psi_new_i] = psi_exp
            psi_by_order[order].append(psi_exp)
    # sort psis for points of same order:
    ordered_sorted_psis = {}
    i = 0
    assert len(sig) == len(sorted_sig)
    while i < len(sig):
        order = sorted_sig[i]
        for j, psi_exp in enumerate(sorted(psi_by_order[order])):
            assert sorted_sig[i + j] == order
            ordered_sorted_psis[i + j + 1] = psi_exp
        while i < len(sig) and sorted_sig[i] == order:
            i += 1
    return (tuple(sorted_sig), tuple(sorted(ordered_sorted_psis.items())))
