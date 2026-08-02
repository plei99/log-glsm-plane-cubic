r"""
Spin Stratum: A stratum with signature of even type will has spin structure, i.e.
its connected components can be partitioned to two non-empty sets.

cf. [KZ03]_ and [Boi15]_

"""

from copy import deepcopy
import sage.misc.persist

from sage.rings.all import QQ
from sage.matrix.constructor import matrix, vector
from sage.misc.cachefunc import cached_function
from sage.combinat.all import SetPartitions

import admcycles.admcycles
import admcycles.stratarecursion

import admcycles.diffstrata.levelgraph
import admcycles.diffstrata.bic
import admcycles.diffstrata.stratatautring
import admcycles.diffstrata.embeddedlevelgraph

from admcycles.stable_graph import StableGraph
from admcycles.diffstrata.generalisedstratum import GeneralisedStratum
from admcycles.diffstrata.additivegenerator import AdditiveGenerator
from admcycles.diffstrata.elgtautclass import ELGTautClass
from admcycles.diffstrata.sig import Signature

from admcycles.identify_classes import SolveMultLin


def Spin_strataclass(sig, res_cond=[]):
    r"""
    Compute the spin stratum class of given signature of even type (may
    have paired simple poles) and residue conditions.

    INPUT:

    - sig (tuple)
    - res_cond (list): It is a list of lists corresponding to residue
      conditions. For example, sig=(4,-2,-1,-1) with the
      simple poles being paired. Then
      res_cond=[[(0,2), (0,3)]]. Here the zero indicate
      the components of a generalized stratum.

    EXAMPLES::

        sage: from admcycles.diffstrata.spinstratum import Spin_strataclass
        sage: Spin_strataclass((2,2)).basis_vector()
        (159/4, -179/48, -7/24, -7/24, -131/48, -23/24, -131/48, 149/24, -83/16, 271/48, 77/48, 77/48, -193/8, -185/48, 127/48, 395/48, 221/48, -73/8, -185/48, 127/48, 395/48, 221/48, -73/8, -185/48, -41/48, 389/48, 389/48, -23/8, 51/8, -139/16, -323/16, 1/8, -25/4, -11/4, -11/4, -23/4, 973/48, 37/96, -5/2, 23/96, 23/96, -25/32)


    """
    X = SpinStratum([Signature(sig)], res_cond=res_cond)
    A = AG_with_spin(X, ((), 0))
    E = ELGT_with_spin(X, [(1, A)])
    return E.to_prodtautclass_spin().pushforward()


#######################################################################
#######################################################################
# Define the class of spin strata #####################################
#######################################################################

class SpinStratum(GeneralisedStratum):

    def __init__(self, sig_list, res_cond=None, pair_res_cond=None):
        r"""
        A stratum can have spin structure if there are pairs of simple poles,
        such that their residues r_i + r_i+1 = 0. In this code, we only handle the
        cases of at most 1 pair of simple poles.

        INPUT:

        - sig_list (list): a list of tuples indicating the signatures
          of the components of a generalised stratum

        - res_cond (list): a list whose entries are lists corresponding
          to residue conditions

        EXAMPLES::

            sage: from admcycles.diffstrata.spinstratum import SpinStratum
            sage: from admcycles.diffstrata.sig import Signature
            sage: X = SpinStratum([Signature((4,-2,-1,-1))],
            ....: res_cond=[[(0,2),(0,3)]])

        """

        super().__init__(sig_list, res_cond)

        if res_cond:
            self._res_cond = res_cond
        else:
            self._res_cond = []
        if pair_res_cond:
            self.pair_res_cond = pair_res_cond
        else:
            self.pair_res_cond = []
        self.odd_sing_list = []
        for i in range(0, len(self._sig_list)):
            for j in range(0, len(sig_list[i].sig)):
                if sig_list[i].sig[j] % 2 != 0:
                    self.odd_sing_list.append((i, j))

        self.sim_pole_list = list(self.simple_poles())
        self.sim_pole_list.sort()
        self.hor_edge = []
        self.eff_hor_edge = []
        self.eff_pair_res = []

        if self.odd_sing_list != [] and pair_res_cond is None:
            self.is_legal()

    def __repr__(self):
        return "SpinStratum(sig_list=%r,res_cond=%r)" % (
            self._sig_list, self._res_cond)

    def __str__(self):
        rep = ''
        if self._h0 > 1:
            rep += 'Product of Strata with spin:\n'
        else:
            rep += 'Spin Stratum: '
        for sig in self._sig_list:
            rep += repr(sig.sig) + '\n'
        rep += 'with residue conditions: '
        if not self._res_cond:
            rep += repr([]) + '\n'
        for res in self._res_cond:
            rep += repr(res) + '\n'
        return rep

    def is_legal(self):
        r"""
        This is just a method for initializing a spin stratum, checking whether
        the stratum really has spin structure.

        EXAMPLES::

            sage: from admcycles.diffstrata.spinstratum import SpinStratum
            sage: from admcycles.diffstrata.sig import Signature
            sage: X = SpinStratum([Signature((4,-2,-2))])
            sage: X = SpinStratum([Signature((4,-3,-1))])
            Traceback (most recent call last):
            ...
            ValueError: The Stratum does not have spin components
            sage: X = SpinStratum([Signature((4,-2,-1,-1))])
            sage: X = SpinStratum([Signature((4,-1,-1,-1,-1))])
            Traceback (most recent call last):
            ...
            NotImplementedError: We can only handle one pair of simple poles

        """

        if len(self.odd_sing_list) - len(self.sim_pole_list) > 0:
            raise ValueError("The Stratum does not have spin components")
        elif len(self.sim_pole_list) > 0:
            t = len(self.sim_pole_list)
            if t > 2:
                raise NotImplementedError("We can only handle " +
                                          "one pair of simple poles")
            if not t % 2 == 0:
                raise ValueError("The Stratum does not have spin components")

            M1 = self.smooth_LG.residue_matrix_from_RT
            if self.res_cond != []:
                M2 = M1.stack(self.matrix_from_res_conditions(self.res_cond))
            else:
                M2 = M1
            a = 0
            b = M2.rank()
            c = M1.rank()

            for q in range(0, t - 1):
                for q1 in range(q + 1, t):
                    if self.sim_pole_list[q][0] != self.sim_pole_list[q1][0]:
                        res = [self.sim_pole_list[q], self.sim_pole_list[q1]]
                        T = self.matrix_from_res_conditions(res)
                        M0 = M2.stack(T)
                        a = M0.rank()
                        if a > b:
                            continue
                        else:
                            self.hor_edge.append([self.sim_pole_list[q],
                                                  self.sim_pole_list[q1]])
            sim_pole_list_partitons = SetPartitions(t, [2 for i in
                                                    range(0, int(t / 2))]
                                                    ).list()
            rc = []
            for part in sim_pole_list_partitons:
                pair_res_cond_temp = []
                for q in part:
                    w = list(q)
                    pair_res_cond_temp.append([self.sim_pole_list[i - 1]
                                               for i in w])

                M0 = M2.stack(self.matrix_from_res_conditions(
                              pair_res_cond_temp))
                a = M0.rank()
                if a > b:
                    continue
                else:
                    self.pair_res_cond = pair_res_cond_temp
                    rc = pair_res_cond_temp + deepcopy(self._res_cond)
                    break

            if self.pair_res_cond is None:
                raise ValueError("The Stratum does not have spin components")
            if t == 1 and a == 1:
                self._res_cond = []
            elif not self.is_empty():
                new_rc = []
                for j, q in enumerate(rc):
                    rk = M1.stack(self.matrix_from_res_conditions([q])).rank()
                    if rk > c:
                        M1 = M1.stack(self.matrix_from_res_conditions([q]))
                        c += 1
                        if q in self.pair_res_cond:
                            self.eff_pair_res.append(j)
                        new_rc.append(q)
                new_rc.reverse()
                self._res_cond = new_rc
        return

    def remove_pair_res(self, psis=None):
        r"""
        This method will return the psi polynomial of the stratum in the
        tautological ring of ambient stratum (without a residue condition
        on the pair of simple poles. The formula follows Sauvaget-19 and we pick
        the leg of a simple pole for the elimination of psi class.

        Only the BICs, which have some odd enhancement (so odd even are half
        half), or still preserve the residue condition of the pair of
        simple poles, are left. Hence, one can
        still recursively compute the spin classes.

        INPUT:

        - psis (list): a list of dicts indicating the psi polynomial

        EXAMPLES::


            sage: from admcycles.diffstrata.spinstratum import SpinStratum
            sage: from admcycles.diffstrata.sig import Signature
            sage: from admcycles.admcycles import psiclass
            sage: X = SpinStratum([Signature((4,-2,-1,-1))],res_cond=[[(0,2),(0,3)]])
            sage: cl1 = X.remove_pair_res().to_prodtautclass_spin().pushforward()
            sage: (psiclass(1,1,4)**2*cl1).evaluate()
            5/24

        """
        if not self.eff_pair_res == []:
            # print("removing pair res")
            res_cond = deepcopy(self._res_cond)
            lost_rc = res_cond.pop()
            leg = lost_rc[0]
            new_stratum = SpinStratum(self._sig_list, res_cond,
                                      pair_res_cond=self.pair_res_cond)
            new_AG = new_stratum.additive_generator(((), 0), psis)
            elgt = new_AG.as_taut() * (new_stratum.res_stratum_class(lost_rc) +
                                       new_stratum.xi_with_leg() -
                                       new_stratum.xi_with_leg(leg))
            temp_psi_list = []
            for a in elgt._psi_list:
                if all(kappa % 2 != 0
                       for e, kappa in a[1]._G.LG.prongs_list()):
                    temp_psi_list.append(a)
            psi_list = [(c, AG_with_spin(new_stratum,
                        AG._enh_profile, AG._leg_dict))
                        for c, AG in temp_psi_list]
            elgt_odd_enh = ELGT_with_spin(new_stratum, psi_list)
            # print(str(elgt_odd_enh))
            result = elgt_odd_enh
            return result
        else:
            raise NotImplementedError

    def smooth_con_stratum_prodtautclass(self):
        r"""
        To compute the spin stratum class of a connected stratum. The output will be a
        product tautological class. By making use of standised signature, we will do
        the core computation in the function smooth_stratum_standard().

        EXAMPLES::

            sage: from admcycles.diffstrata.spinstratum import SpinStratum
            sage: from admcycles.diffstrata.sig import Signature
            sage: X=SpinStratum([Signature((4,0,-2))])
            sage: X.smooth_con_stratum_prodtautclass().pushforward().basis_vector()
            (-63/2, 21/2, -20, -63/2, -37/2, 38, 61/2, 61/2, 147/2, 29, 57/2, 29, -41, -60, 21, -63/2, -42, 61/2, -41, -139/2, 147/2, 63/2, -147/2, 21/2, 0, -21/2, -7, -21/2, -19/2, 19/2, 21/2, -21/2, 1, -21/2, 0, 0, 21/2, 21/2, -7, 1, -1, -1, 0, 7)

        """
        res_cond = self._res_cond
        sig = self._sig_list[0].sig
        g = self._sig_list[0].g
        LG = self.smooth_LG.LG
        stgraph0 = LG.stgraph

        if all(x == 0 for x in sig):  # it is just g=1 and the stratum class is just fundclass of M_1,n
            t = deepcopy(smooth_stratum_standard(sig))
            if all(w == 0 for w in t[3]):
                result = admcycles.admcycles.prodtautclass(stgraph0, terms=[])
            else:
                ind = admcycles.admcycles.generating_indices(t[0], t[1], t[2])
                gen = admcycles.admcycles.tautgens(t[0], t[1], t[2])
                result = sum([t[3][i] * gen[ind[i]]
                              for i in range(0, len(t[3]))]).toprodtautclass()
            return result

        L = standardise_sig(sig, res_cond)  # We use the standardised signature
        simplified = L[0]
        num_of_reg = L[1]
        legs = {x2: x1 for x1, x2 in L[2].items()}
        new_rc = L[3]

        # The stratum classes of signature with 0s are essentially just the
        # forgetful pullback of stratum class of signature without 0s.
        # We use the function smooth_stratum_standard to do the pullbacks solving.

        t = deepcopy(smooth_stratum_standard(simplified, tuple(new_rc)))
        if all(w == 0 for w in t[3]):
            cl = admcycles.admcycles.fundclass(g, len(simplified))
        else:
            ind = admcycles.admcycles.generating_indices(t[0], t[1], t[2])
            gen = admcycles.admcycles.tautgens(t[0], t[1], t[2])
            cl = sum([t[3][i] * gen[ind[i]] for i in range(0, len(t[3]))])
        if num_of_reg != 0:
            cl = cl.forgetful_pullback([len(simplified) + i
                                        for i in range(1, L[1] + 1)])
        result = cl.rename_legs(legs, inplace=False).toprodtautclass()
        # print(str(legs)+ "\n"+ str(cl) + "\n" + str(result))
        return result

#######################################################################
#######################################################################
# Define the spin version additive generators #########################
#########################################
# * the spin taut class of those enhanced level graphs with some even
#   enhancement are just considered as zero in this algorithm, because
#   the pushforward of them to \bar M_g,n will be zero.
#######################################################################


class AG_with_spin (AdditiveGenerator):

    def __init__(self, X, enh_profile, leg_dict=None):
        assert isinstance(X, SpinStratum)
        super().__init__(X, enh_profile, leg_dict)
        self.X = X

    def __repr__(self):
        return "AG_with_spin(X=%r,enh_profile=%r,leg_dict=%r)"\
            % (self._X, self._enh_profile, self._leg_dict)

    def __str__(self):
        str = ""
        if self._leg_dict is not None:
            for l in self._leg_dict:
                str += "Spin Psi class %r with exponent %r on level %r * "\
                    % (l, self._leg_dict[l], self._level_dict[l])
        str += "Graph %r" % (self._enh_profile,)
        return str

    def to_prodtautclass_spin(self, rearrange_markings=True):
        r"""
        This function is to compute the product tautological class of an
        additive generator with spin. It is basically the same as
        toprodtautclass, but just includes more case divisions.

        EXAMPLES::

            sage: from admcycles.diffstrata.spinstratum import SpinStratum,AG_with_spin
            sage: from admcycles.diffstrata.sig import Signature
            sage: X = SpinStratum([Signature((4,-2))])
            sage: AG_with_spin(X, ((1,),0)).to_prodtautclass_spin()
            Outer graph : [1, 0] [[3, 4], [1, 2, 5, 6]] [(3, 5), (4, 6)]
            Vertex 0 :
            Graph :      [1] [[1, 2]] []
            Polynomial : -1/2
            Vertex 1 :
            Graph :      [0] [[1, 2, 3, 4]] []
            Polynomial : psi_4
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1] [[1, 2]] []
            Polynomial : -1/2
            Vertex 1 :
            Graph :      [0, 0] [[2, 3, 11], [1, 4, 12]] [(11, 12)]
            Polynomial : 3
            <BLANKLINE>
            <BLANKLINE>
            Vertex 0 :
            Graph :      [1] [[1, 2]] []
            Polynomial : -1/2
            Vertex 1 :
            Graph :      [0, 0] [[3, 4, 11], [1, 2, 12]] [(11, 12)]
            Polynomial : -3

        """
        LG = self._G.LG
        stgraph = LG.stgraph

        # If the level graph has any edge with positive even enhancement,
        # half of the prong matchings give an even spin curve and half of
        # them give an odd spin curve
        if any(kappa % 2 == 0 for e, kappa in LG.prongs_list()):
            return admcycles.admcycles.prodtautclass(stgraph, terms=[])

        # Check whether the prodtautclass will be zero
        # If any level has extra freedom of scaling, the class will be zero on
        # moduli space of curves.
        if any(self.level(l).zeroStratumClass()
               for l in range(self.codim + 1)):
            if rearrange_markings:
                stgraph = self._G.relabel(self._G.standard_markings(), tidyup=False).LG.stgraph
            return admcycles.admcycles.prodtautclass(stgraph, terms=[])

        # To handle the case that the stratum is connected
        if len(LG.genera) == 1:
            adm_psis = admcycles.admcycles.decstratum(stgraph,
                                                      psi=self.leg_dict)
            adm_psis_taut = admcycles.admcycles.tautclass([adm_psis])

            # To handle the case that the genus is 0 and there are no
            # simple poles.
            if LG.genera[0] == 0 and self.X.sim_pole_list == []:
                stratum_class = AdditiveGenerator(self.X, ((), 0)).as_taut()
                stratum_class = stratum_class.to_prodtautclass().pushforward()
                protaut = [adm_psis_taut * stratum_class]
                result = admcycles.admcycles.prodtautclass(stgraph,
                                                           protaut=protaut)
                return result

            # To handle the case that there are only two poles which are a pair
            # of simple poles.
            if (self.X.res_cond == [] or
               len(self.X.eff_pair_res) < len(self.X.pair_res_cond)):
                stratum_class = self.X.smooth_con_stratum_prodtautclass()
                stratum_class = stratum_class.pushforward()
                protaut = [stratum_class * adm_psis_taut]
                result = admcycles.admcycles.prodtautclass(stgraph,
                                                           protaut=protaut)
                return result

            # To resolve residue condition of a pair of simple poles if the
            # residue condition is extra.
            elif (self.X.eff_pair_res != [] and
                  len(self.X.eff_pair_res) == len(self.X.pair_res_cond)):
                Q = self.X.remove_pair_res()
                total = Q.to_prodtautclass_spin().pushforward()
                protaut = [adm_psis_taut * total]
                result = admcycles.admcycles.prodtautclass(stgraph,
                                                           protaut=protaut)
                return result

        # Now it processes to break down the level graph into level
        # stratum and resolve the residue conditions (non simple poles).
        # Then it computes the result by level clutching.
        # The level graph that can come to this stage will have only odd
        # enhancements or the stratum is not connected, i.e. no vertical edges.
        # In the both cases, the spin class depends on the spin class on
        # each level.

        alpha = []
        vertices = []
        for l in range(self.codim + 1):
            psis = self.psis_on_level(l)
            level = self.level(l)
            level = addspin(level)  # level stratum with spin structure

            if (level.eff_pair_res != [] and
               len(level.eff_pair_res) == len(level.pair_res_cond)):
                # To remove residue condition on a pair of simple poles on
                # a disconnected stratum.
                ptc = level.remove_pair_res(psis).to_prodtautclass_spin()
            else:
                T = level.remove_res_cond(psis)
                ptc = ELGT_addspin(T).to_prodtautclass_spin()
            # print(str(ptc0))
            # print(str(ptc1))
            alpha.append(ptc)
            vertices.append(LG.verticesonlevel(LG.internal_level_number(l)))

        if rearrange_markings:  # if we want the labels of stgraph correspond to the markings
            prodtautst = self._G.relabel(self._G.standard_markings(), tidyup=False).LG.stgraph
            prod = admcycles.admcycles.prodtautclass(prodtautst)

        else:
            prod = admcycles.admcycles.prodtautclass(stgraph)

        # To clutch the level stratum classes together
        for l, ptc in enumerate(alpha):
            prod = prod.factor_pullback(vertices[l], ptc)

        return self.stack_factor * prod


def standardise_sig(sig, res_cond=[]):
    r"""
    To standardise the arrangement of entries in the signature, such that it go from low to high and all the 0s
    are placed at the right end. This is just to reduce the required memory and time to compute the stratum classes.
    It returns

            1. the simplified sig (0s removed);
            2. number of 0s;
            3. the leg_dict of rearrangement;
            4. new res_cond

    INPUT:

        - sig (tuple): a tuple of integers
        - res_cond (list, default=[]): the list of res_cond

    EXAMPLES::

        sage: from admcycles.diffstrata.spinstratum import standardise_sig
        sage: standardise_sig((2,4,-1,0,-1), res_cond=[[(0,2),(0,4)]])
        [(-1, -1, 2, 4), 1, {1: 3, 2: 4, 3: 1, 4: 5, 5: 2}, [((0, 0), (0, 1))]]

    """
    res = deepcopy(res_cond)
    new_rc = []
    l = list(set(sig))
    g = Signature(sig).g
    l.sort()
    t = list(sig)
    freq = {x: t.count(x) for x in l}
    counter = deepcopy(freq)
    standard_sig_without_regular = []

    for x in l:
        if x != 0:
            standard_sig_without_regular += freq[x] * [x]
    m = len(standard_sig_without_regular)
    if g == 0 and m < 3:
        standard_sig_without_regular += (3 - m) * [0]

    leg_dict = {}
    for i in range(0, len(sig)):
        if sig[i] == 0:
            p = sum([f for j, f in freq.items() if j != 0])
            leg_dict[i + 1] = p + freq[0] - counter[0] + 1
            counter[0] -= 1
        if sig[i] != 0:
            p = sum([f for j, f in freq.items() if (j != 0 and j < sig[i])])
            leg_dict[i + 1] = p + freq[sig[i]] - counter[sig[i]] + 1
            counter[sig[i]] -= 1
    output = []
    output.append(tuple(standard_sig_without_regular))
    if 0 in l:
        if g == 0 and m < 3:
            output.append(freq[0] - 3 + m)
        else:
            output.append(freq[0])
    else:
        output.append(0)
    output.append(leg_dict)

    for r in res:
        w = []
        for q in r:
            a = (0, leg_dict[q[1] + 1] - 1)
            w.append(a)
        new_rc.append(tuple(w))
    output.append(new_rc)
    return output


@cached_function
def smooth_stratum_standard(sig, res_cond=()):
    r"""
    This function is to compute the spin stratum class of a standardised
    signature i.e. it is ordered and has no 0s (or all are 0s), also the
    stratum has only one component. The results are cached.

    EXAMPLES::

        sage: from admcycles.diffstrata.spinstratum import smooth_stratum_standard
        sage: smooth_stratum_standard((-1,-1,4))
        [2,
         3,
         2,
         (0, 0, 4, 0, -11/2, -12, 0, -7/2, 0, 1/2, 35/2, 4, 0, -12, 0, 0, 0, -2, 19/2, 0, 7/2, 0, 0, 0, 2, 0, -7/2, 0, 2, -2, 5/2, -13/2, 9/2, 8, 0, 0, 0, 17/2, -11, 2, -2, -9/2, 0, 2)]


    """
    assert (sum([x for x in sig]) == -2 and len(sig) == 3) or all(j == 0 for j in sig) or all(j != 0 for j in sig)

    if res_cond != ():
        rc = []
        for q in res_cond:
            res_single_cond = [t for t in q]
            rc.append(res_single_cond)
        X = SpinStratum([Signature(sig)], rc)
    else:
        X = SpinStratum([Signature(sig)])

    # LG = X.smooth_LG.LG
    # stgraph = LG.stgraph
    g = X._sig_list[0].g
    n = X._sig_list[0].n
    r = 0
    if any(x < 0 for x in sig):
        r = len(res_cond) + 1
    k = g + r - 1
    # The codim we are dealing with for computation of stratum class.

    # length=len(admcycles.admcycles.generating_indices(g,n,k))
    check_basic = is_base_case(sig, g, n, k)

    # To check whether the stratum class is the base cases
    if check_basic != []:
        return check_basic

    # Solve the pullback equations
    ListM = Pullback_Matrices(g, n, r)
    prodtaut_list = bdry_pullback_classify(X)
    tensorTaut_list = []
    for i, u in enumerate(prodtaut_list):
        if u != []:
            t = sum(u)
            v = vector(t.totensorTautbasis(k, vecout=True))
            tensorTaut_list.append(v)
        else:
            length = len(ListM[i][0])
            v = vector(QQ, length)
            tensorTaut_list.append(v)

    sol_inTautbasis = SolveMultLin(ListM, tensorTaut_list)
    return [g, n, k, sol_inTautbasis]

# To give the stratum classes with spin for the easy cases.


def is_base_case(sig, g, n, k):
    if g == 0:
        if all(q % 2 == 0 for q in sig):
            return [g, n, k, vector([QQ(1)])]

        # case of (0,-1,-1)
        if n == 3:
            return [g, n, k, vector([-QQ(1)])]

        # case of (2k,-2k,-1,-1)
        if n == 4:
            result = admcycles.admcycles.psiclass(1, 0, 4).basis_vector()
            return [g, n, k, result]

    if g == 1:
        if all(t == 0 for t in sig):
            return [g, n, k, vector([-QQ(1)])]

    if g == 2 and n == 1:
        return [g, n, k, vector([QQ(1) / 2, -QQ(7) / 2, QQ(1) / 2])]

    return []


def bdry_pullback_classify(stratum):
    X = stratum
    sig = X._sig_list[0].sig
    g = X._sig_list[0].g
    n = X._sig_list[0].n
    # r = 0
    # if any(m < 0 for m in sig):
    #    r = len(X.res_cond) + 1
    # k = g + r - 1
    Agraphs = admcycles.admcycles.list_strata(g, n, 1)
    if g != 0 and n >= 3:
        Agraphs = [G for i, G in enumerate(Agraphs) if i != 0]
    else:
        G = Agraphs[0].copy()
        G.tidy_up()
        Agraphs[0] = G

    bddiv_codim1_LG_prodtautclass_list = [[] for G in Agraphs]
    bic_stgraph = [bic.relabel(bic.standard_markings(), tidyup=False) for bic in X.bics]

    for j, G in enumerate(Agraphs):
        sig0 = [sig[i - 1] for i in G._legs[0] if i < n + 1]
        node_sing = 2 * G._genera[0] - 2 - sum(sig0)

        if node_sing == -1:
            sig0.append(-1)
            sig1 = [sig[i - 1] for i in G._legs[1] if i < n + 1] + [-1]
            new_sig_list = [sig0, sig1]
            splt_st0 = GeneralisedStratum([Signature(tuple(new_sig_list[0]))])
            splt_st1 = GeneralisedStratum([Signature(tuple(new_sig_list[1]))])
            elgt0 = splt_st0.additive_generator(((), 0)).as_taut()
            elgt1 = splt_st1.additive_generator(((), 0)).as_taut()
            st0_spin = ELGT_addspin(elgt0).to_prodtautclass_spin()
            st1_spin = ELGT_addspin(elgt1).to_prodtautclass_spin()
            prodtaut = admcycles.admcycles.prodtautclass(G)

            # Note that for a horizontal one edge graph, the splitting of the
            # two vertices will give rise to extra symplectic basis and thus
            # if sum of spins are odd (even), then the horizontal graph
            # has spin even(odd).
            clutch = prodtaut.factor_pullback([0], st0_spin)
            clutch = clutch.factor_pullback([1], -st1_spin)
            bddiv_codim1_LG_prodtautclass_list[j].append(clutch)

        else:
            if len(G._genera) == 1:
                sig_new = ([sig[i - 1] for i in G._legs[0] if i < n + 1] +
                           [-1, -1])
                Y = SpinStratum([Signature(tuple(sig_new))],
                                res_cond=[[(0, n), (0, n + 1)]])
                elgt = Y.additive_generator(((), 0)).as_taut()
                prodtaut = admcycles.admcycles.prodtautclass(G)

                stratum_class = ELGT_addspin(elgt).to_prodtautclass_spin()
                Intersect = prodtaut.factor_pullback([0], stratum_class)
                # print(str(Intersect))
                bddiv_codim1_LG_prodtautclass_list[j].append(Intersect)
            for num, bic in enumerate(bic_stgraph):
                stgraph = bic.LG.stgraph
                a = G._edges[0][0]
                b = G._edges[0][1]

                Astructure_list = admcycles.admcycles.Astructures(stgraph, G)

                if Astructure_list != []:
                    # print(str(Astructure_list))
                    AG = AG_addspin(X.additive_generator(((num,), 0)))
                    prodtaut_temp = AG.to_prodtautclass_spin()
                    # print(str(prodtaut_temp))
                    ell = bic.ell
                    for w in Astructure_list:
                        G_no_edge = StableGraph(G._genera, G._legs, [])
                        image_half_edges = [(w[1][e[0]], w[1][e[1]])
                                            for e in G._edges]
                        st_edges_new = [e for e in stgraph._edges if
                                        (e not in image_half_edges) and
                                        ((e[1], e[0]) not in image_half_edges)]
                        st_new = StableGraph(stgraph._genera,
                                             stgraph._legs, st_edges_new)
                        prodtaut = admcycles.admcycles.prodtautclass(st_new, prodtaut_temp.terms)
                        # print(str(G_no_edge)+"\n"+str(st_new))
                        E_temp = prodtaut.partial_pushforward(G_no_edge,
                                                              w[0], w[1])
                        E = admcycles.admcycles.prodtautclass(G, E_temp.terms)
                        # print(str(E_temp)+"\n" +str(E))
                        kappa_e = bic.LG.prong((w[1][a], w[1][b]))
                        Intersect = QQ(ell / kappa_e) * E
                        bddiv_codim1_LG_prodtautclass_list[j].append(Intersect)

    return bddiv_codim1_LG_prodtautclass_list


# To generate a list of matrices of the clutching pullback maps.
@cached_function
def Pullback_Matrices(g, n, r):
    k = g + r - 1
    L1 = admcycles.admcycles.tautgens(g, n, k)
    Listgen1 = admcycles.admcycles.generating_indices(g, n, k)

    graph_list = admcycles.admcycles.list_strata(g, n, 1)
    if g != 0 and n >= 3:
        graph_list = [G for i, G in enumerate(graph_list) if i != 0]

    M_list = []
    for G in graph_list:
        M0 = []
        if len(G._genera) == 1:
            G = G.copy()
            G.tidy_up()
        for i in Listgen1:
            T = G.boundary_pullback(L1[i])
            M0 = M0 + [T.totensorTautbasis(k, vecout=True)]
        M = matrix(M0)
        M_list.append(M)
    return M_list


#######################################################################
#######################################################################
# Define the spin version ELGTautClass ################################
#######################################################################


class ELGT_with_spin (ELGTautClass):

    def __init__(self, X, psi_list):
        assert isinstance(X, SpinStratum)
        super().__init__(X, psi_list)
        self.X = X

    def __repr__(self):
        return "ELGT_SpinClass(X=%r,psi_list=%r)"\
            % (self._X, self._psi_list)

    def __str__(self):
        str = "ELGT_SpinClass on %s\n" % self._X
        for coeff, psi in self._psi_list:
            str += "%s * %s + \n" % (coeff, psi)
        return str

    def to_prodtautclass_spin(self):
        G = self.X.smooth_LG
        # print(str(G))
        stgraph = G.LG.stgraph
        total = admcycles.admcycles.prodtautclass(stgraph, terms=[])
        for c, AG in self._psi_list:
            legdictlen = len(AG._G.dmp)
            ptc = AG.to_prodtautclass_spin()
            vertex_map = {}
            for v, _ in enumerate(G.LG.genera):
                mp_on_stratum = G.dmp[G.LG.legs[v][0]]
                l_AG = AG._G.dmp_inv[mp_on_stratum]
                LG = AG._G.LG
                v_AG = LG.vertex(l_AG)
                UG_v = LG.UG_vertex(v_AG)
                T = LG.UG_without_infinity()
                T = T.connected_component_containing_vertex(UG_v)
                for w, g, kind in T:
                    if kind != 'LG':
                        continue
                    vertex_map[w] = v
            leg_map = {i: i for i in range(1, legdictlen + 1)}
            # print(str(stgraph)+"\n"+ str(leg_map)+ "\n"+ str(vertex_map)+ "\n" + str(self)+
            # "\n" + str(ptc))
            pf = ptc.partial_pushforward(stgraph, vertex_map, leg_map)
            # print(str(pf) + "\n" + str(ptc)+ "\n" + str(c))
            total += c * pf
        return total


#############################################################################
# Some useful functions #####################################################
#############################################################################


# turn a ELGTautClass to spin version

def ELGT_addspin(elgt):
    X = SpinStratum(elgt._X._sig_list, elgt._X._res_cond)
    psi_list = [(c, AG_addspin(AG)) for c, AG in elgt._psi_list]
    return ELGT_with_spin(X, psi_list)


# turn a AdditiveGenerator to spin version

def AG_addspin(ag):
    assert isinstance(ag, AdditiveGenerator)
    X = addspin(ag._X)
    return AG_with_spin(X, ag._enh_profile, ag._leg_dict)


# turn a stratum to spin version

def addspin(stratum):
    assert isinstance(stratum, GeneralisedStratum)
    return SpinStratum(stratum._sig_list, stratum._res_cond)


def save_spin_stratum():
    cachedict = dict(smooth_stratum_standard.cache)
    sage.misc.persist.save(cachedict, 'spin_stratum_cache')


def load_spin_stratum():
    cachedict = sage.misc.persist.load('spin_stratum_cache.sobj')
    for a, b in cachedict.items():
        smooth_stratum_standard.set_cache(b, *a[0])
