r"""Genus-two/three complete-intersection calculations with labeled marks.

The ambient catalogue and degree-two GKM data are the same 43 presentations
as the genus-four pass. Markings are labeled and psi means a stable-map
descendant. Class searches and numerical reductions have separate outputs."""

from sage.all import *

from itertools import combinations_with_replacement
from itertools import permutations
from itertools import product
from log_glsm.core.hodge import (
    HodgeIntegralBackend,
)
from log_glsm.cycles.marked import (
    genus_four_backend,
)
from log_glsm.targets.catalogue import (
    genus_four_degrees,
    genus_four_kind,
    known_genus_four_ordinary,
)
from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersection,
    grassmannian_known_ordinary_gw,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersection,
    toric_known_ordinary_gw,
)

__all__ = ['higher_marked_backend', 'higher_marked_ordinary_seed', 'higher_marked_profiles', 'marked_profile_permutations']


def higher_marked_profiles(canonical=False,max_markings=3):
    if max_markings not in ZZ or max_markings<0:raise ValueError('nonnegative marking bound required')
    alphabet=('1','psi','H')
    if canonical:
        return tuple(p for n in range(max_markings+1) for p in combinations_with_replacement(alphabet,n))
    return tuple(p for n in range(max_markings+1) for p in product(alphabet,repeat=n))

def higher_marked_backend(target,genus,base=4,hodge=None):
    if genus not in (2,3):raise ValueError('this pass is for genera two and three')
    return genus_four_backend(target,base,hodge,genus)

def higher_marked_ordinary_seed(target,genus,beta):
    if genus not in (2,3):raise ValueError('this pass is for genera two and three')
    beta=tuple(beta)
    if beta not in genus_four_degrees(target):raise ValueError('outside the bounded curve catalogue')
    kind=genus_four_kind(target)
    if kind=='K3':return QQ.zero()
    c1degree=(target.index*beta[0] if isinstance(target,ProjectiveCompleteIntersection)
              else target.ambient.pairing(target.c1,beta))
    if (1-genus)*(target.dimension-3)+c1degree:return QQ.zero()
    if not any(beta):
        if isinstance(target,ProjectiveCompleteIntersection):
            ch=target.chern_coefficients();C=target.degree*(ch[3]-ch[1]*ch[2])
        else:C=target.chern_numbers()['C']
        return (-1)**genus*C*HodgeIntegralBackend.top_hodge_triple_constant(genus)
    if kind!='CY3':return None
    if isinstance(target,GrassmannianCompleteIntersection):
        return grassmannian_known_ordinary_gw(target,genus,beta)
    if isinstance(target,ToricCompleteIntersection):return toric_known_ordinary_gw(target,genus,beta)
    if beta==(1,):
        return known_genus_four_ordinary(target,beta)*172800/(240 if genus==2 else 6048)
    return None

def marked_profile_permutations(profile):
    """One old-label -> new-label map per distinct ordered profile."""
    profile=tuple(profile);seen=set();result=[]
    for order in permutations(range(len(profile))):
        renamed=tuple(profile[old] for old in order)
        if renamed in seen:continue
        seen.add(renamed)
        mapping=tuple(order.index(old)+1 for old in range(len(profile)))
        result.append((renamed,mapping))
    return tuple(result)
