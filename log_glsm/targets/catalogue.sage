r"""targets / catalogue."""

from sage.all import *


from log_glsm.core.hodge import (
    HodgeIntegralBackend,
)
from log_glsm.targets.grassmannian import (
    GrassmannianCompleteIntersection,
    grassmannian_line_count,
    grassmannian_sample,
    grassmannian_seven_hyperplane_cy3,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
    requested_complete_intersections,
)
from log_glsm.targets.toric import (
    ToricCompleteIntersection,
    toric_degree_tuples,
    toric_known_ordinary_gw,
    toric_primitive_curve_counts,
    toric_sample,
)

__all__ = ['genus_four_degrees', 'genus_four_kind', 'genus_four_targets', 'known_genus_four_ordinary']


def genus_four_targets():
    return ((ProjectiveCompleteIntersection(2,(3,),'elliptic cubic'),
             ProjectiveCompleteIntersection(3,(4,),'quartic K3'))
            +requested_complete_intersections()+toric_sample()+grassmannian_sample()
            +(grassmannian_seven_hyperplane_cy3(),))

def genus_four_degrees(target):
    return ((0,),(1,),(2,)) if isinstance(target,ProjectiveCompleteIntersection) else toric_degree_tuples(target.ambient)

def genus_four_kind(target):
    if not isinstance(target,ProjectiveCompleteIntersection):return target.kind
    if target.index:return 'Fano'
    return {1:'elliptic',2:'K3',3:'CY3'}[target.dimension]

def known_genus_four_ordinary(target,beta):
    """Unmarked ordinary numbers, not zero-vertex twisted numbers."""
    kind=genus_four_kind(target);beta=tuple(beta)
    if kind=='K3':return QQ.zero()
    c1degree=(target.index*beta[0] if isinstance(target,ProjectiveCompleteIntersection)
              else target.ambient.pairing(target.c1,beta))
    vdim=-3*(target.dimension-3)+c1degree
    if vdim:return QQ.zero()
    if not any(beta):
        if isinstance(target,ProjectiveCompleteIntersection):
            ch=target.chern_coefficients();C=target.degree*(ch[3]-ch[1]*ch[2])
        else:C=target.chern_numbers()['C']
        return QQ(C)*HodgeIntegralBackend.top_hodge_triple_constant(4)
    if kind!='CY3':return None
    if isinstance(target,ProjectiveCompleteIntersection) and beta==(1,):
        P=PolynomialRing(QQ,2,('x','y'));x,y=P.gens()
        integrand=(x-y)*prod(prod(a*x+(d-a)*y for a in range(d+1)) for d in target.degrees)
        lines=integrand.monomial_coefficient(x**target.N*y**(target.N-1))
        return QQ(lines)/172800
    if isinstance(target,GrassmannianCompleteIntersection):
        if beta==(1,):return QQ(grassmannian_line_count(target))/172800
        if target.ambient.n==7 and target.degrees==(1,)*7 and beta==(2,):return QQ(1225+32*196)/172800
    elif isinstance(target,ToricCompleteIntersection):
        counts=toric_primitive_curve_counts(target)
        if beta in counts:return QQ(counts[beta])/172800
        old=toric_known_ordinary_gw(target,3,beta)
        if old is not None:return old*QQ(4*6048)/172800
    return None
