"""Independent geometry and numerical checks for the Grassmannian extension."""
from itertools import combinations
from log_glsm.targets.grassmannian import (
    GrassmannianAmbient,
    GrassmannianCompleteIntersection,
    grassmannian_known_ordinary_gw,
    grassmannian_line_count,
    grassmannian_sample,
)
from log_glsm.targets.projective import (
    ProjectiveCompleteIntersection,
    cy3_line_count_localization,
)
from log_glsm.targets.toric import (
    toric_genus_zero_number,
)


def test_grassmannian_geometry():
    for n,degree,points,curves in ((3,1,3,3),(4,2,6,12),(5,5,10,30)):
        X=GrassmannianAmbient(n);H=X.divisor((1,));D=X.dimension
        assert X.integral(H**D)==degree
        assert X.integral(X.homogeneous(X.total_chern,D))==points
        assert X.integral(X.homogeneous(X.total_chern,1)*H**(D-1))==n*degree
        assert X.record()['fixed_points']==points and X.record()['invariant_curves']==curves
        # Independent fixed-point integration of H^D and every top Chern
        # number used by surface/threefold targets, versus the root algebra.
        for base in (4,5):
            a=X.ray(base)
            for parts in ((1,)*D,(D,),(D-1,1),(D-2,2)):
                ordinary=X.integral(prod(X.homogeneous(X.total_chern,k) for k in parts))
                localized=QQ.zero()
                for i in range(len(X.points)):
                    ws=[X.character_value(w,a) for w in X.tangent[i]]
                    elementary=lambda k:sum(prod(v) for v in combinations(ws,k))
                    localized+=prod(elementary(k) for k in parts)/prod(ws)
                assert ordinary==localized,(n,parts,base)


def test_quadric_geometry():
    # Gr(2,4) is Q4 in P5. Adjunction must reproduce every Chern number
    # of the same target presented as a projective complete intersection.
    for t in grassmannian_sample()[:8]:
        p=ProjectiveCompleteIntersection(5,(2,)+t.degrees)
        c=p.chern_coefficients();actual=t.chern_numbers()
        assert actual['euler']==p.euler_characteristic()
        if t.dimension==2:assert actual['c1_squared']==p.degree*c[1]**2
        else:
            assert actual['c1_cubed']==p.degree*c[1]**3
            assert actual['c1_c2']==p.degree*c[1]*c[2]


def test_cy_curve_counts():
    # Independent enumerative counts; genus-zero degree two includes covers.
    expected={'G4':(1280,92448),'G11':(540,QQ(25245)/2),'G12':(400,5590)}
    for t in grassmannian_sample():
        if t.kind!='CY3':continue
        line,degree_two=expected[t.slug]
        assert grassmannian_line_count(t)==line
        assert grassmannian_line_count(t,[2**i for i in range(t.ambient.n)])==line
        for base in (4,5):
            assert toric_genus_zero_number(t,(1,),base)==line
            assert toric_genus_zero_number(t,(2,),base)==degree_two
        for g,den in ((2,240),(3,6048)):
            assert grassmannian_known_ordinary_gw(t,g,(1,))==QQ(line)/den
            assert grassmannian_known_ordinary_gw(t,g,(2,)) is None
    t=grassmannian_sample()[3]
    assert grassmannian_line_count(t)==cy3_line_count_localization(ProjectiveCompleteIntersection(5,(2,4)))


def test_theory_boundaries():
    for t in grassmannian_sample():
        if t.kind=='K3':
            assert t.chern_numbers()=={'c1_squared':0,'euler':24}
            assert grassmannian_known_ordinary_gw(t,3,(2,))==0
        if t.kind=='Fano':assert grassmannian_known_ordinary_gw(t,3,(1,)) is None
    for n,degrees in ((6,(1,)*5),(5,(0,1,1)),(5,(1,1,4)),(4,(1,1,1))):
        try:GrassmannianCompleteIntersection(n,degrees)
        except ValueError:pass
        else:raise AssertionError('unsupported target accepted')


test_grassmannian_geometry()
test_quadric_geometry()
test_cy_curve_counts()
test_theory_boundaries()
print('Grassmannian geometry, quadric comparison, CY curve counts, and theory-boundary checks passed')
