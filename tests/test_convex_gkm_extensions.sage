"""Independent incidence geometry, equivariant K classes, and old F1 behavior."""

from log_glsm.cycles.marked import (
    genus_four_backend,
)
from log_glsm.targets.convex_gkm import (
    convex_gkm_ordinary_seed,
    convex_gkm_targets,
    partial_flag_fibre_lines,
    partial_flag_primitive_counts,
)
from log_glsm.targets.toric import (
    toric_genus_zero_number,
    toric_sample,
)
targets=convex_gkm_targets()
for t in targets:
    A=t.ambient
    assert len(A.tangent)==len(A.points)
    for (i,j),edge in A.edges.items():
        w=vector(edge['weight'])
        end=[tuple(-w)]+[tuple(vector(ell)-d*w) for d,ell in edge['normals']]
        assert sorted(end)==sorted(A.tangent[j]),(A.slug,i,j)
        assert 2+sum(d for d,_ in edge['normals'])==A.pairing(A.anticanonical,edge['class'])
    for base in (4,5):A.ray(base)
    print(t.slug,t.record(),flush=True)
by={t.slug:t for t in targets}
assert by['NF1'].chern_numbers()=={'c1_squared':8,'euler':4}
assert by['NF4'].chern_numbers()=={'c1_squared':8,'euler':4}
assert by['NF5'].chern_numbers()=={'c1_cubed':54,'c1_c2':24,'euler':6,'C':-18}
assert by['NF2'].chern_numbers()=={'c1_squared':0,'euler':24}
assert by['PF1'].chern_numbers()=={'c1_squared':6,'euler':6}
assert by['PF2'].chern_numbers()=={'c1_squared':0,'euler':24}
# Incidence descriptions give an independent product-projective calculation.
for key,N,degrees in (('PF1',2,((1,1),(1,1))),('PF2',2,((1,1),(2,2))),
                       ('PF5',3,((1,1),(1,1),(1,1))),('PF6',3,((1,1),(1,1),(2,2)))):
    P=PolynomialRing(QQ,2,('u','v'));u,v=P.gens();X=P.quotient(P.ideal(u**(N+1),v**(N+1)))
    c=X((1+u)**(N+1)*(1+v)**(N+1));fund=X.one()
    for a,b in degrees:
        L=X(a*u+b*v);fund*=L;c*=sum((-L)**k for k in range(2*N+1))
    ch=lambda k:X(P({m:co for m,co in c.lift().dict().items() if sum(m)==k}))
    integrate=lambda p:QQ(X(fund*p).lift().monomial_coefficient(u**N*v**N))
    if by[key].dimension==2:expected={'c1_squared':integrate(ch(1)**2),'euler':integrate(ch(2))}
    else:
        c12=integrate(ch(1)*ch(2));c3=integrate(ch(3))
        expected={'c1_cubed':integrate(ch(1)**3),'c1_c2':c12,'euler':c3,'C':c3-c12}
    assert by[key].chern_numbers()==expected,(key,expected)
for t in toric_sample():
    if 'F1' in t.ambient.factors:
        for g in (2,3):
            import json
            with open('tests/fixtures/f1-numbers.json') as fixture:
                old=json.load(fixture)[t.slug+'-g%s-numbers.json'%g]
            assert old['target']['chern_numbers']=={k:str(v) for k,v in t.chern_numbers().items()}
            b=genus_four_backend(t,genus=g)
            for row in old['rows'][:2]:
                value=b.numerical_localization(tuple(row['beta']))
                assert not value['poles'] and value['finite']==QQ(row['twisted_unmarked'])
print('Incidence, endpoint K-theory, and existing F1 regressions passed',flush=True)
for key in ('NF2','NF3','NF4','PF1','PF4'):
    t=by[key];beta=(1,)+(0,)*(t.ambient.picard_rank-1)
    v=[genus_four_backend(t,base,genus=2).numerical_localization(beta) for base in (4,5)]
    assert all(not x['poles'] for x in v) and v[0]['finite']==v[1]['finite'],(key,v)
    print('Independent rays',key,beta,v[0]['finite'],flush=True)
print('Convex GKM extension tests passed',flush=True)
assert partial_flag_fibre_lines()==33
assert convex_gkm_ordinary_seed(by['PF4'],2,(1,0))==QQ(11)/80
assert convex_gkm_ordinary_seed(by['PF4'],3,(1,0))==QQ(11)/2016
assert convex_gkm_ordinary_seed(by['PF4'],2,(2,0))==QQ(11)/40
assert convex_gkm_ordinary_seed(by['PF4'],3,(2,0))==QQ(11)/252
print('Flag fibre-line count and multiple-cover values passed',flush=True)
for k in (1,2):
    for base in (4,5):assert toric_genus_zero_number(by['PF4'],(k,0),base)==QQ(33)/k**3
print('Independent genus-zero localization verifies 33 and 33/8 on two rays',flush=True)
assert partial_flag_primitive_counts()=={(1,0):33,(0,1):88,(0,2):5,(1,1):634}
for beta,expected in (((0,1),88),((0,2),16),((1,1),634)):
    for base in (4,5):assert toric_genus_zero_number(by['PF4'],beta,base)==expected
assert convex_gkm_ordinary_seed(by['PF4'],2,(0,2))==QQ(181)/240
assert convex_gkm_ordinary_seed(by['PF4'],3,(0,2))==QQ(709)/6048
assert convex_gkm_ordinary_seed(by['PF4'],2,(1,1))==QQ(317)/120
assert convex_gkm_ordinary_seed(by['PF4'],3,(1,1))==QQ(317)/3024
print('All PF4 bounded primitive counts, conics, and higher-genus formulas passed',flush=True)
