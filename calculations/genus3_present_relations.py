"""Verify simplified displayed relations and export exact readable formulae."""
import json
from pathlib import Path
from sage.all import *

from log_glsm.cycles.plane import (
    O3UnmarkedCycleBackend,
)
out=Path('results/genus3/lower-cycles-2026-11-05')
b=O3UnmarkedCycleBackend()

def I(k=None,l=()):return ('interior',k,tuple(l))
def J(p=0,l=(),a=0):return ('boundary12',0,(() if a==0 else (1,)),p,tuple(l))

# These are simplified using kappa_0=4 and psi=lambda_1 on Mbar_(1,1).
displayed={
 'R3a':[(3,I(3)),(-5,I(2,(1,))),(7,I(1,(2,))),(-8,I(None,(1,2))),(-12,I(None,(3,))),
        (-5,J(0,(2,))),(3,J(1,(1,))),(-1,J(2)),(-2,J(0,(1,),1)),(2,J(1,(),1))],
 'R3b':[(2,I(1,(2,))),(-1,I(2,(1,))),(-4,I(None,(1,2))),
        (-2,J(0,(2,))),(1,J(1,(1,))),(1,J(1,(),1)),(-1,J(0,(1,),1))],
 'R4': [(-11,I(4)),(17,I(3,(1,))),(-24,I(2,(2,))),(7,I(1,(1,2))),
        (11,I(1,(3,))),(-32,I(None,(1,3))),
        (2,J(3)),(-5,J(2,(1,))),(9,J(1,(2,))),(-4,J(0,(1,2))),
        (1,J(1,(1,),1)),(-2,J(0,(2,),1))],
 'R5': [(36,I(5)),(-116,I(4,(1,))),(280,I(3,(2,))),(-173,I(2,(1,2))),
        (-31,I(2,(3,))),(279,I(1,(1,3))),(-460,I(None,(2,3))),
        (6,J(4)),(-7,J(3,(1,))),(-28,J(2,(2,))),(44,J(1,(1,2))),
        (8,J(3,(),1)),(-32,J(2,(1,),1)),(78,J(1,(2,),1)),(-46,J(0,(1,2),1))],
}
certificates=[]
for name,terms in displayed.items():
    r=int(name[1])
    cycle=sum((co*b.monomial_class(key) for co,key in terms),b.R.zero())
    v=vector(QQ,cycle.vector(r))
    M=load(str(out/('pixton-codimension-%s.sobj'%r)))
    w=M.transpose().solve_right(v)
    assert w*M==v
    certificates.append({'name':name,'codimension':r,
        'terms':[{'coefficient':str(co),'key':key} for co,key in terms],
        'strata_vector':[str(x) for x in v],
        'pixton_weights':[[int(i),str(x)] for i,x in enumerate(w) if x]})
    print(name,'verified',flush=True)

def simplified(terms):
    answer={}
    for term in terms:
        key=term['key'];co=QQ(term['coefficient'])
        if key[0]=='interior':
            k=key[1]
            if k==0:k=None;co*=4
            key=I(k,key[2])
        else:
            left_degree=key[1]+sum(key[2])
            assert left_degree<=1
            key=J(key[3],key[4],left_degree)
        answer[key]=answer.get(key,0)+co
    return {k:c for k,c in answer.items() if c}

# Check that each displayed equation actually comes from the computed poles,
# in addition to checking its independent Pixton certificate above.
origins={
    'R3a': [(-2,QQ(36)/49)],
    'R3b': [(-2,QQ(216)/245),(-1,QQ(9)/70)],
    'R4': [(-1,QQ(9)/35)],
    'R5': [(-1,QQ(9)/7)],
}
for cert in certificates:
    name=cert['name']; r=cert['codimension']; combined={}
    poles=json.loads((out/('pixton-certificates-codimension-%s.json'%r)).read_text())
    for power,scale in origins[name]:
        pole=next(p for p in poles if p['epsilon_power']==power)
        for key,co in simplified(pole['formal_terms']).items():
            combined[key]=combined.get(key,0)+scale*co
    combined={k:c for k,c in combined.items() if c}
    assert combined=={k:QQ(c) for c,k in displayed[name]},('incorrect pole origin',name)
    cert['pole_origin']=[{'epsilon_power':power,'scale':str(scale)}
                         for power,scale in origins[name]]
(out/'displayed-relations.json').write_text(json.dumps(certificates,indent=2,default=int)+'\n')

def lambdas(indices,symbol='\\lambda'):
    return ''.join('%s_{%s}'%(symbol,i) for i in indices)
def monomial(key):
    if key[0]=='interior':
        return (('\\kappa_{%s}'%key[1]) if key[1] is not None else '')+lambdas(key[2]) or '1'
    left=lambdas(key[2],'a') or '1'
    right=(('\\psi^{%s}'%key[3]) if key[3] else '')+lambdas(key[4],'b') or '1'
    return '\\xi_*\\bigl(%s\\boxtimes %s\\bigr)'%(left,right)
def equation(terms):
    pieces=[]
    for key,co in terms.items():
        sign='+' if co>0 else '-'
        scalar='' if abs(co)==1 else str(latex(abs(co)))+' '
        pieces.append(sign+scalar+monomial(key))
    return ' '.join(pieces).lstrip('+')

lines=['Exact lower-codimension genus-three formulas.','',
 'Write a_1=lambda_1 on Mbar_(1,1), b_i=lambda_i on Mbar_(2,1), and psi=psi_1 on Mbar_(2,1).',
 'The map xi glues their marked points. Global lambda and kappa classes lie on Mbar_3.',
 'These Z classes use the dual Euler convention; the CJR infinity aggregate is B=-Z.','']
for rec in json.loads((out/'cycle-formulas.json').read_text()):
    if rec['codimension']>5:continue
    terms=simplified(rec['terms'])
    lines+=['\\[ Z_{%s,%s} = %s. \\]'%(rec['ambient_degree'],rec['codimension'],equation(terms) if terms else '0'),'']
for name,terms in displayed.items():
    lines+=[r'\[ %s:\quad %s =0. \]'%(name,equation({k:QQ(c) for c,k in terms})),'']
(out/'explicit-formulas.md').write_text('\n'.join(lines))

lines=[
    'All nine auxiliary-pole relations on Mbar_3, with primitive integer coefficients.',
    '',
    'Notation: a_1=lambda_1 on Mbar_(1,1); b_i=lambda_i and psi=psi_1 on Mbar_(2,1).',
    'The gluing map xi maps their product to Mbar_3. Global lambda and kappa classes lie on Mbar_3.',
    'Each expression is the indicated pole after kappa_0=4 and psi_1=lambda_1 on Mbar_(1,1), multiplied by a rational scalar.',
    '',
]
for r in (3,4,5):
    poles=json.loads((out/('pixton-certificates-codimension-%s.json'%r)).read_text())
    for pole in poles:
        terms=simplified(pole['formal_terms'])
        denominator=lcm([c.denominator() for c in terms.values()])
        divisor=gcd([ZZ(denominator*c) for c in terms.values()])
        scale=denominator/divisor
        if next(iter(terms.values()))<0:scale=-scale
        terms={k:scale*c for k,c in terms.items()}
        lines += ['Codimension %s, auxiliary power %s, scale %s.'%
                  (r,pole['epsilon_power'],scale),'',
                  r'\[ %s = 0. \]'%equation(terms),'']
(out/'all-pole-relations.md').write_text('\n'.join(lines))
