"""Export readable relation bases and replay exact rational certificates."""
import gzip,json
from pathlib import Path
from sage.all import *
out=Path('results/quartic-k3/2026-11-05')


def freeze(value):
    return tuple(freeze(x) for x in value) if isinstance(value,list) else value


def normalized(terms,genus):
    result={}
    for term in terms:
        key=freeze(term['key']);co=QQ(term['coefficient'])
        if key[0]=='interior' and key[1]==0:
            key=('interior',None,key[2]);co*=2*genus-2
        result[key]=result.get(key,0)+co
    return {k:c for k,c in result.items() if c}


graph_labels={}
def monomial(key,genus):
    def lambdas(indices,symbol=r'\lambda'):
        return ''.join(r'%s_{%s}'%(symbol,i) for i in indices)
    if key[0]=='interior':
        return ((r'\kappa_{%s}'%key[1] if key[1] is not None else '')+lambdas(key[2])) or '1'
    if key[0]=='boundary':
        h,p,l,q,m=key[1:]
        left=(r'\psi_1^{%s}'%p if p else '')+lambdas(l,'a') or '1'
        right=(r'\psi_2^{%s}'%q if q else '')+lambdas(m,'b') or '1'
        return r'\xi_{%s,%s,*}(%s\boxtimes %s)'%(h,genus-h,left,right)
    data,powers,indices=key[1:]
    if data not in graph_labels:graph_labels[data]=len(graph_labels)+1
    factors=[]
    for v,(ps,ls) in enumerate(zip(powers,indices),1):
        part=''.join(r'\psi_{%s,%s}^{%s}'%(v,i,p) for i,p in enumerate(ps,1) if p)
        part+=''.join(r'\lambda_{%s,%s}'%(v,i) for i in ls)
        factors.append(part or '1')
    return r'\Theta_{%s}\left[%s\right]'%(graph_labels[data],r'\boxtimes '.join(factors))


def equation(terms,genus):
    parts=[]
    for key,co in terms.items():
        parts.append(('-' if co<0 else '+')+
                     ('' if abs(co)==1 else latex(abs(co))+' ')+monomial(key,genus))
    return ' '.join(parts).lstrip('+') or '0'


def primitive(terms):
    if not terms:return terms
    denominator=lcm([c.denominator() for c in terms.values()])
    divisor=gcd([ZZ(denominator*c) for c in terms.values()])
    scale=denominator/divisor
    if next(iter(terms.values()))<0:scale=-scale
    return {k:scale*c for k,c in terms.items()}


intro=[
    'All classes have rational coefficients. Global kappa and lambda classes lie on Mbar_g.',
    'The map xi_(h,k) glues Mbar_(h,1) times Mbar_(k,1); a_i and b_i are Hodge classes on its two factors.',
    'Theta_j means: glue the decorated graph listed in the graph dictionary below, then forget all its temporary tail markings.',
    'In Theta_j, psi_(v,i) is the i-th local flag in the listed order at vertex v; lambda_(v,i) is its i-th Hodge class.',
    'The gluing maps are actual pushforwards, without extra automorphism divisors. Fixed-graph and deck factors are already in the coefficients.',
    '',
]
lines=['# Independent quartic localization relations','',
       'These are a row-reduced basis of the computed relation space. Each S_(g,r,j) is the decorated gluing pushforward specified in the strata dictionary for that genus and codimension.',
       'Graph data lists vertex genera, ordered flags at each vertex, and paired node flags. Polynomial subscripts index vertices starting at zero. No extra graph automorphism factor is inserted.',
       'The original localization expressions and their origins are retained losslessly in final-relations-g*-r*.json.gz.','']
coverage=[];certificate_count=0
for g in (2,3):
    for r in range(3*g-2):
        path=out/('final-relations-g%s-r%s.json'%(g,r))
        records=json.loads(path.read_text() if path.exists() else gzip.decompress(
            Path(str(path)+'.gz').read_bytes()).decode())
        M=load(str(out/('pixton-g%s-r%s.sobj'%(g,r))))
        rows=[];weights=[]
        for index,record in enumerate(records,1):
            v=vector(QQ,record['strata_vector']);u=vector(QQ,M.nrows())
            for i,co in record['pixton_weights']:u[i]=QQ(co)
            rows.append(v);weights.append(u)
        if rows:
            V=matrix(QQ,rows);U=matrix(QQ,weights)
            assert U*M==V
            assert V.rank()==len(records)
            B=V.row_space().basis_matrix()
            W=M.transpose().solve_right(B.transpose()).transpose()
            assert W*M==B
            reduced=[]
            for index,(v,u) in enumerate(zip(B.rows(),W.rows()),1):
                reduced.append({'strata_vector':[str(x) for x in v],
                    'pixton_weights':[[int(i),str(x)] for i,x in enumerate(u) if x]})
                co=primitive({i:x for i,x in enumerate(v) if x})
                pieces=[('-' if x<0 else '+')+('' if abs(x)==1 else latex(abs(x))+' ')+
                        r'S_{%s,%s,%s}'%(g,r,i) for i,x in co.items()]
                lines+=['## Genus %s, codimension %s, relation %s'%(g,r,index),'',
                        r'\[ %s=0. \]'%' '.join(pieces).lstrip('+'),'']
            (out/('basis-relations-g%s-r%s.json'%(g,r))).write_text(json.dumps(reduced,indent=2)+'\n')
            strata=json.loads((out/('strata-g%s-r%s.json'%(g,r))).read_text())
            lines+=['Strata dictionary for (g,r)=(%s,%s):'%(g,r),'']
            lines+=['- S_(%s,%s,%s): `%s`.'%(g,r,i,s.replace('\n','; '))
                    for i,s in enumerate(strata)]
            lines+=['']
        coverage.append({'genus':g,'codimension':r,'relation_rank':len(records),
                         'pixton_rank':int(M.rank()),'strata_columns':int(M.ncols())})
        certificate_count+=len(records)

finite_lines=['# Quartic twisted classes and CJR aggregate identities','']+intro+[
    'For each entry with r >= 2g-2, the entire CJR infinity-containing aggregate is B_(D,r)=-Z_(D,r).',
    'The t power is 2g-2-r. Positive t powers are recorded as twisted classes only, not as CJR vanishing relations.','']
for g in (2,3):
    records=json.loads((out/('finite-classes-g%s.json'%g)).read_text())
    for d in (2,3):records+=json.loads((out/('finite-classes-degree%s-g%s.json'%(d,g))).read_text())
    for rec in records:
        if not rec['is_cjr_relation']:continue
        finite_lines+=['## Genus %s, ambient degree %s, codimension %s'%(g,rec['ambient_degree'],rec['codimension']),'',
                       r'\[ Z^{(%s)}_{%s,%s}=%s. \]'%(g,rec['ambient_degree'],rec['codimension'],
                       equation(normalized(rec['terms'],g),g)),'']

dictionary=['## Graph dictionary','',
    'Each row lists vertex genera, ordered flags at each vertex, paired node flags, and the number n of temporary markings.',
    'Flags labelled 1 through n are forgotten after gluing; all remaining flags are node branches.','']
for data,label in graph_labels.items():
    dictionary+=['- Theta_%s: `%s`.'%(label,json.dumps(data))]
(out/'independent-relations.md').write_text('\n'.join(lines)+'\n')
(out/'cjr-aggregate-classes.md').write_text('\n'.join(finite_lines+dictionary)+'\n')
(out/'coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
print('Replayed',certificate_count,'exact Pixton certificates and verified all saved basis ranks.')

# Check the two short formulas printed in the explanatory report separately.

from log_glsm.cycles.quartic import (
    QuarticK3GraphCycles,
)
from log_glsm.targets.quartic import (
    QuarticK3Cycles,
)
displayed=[]
for g,r,terms in (
        (2,1,{('interior',1,()):1,('interior',None,(1,)):-2,
              ('boundary',1,0,(),0,()):-QQ(1)/2}),
        (3,2,{('interior',2,()):1,('interior',1,(1,)):-1,
              ('interior',None,(2,)):4,('boundary',1,0,(),1,()):-1,
              ('boundary',1,0,(),0,(1,)):1})):
    b=QuarticK3Cycles(g)
    v=vector(QQ,b.class_from_terms(terms).vector(r))
    M=load(str(out/('pixton-g%s-r%s.sobj'%(g,r))))
    u=M.transpose().solve_right(v)
    assert u*M==v
    displayed.append({'genus':g,'codimension':r,
        'terms':[{'key':k,'coefficient':str(c)} for k,c in terms.items()],
        'strata_vector':[str(c) for c in v],
        'pixton_weights':[[int(i),str(c)] for i,c in enumerate(u) if c]})
(out/'displayed-relations.json').write_text(json.dumps(displayed,indent=2,default=int)+'\n')
print('Verified the genus-two Picard and genus-three codimension-two displayed formulas.')

# Give the genus-two t^0 classes in the small geometric basis (kappa_2,lambda_2).

b=QuarticK3GraphCycles(2)
M=load(str(out/'pixton-g2-r2.sobj'))
A=M.stack(matrix(QQ,[b.R.kappa(2).vector(2),b.R.lambdaclass(2).vector(2)]))
assert A.rank()==A.ncols()
records=json.loads((out/'finite-classes-g2.json').read_text())
for d in (2,3):records+=json.loads((out/('finite-classes-degree%s-g2.json'%d)).read_text())
small=[]
for record in records:
    if record['codimension']!=2:continue
    v=vector(QQ,M.ncols())
    for term in record['terms']:
        v+=QQ(term['coefficient'])*vector(QQ,b.monomial_class(freeze(term['key'])).vector(2))
    u=A.transpose().solve_right(v)
    assert u*A==v
    small.append({'degree':record['ambient_degree'],'kappa2':str(u[-2]),'lambda2':str(u[-1]),
                  'weights':[str(c) for c in u],'strata_vector':[str(c) for c in v]})
(out/'genus-two-constant-term.json').write_text(json.dumps(small,indent=2)+'\n')
print('Verified all four genus-two constant terms in the basis kappa_2, lambda_2.')

# Keep the full formal localization expressions, but avoid hundreds of MB of
# repeated graph keys and JSON indentation in the delivered research data.
for pattern in ('final-relations-*.json','combined-relations-*.json','basis-relations-*.json'):
    for path in out.glob(pattern):
        original=path.read_bytes()
        packed=gzip.compress(original,compresslevel=9,mtime=0)
        assert gzip.decompress(packed)==original
        Path(str(path)+'.gz').write_bytes(packed)
        path.unlink()
