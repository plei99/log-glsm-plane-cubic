"""Independently row-reduce every saved scalar equation over QQ.

Unlike selecting a square subsystem, this checks the augmented matrix of
all rows, including dependent rows, for consistency. Zero products are
removed before deciding whether a relation is linear.
"""
import json,os
from pathlib import Path
from sage.all import *

from log_glsm.infinity.probes import (
    CompiledLocalizationRelation,
)
from log_glsm.infinity.solver import (
    EffectiveVertex,
)

out=Path(os.environ.get('GENUS3_OUTPUT','results/genus3/corrected-2026-11-05/d1-extension'))
checkpoint=json.loads((out/'checkpoint.json').read_text())
assert checkpoint['version']==10
relations=[CompiledLocalizationRelation.from_record(r,QQ) for r in checkpoint['relations']]
values={EffectiveVertex.from_record(r['vertex']):QQ(r['value'])
        for r in checkpoint['solver']['values']}
original_count=len(values)
roots=[EffectiveVertex.from_record(v) for v in checkpoint['roots']]


def reduced_relation(relation):
    rhs=QQ(relation.known_gw-relation.twisted_zero_level)
    terms={}
    for coefficient,factors in relation.terms:
        coefficient=QQ(coefficient)
        unknown=[]
        for factor in factors:
            if factor in values:
                coefficient*=values[factor]
            else:
                unknown.append(factor)
        if not coefficient:
            continue
        key=tuple(sorted(unknown,key=lambda v:v.order_key()))
        if not key:
            rhs-=coefficient
        else:
            terms[key]=terms.get(key,QQ.zero())+coefficient
    return {k:v for k,v in terms.items() if v},rhs


for iteration in range(20):
    rows=[];nonlinear=[];closed=0
    for relation in relations:
        terms,rhs=reduced_relation(relation)
        if not terms:
            assert rhs==0, ('inconsistent constant relation',relation.to_record(),str(rhs))
            closed+=1
        elif all(len(k)==1 for k in terms):
            rows.append((relation,terms,rhs))
        else:
            nonlinear.append((relation,terms,rhs))
    columns=tuple(sorted({k[0] for _,terms,_ in rows for k in terms},key=lambda v:v.order_key()))
    index={v:i for i,v in enumerate(columns)}
    A=matrix(QQ,len(rows),len(columns),sparse=True)
    b=vector(QQ,[rhs for _,_,rhs in rows])
    for i,(_,terms,_) in enumerate(rows):
        for (v,),coefficient in terms.items():
            A[i,index[v]]=coefficient
    R=A.augment(matrix(QQ,len(rows),1,b)).dense_matrix().rref()
    pivots=R.pivots()
    assert len(columns) not in pivots, 'the full augmented system is inconsistent'
    free=tuple(i for i in range(len(columns)) if i not in set(pivots))
    new={columns[p]:R[i,-1] for i,p in enumerate(pivots)
         if all(not R[i,j] for j in free)}
    print('iteration',iteration,'rows',len(rows),'columns',len(columns),'rank',len(pivots),
          'nullity',len(free),'new values',len(new),'nonlinear',len(nonlinear),flush=True)
    if not new:
        break
    values.update(new)
else:
    raise RuntimeError('linear substitution did not stabilize')

# Check an explicit affine parameterization against every linear row.
particular=vector(QQ,len(columns))
K=matrix(QQ,len(columns),len(free))
for j,f in enumerate(free):
    K[f,j]=1
for i,p in enumerate(pivots):
    particular[p]=R[i,-1]
    for j,f in enumerate(free):
        K[p,j]=-R[i,f]
assert A*particular==b
assert not A*K


def affine_record(v):
    if v in values:
        return {'vertex':v.to_record(),'constant':str(values[v]),'terms':[],'determined':True}
    if v not in index:
        return {'vertex':v.to_record(),'determined':False,'outside_linear_system':True}
    i=index[v]
    return {'vertex':v.to_record(),'constant':str(particular[i]),'determined':False,
            'terms':[{'coefficient':str(K[i,j]),'free_vertex':columns[f].to_record()}
                     for j,f in enumerate(free) if K[i,j]]}

report={'checkpoint_version':checkpoint['version'],'relations':len(relations),
        'closed_relations':closed,'linear_rows':len(rows),'columns':len(columns),
        'rank':len(pivots),'nullity':len(free),'nonlinear_rows':len(nonlinear),
        'all_linear_rows_verified':True,'additional_solved_values':len(values)-original_count,
        'total_solved_values':len(values),'roots':[affine_record(v) for v in roots],
        'genus_three_affine_relations':[affine_record(v) for v in columns if v.genus==3]}
(out/'linear-analysis.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
(out/'verified-values.json').write_text(json.dumps([
    {'vertex':v.to_record(),'value':str(value)}
    for v,value in sorted(values.items(),key=lambda item:item[0].order_key())],indent=2)+'\n')
# A portable exact certificate, with the meanings of every row and column.
certificate={'columns':[v.to_record() for v in columns],
 'row_probes':[r.probe.to_record() for r,_,_ in rows],
 'row_t_powers':[int(r.t_power) for r,_,_ in rows],
 'matrix_entries':[[int(i),int(j),str(value)] for (i,j),value in A.dict().items()],
 'rhs':[str(x) for x in b], 'particular_solution':[str(x) for x in particular],
 'free_columns':[int(i) for i in free],
 'kernel_entries':[[int(i),int(j),str(K[i,j])] for i in range(K.nrows())
                   for j in range(K.ncols()) if K[i,j]]}
(out/'linear-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
print('ROOTS',json.dumps(report['roots'],indent=2),flush=True)

# Eliminate larger-contact vertices first to expose shorter consequences.
order=sorted(range(len(columns)),key=lambda i:(-len(columns[i].contacts),
    -columns[i].genus,-columns[i].ambient_degree,-columns[i].psi_min))
reverse_reduced=A.matrix_from_columns(order).augment(
    matrix(QQ,len(rows),1,b)).dense_matrix().rref()
short=[]
for row in reverse_reduced.rows():
    terms=[{'vertex':columns[order[j]].to_record(),'coefficient':str(row[j])}
           for j in range(len(columns)) if row[j]]
    if row[-1] and any(x['vertex']['genus']==3 for x in terms):
        short.append({'rhs':str(row[-1]),'terms':terms})
short.sort(key=lambda x:len(x['terms']))
(out/'short-genus3-relations.json').write_text(json.dumps(short[:10],indent=2)+'\n')
root_rows=[index[v] for v in roots if v in index]
report['one_contact_freedom_rank']=int(K.matrix_from_rows(root_rows).rank())
(out/'linear-analysis.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print('Shortest genus-three consequence has',len(short[0]['terms']) if short else 0,'terms')

# Eliminate lower-genus unknowns as well, to give relations solely in g=3.
order=sorted(range(len(columns)),key=lambda i:(columns[i].genus==3,
    -len(columns[i].contacts),-columns[i].ambient_degree,-columns[i].psi_min))
R3=A.matrix_from_columns(order).augment(matrix(QQ,len(rows),1,b)).dense_matrix().rref()
pure=[]
for row in R3.rows():
    terms=[{'vertex':columns[order[j]].to_record(),'coefficient':str(row[j])}
           for j in range(len(columns)) if row[j]]
    if row[-1] and terms and all(x['vertex']['genus']==3 for x in terms):
        pure.append({'rhs':str(row[-1]),'terms':terms})
pure.sort(key=lambda x:len(x['terms']))
(out/'pure-genus3-relations.json').write_text(json.dumps(pure[:10],indent=2)+'\n')
print('Shortest pure genus-three inhomogeneous relation has',len(pure[0]['terms']) if pure else 0,'terms')
