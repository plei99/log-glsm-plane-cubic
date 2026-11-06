"""Trace the six-term genus-three equation to unmarked localization probes.

All rows have the determined values from the completed rerun substituted.
The last coordinate is the right-hand side, so checking row equality also
checks the numerical constant. Punctured axioms are retained explicitly.
"""
import json
from pathlib import Path
from sage.all import *
p=Path('results/genus3/corrected-2026-11-05/d1-extension')
c=json.loads((p/'linear-certificate.json').read_text())
n=len(c['columns']);m=len(c['rhs'])
A=matrix(QQ,m,n,sparse=True)
for i,j,x in c['matrix_entries']:A[i,j]=QQ(x)
b=vector(QQ,c['rhs'])
def key(v):return(v['genus'],v['ambient_degree'],tuple(v['contacts']),v['psi_min'],tuple(v['insertions']))
index={key(v):i for i,v in enumerate(c['columns'])}
target=vector(QQ,n+1)
for coeff,v in [(16,(3,0,(-4,-2),2,(0,0))),(96,(3,0,(-4,-2),1,(1,0))),
 (240,(3,0,(-4,-2),0,(2,0))),(1,(3,1,(-2,),4,(0,))),
 (-125,(3,0,(-5,),1,(0,))),(-750,(3,0,(-5,),0,(1,)))]:target[index[v]]=coeff
target[-1]=QQ(113)/144
M=A.augment(matrix(QQ,m,1,b))
axioms=[i for i,r in enumerate(c['row_probes']) if r['label'].startswith('CJR ') and r['label']!='CJR Chow probe']
print('axioms',len(axioms),flush=True)
R=M.matrix_from_rows(axioms).dense_matrix().rref()
piv=R.pivots()
def reduce(v):
 v=vector(QQ,v)
 for i,p in enumerate(piv):
  if v[p]:v-=v[p]*R.row(i)
 return v
t=reduce(target)
rows=[i for i in range(m) if i not in set(axioms)]
rows.sort(key=lambda i:(c['row_probes'][i]['ambient_degree'],len(c['row_probes'][i]['insertions']),c['row_probes'][i]['label']!='CJR Chow probe'))
B=matrix(QQ,[reduce(M.row(i)) for i in rows])
print('residual rank',B.rank(),'rows',len(rows),flush=True)
for degree in (0,1,2):
 selected=[j for j,i in enumerate(rows) if c['row_probes'][i]['ambient_degree']<=degree]
 Bsub=B.matrix_from_rows(selected)
 try:w=Bsub.transpose().solve_right(t)
 except ValueError:
  print('not in degree',degree,flush=True);continue
 print('SOURCE DEGREE',degree,flush=True)
 for j,coeff in enumerate(w):
  if coeff:
   i=rows[selected[j]]
   print('WEIGHT',coeff,'t',c['row_t_powers'][i],'PROBE',c['row_probes'][i],flush=True)
 assert w*Bsub==t
 source_weights=[]
 combination=vector(QQ,n+1)
 for j,coeff in enumerate(w):
  if coeff:
   i=rows[selected[j]]
   source_weights.append({'certificate_row':int(i),'weight':str(coeff),
      'probe':c['row_probes'][i],'t_power':c['row_t_powers'][i]})
   combination+=coeff*M.row(i)
 ax_weights=M.matrix_from_rows(axioms).transpose().solve_right(target-combination)
 assert ax_weights*M.matrix_from_rows(axioms)+combination==target
 certificate={'statement':'boxed equation = -24 R(g=3,D=0,t^-4) - R(g=3,D=1,t^-4), modulo punctured axioms after substituting determined values',
   'localization_rows':source_weights,
   'axiom_rows':[{'certificate_row':int(axioms[j]),'weight':str(coeff),
                  'probe':c['row_probes'][axioms[j]]}
                 for j,coeff in enumerate(ax_weights) if coeff],
   'target_coefficients':[str(x) for x in target[:-1]],
   'target_rhs':str(target[-1]),'all_coordinates_verified':True}
 (p/'boxed-equation-origin.json').write_text(json.dumps(certificate,indent=2)+'\n')
 print('Verified original-row combination with',len(certificate['axiom_rows']),'punctured axiom rows.',flush=True)
 break
