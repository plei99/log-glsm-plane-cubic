"""Check the displayed genus-three consequences against the exact certificate."""
import json
import os
from pathlib import Path
from sage.all import *

out = Path(os.environ.get('GENUS3_OUTPUT',
    'results/genus3/corrected-2026-11-05/d1-extension'))
c = json.loads((out/'linear-certificate.json').read_text())
values = json.loads((out/'verified-values.json').read_text())

def key(v):
    return (v['genus'], v['ambient_degree'], tuple(v['contacts']),
            v['psi_min'], tuple(v['insertions']))

columns = {key(v): i for i, v in enumerate(c['columns'])}
known = {key(r['vertex']): QQ(r['value']) for r in values}
p = vector(QQ, c['particular_solution'])
K = matrix(QQ, len(columns), len(c['free_columns']), sparse=True)
for i, j, value in c['kernel_entries']:
    K[i, j] = QQ(value)

def check(terms, rhs=0):
    row = vector(QQ, len(columns))
    constant = QQ.zero()
    for coefficient, vertex in terms:
        coefficient = QQ(coefficient)
        if vertex in known:
            constant += coefficient*known[vertex]
        else:
            row[columns[vertex]] += coefficient
    assert row*p + constant == QQ(rhs)
    assert not row*K

A = (3, 0, (-5,), 0, (1,))
B = (3, 0, (-5,), 1, (0,))
check([(1, (3, 0, (-5, -1), 1, (1, 0))), (-5, A)])
check([(1, (3, 0, (-5, -1), 1, (0, 1))), (-5, A)])
check([(1, (3, 0, (-5, -1), 2, (0, 0))), (-5, B), (15, A)])
check([(1, (3, 0, (-5, -1), 0, (1, 1)))])
check([(1, (3, 0, (-5, -1), 0, (2, 0)))])
check([
    (1, (3, 0, (-4, -2), 2, (0, 0))),
    (6, (3, 0, (-4, -2), 1, (1, 0))),
    (15, (3, 0, (-4, -2), 0, (2, 0))),
    ('1/16', (3, 1, (-2,), 4, (0,))),
    ('-125/16', B), ('-375/8', A),
    ('-28102523/3408768', (2, 0, (-3,), 0, (1,))),
], '12021517/818104320')
if (2, 0, (-3,), 0, (1,)) in known:
    check([(1, (2, 0, (-3,), 0, (1,)))], '1/240')
    check([(1, (2, 0, (-3,), 1, (0,)))], '-1/80')
    check([
        (1, (3, 0, (-4, -2), 2, (0, 0))),
        (6, (3, 0, (-4, -2), 1, (1, 0))),
        (15, (3, 0, (-4, -2), 0, (2, 0))),
        ('1/16', (3, 1, (-2,), 4, (0,))),
        ('-125/16', B), ('-375/8', A),
    ], '113/2304')
print('All displayed genus-three logarithmic consequences verified over QQ.')

# Enumerate primary-contact candidates with no (-1) contact, modulo relabeling.
from itertools import product
profiles = [(0, (-5,)), (0, (-4, -2)), (0, (-3, -3)),
            (0, (-3, -2, -2)), (0, (-2, -2, -2, -2)), (1, (-2,))]
candidates = set()
for degree, contacts in profiles:
    for powers in product(range(3), repeat=len(contacts)):
        psi = len(contacts) + 3*degree - sum(powers)
        if psi >= 0:
            pairs = sorted(zip(contacts, powers))
            candidates.add((3, degree, tuple(x[0] for x in pairs),
                            psi, tuple(x[1] for x in pairs)))
present = sorted(v for v in candidates if v in columns)
missing = sorted(v for v in candidates if v not in columns and v not in known)
summary = {'primary_candidates_without_minus_one_contacts': len(candidates),
           'present_unknowns': len(present),
           'known_candidates': sum(v in known for v in candidates),
           'missing_candidates': missing,
           'projection_rank': int(K.matrix_from_rows([columns[v] for v in present]).rank())}
(out/'primary-candidate-analysis.json').write_text(json.dumps(summary, indent=2)+'\n')
print(summary)
if not missing:
    projected_kernel = K.matrix_from_rows([columns[v] for v in present])
    projected_particular = vector(QQ, [p[columns[v]] for v in present])
    constraints = projected_kernel.left_kernel().basis_matrix()
    assert not constraints*projected_kernel
    rhs = constraints*projected_particular
    primary = {'variables': [c['columns'][columns[v]] for v in present],
               'constraint_rank': int(constraints.nrows()),
               'equations': [{'coefficients': [str(x) for x in row],
                              'rhs': str(rhs[i])}
                             for i, row in enumerate(constraints.rows())]}
    (out/'primary-genus3-constraints.json').write_text(json.dumps(primary, indent=2)+'\n')
