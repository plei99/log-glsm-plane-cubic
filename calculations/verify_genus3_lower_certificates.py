"""Replay saved rational Pixton certificates and check their strata coordinates."""
import json
from pathlib import Path
from sage.all import *


from log_glsm.cycles.plane import (
    O3UnmarkedCycleBackend,
)
backend = O3UnmarkedCycleBackend()
out = Path('results/genus3/lower-cycles-2026-11-05')


def tuple_key(value):
    return tuple(tuple_key(x) for x in value) if isinstance(value, list) else value


def check(record, relation_matrix, codimension):
    v = vector(QQ, record['strata_vector'])
    weights = vector(QQ, relation_matrix.nrows())
    for index, coefficient in record['pixton_weights']:
        weights[index] = QQ(coefficient)
    assert weights * relation_matrix == v
    terms = record.get('formal_terms', record.get('terms'))
    if terms is not None:
        cycle = sum((QQ(term['coefficient']) * backend.monomial_class(tuple_key(term['key']))
                     for term in terms), backend.R.zero())
        assert vector(QQ, cycle.vector(codimension)) == v


count = 0
displayed = json.loads((out / 'displayed-relations.json').read_text())
for r in (3, 4, 5):
    matrix_path = out / ('pixton-codimension-%s.sobj' % r)
    relation_matrix = load(str(matrix_path))
    strata = json.loads((out / ('strata-order-codimension-%s.json' % r)).read_text())
    assert [s['stratum'] for s in strata] == [str(s) for s in backend.R.generators(r)]
    for prefix in ('pixton-certificates', 'path-certificates'):
        for record in json.loads((out / ('%s-codimension-%s.json' % (prefix, r))).read_text()):
            check(record, relation_matrix, r)
            count += 1
    for record in displayed:
        if record['codimension'] == r:
            check(record, relation_matrix, r)
            count += 1
assert count == 19
print('Verified 19 exact certificates, including the formal pole and displayed equations.')
