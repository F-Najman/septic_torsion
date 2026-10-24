#!/usr/bin/env python3
"""Independently replay the exported degree-seven X1(27) class sieve.

This checks integer group arithmetic, not the geometric ClassGroup maps.
The geometric enumeration and all source fibre patterns are separate inputs.
"""
import ast
from collections import Counter
import json
from pathlib import Path
import sys

from sympy import Matrix, diag
from sympy.matrices.normalforms import hermite_normal_form

source = Path(sys.argv[1] if len(sys.argv) > 1 else
              'data/N54_class_arithmetic_data.py')
data = {}
for statement in ast.parse(source.read_text()).body:
    assert isinstance(statement, ast.Assign) and len(statement.targets) == 1
    data[statement.targets[0].id] = ast.literal_eval(statement.value)
assert data['complete'] is True


def presentation(name):
    rows = data['relations' + name]
    assert rows and len(rows[0]) == len(rows) + 1
    assert all(len(row) == len(rows[0]) and row[-1] == 0 for row in rows)
    finite = Matrix([row[:-1] for row in rows])
    assert finite.det() != 0
    return finite


R2, R5 = presentation('2'), presentation('5')
n2, n5 = R2.cols, R5.cols


def finite(vector, n):
    assert len(vector) == n + 1 and vector[-1] == 0
    return vector[:-1]


C2 = [finite(v, n2) for v in data['cusp2']]
C5 = [finite(v, n5) for v in data['cusp5']]
assert len(C2) == len(C5)
M5 = hermite_normal_form(Matrix(R5.tolist() + C5).T)
assert M5.shape == (n5, n5)
I5 = M5.inv()
joint_relations = diag(R2, R5).tolist()
joint_gens = [a+b for a, b in zip(C2, C5)]
MJ = hermite_normal_form(Matrix(joint_relations + joint_gens).T)
assert MJ.shape == (n2+n5, n2+n5)
IJ = MJ.inv()


def member(inverse, vector):
    return all(entry.q == 1 for entry in inverse * Matrix(vector))


def residue_key(vector):
    # A canonical class representative in Z^n / row-lattice(R2).
    coordinates = R2.T.inv() * Matrix(vector)
    return tuple(entry % 1 for entry in coordinates)


candidates2 = {residue_key(finite(v, n2)) for v in data['candidate2']}
assert len(candidates2) == 90
members = passing = 0
by_degree = Counter()
for degree, q5, was_member, was_pass, im2 in data['tests']:
    assert degree in range(3, 8)
    q5 = finite(q5, n5)
    is_member = member(I5, q5)
    assert int(is_member) == was_member
    is_pass = False
    if is_member:
        members += 1
        im2 = finite(im2, n2)
        assert member(IJ, im2+q5), 'Incompatible reduction map'
        is_pass = residue_key(im2) in candidates2
    else:
        assert im2 == []
    assert int(is_pass) == was_pass
    passing += is_pass
    by_degree[degree] += 1

assert by_degree[3] == 4*636
assert by_degree[4] == 19*192
assert by_degree[5] == 44*48
assert by_degree[7] == 828+4*171
assert by_degree[6] >= 212*9+4*30*9
assert (by_degree[6]-212*9-4*30*9) % 9 == 0
print(json.dumps({
    'tests': len(data['tests']), 'by_noncuspidal_degree': dict(by_degree),
    'members_at_5': members, 'passing_both_primes': passing,
    'rational_cusp_image_index_at_5': int(abs(M5.det())),
    'inert_quadratic_diamond_representatives':
        (by_degree[6]-212*9-4*30*9)//9,
    'complete': True,
}, indent=2))
assert passing == 0
print('N54_INTEGER_CLASS_REPLAY_DONE complete=true passing=0')
