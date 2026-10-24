#!/usr/bin/env python3
"""Independent exact row-lattice replay of the N54 two-prime certificate."""
import ast
from collections import Counter
from pathlib import Path
from sympy import Matrix, diag
from sympy.matrices.normalforms import hermite_normal_form

BASE = Path(__file__).resolve().parent / "data"
data = {}
for s in ast.parse((BASE / "N54_class_arithmetic_data.py").read_text()).body:
    assert isinstance(s, ast.Assign) and len(s.targets) == 1
    data[s.targets[0].id] = ast.literal_eval(s.value)
assert data["complete"] is True

def relation(name):
    raw = data["relations"+name]
    n = len(raw)
    assert all(len(v) == n+1 and v[-1] == 0 for v in raw)
    R = Matrix([v[:-1] for v in raw])
    assert R.det() != 0
    return R

R2, R5 = relation("2"), relation("5")
n2, n5 = R2.cols, R5.cols

def finite(v, n):
    assert len(v) == n+1 and v[-1] == 0
    return v[:-1]

C2 = [finite(v, n2) for v in data["cusp2"]]
C5 = [finite(v, n5) for v in data["cusp5"]]
assert len(C2) == len(C5) == 14
L2 = hermite_normal_form(Matrix(R2.tolist()+C2).T)
L5 = hermite_normal_form(Matrix(R5.tolist()+C5).T)
LJ = hermite_normal_form(Matrix(diag(R2, R5).tolist()+[a+b for a,b in zip(C2,C5)]).T)
assert abs(L2.det()) == 1                     # H2 equals the whole finite Pic0.
order2 = abs(R2.det())
order5 = abs(R5.det())
orderH5 = order5 / abs(L5.det())
orderHJ = order2*order5 / abs(LJ.det())
assert order2 == orderH5 == orderHJ == 3*3*52497
# HJ -> H5 is onto by construction, and equal orders make it injective.
# Thus these actual cusp vectors define a UNIQUE cross-prime map.

I5, IJ = L5.inv(), LJ.inv()

def member(I,v):
    return all(x.q == 1 for x in I*Matrix(v))

R2inv = R2.T.inv()
def key(v):
    return tuple(x % 1 for x in R2inv*Matrix(v))

candidates = {key(finite(v,n2)) for v in data["candidate2"]}
assert len(candidates) == 90
counts = Counter()
members = passing = 0
for degree,v,reported_member,reported_pass,im2 in data["tests"]:
    counts[degree] += 1
    v = finite(v,n5)
    inH = member(I5,v)
    assert bool(reported_member) == inH
    passed = False
    if inH:
        members += 1
        w = finite(im2,n2)
        assert member(IJ,w+v)
        passed = key(w) in candidates
    else:
        assert im2 == []
    assert bool(reported_pass) == passed
    passing += passed
assert counts == {3:2544,4:3648,5:2112,6:2997,7:1512}
assert sum(counts.values()) == 12813 and members == 28 and passing == 0
print("N54_SMALL_CYCLIC_INTEGER_AUDIT_DONE actual_relations=true cross_prime_map_unique=true tests=12813 members=28 passing=0 complete=true")
