#!/usr/bin/env python3
"""Replay every paired-lift C57 operator in the actual Picard presentation."""
import ast
from collections import Counter
import hashlib
import itertools
from pathlib import Path
from sympy import Matrix

HERE = Path(__file__).resolve().parent
BASE = HERE / "data"
path = BASE / "C57_good7_class_arithmetic.txt"
raw = path.read_text().replace("[*", "[").replace("*]", "]")
invariants, relations, cases = ast.literal_eval(raw)
assert invariants == [105,210,0] and len(cases) == 2
assert len(relations) == 2 and all(len(r) == 3 and r[2] == 0 for r in relations)
R = Matrix([r[:-1] for r in relations])
assert abs(R.det()) == 105*210
RI = R.T.inv()

def signature(v):
    assert len(v) == 3 and v[2] == 0
    return tuple(x % 1 for x in RI*Matrix(v[:-1]))

def add(a,b):
    return [x+y for x,y in zip(a,b)]

total = 0
for expected_index, record in enumerate(cases,1):
    index, key, P, sP, targets, exported = record
    assert index == expected_index and len(key) == 21
    assert len(P) == len(sP) == 3 and P[-1] == sP[-1] == 7
    assert len(targets) == 3
    assert sum(degree*mult for degree,mult,U,V in targets) == 42
    for degree,mult,U,V in targets:
        assert degree == 14 and mult == 1 and U[-1] == V[-1] == degree
    expected = []
    # All possible choices of target lifts. Order is immaterial; use a
    # multiset comparison, retaining duplicated classes and multiplicities.
    for choices in itertools.product(*(range(t[1]+1) for t in targets)):
        value = [4*(x-y) for x,y in zip(sP,P)]
        for a,(degree,mult,U,V) in zip(choices,targets):
            # T=aU+(m-a)V, sT=aV+(m-a)U.
            value = add(value,[(mult-2*a)*(u-v) for u,v in zip(U,V)])
        expected.append(signature(value))
    got = [signature(v) for v in exported]
    assert len(expected) == len(got) == 8
    assert Counter(expected) == Counter(got)
    assert all(any(x != 0 for x in sig) for sig in expected)
    print(f"C57_INDEPENDENT_CLASS_CASE index={index} choices=8 nonzero=8")
    total += len(expected)
assert total == 16
print("C57_INDEPENDENT_CLASS_REPLAY_DONE actual_relations=true cases=2 tests=16 zero_classes=0 complete=true")
print("input_sha256="+hashlib.sha256(path.read_bytes()).hexdigest())
