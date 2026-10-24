#!/usr/bin/env python3
"""Replay all C57 lift choices in the exact exported Picard presentation.

The field/class-group construction and place classes are Magma inputs.
This independent replay checks their integer arithmetic, complete lift
coverage, nonvanishing, and invariance of the outcome under deck choice.
"""
import ast
from itertools import product
from pathlib import Path
import re
import sys

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "data" / (
    "C57_good7_class_arithmetic.txt")
raw = path.read_text().replace("[*", "[").replace("*]", "]")
invariants, relations, records = ast.literal_eval(raw)
assert invariants == [105, 210, 0]
h = hermite_normal_form(Matrix(relations).T)
assert h.cols == 2 and h.rows == 3


def is_zero(v):
    try:
        solution, parameters = h.gauss_jordan_solve(Matrix(v))
    except ValueError:
        return False
    assert not parameters.rows
    return all(x.q == 1 for x in solution)


def same(a, b):
    return is_zero([x - y for x, y in zip(a, b)])


assert len(records) == 2
keys, total = set(), 0
for index, key, p, sp, targets, recorded in records:
    assert index in (1, 2) and tuple(key) not in keys
    keys.add(tuple(key))
    assert len(targets) == 3 and sum(d * m for d, m, _, _ in targets) == 42
    assert all(d == 14 and m == 1 for d, m, _, _ in targets)
    assert p[2] == sp[2] == 7
    assert all(c[2] == sc[2] == d for d, _, c, sc in targets)
    # Magma's successive multi-index sequence comprehensions vary the
    # first target coordinate fastest; Python product varies it slowest.
    choices = [tuple(reversed(t)) for t in product(
        *(range(m + 1) for _, m, _, _ in reversed(targets)))]
    assert len(choices) == len(recorded) == 8
    for choice, observed in zip(choices, recorded):
        expected = [4 * (sp[j] - p[j]) + sum(
            (m - 2 * a) * (c[j] - sc[j])
            for a, (_, m, c, sc) in zip(choice, targets)) for j in range(3)]
        assert same(expected, observed), (index, choice, expected, observed)
        assert expected[2] == 0 and not is_zero(expected)
        # Interchanging source and all target lifts negates the divisor.
        other_deck_choice = [-v for v in expected]
        assert not is_zero(other_deck_choice)
        total += 1
assert total == 16
print({"level": 57, "picard_invariants": invariants,
       "source_cases": len(records), "operator_classes": total,
       "nonzero_classes": total, "all_lift_choices_covered": True,
       "independent_integer_HNF_replay": True})
