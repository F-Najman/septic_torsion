#!/usr/bin/env python3
"""Independent, order-free reconstruction of the sixteen C57 lift tests."""
import ast
from collections import Counter
from itertools import product
import json
from pathlib import Path

p = Path(__file__).resolve().parent / 'data/C57_good7_class_arithmetic.txt'
inv, relations, cases = ast.literal_eval(
    p.read_text().replace('[*', '[').replace('*]', ']'))
assert inv == [105, 210, 0]
# Check the actual coordinate relations, rather than assuming invariant
# factors describe the exported coordinates.
assert relations == [[105, 0, 0], [0, 210, 0]]
assert len(cases) == 2


def canonical(v):
    assert len(v) == 3
    return v[0] % 105, v[1] % 210, v[2]


total = 0
for index, key, point, deck_point, targets, observed in cases:
    assert index in (1, 2) and len(targets) == 3
    assert point[2] == deck_point[2] == 7
    assert all(d == 14 and m == 1 and a[2] == b[2] == 14
               for d, m, a, b in targets)
    reconstructed = []
    for choices in product((0, 1), repeat=3):
        image = [0, 0, 0]
        deck_image = [0, 0, 0]
        for choose, (_, _, a, b) in zip(choices, targets):
            selected, other = (a, b) if choose else (b, a)
            image = [x+y for x, y in zip(image, selected)]
            deck_image = [x+y for x, y in zip(deck_image, other)]
        assert image[2] == deck_image[2] == 42
        v = [deck_image[j]-image[j]+4*(deck_point[j]-point[j])
             for j in range(3)]
        assert v[2] == 0
        c = canonical(v)
        assert c != (0, 0, 0)
        reconstructed.append(c)
    assert len(observed) == 8
    assert Counter(reconstructed) == Counter(map(canonical, observed))
    total += len(reconstructed)
assert total == 16
print(json.dumps({'level': 57, 'source_cases': 2, 'tests': total,
                  'nonzero': total, 'actual_coordinate_relations_checked': True,
                  'all_lift_multisets_reconstructed': True, 'complete': True}))
print('C57_ROOT_LIFT_MULTISET_REPLAY_DONE tests=16 nonzero=16 complete=true')
