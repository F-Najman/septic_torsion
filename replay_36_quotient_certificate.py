"""Independent integer-lattice/combinatorial replay of X1(36)/F5.

Magma supplies exact first-support diamond representatives and class maps;
this script verifies every necessary
finite-group test with HNF over Z, without using Magma's subgroup or
quotient operations. It also checks place counts from Newton sums.
"""
import ast
from collections import defaultdict
from itertools import combinations_with_replacement
from math import lcm, prod
from pathlib import Path
from sympy import Matrix, diag
from sympy.matrices.normalforms import hermite_normal_form

source = Path(__file__).resolve().parent / "data/certificate_36_quotient_data.py"
tree = ast.parse(source.read_text())
data = {}
for statement in tree.body:
    assert isinstance(statement, ast.Assign) and len(statement.targets) == 1
    assert isinstance(statement.targets[0], ast.Name)
    data[statement.targets[0].id] = ast.literal_eval(statement.value)
orders = tuple(data["invariants"])
rank = len(orders)
zero = (0,) * rank

def norm(vector):
    assert len(vector) == rank
    return tuple(int(a) % n for a, n in zip(vector, orders))

def add(a, b):
    return tuple((x + y) % n for x, y, n in zip(a, b, orders))

cc = [norm(x) for x in data["cusp_classes"]]
qc = [norm(x) for x in data["rational_cusp_classes"]]
points = [[norm(x) for x in part] for part in data["low_representative_classes"]]
smallpoints = [[norm(x) for x in part] for part in data["small_point_classes"]]
assert [len(x) for x in points] == [0,0,8,24,110]
assert data["representative_counts"] == [0,0,8,24,110,480,1788]
assert [len(x) for x in smallpoints] == [48,144]
relations = diag(*orders).row_join(Matrix(qc).T)
lattice = hermite_normal_form(relations)
assert lattice.shape == (rank, rank)
assert prod(orders) == 5061852020736
assert prod(orders) // int(lattice.det()) == 12 * 252 * 7812
inverse = lattice.inv()
constraints = []
for row in inverse.tolist():
    denominator = lcm(*(x.q for x in row))
    if denominator > 1:
        constraints.append((tuple(int(x * denominator) % denominator for x in row), denominator))

def coset(a):
    return tuple(sum(x * y for x, y in zip(row, a)) % modulus
                 for row, modulus in constraints)

identity = coset(zero)
assert all(coset(x) == identity for x in qc)

def effective(weights, classes, degree):
    result = set()
    number = 0
    def visit(remaining, start, current):
        nonlocal number
        if remaining == 0:
            result.add(current)
            number += 1
            return
        for i in range(start, len(weights)):
            if weights[i] <= remaining:
                visit(remaining - weights[i], i, add(current, classes[i]))
    visit(degree, 0, zero)
    return result, number

qeffective, qnumber = effective(data["rational_cusp_degrees"], qc, 7)
assert qnumber == len(qeffective) == 19159
allcusps, number = effective(data["cusp_degrees"], cc, 7)
assert number == 188865
candidates = {a for a in allcusps if coset(a) == identity and a not in qeffective}
print("C36_QUOTIENT_PYTHON_ALL_CUSPIDAL", number, len(allcusps), len(candidates), flush=True)
complements = {}
for k in range(5):
    out, count = effective(data["cusp_degrees"], cc, k)
    bycoset = defaultdict(set)
    for a in out:
        bycoset[coset(a)].add(a)
    complements[k] = bycoset

def test_parts(parts, degree):
    passes = effective_count = total = 0
    for a in parts:
        total += 1
        negative_key = coset(tuple(-x for x in a))
        for b in complements[7 - degree].get(negative_key, ()):
            gc = add(a, b)
            assert coset(gc) == identity
            passes += 1
            if gc in qeffective:
                effective_count += 1
            else:
                candidates.add(gc)
    return total, passes, effective_count, len(candidates)

for d in range(3, 6):
    print("C36_QUOTIENT_PYTHON_TYPE", [d], test_parts(points[d - 1], d), flush=True)
print("C36_QUOTIENT_PYTHON_TYPE", [3, 3],
      test_parts((add(a, b) for a in points[2] for b in smallpoints[0]), 6), flush=True)
print("C36_QUOTIENT_PYTHON_TYPE", [3, 4],
      test_parts((add(a, b) for a in points[2] for b in smallpoints[1]), 7), flush=True)
# Independently reconstruct the genus7 rational-cusp subgroup over Z.
yorders = tuple(data["target_invariants"])
assert yorders == (6,84,252)
ygens = data["target_cusp_classes"]
ylattice = hermite_normal_form(diag(*yorders).row_join(Matrix(ygens).T))
assert prod(yorders) // int(ylattice.det()) == 2*42*126
yinverse = ylattice.inv()
yconstraints = []
for row in yinverse.tolist():
    den = lcm(*(x.q for x in row))
    if den > 1:
        yconstraints.append((tuple(int(x*den)%den for x in row),den))
def ymember(a):
    assert len(a) == 3
    return all(sum(x*y for x,y in zip(row,a))%den == 0 for row,den in yconstraints)
assert all(ymember(a) for a in ygens)
records = list(zip(data["high_degrees"],data["high_indices"],
                   data["high_target_classes"],data["high_source_classes"],strict=True))
assert len(records) == 2268
seen = set()
passed = {6:[],7:[]}
for degree,index,yclass,xclass in records:
    key = (degree,index)
    assert key not in seen
    seen.add(key)
    assert degree in (6,7) and 1 <= index <= data["representative_counts"][degree-1]
    assert bool(xclass) == ymember(yclass)
    if xclass:
        passed[degree].append(norm(xclass))
assert seen == {(d,i) for d in (6,7) for i in range(1,data["representative_counts"][d-1]+1)}
for d in (6,7):
    print("C36_QUOTIENT_PYTHON_HIGH",[d],"y_passes",len(passed[d]),
          "fine",test_parts(passed[d],d),flush=True)
assert candidates == {norm(a) for a in data["magma_candidates"]}
known = {norm(a) for a in data["known_classes"]}
assert len(known) == 12
unknown = candidates - known

fp = data["frobenius_coefficients"]
assert len(fp) == 35 and fp[-1] == 1
traces, closed = [], []
for n in range(1, 8):
    trace = -n * fp[34 - n]
    trace -= sum(fp[34 - j] * traces[n - j - 1] for j in range(1, n))
    traces.append(trace)
    count = 5 ** n + 1 - trace
    count -= sum(d * closed[d - 1] for d in range(1, n) if n % d == 0)
    assert count % n == 0
    closed.append(count // n)
    assert closed[-1] == data["point_counts"][n - 1] + data["cusp_degrees"].count(n)
assert data["point_counts"] == [0,0,48,144,660,2856,10728]
print("C36_QUOTIENT_PYTHON_FROBENIUS_COUNTS", data["point_counts"], closed, flush=True)
print("C36_QUOTIENT_PYTHON_REPLAY_DONE candidates", len(candidates), "known", len(known),
      "unknown", len(unknown), "complete=true", flush=True)
assert not unknown
