#!/usr/bin/env python3
"""Replay exact finite-group arithmetic and coverage for the N26/N28 sieve.

This checks exported relations, rather than identifying Magma group
coordinates with invariant-factor coordinates. Geometric enumeration and
class-map correctness remain the independently audited Magma stage.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent
EXPECTED = {
    26: {"reps": [6, 26, 48], "points": [36, 122, 288],
         "H": [2, 14, 266, 3990, 11970]},
    28: {"reps": [8, 16, 52], "points": [48, 64, 312],
         "H": [2, 4, 4, 4, 8, 24, 936, 936]},
}


def read_one(path):
    s = path.read_text().replace("[*", "[").replace("*]", "]")
    # json.loads rejects trailing expressions, unlike a first-record reader.
    return json.loads(s)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int, choices=EXPECTED)
    args = parser.parse_args()
    N = args.N
    expected = EXPECTED[N]
    raw = HERE / f"data/N{N}_D7_full_class_arithmetic.txt"
    logpath = HERE / f"logs/m2/X1_2_{N}_classes.log"
    inv, relations, hgens, records = read_one(raw)
    assert inv[-1] == 0 and all(a > 0 for a in inv[:-1])
    n = len(inv) - 1
    vectors = relations + hgens + [r[2] for r in records]
    assert all(len(v) == n + 1 and v[-1] == 0 for v in vectors)
    R = hermite_normal_form(Matrix([v[:-1] for v in relations]).T)
    H = hermite_normal_form(Matrix([v[:-1] for v in relations + hgens]).T)
    assert R.shape == H.shape == (n, n)
    gorder = abs(int(R.det()))
    hindex = abs(int(H.det()))
    horder = gorder // hindex
    assert gorder == math.prod(inv[:-1])
    assert gorder % hindex == 0
    assert horder == math.prod(expected["H"])
    assert gorder == {26: 124501550040000, 28: 19291718025216}[N]
    hb = [[int(H[i, j]) for j in range(n)] for i in range(n)]
    assert all(hb[i][j] == 0 for i in range(n) for j in range(i))

    def key(v):
        v = list(v[:-1])
        for i in range(n - 1, -1, -1):
            q = v[i] // hb[i][i]
            for j in range(i + 1):
                v[j] -= q * hb[j][i]
        return tuple(v)

    zero = (0,) * n
    assert all(key(v) == zero for v in relations + hgens)
    wanted = [(d, i) for d, count in zip(range(5, 8), expected["reps"])
              for i in range(1, count + 1)]
    assert [(d, i) for d, i, _, inside in records] == wanted
    remaining = []
    for d, i, vector, inside in records:
        assert isinstance(inside, bool)
        assert (key(vector) == zero) == inside
        if inside:
            remaining.append([d, i])
    log = logpath.read_text()
    assert not re.search(r"error|interrupt|segmentation|failed", log, re.I)
    marker = (f"M2_D7_FULL_CLASS_DONE N={N} tested={len(wanted)} "
              f"remaining={len(remaining)} complete=true")
    assert marker in log
    for d, count in zip(range(5, 8), expected["reps"]):
        matches = re.findall(rf"M2_D7_FULL_CLASS_TYPE\s+{N}\s+{d}\s+"
                             rf"{count}\s+(\d+)\s+", log)
        assert len(matches) == 1
        assert int(matches[0]) == sum(a == d for a, _ in remaining)
    result = {
        "status": "passed", "N": N, "p": 3,
        "place_counts_degrees_5_6_7": expected["points"],
        "diamond_representatives_degrees_5_6_7": expected["reps"],
        "tested": len(wanted), "finite_J_order": gorder,
        "rational_cuspidal_image_order": horder,
        "image_index": hindex, "remaining": remaining,
        "arithmetic_sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(
            (HERE / "class_mixed_degree7_full.m").read_bytes()).hexdigest(),
    }
    (HERE / f"data/N{N}_D7_full_class_arithmetic.json").write_text(
        json.dumps([inv, relations, hgens, records], indent=2) + "\n")
    (HERE / f"data/N{N}_D7_full_class_replay.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
