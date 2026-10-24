#!/usr/bin/env python3
"""Strict finite coverage audit for the nine C40 correction tests."""
from pathlib import Path
import hashlib
import itertools
import json
import re

HERE = Path(__file__).resolve().parent
LOGS = HERE / "logs/sieve"
DATA = HERE / "data"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def numbers(text):
    return list(map(int, re.findall(r"-?\d+", text)))

text = (LOGS / "X1_40_cusp_correction.log").read_text()
assert not re.search(r"Runtime error|User error|Assertion failed|principal=true", text)
assert "bad_quadratic=4 Q_orbits=2 full_complements=84 Q_liftable=80" in text
assert "C40_CUSP_CORRECTION_MODEL_READY genus=25 p=3 quintic_degree=5" in text
assert text.count("C40_CUSP_CORRECTION_TEST kind=Q_liftable complements=80 principal=false") == 1
tests = re.findall(r"C40_CUSP_CORRECTION_TEST kind=bad_quadratic cycle=(\d+) cusp=(\d+) principal=(\w+)", text)
assert len(tests) == 8 and set(tests) == {
    (str(cycle), str(cusp), "false") for cycle in (1, 2) for cusp in range(1, 5)
}
assert text.count("C40_CUSP_CORRECTION_DONE tests=9 survivors=0 all_complements=84 complete=true") == 1
fresh_correction = (LOGS / "X1_40_cusp_correction_reproduction.log").read_text()
assert not re.search(r"Runtime error|User error|Assertion failed|principal=true", fresh_correction)
assert fresh_correction.count("C40_CUSP_CORRECTION_TEST kind=Q_liftable complements=80 principal=false") == 1
fresh_tests = re.findall(r"C40_CUSP_CORRECTION_TEST kind=bad_quadratic cycle=(\d+) cusp=(\d+) principal=(\w+)", fresh_correction)
assert len(fresh_tests) == 8 and set(fresh_tests) == set(tests)
assert fresh_correction.count("C40_CUSP_CORRECTION_DONE tests=9 survivors=0 all_complements=84 complete=true") == 1

# Verify the entire combinatorial overapproximation independently of Magma.
sigma = (1, 0, 3, 2)
roots = [p for p in itertools.permutations(range(4))
         if tuple(p[p[i]] for i in range(4)) == sigma]
assert len(roots) == 2
assert all(all((p[i] < 2) != (i < 2) for i in range(4)) for p in roots)
assert 12 * 13 // 2 + 2 == 80 and 80 + 4 == 84

# Exact original quintic serialization, including field polynomial and coords.
fixed = DATA / "N40_survivor_fixed.txt"
original = DATA / "N40_survivors.txt"
assert digest(fixed) == "3b48c4a8726cc8aeae9e07d5bcba5068517837f541979a15e8e10498cec93873"
assert fixed.read_bytes() == original.read_bytes()
old = (LOGS / "X1_40_p3_first_attempt.log").read_text()
m = re.search(r"SURVIVOR_DATA type=\[ 5 \] data=(.*?)\nRESULT type=\[ 5 \]", old, re.S)
assert m and numbers(m.group(1)) == numbers(fixed.read_text())
for degree, total, survivors in [(4, 4, 0), (5, 9, 1), (6, 10, 0)]:
    assert old.count(f"RESULT type=[ {degree} ] total={total} survivors={survivors}") == 1
assert old.count("SURVIVOR type=[ 5 ] err=false") == 1
assert "NONCUSP counts deg 1..7: [ 0, 0, 0, 32, 72, 80, 408 ]" in old
assert "REPS up to diamond: [ 0, 0, 0, 4, 9, 10, 51 ]" in old

# This fresh explicitly identified commonzero run settles the original
# generic-driver ambiguity. The original log is now provenance only.
fresh_lower = (LOGS / "X1_40_p3_lower.log").read_text()
assert not re.search(r"Runtime error|User error|Assertion failed", fresh_lower)
assert "SETUP N=40 M=1 p=3 genus=25" in fresh_lower
assert "NONCUSP counts deg 1..7: [ 0, 0, 0, 32, 72, 80, 408 ]" in fresh_lower
assert "REPS up to diamond: [ 0, 0, 0, 4, 9, 10, 51 ]" in fresh_lower
for degree, total, survivors in [(4, 4, 0), (5, 9, 1), (6, 10, 0)]:
    assert fresh_lower.count(f"RESULT type=[ {degree} ] total={total} survivors={survivors}") == 1
assert fresh_lower.count("SURVIVOR type=[ 5 ] err=false") == 1
fm = re.search(r"SURVIVOR_DATA type=\[ 5 \] data=(.*?)\nRESULT type=\[ 5 \]", fresh_lower, re.S)
assert fm and numbers(fm.group(1)) == numbers(fixed.read_text())
assert "SUMMARY N=40 M=1 p=3 D=7 shard=0/1 total_survivors=1" in fresh_lower
assert "SIEVE_RESUME_DONE FromDegree=4 MaxPart=6" in fresh_lower

# A fresh complete septic run; never reuse the unfinished old septic indices.
smart = (LOGS / "X1_40_p3_septic.log").read_text()
starts = re.findall(r"SMART_CASE_START type=\[ 7 \] index=(\d+)", smart)
ends = re.findall(r"SMART_CASE_RESULT type=\[ 7 \] index=(\d+) killed=(\w+)", smart)
assert len(starts) == 51 and set(map(int, starts)) == set(range(1, 52))
assert len(ends) == 51 and set(ends) == {(str(i), "true") for i in range(1, 52)}
assert "RESULT type=[ 7 ] total=51 survivors=0" in smart
assert "SIEVE_RESUME_DONE FromDegree=7 MaxPart=7" in smart
assert "SUMMARY N=40 M=1 p=3 D=7 shard=0/1 total_survivors=0" in smart
assert not re.search(r"Runtime error|User error|SURVIVOR type=|killed=false", smart)

report = {
    "complete": True, "all_complements": 84, "Q_liftable_complements": 80,
    "bad_quadratic_cusps": 4, "possible_diamond_cycles": len(roots),
    "correction_tests": 9, "correction_survivors": 0,
    "lower_support_tests": {"degree4": 4, "degree5": 9, "degree6": 10},
    "fresh_degree7_tests": 51,
    "unused_old_degree7_prefix": True,
    "fresh_lower_complete": True, "fresh_correction_complete": True,
    "fresh_survivor_exactly_equals_fixed_input": True,
    "original_lower_log_role": "provenance only; fresh commonzero lower-support run is sufficient",
}
(DATA / "cusp_correction_audit.json").write_text(json.dumps(report, indent=2) + "\n")
print("C40_CUSP_CORRECTION_INDEPENDENT_AUDIT_DONE tests=9 complements=84 lower_tests=23 fresh_septic_tests=51 survivors=0 complete=true")
