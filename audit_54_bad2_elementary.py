#!/usr/bin/env python3
"""Exact elementary checks for the septic X1(54) all-bad-at2 lemma.

The DR semistable-model theorem and published rational-Jacobian bound are
mathematical inputs, not conclusions of this short arithmetic check.
"""
from collections import Counter
from fractions import Fraction
from math import gcd


def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def cusp_orbits(n):
    units = {a for a in range(n) if gcd(a, n) == 1}
    result = []
    for d in range(1, n + 1):
        if n % d:
            continue
        m = n // d
        h = {a for a in units if a % m == 1 % m} if d > 2 else {
            a for a in units if a % m in {1 % m, -1 % m}}
        size = len(units) // len(h)
        geometric = phi(d) * phi(m) // 2
        assert geometric % size == 0
        result.extend([(d, size, h)] * (geometric // size))
    return result


assert {pow(2, j, 27) for j in range(18)} == {a for a in range(27) if gcd(a, 27) == 1}
assert pow(2, 9, 27) == 26
print('N54_BAD2_FROBENIUS primitive_mod27=true order=18')

# Gal(Q(zeta54)/Q)=Gal(Q(zeta27)/Q). Frobenius at 2 is represented
# by 29 mod54, the odd unit congruent to 2 mod27, not by 2 mod54.
sigma = 29
assert sigma % 27 == 2 and gcd(sigma, 54) == 1
orbits54 = cusp_orbits(54)
global_degrees = Counter(size for d, size, h in orbits54)
residue_degrees = Counter()
for d, size, h in orbits54:
    f = next(j for j in range(1, 19) if pow(sigma, j, 54) in h)
    assert f == size
    residue_degrees[f] += size // f
assert global_degrees == residue_degrees == {1: 18, 2: 6, 6: 2, 9: 2}
print('N54_BAD2_CUSP_ORBITS', dict(sorted(global_degrees.items())), 'inert_in_every_cusp_field=true')

cusps27 = sum(size for d, size, h in cusp_orbits(27))
cusps54 = sum(size for d, size, h in orbits54)
assert (cusps27, cusps54) == (30, 60)
i27 = Fraction(27**2, 2) * Fraction(8, 9)
i54 = Fraction(54**2, 2) * Fraction(3, 4) * Fraction(8, 9)
assert (i27, i54) == (324, 972)
g27 = 1 + i27 / 12 - Fraction(cusps27, 2)
g54 = 1 + i54 / 12 - Fraction(cusps54, 2)
assert (g27, g54) == (13, 52)
# Two genus13 components and27 transverse nodes give arithmetic genus52.
assert 2 * g27 + 27 - 1 == g54
print('N54_BAD2_MODEL_COUNTS geometric_cusps=60 component_cusps=30 components=2 genus=52 nodes=27')

bound = Fraction(325, 32768) * i54
assert bound > 9 > 7
print('N54_BAD2_GONALITY_BOUND', str(bound), '>9>7')

cusp_invariants = [3, 3, 3, 9, 9, 1102437, 1529080119]
assert all(a % 2 == 1 for a in cusp_invariants)
assert all(index % 2 == 1 for index in [1, 3])
print('N54_BAD2_PUBLISHED_TORSION odd_cuspidal_invariants=true possible_indices=1,3 odd_JQ=true')
print('N54_BAD2_ELEMENTARY_AUDIT_DONE complete=true')
