"""Dirichlet characters in the Conrey labelling (copied from lmfdb_ranks.py of the sextic paper)."""
from fractions import Fraction as Fr
from math import gcd


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def primitive_root(q, p):
    phi = q - q // p
    fac = [r for r in range(2, phi + 1) if phi % r == 0 and all(r % s for s in range(2, r))]
    for g in range(2, q):
        if gcd(g, q) == 1 and all(pow(g, phi // r, q) != 1 for r in fac):
            return g


def dlog(g, x, q):
    y, k = 1, 0
    while y != x % q:
        y = y * g % q; k += 1
    return k


def factor(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def conrey_value(n, c, a):
    """chi_n(c, a) = exp(2 pi i t), returns t in [0,1)."""
    t = Fr(0)
    for p, e in factor(n).items():
        q = p ** e; phi = q - q // p
        g = primitive_root(q, p)
        t += Fr(dlog(g, c, q) * dlog(g, a, q), phi)
    return t % 1
