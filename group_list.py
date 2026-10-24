# Independent enumeration of minimal obstructions for degree 7.
# Groups: C_n (m=1) and C_2 x C_2k (m=2); primes dividing the order in S(7)={2,..,23} (DKSS).
from sympy import primefactors
S=[2,3,5,7,11,13,17,19,23]
G=set()
for n in list(range(1,29))+[30,33,36]: G.add((1,n))
for k in range(1,11): G.add((2,2*k))
assert len(G)==41
def ok(n): return all(p in S for p in primefactors(n))
def maxsubs(g):
    m,n=g; out=set()
    if m==1:
        for p in primefactors(n): out.add((1,n//p))
    else:
        out.add((1,n))                                  # cyclic subgroup of index 2
        for p in primefactors(n):
            if p==2:
                if (n//2)%2==0: out.add((2,n//2))
                # C2 x C_{n/2} with n/2 odd is cyclic C_n... handled: subgroup C_2 x C_{n/2} ~ C_{n} only if gcd; use (1,n) above
                else: out.add((1,n))
            else: out.add((2,n//p))
    return out
# Every proper subgroup lies in a maximal subgroup; a group is a minimal obstruction iff not in G and all maximal subgroups in G.
B=max(n for m,n in G if m==1)
obs=[]
for n in range(2,23*23*B+1):
    if ok(n) and (1,n) not in G and all(s in G for s in maxsubs((1,n))): obs.append((1,n))
for k in range(1,23*B):
    n=2*k
    if ok(n) and (2,n) not in G and all(s in G for s in maxsubs((2,n))): obs.append((2,n))
print(len(obs)); print(sorted(n for m,n in obs if m==1)); print(sorted(n for m,n in obs if m==2))
