# Independent LMFDB check of the rank/character inputs of the degree-7 draft.
import json, sys
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, __import__("os").path.dirname(__file__))
from conrey import conrey_value, divisors
import psycopg2
# kernel conditions: every positive-rank newform of level dividing n has chi(a)=1 for a in list
KER = {119:[3],121:[2],125:[3],133:[64],143:[67],169:[3],187:[122],209:[34],221:[3],247:[144],
       289:[3],323:[101],361:[3], 77:[3],85:[2],91:[3],95:[11],
       81:[2],99:[2],115:[2],161:[2],253:[2],299:[2],391:[2],437:[2],529:[2],
       57:[2,8],65:[3,6,8,4]}
# rank zero claims for J_1(n) (no positive-rank newform of level dividing n)
RANK0 = [25,27,32,33,34,35,36,38,39,40,42,45,46,49,51,54,55,69]
# quotient claims: J(X_1(n)/<h>) rank 0 <=> every positive-rank form of level dividing n has chi(h) != 1
QUOT = {(63,8):'Y63', (44,23):'X1(2,22)', (48,25):'X1(2,24)', (52,27):'X1(2,26)', (56,29):'X1(2,28)', (60,31):'X1(2,30)', (36,19):'X1(2,18)'}
levels = sorted({d for n in list(KER)+RANK0+[n for n,h in QUOT] for d in divisors(n) if d >= 11})
conn = psycopg2.connect(host="devmirror.lmfdb.xyz", dbname="lmfdb", user="lmfdb", password="lmfdb", connect_timeout=60)
cur = conn.cursor()
cur.execute("""select label, level, conrey_index, char_order, dim, analytic_rank, analytic_rank_proved
               from mf_newforms where weight=2 and level = ANY(%s)""", (levels,))
rows = cur.fetchall()
pos = [r for r in rows if r[5] is not None and r[5] > 0]
unknown = [r for r in rows if r[5] is None]
print("newforms fetched:", len(rows), " positive rank:", len(pos), " analytic rank unknown:", len(unknown))
for r in unknown[:10]: print("  UNKNOWN RANK", r)
ok = True
def posforms(n): return [r for r in pos if n % r[1] == 0]
for n, As in sorted(KER.items()):
    for r in posforms(n):
        for a in As:
            v = conrey_value(r[1], int(r[2]), a % r[1])
            if v != 0: ok=False; print("FAIL kernel", n, a, r, "chi(a)=exp(2pi i", v, ")")
        if not r[6]: print("  rank not proved:", r)
    print("KER", n, As, "posforms:", [(r[0], r[5]) for r in posforms(n)])
for n in RANK0:
    pf = posforms(n)
    if pf: ok=False
    print("RANK0", n, "OK" if not pf else ("FAIL", [(r[0],r[5]) for r in pf]))
for (n,h),name in sorted(QUOT.items()):
    bad=[r for r in posforms(n) if conrey_value(r[1], int(r[2]), h % r[1]) == 0]
    if bad: ok=False
    print("QUOT", name, "level", n, "h", h, "posforms", [(r[0],r[5]) for r in posforms(n)], "OK" if not bad else ("FAIL", bad))
print("RANKS7_DONE", "OK" if ok else "FAILURES")
