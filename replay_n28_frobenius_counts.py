#!/usr/bin/env python3
"""Independent integer Newton replay of the modular-symbol count export."""
import collections
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
source = HERE / "data/N28_D7_frobenius_counts.txt"
co, saved_traces, saved_closed, cusps, saved_nc = json.loads(
    source.read_text().replace("[*", "[").replace("*]", "]"))
assert len(co) == 51 and co[-1] == 1 and co[0] == 3**25
assert all(co[i] == 3**(25-i)*co[50-i] for i in range(26))
traces = []
closed = []
for n in range(1, 8):
    traces.append(-n*co[50-n] - sum(
        co[50-i]*traces[n-i-1] for i in range(1, n)))
    count = 3**n + 1 - traces[-1] - sum(
        d*closed[d-1] for d in range(1, n) if n % d == 0)
    assert count % n == 0
    closed.append(count//n)
record = next(r for r in json.loads((HERE / "data/cusp_conditions_m2_degree7.json").read_text())
              if r["n"] == 28)
degrees = collections.Counter(record["frobenius_degrees"])
assert cusps == [degrees[i] for i in range(1, 8)]
nc = [a-b for a, b in zip(closed, cusps)]
assert traces == saved_traces and closed == saved_closed
assert nc == saved_nc == [0, 0, 0, 0, 48, 64, 312]
print(json.dumps({"status": "passed", "NC_degrees_1_to_7": nc,
                  "traces": traces, "closed_places": closed,
                  "frobenius_export_sha256": hashlib.sha256(
                      source.read_bytes()).hexdigest()}, sort_keys=True))
