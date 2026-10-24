#!/usr/bin/env python3
"""Extract SURVIVOR_DATA entries (possibly wrapped over several lines by Magma) from a sieve log.
Usage: extract_survivors.py LOG OUT  -> one Magma list expression per line in OUT."""
import sys, re
markers = ("TYPE ", "PROGRESS ", "RESULT ", "SURVIVOR ", "SURVIVOR_DATA ", "SUMMARY ", "SIEVE_DONE", "ERROR ")
lines = open(sys.argv[1]).read().split("\n")
out = []; cur = None
for l in lines:
    if l.startswith("SURVIVOR_DATA "):
        if cur is not None: out.append(cur)
        cur = l.split("data=", 1)[1]
    elif cur is not None:
        if l.startswith(markers) or l.strip() == "":
            out.append(cur); cur = None
        else:
            cur += " " + l.strip()
if cur is not None: out.append(cur)
open(sys.argv[2], "w").write("\n".join(out) + "\n")
print("%d survivors extracted" % len(out))
