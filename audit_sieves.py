#!/usr/bin/env python3
"""Audit all-partition degree-seven sieve coverage and independent counts."""
from pathlib import Path
import json
import re

HERE=Path(__file__).resolve().parent


def partitions(n, k=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n,k or n),0,-1):
        for rest in partitions(n-first,first):
            yield (first,)+rest


def audit(n,p):
    path=HERE/f'logs/sieve/X1_{n}_p{p}.log'
    text=path.read_text()
    assert not re.search(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',text,re.I|re.M)
    assert f'SUMMARY N={n} M=1 p={p} D=7 shard=0/1 total_survivors=0' in text
    assert 'SIEVE_DONE' in text
    done={}
    for m in re.finditer(r'TYPE \[([^\]]+)\]: no divisors|RESULT type=\[([^\]]+)\] total=(\d+) survivors=(\d+)',text):
        part=tuple(map(int,(m.group(1) or m.group(2)).replace(',',' ').split()))
        assert part not in done
        done[part]=int(m.group(4) or 0)
    expected={part for degree in range(1,8) for part in partitions(degree)}
    assert set(done)==expected and not any(done.values())
    countline=re.search(r'NONCUSP counts deg 1\.\.7: \[([^\]]+)\]',text)
    counts=[int(x) for x in countline.group(1).replace(',',' ').split()]
    independent=(HERE/'logs/counts/independent_counts.log').read_text()
    m=re.search(rf'DEGREE7_INDEPENDENT_COUNT m=1 n={n} p={p} .*?closed_place_counts=\[([^\]]+)\]',independent)
    assert m
    check=[int(x) for x in m.group(1).replace(',',' ').split()]
    assert counts==check
    return {'degree':7,'level':n,'prime':p,'partition_types':len(expected),'closed_place_counts':counts,'survivors':0,'independent_counts_agree':True,'complete':True}


if __name__=='__main__':
    result=[audit(n,p) for n,p in [(34,3),(38,3)]]
    (HERE/'data/sieve_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
