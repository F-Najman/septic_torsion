#!/usr/bin/env python3
"""Audit individual C42/F5 nonprincipality tests for degree seven.

The inherited degree-six exclusion is not an input. Only its completed
individual tests with still-admissible Hecke primes are reused. In particular,
the old quartic tests and interrupted initial sextic loop are not used.
The default mode requires every new computation to finish. --partial reports
honest progress and writes a separate, explicitly incomplete coverage file.
This is a coverage/provenance audit, not a new computation of principality.
"""
from pathlib import Path
from collections import Counter
import argparse
import ast
import hashlib
import json
import math
import re

HERE = Path(__file__).resolve().parent
LOGS = HERE / 'logs'
DATA = HERE / 'data'
OLD_SHA = 'fa129d97164e224e5366e335594c412bd5835e7b7ebc40cc9828478caf30c6cf'
NEW_SHA = '24799ed61dbb19f9ac640c5af4449704cff137d744dff95e1e24c9452355b0fc'
MIXED_SHA = 'c7bf5b2975b20bbc3220714ece6b6170a2458f4bf7f08eab3c5aa26b898114d7'
ERROR = re.compile(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR', re.I | re.M)


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    s = p.read_text(errors='strict')
    assert not ERROR.search(s), (str(p), 'Magma error')
    return s


def partitions(n, least=1):
    if n == 0:
        yield ()
    for a in range(least, n + 1):
        for rest in partitions(n - a, a):
            yield (a,) + rest


def phi(n):
    return sum(math.gcd(a, n) == 1 for a in range(1, n + 1))


def cusp_audit():
    """DEHMZ Lemma 2.7, independently recomputed with integer arithmetic."""
    units = {a for a in range(42) if math.gcd(a, 42) == 1}
    hp = {pow(5, j, 42) for j in range(12)}
    assert hp == {1, 5, 17, 25, 37, 41}
    orbits = []
    for d in range(1, 43):
        if 42 % d:
            continue
        m = 42 // d
        if d > 2:
            h = {a for a in units if a % m == 1 % m}
        else:
            h = {a for a in units if a % m in {1 % m, -1 % m}}
        size = len(units) // len(h)
        geometric = phi(d) * phi(m) // 2
        assert geometric % size == 0
        f = next(j for j in range(1, 13) if pow(5, j, 42) in h)
        for _ in range(geometric // size):
            orbits.append((size, f, h))
    q_degrees = Counter(size for size, f, h in orbits)
    finite_degrees = Counter()
    for size, f, h in orbits:
        finite_degrees[f] += size // f
    assert q_degrees == {1: 12, 2: 6, 6: 4}
    assert finite_degrees == {1: 12, 2: 6, 3: 4, 6: 2}
    allowed = {}
    for k in range(8):
        allowed[k] = set(units)
        for size, f, h in orbits:
            if f <= k:
                allowed[k] &= {a * b % 42 for a in hp for b in h}
        assert allowed[k] == (units if k <= 2 else hp)
    return allowed, {'Q_orbit_degree_counts': dict(q_degrees),
                     'F5_closed_cusp_degree_counts': dict(finite_degrees),
                     'admissible_units_by_complement_bound': {k: sorted(v) for k, v in allowed.items()}}


def fixed_table(p, expected_sha, n, d, component=None):
    """Parse the exact table and check distinct closed Frobenius orbits.

    All points are normalized by their last coordinate 1. Frobenius is an
    explicitly computed GF(5)-linear map in the recorded polynomial basis.
    The original Magma exporter separately checked equality to every original
    normalized place, and the hash orbit routine checked whole-input coverage.
    """
    assert digest(p) == expected_sha, (str(p), 'fixed input changed')
    s = p.read_text()
    assert s.count(':=') == 1
    rows = ast.literal_eval(s.split(':=')[1].strip().removesuffix(';').replace('[*', '[').replace('*]', ']'))
    if component is not None:
        assert len(rows) == 2
        rows = rows[component]
    assert len(rows) == n
    polys = {tuple(r[1]) for r in rows}
    assert len(polys) == 1
    poly = next(iter(polys))
    assert len(poly) == d + 1 and poly[-1] == 1

    def mul(a, b):
        c = [0] * (2 * d - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                c[i + j] = (c[i + j] + ai * bj) % 5
        for j in range(2 * d - 2, d - 1, -1):
            for i in range(d):
                c[j - d + i] = (c[j - d + i] - c[j] * poly[i]) % 5
        return tuple(c[:d])

    basis_fifth = []
    for j in range(d):
        v = tuple(int(i == j) for i in range(d))
        v2 = mul(v, v)
        basis_fifth.append(mul(mul(v2, v2), v))

    def frob(v):
        return tuple(sum(v[j] * basis_fifth[j][i] for j in range(d)) % 5 for i in range(d))

    keys = set()
    for degree, defpoly, coords in rows:
        assert degree == d and tuple(defpoly) == poly
        assert len(coords) == 3 and all(len(v) == d for v in coords)
        assert all(type(a) is int and 0 <= a < 5 for v in coords for a in v)
        assert coords[-1] == [1] + [0] * (d - 1)
        pt = tuple(map(tuple, coords))
        orbit = []
        x = pt
        for _ in range(d):
            orbit.append(x)
            x = tuple(map(frob, x))
        assert x == pt and len(set(orbit)) == d
        key = min(orbit)
        assert key not in keys, 'duplicate closed place in fixed table'
        keys.add(key)
    return {'count': n, 'degree': d, 'sha256': expected_sha,
            'independent_distinct_Frobenius_orbits': len(keys), 'field_polynomial': poly}


def inherited_audit(allowed):
    s = read(LOGS / 'c42/degree6_prefix_session.log')
    assert 'NONCUSP counts deg 1..6: [ 0, 0, 4, 183, 648, 2286 ]' in s
    assert 'REPS up to diamond: [ 0, 0, 2, 33, 108, 392 ]' in s
    assert 'SURVIVOR type=' not in s and 'C42_CUBIC_PAIR_SURVIVOR' not in s
    # This historical log is intentionally interrupted in the UNUSED sextic
    # prefix. Only the independently completed branches below are used.
    assert 'SIEVE_DONE' not in s
    for d, n, qs in [(3, 2, [17]), (5, 108, [11, 13, 17])]:
        m = re.findall(rf'TYPE \[ {d} \]: {n} divisors \(up to diamond\), qlist=\[ ([\d, ]+) \]', s)
        assert len(m) == 1 and [int(a) for a in m[0].split(',')] == qs
        assert len(re.findall(rf'RESULT type=\[ {d} \] total={n} survivors=0 time=', s)) == 1
        assert set(qs) <= allowed[7 - d]
    assert [int(i) for i in re.findall(r'C42_EXPORT_INPUT (\d+)', s)] == list(range(1, 393))
    assert len(re.findall(r'C42_EXPORT_STEP 4 ', s)) == 392
    assert s.count('C42_FIXED_REPRESENTATIVES_EXPORTED count=392 roundtrip=true') == 1
    assert re.findall(r'C42_CUBIC_PAIR_TEST (\d+) (true|false) ', s) == [(str(i), 'true') for i in range(1, 9)]
    assert s.count('C42_CUBIC_PAIR_DONE total=8 survivors=0 complete=true') == 1
    pair_source = (HERE / 'export_42_fixed_reps.m').read_text()
    assert 'qlist:=[11,13,17]' in pair_source
    assert 'multiplicities[si]*HeckeCache[key]' in pair_source
    assert {11, 13, 17} <= allowed[1]
    table = fixed_table(DATA / 'c42_fixed_representatives.m', OLD_SHA, 392, 6)
    witnesses = []
    for shard in range(4):
        t = read(LOGS / f'c42/degree6_shard{shard}.log')
        assert t.count(f'C42_FIXED_SHARD_SETUP shard={shard}/4 count=98 total=392 genus=25 p=5 qlist=[11,13,17]') == 1
        records = [tuple(map(int, x)) for x in re.findall(r'C42_FIXED_CASE shard=(\d+) index=(\d+) nonprincipal_q=(\d+)', t)]
        expected = [i for i in range(1, 393) if i % 4 == shard]
        assert [i for sh, i, q in records] == expected
        assert all(sh == shard and q in allowed[1] and q == 11 for sh, i, q in records)
        assert 'C42_FIXED_SURVIVOR' not in t
        assert len(re.findall(rf'C42_FIXED_SHARD_DONE shard={shard}/4 tested=98 survivors=0 indices=\[\s*\] complete=true', t)) == 1
        witnesses.extend({'index': i, 'q': q, 'shard': sh} for sh, i, q in records)
    assert sorted(r['index'] for r in witnesses) == list(range(1, 393))
    return {'complete': True, 'used_counts': {'3': 2, '5': 108, '6': 392, '3,3': 8},
            'discarded_branches': ['old quartic tests: complement degree changed', 'interrupted original sextic loop'],
            'prime_information': {'3': {'actual_q': 17}, '5': {'tried_q': [11, 13, 17], 'individual_winner_not_logged': True},
                                  '6': {'all_actual_q': 11}, '3,3': {'tried_q': [11, 13, 17], 'individual_winner_not_logged': True}},
            'fixed_table': table, 'sextic_witnesses': witnesses}


def current_audit(logdir, allowed):
    pending = []
    validation = read(LOGS / 'tests/test_commonzero_cyclic.log')
    assert 'COMMONZERO_CYCLIC_REFERENCE_EQUAL' in validation
    assert 'COMMONZERO_CYCLIC_TEST_DONE' in validation
    validation_bc = read(LOGS / 'tests/test_cyclic_bc_hash_36.log')
    assert 'CYCLIC_BC_HASH36_VALIDATED places=48 reps=8 complete=true' in validation_bc
    export = read(logdir / 'export7.log')
    assert export.count('N42_D7_EXPORT_DONE places=11160 representatives=1860 roundtrip=true complete=true') == 1
    table = fixed_table(DATA / 'N42_degree7_representatives.m', NEW_SHA, 1860, 7)
    counts = read(LOGS / 'counts/additional.log')
    line = next(line for line in counts.splitlines() if 'DEGREE7_FAST_COUNT m=1 n=42 p=5 ' in line)
    rat = ast.literal_eval(re.search(r'rational_point_counts=(\[[^\]]+\])', line).group(1))
    closed = ast.literal_eval(re.search(r'closed_place_counts=(\[[^\]]+\])', line).group(1))
    assert closed == [0, 0, 4, 183, 648, 2286, 11160]
    for d in range(1, 8):
        assert rat[d - 1] == sum(e * closed[e - 1] for e in range(1, d + 1) if d % e == 0)
    assert 'DEGREE7_FAST_COUNTS_DONE complete=true' in counts
    lower = read(logdir / 'lower7.log')
    grade4 = [(int(i), b) for i, b in re.findall(r'N42_D7_GRADE4 (\d+) (true|false) ', lower)]
    assert grade4 == [(i, 'false') for i in range(1, len(grade4) + 1)] and len(grade4) <= 33
    assert 17 in allowed[3]
    # The stopped original mixed34 prefix is not used. All mixed pairs are
    # independently supplied by four new shards of one round-tripped table.
    unused_mixed = re.findall(r'N42_D7_MIXED34 (\d+) (\d+) (true|false) ', lower)
    assert {11, 13, 17} <= allowed[0]
    if len(grade4) != 33:
        pending.append('new quartic completion')
    mixed_export = read(logdir / 'export_mixed34.log')
    assert mixed_export.count('N42_D7_MIXED_EXPORT_DONE cubic_reps=2 quartic_places=183 roundtrip=true complete=true') == 1
    mixed_tables = [fixed_table(DATA / 'N42_mixed34_fixed.m', MIXED_SHA, 2, 3, component=0),
                    fixed_table(DATA / 'N42_mixed34_fixed.m', MIXED_SHA, 183, 4, component=1)]
    mixed, mixed_shards = [], []
    for shard in range(4):
        name = f'mixed34_shard{shard}.log'
        s = read(logdir / name)
        expected = [(i, j) for i in range(1, 3) for j in range(1, 184) if (j - 1) % 4 == shard]
        setup_count = len(re.findall(rf'N42_D7_MIXED_SHARD_SETUP {shard} ', s))
        if setup_count == 0:
            assert 'N42_D7_MIXED_SHARD_CASE' not in s and 'N42_D7_MIXED_SHARD_DONE' not in s
            pending.append(f'mixed34 shard {shard}: initializing')
            mixed_shards.append({'shard': shard, 'tested': 0, 'expected': len(expected), 'complete': False,
                                 'log_sha256': None})
            continue
        assert setup_count == 1
        records = [(int(sh), int(i), int(j), b) for sh, i, j, b in re.findall(r'N42_D7_MIXED_SHARD_CASE (\d+) (\d+) (\d+) (true|false) ', s)]
        assert [(i, j) for sh, i, j, b in records] == expected[:len(records)] and len(records) <= len(expected)
        assert all(sh == shard and b == 'true' for sh, i, j, b in records)
        marker = f'N42_D7_MIXED_SHARD_DONE shard={shard}/4 tested={len(expected)} passing=0 complete=true'
        done = marker in s
        if done:
            assert len(records) == len(expected) and s.count(marker) == 1
        else:
            pending.append(f'mixed34 shard {shard}: {len(records)}/{len(expected)}')
        mixed.extend({'cubic_index': i, 'quartic_index': j, 'shard': sh} for sh, i, j, b in records)
        mixed_shards.append({'shard': shard, 'tested': len(records), 'expected': len(expected), 'complete': done,
                             'log_sha256': digest(logdir / name) if done else None})
    assert len({(r['cubic_index'], r['quartic_index']) for r in mixed}) == len(mixed)
    witnesses, shards = [], []
    for shard in range(8):
        name = f'test7_shard{shard}.log'
        s = read(logdir / name)
        expected = [i for i in range(1, 1861) if i % 8 == shard]
        setup = re.findall(r'N42_D7_SHARD_SETUP (\d+) (\d+) (\d+) ', s)
        assert setup == [(str(shard), '8', str(len(expected)))]
        records = [tuple(map(int, m)) for m in re.findall(r'N42_D7_SHARD_CASE index=(\d+) q=(\d+) nonprincipal=true', s)]
        assert [i for i, q in records] == expected[:len(records)] and len(records) <= len(expected)
        assert all(q in {11, 13, 17} and q in allowed[0] for i, q in records)
        assert 'N42_D7_SHARD_SURVIVOR' not in s
        marker = f'N42_D7_SHARD_DONE shard={shard}/8 tested={len(expected)} total={len(expected)} passing=0 complete=true'
        done = marker in s
        if done:
            assert len(records) == len(expected) and s.count(marker) == 1
        else:
            pending.append(f'septic shard {shard}: {len(records)}/{len(expected)}')
        witnesses.extend({'index': i, 'q': q, 'shard': shard} for i, q in records)
        shards.append({'shard': shard, 'tested': len(records), 'expected': len(expected), 'complete': done,
                       'log_sha256': digest(logdir / name) if done else None})
    assert len({r['index'] for r in witnesses}) == len(witnesses)
    if not pending:
        assert sorted(r['index'] for r in witnesses) == list(range(1, 1861))
        assert {(r['cubic_index'], r['quartic_index']) for r in mixed} == {(i, j) for i in range(1, 3) for j in range(1, 184)}
    return {'complete': not pending, 'pending': pending, 'fixed_table': table,
            'closed_noncuspidal_place_counts': closed,
            'quartic_test_count': len(grade4), 'mixed34_test_count': len(mixed),
            'ignored_interrupted_mixed_prefix_count': len(unused_mixed),
            'mixed_fixed_tables': mixed_tables, 'mixed_shards': mixed_shards,
            'mixed34_witnesses': mixed,
            'shards': shards, 'septic_witnesses': sorted(witnesses, key=lambda r: r['index'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partial', action='store_true')
    args = parser.parse_args()
    allowed, cusp = cusp_audit()
    inherited = inherited_audit(allowed)
    logdir = LOGS / 'c42'
    current = current_audit(logdir, allowed)
    partition_set = [part for d in range(1, 8) for part in partitions(d)]
    active = [part for part in partition_set if min(part) >= 3]
    assert len(partition_set) == 44 and set(active) == {(3,), (4,), (5,), (6,), (7,), (3, 3), (3, 4)}
    expected = {'3': 2, '4': 33, '5': 108, '6': 392, '7': 1860, '3,3': 8, '3,4': 366}
    assert sum(expected.values()) == 2769
    report = {'complete': current['complete'], 'degree': 7, 'level': 42, 'p': 5,
              'positive_noncuspidal_degree_partitions': 44, 'nonempty_types': expected,
              'all_cuspidal_case': 'separate mathematical rank-zero / torsion reduction argument',
              'cusp_admissibility': cusp, 'inherited_individual_tests': inherited, 'new_tests': current}
    if current['pending'] and not args.partial:
        raise AssertionError('INCOMPLETE: ' + '; '.join(current['pending']))
    name = 'coverage42_degree7_partial.json' if args.partial else 'coverage42_degree7.json'
    (DATA / name).write_text(json.dumps(report, indent=2) + '\n')
    print('C42_DEGREE7_COVERAGE_AUDIT complete=' + str(current['complete']).lower(),
          'inherited_tests=510', f'new_grade4={current["quartic_test_count"]}/33',
          f'new_mixed34={current["mixed34_test_count"]}/366',
          f'new_septic={len(current["septic_witnesses"])}/1860',
          'positive_partitions=44 all_cuspidal_case=separate')
    if current['pending']:
        print('PENDING: ' + '; '.join(current['pending']))


if __name__ == '__main__':
    main()
