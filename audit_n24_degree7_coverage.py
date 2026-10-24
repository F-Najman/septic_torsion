#!/usr/bin/env python3
"""Check complete N24 support coverage and the exported residual classes.

This is a coverage/provenance replay. The exact divisor class maps and
comparison with known characteristic-zero witnesses are the Magma stage.
"""
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPECTED = {(3,): 4, (4,): 32, (5,): 168, (6,): 684, (7,): 2796,
            (3, 3): 64, (3, 4): 384}


def read_one(path):
    return json.loads(path.read_text().replace("[*", "[").replace("*]", "]"))


def main():
    log = (HERE / "logs/m2/X1_2_24_classes.log").read_text()
    assert not re.search(r"error|interrupt|segmentation|failed", log, re.I)
    matches = re.findall(r"N24_D7_TYPE\s+\[\s*([\d,\s]+?)\s*\]"
                         r"\s+(\d+)\s+(\d+)\s+([\d.]+)", log)
    counts = {}
    passes = {}
    for pat, count, passed, seconds in matches:
        pat = tuple(map(int, re.findall(r"\d+", pat)))
        assert pat not in counts
        counts[pat] = int(count)
        passes[pat] = int(passed)
    assert counts == EXPECTED
    assert "N24_D7_KNOWN_REDUCIBLE 64" in log
    statuses = re.findall(r"N24_D7_CLASS_STATUS\s+(\d+)\s+(\d+)\s+(\d+)", log)
    assert len(statuses) == 1
    classes, known, unknown = map(int, statuses[0])
    assert classes == known + unknown
    records = read_one(HERE / "data/N24_D7_survivor_classes.txt")
    N, D, signatures, degrees, base, vectors = read_one(
        HERE / "data/N24_D7_cusp_lifts.txt")
    assert (N, D) == (24, 7)
    assert len(signatures) == len(degrees) == 19
    assert sorted(degrees) == [1] * 8 + [2] * 6 + [4] * 5
    assert degrees[base - 1] == 1
    assert len(vectors) == unknown
    assert all(len(v) == 19 for v in vectors)
    support_keys = set()
    full_classes = set()
    for pat, indices, complement, full_class in records:
        pat = tuple(pat)
        assert pat in EXPECTED and len(indices) == len(pat)
        assert all(len(v) == 9 and v[-1] == 0 for v in [complement, full_class])
        limits = {(3,): [4], (4,): [32], (5,): [168], (6,): [684],
                  (7,): [2796], (3, 3): [4, 16], (3, 4): [4, 96]}[pat]
        assert all(1 <= i <= n for i, n in zip(indices, limits))
        support_keys.add((pat, tuple(indices)))
        full_classes.add(tuple(full_class))
    assert len(full_classes) == classes
    assert Counter(pat for pat, ix in support_keys) == Counter(
        {pat: n for pat, n in passes.items() if n})
    marker = (f"N24_D7_DONE tested={sum(EXPECTED.values())} "
              f"retained={sum(passes.values())} rational_classes={classes} "
              f"unresolved_classes={unknown} complete=true")
    assert marker in log
    result = {
        "status": "passed" if unknown == 0 else "complete-with-unresolved-classes",
        "tested": sum(EXPECTED.values()),
        "coverage": {str(p): n for p, n in counts.items()},
        "support_passes": {str(p): n for p, n in passes.items()},
        "distinct_rational_classes": classes,
        "known_reducible_classes": known, "unresolved_classes": unknown,
        "source_sha256": hashlib.sha256(
            (HERE / "fine_2_24_degree7.m").read_bytes()).hexdigest(),
    }
    (HERE / "data/N24_D7_coverage_replay.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
