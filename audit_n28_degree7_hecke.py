#!/usr/bin/env python3
"""Combine completed lower-degree prefix with two fresh complete seven shards.

The interrupted tail of the old run contributes no degree-seven evidence.
"""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parse_cases(text, pattern):
    data = {}
    for m in re.finditer(pattern, text):
        key = tuple(map(int, m.group(1, 2)))
        data.setdefault(key, []).append((int(m.group(3)), m.group(4) == "true"))
    return data


def check_negative(tests):
    assert tests in [[(5, False)], [(5, True), (11, False)]]


def main():
    old = (HERE / "logs/m2/X1_2_28_sieve_interrupted.log").read_text()
    endpoint = re.search(r"SEPTIC_M2_28_TYPE\s+6\s+16\s+0\s+[\d.]+", old)
    assert endpoint is not None
    prefix = old[:endpoint.end()]
    assert not re.search(r"error|interrupt|segmentation|failed", prefix, re.I)
    lower = parse_cases(prefix, r"SEPTIC_M2_28_CASE\s+(\d+)\s+(\d+)"
                        r"\s+(\d+)\s+(true|false)")
    assert set(lower) == {(d, i) for d, n in [(5, 8), (6, 16)]
                          for i in range(1, n + 1)}
    for tests in lower.values():
        check_negative(tests)
    lower_text = prefix + "\n"
    definition, records = json.loads((HERE / "data/N28_D7_hecke_input.txt").read_text()
                                     .replace("[*", "[").replace("*]", "]"))
    assert len(definition) == 8 and definition[-1] == 1 and len(records) == 52
    assert all(len(pair) == 2 and all(len(v) == 7 for v in pair)
               for pair in records)
    exportlog = (HERE / "logs/m2/X1_2_28_export.log").read_text()
    assert "N28_D7_FIXED_INPUT_DONE records=52 complete=true" in exportlog
    assert list(map(int, re.findall(r"N28_D7_FIXED_INPUT\s+(\d+)\s+roundtrip=true",
                                   exportlog))) == list(range(1, 53))
    allindices = set()
    fresh_tests = []
    for shard in [0, 1]:
        text = (HERE / f"logs/m2/X1_2_28_hecke_shard{shard}.log").read_text()
        assert not re.search(r"error|interrupt|segmentation|failed", text, re.I)
        cases = parse_cases(text, r"N28_D7_FRESH_HECKE_CASE\s+(\d+)\s+(\d+)"
                            r"\s+(\d+)\s+(true|false)")
        wanted = {(shard, i) for i in range(1, 53) if (i - 1) % 2 == shard}
        assert set(cases) == wanted
        for tests in cases.values():
            check_negative(tests)
            fresh_tests.append(tests)
        assert f"N28_D7_FRESH_HECKE_DONE shard={shard} tested=26 remaining=0 complete=true" in text
        assert not (allindices & {i for s, i in cases})
        allindices |= {i for s, i in cases}
    assert allindices == set(range(1, 53))
    # Preserve only the completed degree-five/six evidence as a separate file.
    # The original log remains bound, including its unused interrupted tail.
    (HERE / "data/N28_D7_lower_completed_prefix.log").write_text(lower_text)
    result = {"status": "passed", "prime": 3, "level": 28, "degree": 7,
              "complete_lower_cases": {"5": 8, "6": 16},
              "fresh_degree_seven_cases": 52, "total": 76,
              "remaining": 0,
              "degree_five_A11_tests": sum(len(v) == 2 for (d, i), v
                                           in lower.items() if d == 5),
              "degree_six_A11_tests": sum(len(v) == 2 for (d, i), v
                                          in lower.items() if d == 6),
              "fresh_degree_seven_A11_tests": sum(len(v) == 2
                                                  for v in fresh_tests),
              "old_interrupted_degree_seven_tail_used": False,
              "completed_lower_prefix_sha256": hashlib.sha256(
                  lower_text.encode()).hexdigest(),
              "fixed_input_sha256": hashlib.sha256(
                  (HERE / "data/N28_D7_hecke_input.txt").read_bytes()).hexdigest(),
              "shard_source_sha256": hashlib.sha256(
                  (HERE / "n28_degree7_hecke_shard.m").read_bytes()).hexdigest()}
    (HERE / "data/N28_D7_fresh_hecke_coverage.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
