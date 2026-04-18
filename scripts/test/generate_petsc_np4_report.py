#!/usr/bin/env python3

"""Generate PETSc np=1 vs np=4 comparison artefacts from ctest JUnit XML."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
import xml.etree.ElementTree as ET


@dataclass(frozen=True)
class TestCase:
    name: str
    status: str
    time_s: float
    system_out: str


def parse_junit(path: Path) -> dict[str, TestCase]:
    tree = ET.parse(path)
    root = tree.getroot()
    cases: dict[str, TestCase] = {}
    for testcase in root.findall("testcase"):
        name = testcase.attrib["name"]
        time_s = float(testcase.attrib.get("time", "0") or "0")
        status_attr = testcase.attrib.get("status", "run").lower()
        has_failure = testcase.find("failure") is not None
        has_skipped = testcase.find("skipped") is not None
        if has_failure or status_attr == "fail":
            status = "FAIL"
        elif has_skipped or status_attr in {"notrun", "skip", "skipped"}:
            status = "SKIP"
        else:
            status = "PASS"
        system_out = (testcase.findtext("system-out") or "").strip()
        cases[name] = TestCase(name=name, status=status, time_s=time_s, system_out=system_out)
    return cases


def normalise_np4_name(np4_name: str) -> str:
    """Map generated -np4 tests back to baseline names."""
    return np4_name.replace("-np4", "", 1)


def classify_failure_reason(system_out: str, status: str) -> str:
    if status != "FAIL":
        return ""
    lower = system_out.lower()
    if "partitioned_msh_cfg4.bin does not exist" in lower:
        return "missing cfg4 partitioned mesh"
    if "error reading file" in lower and ".vtu" in lower:
        return "missing np4 result file for vtkdiff"
    if "no such file or directory" in lower and "vtkdiff" in lower:
        return "missing vtkdiff tool"
    if "timeout" in lower:
        return "timeout"
    if system_out:
        return system_out.splitlines()[0][:160]
    return "failed"


def latex_escape(text: str) -> str:
    replacements = {
        "\\": "\\textbackslash{}",
        "&": "\\&",
        "%": "\\%",
        "$": "\\$",
        "#": "\\#",
        "_": "\\_",
        "{": "\\{",
        "}": "\\}",
        "~": "\\textasciitilde{}",
        "^": "\\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def fmt_float(value: float | None, digits: int = 3) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-junit", type=Path, required=True)
    parser.add_argument("--np4-junit", type=Path, required=True)
    parser.add_argument("--csv-out", type=Path, required=True)
    parser.add_argument("--tex-rows-out", type=Path, required=True)
    parser.add_argument("--summary-json-out", type=Path, required=True)
    args = parser.parse_args()

    baseline_cases = parse_junit(args.baseline_junit)
    np4_cases = parse_junit(args.np4_junit)

    rows: list[dict[str, str | float | None]] = []
    mapped_baseline_names: set[str] = set()

    for np4_name in sorted(np4_cases):
        np4_case = np4_cases[np4_name]
        baseline_name = normalise_np4_name(np4_name)
        baseline_case = baseline_cases.get(baseline_name)
        mapped_baseline_names.add(baseline_name)

        baseline_status = baseline_case.status if baseline_case else "MISSING"
        baseline_time = baseline_case.time_s if baseline_case else None
        np4_status = np4_case.status
        np4_time = np4_case.time_s
        delta = np4_time - baseline_time if baseline_time is not None else None
        ratio = np4_time / baseline_time if baseline_time not in (None, 0.0) else None
        reason = classify_failure_reason(np4_case.system_out, np4_status)

        rows.append(
            {
                "baseline_test": baseline_name,
                "np4_test": np4_name,
                "baseline_status": baseline_status,
                "baseline_time_s": baseline_time,
                "np4_status": np4_status,
                "np4_time_s": np4_time,
                "delta_s": delta,
                "ratio_np4_over_np1": ratio,
                "np4_failure_reason": reason,
            }
        )

    baseline_only = sorted(set(baseline_cases.keys()) - mapped_baseline_names)
    np4_unmatched = sorted(
        np4_name for np4_name in np4_cases if normalise_np4_name(np4_name) not in baseline_cases
    )

    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "baseline_test",
                "np4_test",
                "baseline_status",
                "baseline_time_s",
                "np4_status",
                "np4_time_s",
                "delta_s",
                "ratio_np4_over_np1",
                "np4_failure_reason",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "baseline_test": row["baseline_test"],
                    "np4_test": row["np4_test"],
                    "baseline_status": row["baseline_status"],
                    "baseline_time_s": fmt_float(row["baseline_time_s"], 6)
                    if row["baseline_time_s"] is not None
                    else "",
                    "np4_status": row["np4_status"],
                    "np4_time_s": fmt_float(row["np4_time_s"], 6)
                    if row["np4_time_s"] is not None
                    else "",
                    "delta_s": fmt_float(row["delta_s"], 6) if row["delta_s"] is not None else "",
                    "ratio_np4_over_np1": fmt_float(row["ratio_np4_over_np1"], 6)
                    if row["ratio_np4_over_np1"] is not None
                    else "",
                    "np4_failure_reason": row["np4_failure_reason"],
                }
            )

    with args.tex_rows_out.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                "{} & {} & {} & {} & {} & {} & {} & {} \\\\\n".format(
                    latex_escape(str(row["baseline_test"])),
                    row["baseline_status"],
                    fmt_float(row["baseline_time_s"], 3) if row["baseline_time_s"] is not None else "--",
                    row["np4_status"],
                    fmt_float(row["np4_time_s"], 3) if row["np4_time_s"] is not None else "--",
                    fmt_float(row["delta_s"], 3) if row["delta_s"] is not None else "--",
                    fmt_float(row["ratio_np4_over_np1"], 3)
                    if row["ratio_np4_over_np1"] is not None
                    else "--",
                    latex_escape(str(row["np4_failure_reason"])),
                )
            )

        if baseline_only:
            handle.write("% baseline-only (no np4 mapping)\n")
            for name in baseline_only:
                handle.write(
                    "{} & PASS & -- & MISSING & -- & -- & -- & no generated np4 test \\\\\n".format(
                        latex_escape(name)
                    )
                )

        if np4_unmatched:
            handle.write("% np4-only (no baseline mapping)\n")
            for name in np4_unmatched:
                handle.write(
                    "{} & MISSING & -- & FAIL & -- & -- & -- & np4 has no baseline mapping \\\\\n".format(
                        latex_escape(name)
                    )
                )

    baseline_pass = sum(1 for case in baseline_cases.values() if case.status == "PASS")
    baseline_fail = sum(1 for case in baseline_cases.values() if case.status == "FAIL")
    np4_pass = sum(1 for case in np4_cases.values() if case.status == "PASS")
    np4_fail = sum(1 for case in np4_cases.values() if case.status == "FAIL")
    pair_pass_pass = sum(
        1
        for row in rows
        if row["baseline_status"] == "PASS" and row["np4_status"] == "PASS"
    )
    pair_pass_fail = sum(
        1
        for row in rows
        if row["baseline_status"] == "PASS" and row["np4_status"] == "FAIL"
    )
    pair_fail_fail = sum(
        1
        for row in rows
        if row["baseline_status"] == "FAIL" and row["np4_status"] == "FAIL"
    )

    summary = {
        "baseline_total": len(baseline_cases),
        "baseline_pass": baseline_pass,
        "baseline_fail": baseline_fail,
        "np4_total": len(np4_cases),
        "np4_pass": np4_pass,
        "np4_fail": np4_fail,
        "paired_total": len(rows),
        "pair_pass_pass": pair_pass_pass,
        "pair_pass_fail": pair_pass_fail,
        "pair_fail_fail": pair_fail_fail,
        "baseline_only_count": len(baseline_only),
        "np4_unmatched_count": len(np4_unmatched),
    }

    args.summary_json_out.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

