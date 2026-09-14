#!/usr/bin/env python3
"""Shared QC-flagging primitives.

Pure functions used to flag batch samples based on QUAST and BactInspector
outputs (batch-relative MAD thresholds, excluded/failed sample handling,
and QC-control sample identification). Kept separate from any report/
dashboard generator so scripts that only need flagging logic don't have
to import one.
"""

from __future__ import annotations

import csv
import re
import statistics
from dataclasses import dataclass
from pathlib import Path

MAD_SCALE = 1.4826
MAD_MULTIPLIER = 3

QC_CONTROL_PATTERN = re.compile(r"^NGIVPZEC\d+(?:_S\d+)?$", re.IGNORECASE)


def robust_limits(values):
    median = statistics.median(values)
    raw_mad = statistics.median(abs(value - median) for value in values)
    scaled_mad = MAD_SCALE * raw_mad
    return {
        "median": median,
        "raw_mad": raw_mad,
        "scaled_mad": scaled_mad,
        "lower": median - MAD_MULTIPLIER * scaled_mad,
        "upper": median + MAD_MULTIPLIER * scaled_mad,
    }


@dataclass(frozen=True)
class Flag:
    sample: str
    step: str
    reason: str
    severity: str = "review"


def read_delimited(path: Path, delimiter: str = "\t") -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), list(reader)


def excluded_rows(path: Path | None) -> list[dict[str, str]]:
    return read_delimited(path)[1] if path and path.is_file() else []


def baseline_exclusions(samples, manual: str = "") -> set[str]:
    automatic = {sample for sample in samples if QC_CONTROL_PATTERN.fullmatch(sample)}
    specified = {sample.strip() for sample in manual.split(",") if sample.strip()}
    return automatic | specified


def parse_float(value: str) -> float | None:
    try:
        return float(value.replace(",", ""))
    except ValueError:
        return None


def numeric_columns(headers: list[str], rows: list[dict[str, str]]) -> list[str]:
    columns = []
    for header in headers[1:]:
        values = [row.get(header, "").strip() for row in rows]
        if values and all(value and parse_float(value) is not None for value in values):
            columns.append(header)
    return columns


def compute_bactinspector_flags(
    rows: list[dict[str, str]],
    headers: list[str],
    excluded_data: list[dict[str, str]],
) -> tuple[list[Flag], set[str]]:
    """Return the existing categorical and excluded/failed BactInspector flags."""
    del headers  # Kept in the pure-function interface for symmetry with QUAST input.
    flags = []
    flagged = set()
    for row in rows:
        sample = row.get("file", "")
        reasons = []
        if row.get("result", "").strip().lower() == "uncertain":
            reasons.append("result = uncertain")
        if row.get("species", "").strip().lower() == "no significant matches":
            reasons.append('species = "No significant matches"')
        flags.extend(Flag(sample, "BactInspector", reason) for reason in reasons)
        if reasons:
            flagged.add(sample)
    for row in excluded_data:
        sample = row.get("sample_id", "")
        flags.append(
            Flag(
                sample,
                "BactInspector",
                f"{row.get('status', 'excluded')}: {row.get('detail', '')}",
                "failure",
            )
        )
        flagged.add(sample)
    return flags, flagged


def compute_quast_flags(
    rows: list[dict[str, str]],
    headers: list[str],
    excluded_from_baseline: set[str],
    excluded_data: list[dict[str, str]],
) -> tuple[
    list[Flag],
    set[str],
    dict[str, list[str]],
    dict[str, dict[str, float]],
    dict[str, list[float]],
]:
    """Return existing QUAST MAD and excluded/failed flags plus computed plot data."""
    flags, flagged = [], set()
    reasons_by_sample: dict[str, list[str]] = {}
    baseline_rows = [row for row in rows if row.get(headers[0], "") not in excluded_from_baseline]
    if not baseline_rows:
        raise ValueError("QUAST baseline must contain at least one sample")
    limits_by_column = {}
    values_by_column = {}
    for column in numeric_columns(headers, rows):
        values = [parse_float(row[column]) for row in rows]
        baseline_values = [parse_float(row[column]) for row in baseline_rows]
        limits = robust_limits(baseline_values)
        limits_by_column[column] = limits
        values_by_column[column] = values
        for row, value in zip(rows, values):
            if abs(value - limits["median"]) > 3 * limits["scaled_mad"]:
                sample = row.get(headers[0], "")
                direction = "below" if value < limits["lower"] else "above"
                reason = (
                    f"{column} = {row[column]} ({direction} batch median \u00b1 3 scaled MAD "
                    f"[{limits['lower']:.3g}, {limits['upper']:.3g}])"
                )
                flags.append(Flag(sample, "QUAST", reason))
                flagged.add(sample)
                reasons_by_sample.setdefault(sample, []).append(reason)
    for row in excluded_data:
        sample = row.get("sample_id", "")
        reason = f"{row.get('status', 'excluded')}: {row.get('detail', '')}"
        flags.append(Flag(sample, "QUAST", reason, "failure"))
        flagged.add(sample)
    return flags, flagged, reasons_by_sample, limits_by_column, values_by_column
