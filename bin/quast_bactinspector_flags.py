#!/usr/bin/env python3
"""Emit non-control samples flagged by batch QUAST or BactInspector QC."""

from __future__ import annotations

import argparse
from pathlib import Path

from qc_flagging import (
    QC_CONTROL_PATTERN,
    baseline_exclusions,
    compute_bactinspector_flags,
    compute_quast_flags,
    excluded_rows,
    read_delimited,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quast-summary", required=True, type=Path)
    parser.add_argument("--quast-excluded", required=True, type=Path)
    parser.add_argument("--bactinspector-summary", required=True, type=Path)
    parser.add_argument("--bactinspector-excluded", required=True, type=Path)
    parser.add_argument(
        "--expected-samples",
        required=True,
        help="Comma-separated list of sample IDs expected in this batch.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    expected_sample_ids = {
        sample.strip() for sample in args.expected_samples.split(",") if sample.strip()
    }

    quast_headers, quast_rows = read_delimited(args.quast_summary)
    bactinspector_headers, bactinspector_rows = read_delimited(args.bactinspector_summary)
    excluded_from_baseline = baseline_exclusions(expected_sample_ids)

    _, quast_flagged, _, _, _ = compute_quast_flags(
        quast_rows,
        quast_headers,
        excluded_from_baseline,
        excluded_rows(args.quast_excluded),
    )
    _, bactinspector_flagged = compute_bactinspector_flags(
        bactinspector_rows,
        bactinspector_headers,
        excluded_rows(args.bactinspector_excluded),
    )

    flagged = quast_flagged | bactinspector_flagged
    for sample_id in sorted(flagged):
        if not QC_CONTROL_PATTERN.fullmatch(sample_id):
            print(sample_id)


if __name__ == "__main__":
    main()
