#!/usr/bin/env python3
"""Extract target-taxon percentages from one Kraken2 report."""

import argparse
import csv
from decimal import Decimal
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target-genus-taxid", required=True, type=int)
    parser.add_argument("--target-species-taxid", required=True, type=int)
    return parser.parse_args()


def read_metrics(report, genus_taxid, species_taxid):
    target_taxids = {
        "0": "unclassified",
        str(genus_taxid): "genus",
        str(species_taxid): "species",
    }
    percentages = {name: Decimal("0.00") for name in target_taxids.values()}
    root_clade_reads = None
    unclassified_reads = None
    genus_direct_reads = None
    with report.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(
                    f"{report}:{line_number}: expected 6 tab-separated fields"
                )
            taxid = fields[4].strip()
            if taxid in target_taxids:
                value = fields[0].strip().replace(",", ".")
                percentages[target_taxids[taxid]] = Decimal(value)
            rank = fields[3].strip()
            if rank == "U" and taxid == "0":
                unclassified_reads = int(fields[1].strip())
            elif rank == "R" and taxid == "1":
                root_clade_reads = int(fields[1].strip())
            elif rank == "G" and taxid == str(genus_taxid):
                genus_direct_reads = int(fields[2].strip())

    missing = [
        name
        for name, value in (
            ("unclassified taxid 0", unclassified_reads),
            ("root taxid 1", root_clade_reads),
            (f"genus taxid {genus_taxid}", genus_direct_reads),
        )
        if value is None
    ]
    if missing:
        raise ValueError(f"{report}: missing required row(s): {', '.join(missing)}")
    total_reads = root_clade_reads + unclassified_reads
    if total_reads == 0:
        raise ValueError(f"{report}: total read count is zero")
    genus_only_pct = Decimal(genus_direct_reads) * Decimal("100") / Decimal(total_reads)
    return percentages, genus_only_pct


def main():
    args = parse_args()
    percentages, genus_only_pct = read_metrics(
        args.report, args.target_genus_taxid, args.target_species_taxid
    )
    remainder = (
        Decimal("100")
        - percentages["unclassified"]
        - percentages["species"]
    )

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "unclassified_pct",
                "genus_pct",
                "species_pct",
                "other_pct",
                "genus_only_unresolved_pct",
            ]
        )
        writer.writerow(
            [
                args.sample_id,
                f'{percentages["unclassified"]:.2f}',
                f'{percentages["genus"]:.2f}',
                f'{percentages["species"]:.2f}',
                f"{remainder:.2f}",
                f"{genus_only_pct:.6f}",
            ]
        )


if __name__ == "__main__":
    main()
