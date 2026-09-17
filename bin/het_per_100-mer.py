#!/usr/bin/env python3
"""Summarise heterozygous sites across the legacy fixed FA1090 bins.

Reference-mapped branch only — BIN_SIZE/BIN_COUNT are sized specifically to
the FA1090 reference genome's length, and CHROM is intentionally discarded
(safe only because FA1090 is a single-chromosome reference). Do not reuse
this for assembly-mapped VCFs, which are multi-contig and vary in length
per sample.
"""

import argparse


BIN_SIZE = 21_848
BIN_COUNT = 100


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="CHROM/POS/REF/ALT/AF table from bcftools query")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def create_bins():
    return {index * BIN_SIZE: 0 for index in range(1, BIN_COUNT + 1)}


def fill_bins(input_path, bins):
    with open(input_path, encoding="utf-8") as handle:
        for line in handle:
            if "#" in line or not line.strip():
                continue
            position = int(line.split()[1])
            upper_bound = ((position // BIN_SIZE) + 1) * BIN_SIZE
            bins[upper_bound] += 1
    return bins


def count_bins_by_site_count(bins):
    return {
        site_count: sum(value == site_count for value in bins.values())
        for site_count in set(bins.values())
    }


def main():
    args = parse_args()
    bins = fill_bins(args.input, create_bins())
    distribution = count_bins_by_site_count(bins)

    with open(args.output, "w", encoding="utf-8") as output:
        for site_count, bin_count in sorted(distribution.items()):
            output.write(f"{args.sample_id}\t{site_count}\t{bin_count}\n")


if __name__ == "__main__":
    main()
