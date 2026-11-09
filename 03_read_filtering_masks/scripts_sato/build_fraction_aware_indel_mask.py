#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--analysis-class", required=True)
    parser.add_argument("--region-bed", required=True)
    parser.add_argument("--minimum-support", type=int, default=3)
    parser.add_argument("--minimum-fraction", type=float, default=0.05)
    parser.add_argument("--flank", type=int, default=5)
    parser.add_argument("--mask-bed", required=True)
    parser.add_argument("--event-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def read_regions(path):
    bounds = {}
    total_bp = 0

    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            chrom, start, end = line.rstrip().split("\t")[:3]
            start = int(start)
            end = int(end)

            bounds[chrom] = (start, end)
            total_bp += end - start

    if not bounds:
        raise RuntimeError(f"No regions found in {path}")

    return bounds, total_bp


def parse_pileup_bases(text):
    usable_depth = 0
    events = defaultdict(int)

    i = 0
    n = len(text)

    while i < n:
        char = text[i]

        if char == "^":
            i += 2
            continue

        if char == "$":
            i += 1
            continue

        if char in ".,ACGTNacgtn":
            usable_depth += 1
            i += 1

            if i < n and text[i] in "+-":
                sign = text[i]
                i += 1

                number_start = i
                while i < n and text[i].isdigit():
                    i += 1

                if number_start == i:
                    raise RuntimeError(
                        f"Malformed indel length in pileup string: {text}"
                    )

                length = int(text[number_start:i])
                sequence = text[i:i + length].upper()

                if len(sequence) != length:
                    raise RuntimeError(
                        f"Truncated indel sequence in pileup string: {text}"
                    )

                indel_type = "insertion" if sign == "+" else "deletion"
                events[(indel_type, length, sequence)] += 1
                i += length

            continue

        if char in "*#<>":
            i += 1
            continue

        raise RuntimeError(
            f"Unrecognized pileup character {char!r} in {text!r}"
        )

    return usable_depth, events


def merge_intervals(intervals):
    merged = []

    for chrom in sorted(intervals):
        ordered = sorted(intervals[chrom])

        if not ordered:
            continue

        current_start, current_end = ordered[0]

        for start, end in ordered[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((chrom, current_start, current_end))
                current_start, current_end = start, end

        merged.append((chrom, current_start, current_end))

    return merged


def main():
    args = parse_args()
    bounds, total_region_bp = read_regions(args.region_bed)

    intervals = defaultdict(list)
    event_rows = []

    pileup_sites = 0
    sites_with_indel = 0
    total_indel_alleles = 0
    passing_indel_alleles = 0

    for line in sys.stdin:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 5:
            raise RuntimeError(f"Unexpected mpileup row: {line.rstrip()}")

        chrom = fields[0]
        position_1based = int(fields[1])
        pileup_bases = fields[4]

        if chrom not in bounds:
            continue

        pileup_sites += 1
        usable_depth, events = parse_pileup_bases(pileup_bases)

        if not events:
            continue

        sites_with_indel += 1

        for (indel_type, length, sequence), support in events.items():
            total_indel_alleles += 1

            fraction = (
                support / usable_depth
                if usable_depth > 0
                else 0.0
            )

            passed = (
                usable_depth > 0
                and support >= args.minimum_support
                and fraction >= args.minimum_fraction
            )

            if not passed:
                continue

            passing_indel_alleles += 1

            anchor_0based = position_1based - 1
            region_start, region_end = bounds[chrom]

            if indel_type == "insertion":
                event_start = anchor_0based
                event_end = anchor_0based + 1
            else:
                # Anchor base plus deleted reference interval
                event_start = anchor_0based
                event_end = anchor_0based + length + 1

            mask_start = max(region_start, event_start - args.flank)
            mask_end = min(region_end, event_end + args.flank)

            intervals[chrom].append((mask_start, mask_end))

            event_rows.append({
                "sample": args.sample,
                "species": args.species,
                "analysis_class": args.analysis_class,
                "chromosome": chrom,
                "anchor_position_1based": position_1based,
                "indel_type": indel_type,
                "indel_length_bp": length,
                "indel_sequence": sequence,
                "usable_depth": usable_depth,
                "supporting_reads": support,
                "supporting_fraction": f"{fraction:.8f}",
                "mask_start_0based": mask_start,
                "mask_end_0based": mask_end,
            })

    merged = merge_intervals(intervals)
    masked_bp = sum(end - start for _, start, end in merged)
    masked_percent = 100 * masked_bp / total_region_bp

    with open(args.mask_bed, "w") as handle:
        for chrom, start, end in merged:
            handle.write(f"{chrom}\t{start}\t{end}\n")

    event_fields = [
        "sample",
        "species",
        "analysis_class",
        "chromosome",
        "anchor_position_1based",
        "indel_type",
        "indel_length_bp",
        "indel_sequence",
        "usable_depth",
        "supporting_reads",
        "supporting_fraction",
        "mask_start_0based",
        "mask_end_0based",
    ]

    with open(args.event_output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=event_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(event_rows)

    summary = [
        ("sample", args.sample),
        ("species", args.species),
        ("analysis_class", args.analysis_class),
        ("pileup_sites", pileup_sites),
        ("sites_with_any_indel", sites_with_indel),
        ("all_indel_alleles", total_indel_alleles),
        ("minimum_supporting_reads", args.minimum_support),
        ("minimum_supporting_fraction", args.minimum_fraction),
        ("flank_bp", args.flank),
        ("passing_indel_alleles", passing_indel_alleles),
        ("merged_mask_intervals", len(merged)),
        ("total_region_bp", total_region_bp),
        ("masked_bp", masked_bp),
        ("masked_percent", f"{masked_percent:.8f}"),
    ]

    with open(args.summary_output, "w") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in summary:
            handle.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
