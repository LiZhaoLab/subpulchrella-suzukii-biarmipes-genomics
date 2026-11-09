#!/usr/bin/env python3

import argparse
import bisect
from collections import defaultdict
from pathlib import Path


def parse_attributes(text):
    attributes = {}
    for field in text.rstrip().split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attributes[key] = value
    return attributes


def merge_intervals(intervals):
    merged = []

    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    return [(start, end) for start, end in merged]


def merge_sites_to_bed(sites):
    records = []
    by_contig = defaultdict(list)

    for contig, position in sites:
        by_contig[contig].append(position)

    for contig in sorted(by_contig):
        positions = sorted(set(by_contig[contig]))

        run_start = positions[0]
        previous = positions[0]

        for position in positions[1:]:
            if position != previous + 1:
                # Convert 1-based inclusive coordinates to BED.
                records.append((contig, run_start - 1, previous))
                run_start = position
            previous = position

        records.append((contig, run_start - 1, previous))

    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    mrna_info = {}
    transcript_exons = defaultdict(list)
    all_exons = defaultdict(list)

    with open(args.gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            contig, source, feature, start, end, score, strand, phase, attr_text = fields
            start = int(start)
            end = int(end)
            attributes = parse_attributes(attr_text)

            if feature == "mRNA":
                transcript_id = attributes.get("ID")
                if transcript_id:
                    mrna_info[transcript_id] = {
                        "contig": contig,
                        "strand": strand,
                        "gene": attributes.get("Parent", "NA"),
                    }

            elif feature == "exon":
                all_exons[contig].append((start, end))

                for parent in attributes.get("Parent", "").split(","):
                    if parent:
                        transcript_exons[parent].append((start, end))

    merged_exons = {}
    exon_starts = {}

    for contig, intervals in all_exons.items():
        merged = merge_intervals(intervals)
        merged_exons[contig] = merged
        exon_starts[contig] = [start for start, end in merged]

    def overlaps_any_exon(contig, position):
        intervals = merged_exons.get(contig, [])
        starts = exon_starts.get(contig, [])

        index = bisect.bisect_right(starts, position) - 1
        if index < 0:
            return False

        start, end = intervals[index]
        return start <= position <= end

    unique_introns = {}
    transcripts_with_multiple_exons = 0

    for transcript_id, info in mrna_info.items():
        exons = sorted(set(transcript_exons.get(transcript_id, [])))

        if len(exons) < 2:
            continue

        transcripts_with_multiple_exons += 1

        for left, right in zip(exons, exons[1:]):
            intron_start = left[1] + 1
            intron_end = right[0] - 1

            if intron_end < intron_start:
                continue

            intron_length = intron_end - intron_start + 1
            key = (
                info["contig"],
                intron_start,
                intron_end,
                info["strand"],
            )

            if key not in unique_introns:
                unique_introns[key] = set()

            unique_introns[key].add(transcript_id)

    short_introns = []
    raw_sites = set()

    for key, transcript_ids in unique_introns.items():
        contig, intron_start, intron_end, strand = key
        intron_length = intron_end - intron_start + 1

        if not 40 <= intron_length <= 65:
            continue

        if strand == "+":
            selected_start = intron_start + 7
            selected_end = intron_start + 29
        elif strand == "-":
            selected_start = intron_end - 29
            selected_end = intron_end - 7
        else:
            continue

        short_introns.append((
            contig,
            intron_start,
            intron_end,
            strand,
            intron_length,
            ",".join(sorted(transcript_ids)),
            selected_start,
            selected_end,
        ))

        for position in range(selected_start, selected_end + 1):
            raw_sites.add((contig, position))

    final_sites = {
        (contig, position)
        for contig, position in raw_sites
        if not overlaps_any_exon(contig, position)
    }

    bed_records = merge_sites_to_bed(final_sites)

    intron_path = Path(f"{output_prefix}.short_introns.tsv")
    bed_path = Path(f"{output_prefix}.short_intron_8to30.strict.bed")
    summary_path = Path(f"{output_prefix}.short_intron.summary.tsv")

    with open(intron_path, "w") as output:
        output.write(
            "contig\tintron_start_1based\tintron_end_1based\tstrand\t"
            "intron_length\ttranscript_ids\tselected_start_1based\t"
            "selected_end_1based\n"
        )

        for record in sorted(short_introns):
            output.write("\t".join(map(str, record)) + "\n")

    with open(bed_path, "w") as output:
        for contig, start, end in bed_records:
            output.write(f"{contig}\t{start}\t{end}\n")

    with open(summary_path, "w") as output:
        output.write("metric\tvalue\n")
        output.write(f"species\t{args.species}\n")
        output.write(f"mRNA_records\t{len(mrna_info)}\n")
        output.write(
            f"mRNA_with_multiple_exons\t{transcripts_with_multiple_exons}\n"
        )
        output.write(f"unique_introns\t{len(unique_introns)}\n")
        output.write(f"short_introns_40_to_65bp\t{len(short_introns)}\n")
        output.write(f"candidate_sites_before_exon_filter\t{len(raw_sites)}\n")
        output.write(
            f"sites_overlapping_any_annotated_exon\t"
            f"{len(raw_sites - final_sites)}\n"
        )
        output.write(f"final_short_intron_sites\t{len(final_sites)}\n")
        output.write(f"BED_intervals\t{len(bed_records)}\n")
        output.write("definition\tpositions_8_to_30_of_introns_40_to_65bp\n")
        output.write("coordinate_system\tBED_0_based_half_open\n")
        output.write("status\tPASS\n")

    print(f"species\t{args.species}")
    print(f"final_short_intron_sites\t{len(final_sites)}")
    print(f"bed\t{bed_path}")
    print(f"summary\t{summary_path}")
    print("FINAL: PASS")


if __name__ == "__main__":
    main()
