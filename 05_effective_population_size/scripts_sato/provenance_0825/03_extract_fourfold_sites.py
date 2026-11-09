#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq


FOURFOLD_PREFIXES = {
    "CC": "Pro",
    "AC": "Thr",
    "GT": "Val",
    "GC": "Ala",
    "GG": "Gly",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_attributes(text):
    result = {}

    for field in text.rstrip(";").split(";"):
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = unquote(value)

    return result


def load_fasta(path):
    with open_text(path) as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def parse_gff3(path):
    transcripts = {}
    cds_by_transcript = defaultdict(list)

    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attributes = fields
            attrs = parse_attributes(attributes)

            start = int(start)
            end = int(end)

            if feature == "mRNA":
                transcript_id = attrs.get("ID")

                if not transcript_id:
                    continue

                parents = attrs.get("Parent", "").split(",")
                gene_id = parents[0] if parents and parents[0] else transcript_id

                transcripts[transcript_id] = {
                    "gene_id": gene_id,
                    "chrom": chrom,
                    "strand": strand,
                    "partial": attrs.get("partial", "").lower() == "true",
                    "pseudo": (
                        "pseudo" in attrs or
                        attrs.get("pseudogene", "").lower() == "true"
                    ),
                }

            elif feature == "CDS":
                parents = attrs.get("Parent", "").split(",")

                if not parents or not parents[0]:
                    continue

                try:
                    phase_value = int(phase)
                except ValueError:
                    phase_value = None

                for transcript_id in parents:
                    cds_by_transcript[transcript_id].append({
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": phase_value,
                        "partial": attrs.get("partial", "").lower() == "true",
                    })

    return transcripts, cds_by_transcript


def build_cds(transcript_id, transcript, features, genome):
    chrom = transcript["chrom"]
    strand = transcript["strand"]

    if chrom not in genome:
        return None, "contig_missing"

    if transcript["partial"] or any(x["partial"] for x in features):
        return None, "partial"

    if transcript["pseudo"]:
        return None, "pseudo"

    if strand not in {"+", "-"}:
        return None, "invalid_strand"

    if any(x["chrom"] != chrom or x["strand"] != strand for x in features):
        return None, "inconsistent_location"

    if any(x["phase"] is None for x in features):
        return None, "missing_phase"

    if strand == "+":
        ordered = sorted(features, key=lambda x: (x["start"], x["end"]))
    else:
        ordered = sorted(
            features,
            key=lambda x: (x["start"], x["end"]),
            reverse=True
        )

    sequence_parts = []
    positions = []
    cumulative_length = 0

    for index, feature in enumerate(ordered):
        expected_phase = 0 if index == 0 else (3 - cumulative_length % 3) % 3

        if feature["phase"] != expected_phase:
            return None, "phase_inconsistent"

        start = feature["start"]
        end = feature["end"]

        genomic_sequence = str(genome[chrom].seq[start - 1:end]).upper()

        if strand == "+":
            oriented_sequence = genomic_sequence
            oriented_positions = list(range(start - 1, end))
        else:
            oriented_sequence = str(Seq(genomic_sequence).reverse_complement())
            oriented_positions = list(range(end - 1, start - 2, -1))

        sequence_parts.append(oriented_sequence)
        positions.extend(oriented_positions)
        cumulative_length += end - start + 1

    cds_sequence = "".join(sequence_parts)

    if len(cds_sequence) == 0:
        return None, "empty_CDS"

    if len(cds_sequence) % 3 != 0:
        return None, "length_not_multiple_of_3"

    if len(positions) != len(cds_sequence):
        return None, "position_length_mismatch"

    if any(base not in {"A", "C", "G", "T"} for base in cds_sequence):
        return None, "ambiguous_base"

    codons = [
        cds_sequence[index:index + 3]
        for index in range(0, len(cds_sequence), 3)
    ]

    internal_codons = codons[:-1] if codons[-1] in STOP_CODONS else codons

    if any(codon in STOP_CODONS for codon in internal_codons):
        return None, "internal_stop"

    return {
        "transcript_id": transcript_id,
        "gene_id": transcript["gene_id"],
        "chrom": chrom,
        "strand": strand,
        "sequence": cds_sequence,
        "positions": positions,
        "length": len(cds_sequence),
        "codons": codons,
    }, "PASS"


def select_longest_transcripts(transcripts, cds_by_transcript, genome):
    valid_transcripts = []
    rejection_counts = defaultdict(int)

    for transcript_id, features in cds_by_transcript.items():
        transcript = transcripts.get(transcript_id)

        if transcript is None:
            rejection_counts["mRNA_record_missing"] += 1
            continue

        result, status = build_cds(
            transcript_id,
            transcript,
            features,
            genome
        )

        if result is None:
            rejection_counts[status] += 1
        else:
            valid_transcripts.append(result)

    by_gene = defaultdict(list)

    for transcript in valid_transcripts:
        by_gene[transcript["gene_id"]].append(transcript)

    selected = []

    for gene_id, candidates in by_gene.items():
        candidates.sort(
            key=lambda x: (-x["length"], x["transcript_id"])
        )
        selected.append(candidates[0])

    selected.sort(
        key=lambda x: (
            x["chrom"],
            min(x["positions"]),
            x["gene_id"],
            x["transcript_id"],
        )
    )

    return selected, valid_transcripts, rejection_counts


def identify_fourfold_sites(selected):
    candidate_classes = defaultdict(set)
    candidate_transcripts = defaultdict(set)

    for transcript in selected:
        sequence = transcript["sequence"]
        positions = transcript["positions"]
        chrom = transcript["chrom"]

        for index in range(0, len(sequence), 3):
            codon = sequence[index:index + 3]

            if len(codon) != 3:
                continue

            prefix = codon[:2]

            if prefix not in FOURFOLD_PREFIXES:
                continue

            position = positions[index + 2]
            key = (chrom, position)

            candidate_classes[key].add(FOURFOLD_PREFIXES[prefix])
            candidate_transcripts[key].add(transcript["transcript_id"])

    candidates = set(candidate_classes)
    conflicting_sites = set()

    # A candidate is removed if the same genomic position occurs in another
    # selected CDS but is not a four-fold third-codon position there.
    for transcript in selected:
        sequence = transcript["sequence"]
        positions = transcript["positions"]
        chrom = transcript["chrom"]

        local_fourfold = set()

        for index in range(0, len(sequence), 3):
            codon = sequence[index:index + 3]

            if len(codon) == 3 and codon[:2] in FOURFOLD_PREFIXES:
                local_fourfold.add((chrom, positions[index + 2]))

        for position in positions:
            key = (chrom, position)

            if key in candidates and key not in local_fourfold:
                conflicting_sites.add(key)

    retained = sorted(
        candidates - conflicting_sites,
        key=lambda x: (x[0], x[1])
    )

    return retained, conflicting_sites, candidate_classes, candidate_transcripts


def sites_to_bed_intervals(sites):
    if not sites:
        return []

    intervals = []
    current_chrom, current_start = sites[0]
    current_end = current_start + 1

    for chrom, position in sites[1:]:
        if chrom == current_chrom and position == current_end:
            current_end += 1
        else:
            intervals.append((current_chrom, current_start, current_end))
            current_chrom = chrom
            current_start = position
            current_end = position + 1

    intervals.append((current_chrom, current_start, current_end))
    return intervals


def write_outputs(
    species,
    selected,
    valid_transcripts,
    rejection_counts,
    retained_sites,
    conflicting_sites,
    candidate_classes,
    candidate_transcripts,
    output_bed,
    output_sites,
    output_transcripts,
    output_summary,
):
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    intervals = sites_to_bed_intervals(retained_sites)

    with open(output_bed, "w") as handle:
        for chrom, start, end in intervals:
            handle.write(f"{chrom}\t{start}\t{end}\n")

    with open(output_sites, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "chrom",
            "position_1based",
            "BED_start_0based",
            "TE_codon_family",
            "supporting_transcripts",
        ])

        for chrom, position in retained_sites:
            key = (chrom, position)
            writer.writerow([
                chrom,
                position + 1,
                position,
                ",".join(sorted(candidate_classes[key])),
                ",".join(sorted(candidate_transcripts[key])),
            ])

    with open(output_transcripts, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "gene_id",
            "selected_transcript_id",
            "chrom",
            "strand",
            "CDS_length_bp",
        ])

        for transcript in selected:
            writer.writerow([
                transcript["gene_id"],
                transcript["transcript_id"],
                transcript["chrom"],
                transcript["strand"],
                transcript["length"],
            ])

    summary = [
        ("species", species),
        ("CDS_parent_transcripts_in_GFF3",
         len(valid_transcripts) + sum(rejection_counts.values())),
        ("valid_complete_coding_transcripts", len(valid_transcripts)),
        ("genes_with_valid_transcript", len(selected)),
        ("selected_longest_transcripts", len(selected)),
        ("fourfold_candidates_before_overlap_filter",
         len(retained_sites) + len(conflicting_sites)),
        ("sites_conflicting_with_other_selected_CDS",
         len(conflicting_sites)),
        ("final_fourfold_sites", len(retained_sites)),
        ("BED_intervals", len(intervals)),
        ("codon_families", "CCN,ACN,GTN,GCN,GGN"),
        ("transcript_selection",
         "longest_valid_complete_CDS_per_gene"),
        ("coordinate_system", "BED_0_based_half_open"),
    ]

    for reason in sorted(rejection_counts):
        summary.append((f"rejected_{reason}", rejection_counts[reason]))

    summary.append(("status", "PASS"))

    with open(output_summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)


def main():
    parser = argparse.ArgumentParser(
        description="Extract conservative four-fold degenerate sites."
    )

    parser.add_argument("--species", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--output-prefix", required=True)

    args = parser.parse_args()

    prefix = Path(args.output_prefix)

    print("Loading reference FASTA...", file=sys.stderr)
    genome = load_fasta(args.fasta)

    print("Parsing GFF3...", file=sys.stderr)
    transcripts, cds_by_transcript = parse_gff3(args.gff3)

    print("Selecting longest valid CDS per gene...", file=sys.stderr)
    selected, valid_transcripts, rejection_counts = \
        select_longest_transcripts(
            transcripts,
            cds_by_transcript,
            genome
        )

    print("Identifying four-fold sites...", file=sys.stderr)
    (
        retained_sites,
        conflicting_sites,
        candidate_classes,
        candidate_transcripts,
    ) = identify_fourfold_sites(selected)

    write_outputs(
        species=args.species,
        selected=selected,
        valid_transcripts=valid_transcripts,
        rejection_counts=rejection_counts,
        retained_sites=retained_sites,
        conflicting_sites=conflicting_sites,
        candidate_classes=candidate_classes,
        candidate_transcripts=candidate_transcripts,
        output_bed=Path(f"{prefix}.fourfold_5families.strict.bed"),
        output_sites=Path(f"{prefix}.fourfold_5families.strict.sites.tsv"),
        output_transcripts=Path(f"{prefix}.selected_longest_CDS.tsv"),
        output_summary=Path(f"{prefix}.fourfold_5families.summary.tsv"),
    )

    print(f"species\t{args.species}")
    print(f"selected_longest_transcripts\t{len(selected)}")
    print(f"final_fourfold_sites\t{len(retained_sites)}")
    print("FINAL: PASS")


if __name__ == "__main__":
    main()
