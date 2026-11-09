#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


WINDOW_BP = 200_000

SYNTENY_BASE = Path(
    "/ru-auth/local/home/asato/scratch/"
    "synteny_4species/reanalysis_20260824"
)

OUTPUT_DIR = SYNTENY_BASE / "results" / "TE_gene_tracks_200kb"
OUTPUT_TABLE = OUTPUT_DIR / "Dsuz_Dsub.TE_gene_density.200kb.tsv"
OUTPUT_SUMMARY = OUTPUT_DIR / "Dsuz_Dsub.TE_gene_density.200kb.summary.tsv"

SPECIES = {
    "Dsuz": {
        "fai": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dsuz_IsoJpt1.0/"
            "GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna.fai"
        ),
        "gene_gff": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dsuz_IsoJpt1.0/"
            "GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff"
        ),
        "te_gff": Path(
            "/ru-auth/local/home/asato/scratch/TE_reanalysis/"
            "reanalysis_20260824/results/Dsuz/Dsuz_EarlGrey/"
            "Dsuz_summaryFiles/Dsuz.filteredRepeats.gff"
        ),
        "chromosomes": [
            ("NC_092084.1", "X", "X"),
            ("NC_092080.1", "2L", "2L"),
            ("NC_092081.1", "2R", "2R"),
            ("NC_092082.1", "3", "3"),
            ("NC_092083.1", "4", "4"),
        ],
    },
    "Dsub": {
        "fai": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dsub_v1.1/"
            "GCF_014743375.2_RU_Dsub_v1.1_genomic.fna.fai"
        ),
        "gene_gff": Path(
            "/ru-auth/local/home/asato/scratch/genome/"
            "Dsub_v1.1/"
            "GCF_014743375.2_RU_Dsub_v1.1_genomic.gff"
        ),
        "te_gff": Path(
            "/ru-auth/local/home/asato/scratch/TE_reanalysis/"
            "reanalysis_20260824/results/Dsub/Dsub_EarlGrey/"
            "Dsub_summaryFiles/Dsub.filteredRepeats.gff"
        ),
        "chromosomes": [
            ("NC_050613.1", "X", "X"),
            ("NC_050610.1", "2L", "2L"),
            ("NC_050611.1", "2R", "2R"),
            ("NC_050612.1", "3", "3L"),
            ("NW_023665577.1", "3", "3R-2"),
            ("NC_050609.1", "3", "3R-1"),
            ("NC_050608.1", "4", "4"),
        ],
    },
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")


def read_fai(path: Path) -> dict[str, int]:
    lengths = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            lengths[fields[0]] = int(fields[1])
    return lengths


def is_te_feature(feature: str) -> bool:
    if feature == "Unknown":
        return True
    return feature.startswith(
        ("DNA", "LINE", "LTR", "RC", "SINE", "PLE", "Penelope")
    )


def read_te_intervals(
    path: Path,
    retained_chromosomes: set[str],
) -> tuple[dict[str, list[tuple[int, int]]], int]:
    intervals = defaultdict(list)
    retained_rows = 0

    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            feature = fields[2]
            if chrom not in retained_chromosomes or not is_te_feature(feature):
                continue

            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            if start0 < 0 or end0 <= start0:
                continue

            intervals[chrom].append((start0, end0))
            retained_rows += 1

    return intervals, retained_rows


def read_gene_midpoints(
    path: Path,
    retained_chromosomes: set[str],
) -> tuple[dict[str, list[int]], int]:
    midpoints = defaultdict(list)
    retained_rows = 0

    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            if chrom not in retained_chromosomes or fields[2] != "gene":
                continue

            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            if start0 < 0 or end0 <= start0:
                continue

            midpoints[chrom].append((start0 + end0) // 2)
            retained_rows += 1

    return midpoints, retained_rows


def merge_intervals(
    intervals: list[tuple[int, int]],
    chromosome_length: int,
) -> list[tuple[int, int]]:
    clipped = [
        (max(0, start), min(chromosome_length, end))
        for start, end in intervals
        if start < chromosome_length and end > 0
    ]
    clipped = [(start, end) for start, end in clipped if end > start]
    if not clipped:
        return []

    clipped.sort()
    merged = [list(clipped[0])]
    for start, end in clipped[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def overlap_bp(
    merged_intervals: list[tuple[int, int]],
    window_start: int,
    window_end: int,
) -> int:
    total = 0
    for start, end in merged_intervals:
        if end <= window_start:
            continue
        if start >= window_end:
            break
        total += max(0, min(end, window_end) - max(start, window_start))
    return total


def midpoint_counts(
    midpoints: list[int],
    chromosome_length: int,
) -> list[int]:
    number_of_windows = (chromosome_length + WINDOW_BP - 1) // WINDOW_BP
    counts = [0] * number_of_windows
    for midpoint in midpoints:
        if 0 <= midpoint < chromosome_length:
            counts[midpoint // WINDOW_BP] += 1
    return counts


def build_species_rows(species: str, config: dict) -> tuple[list[dict], dict]:
    for key in ("fai", "gene_gff", "te_gff"):
        require_file(config[key])

    lengths = read_fai(config["fai"])
    retained_chromosomes = {row[0] for row in config["chromosomes"]}
    missing_lengths = sorted(retained_chromosomes - set(lengths))
    if missing_lengths:
        raise ValueError(
            f"{species}: chromosomes absent from FAI: {missing_lengths}"
        )

    te_intervals, te_rows = read_te_intervals(
        config["te_gff"], retained_chromosomes
    )
    gene_midpoints, gene_rows = read_gene_midpoints(
        config["gene_gff"], retained_chromosomes
    )

    output_rows = []
    total_covered_bp = 0
    total_windows = 0

    for chromosome_order, (chrom, arm, label) in enumerate(
        config["chromosomes"], start=1
    ):
        chromosome_length = lengths[chrom]
        merged_te = merge_intervals(te_intervals[chrom], chromosome_length)
        te_midpoint_count = midpoint_counts(
            [
                (start + end) // 2
                for start, end in te_intervals[chrom]
            ],
            chromosome_length,
        )
        gene_count = midpoint_counts(
            gene_midpoints[chrom], chromosome_length
        )
        number_of_windows = len(gene_count)

        for window_index in range(number_of_windows):
            start = window_index * WINDOW_BP
            end = min(start + WINDOW_BP, chromosome_length)
            window_bp = end - start
            covered_bp = overlap_bp(merged_te, start, end)
            total_covered_bp += covered_bp
            total_windows += 1

            output_rows.append({
                "species": species,
                "chromosome": chrom,
                "arm": arm,
                "chromosome_label": label,
                "chromosome_order": chromosome_order,
                "chromosome_length_bp": chromosome_length,
                "window_index_1based": window_index + 1,
                "window_start_0based": start,
                "window_end": end,
                "window_midpoint": (start + end) / 2,
                "window_bp": window_bp,
                "TE_loci": te_midpoint_count[window_index],
                "TE_covered_bp": covered_bp,
                "TE_coverage_percent": 100 * covered_bp / window_bp,
                "gene_count": gene_count[window_index],
                "genes_per_Mb": gene_count[window_index] * 1_000_000 / window_bp,
            })

    summary = {
        "species": species,
        "window_size_bp": WINDOW_BP,
        "major_chromosomes": len(config["chromosomes"]),
        "output_windows": total_windows,
        "retained_TE_GFF_rows": te_rows,
        "retained_gene_GFF_rows": gene_rows,
        "unique_TE_covered_bp": total_covered_bp,
        "TE_GFF": str(config["te_gff"]),
        "gene_GFF": str(config["gene_gff"]),
        "FAI": str(config["fai"]),
    }
    return output_rows, summary


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_rows = []
    summaries = []
    for species, config in SPECIES.items():
        rows, summary = build_species_rows(species, config)
        all_rows.extend(rows)
        summaries.append(summary)

    fieldnames = [
        "species",
        "chromosome",
        "arm",
        "chromosome_label",
        "chromosome_order",
        "chromosome_length_bp",
        "window_index_1based",
        "window_start_0based",
        "window_end",
        "window_midpoint",
        "window_bp",
        "TE_loci",
        "TE_covered_bp",
        "TE_coverage_percent",
        "gene_count",
        "genes_per_Mb",
    ]

    with OUTPUT_TABLE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    summary_fields = list(summaries[0])
    with OUTPUT_SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(summaries)

    print(f"output\t{OUTPUT_TABLE}")
    print(f"summary\t{OUTPUT_SUMMARY}")
    print(f"rows\t{len(all_rows)}")
    print("FINAL: PASS")


if __name__ == "__main__":
    main()
