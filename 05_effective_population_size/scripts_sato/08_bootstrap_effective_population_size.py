#!/usr/bin/env python3
"""Estimate Ne from the 2026-10-23 P98 neutral-site diversity results."""

import csv
import math
import random
from collections import defaultdict
from pathlib import Path


BASE = Path("/ru-auth/local/home/asato/scratch/population_polymorphism/Ne_reanalysis_20261023_repeatmasker_p98")
MUTATION_FILE = BASE / "config/mutation_rates.tsv"
OUTPUT = BASE / "results/effective_population_size/effective_population_size_block_bootstrap.tsv"

SAMPLES = ("Dsub_L5_7", "JP-Kan", "JP-Sap", "JP-Sap2", "JP-Tok")
CLASSES = ("X", "autosome")
SITE_TYPES = ("fourfold", "short_intron")
BOOTSTRAP_REPLICATES = 10000
BLOCK_SIZE_BP = 1_000_000
RANDOM_SEED = 20260826


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def percentile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def input_path(sample, analysis_class, site_type):
    if site_type == "short_intron":
        root = BASE / "grenedalf_neutral_len40to65" / analysis_class / site_type
        label = f"{sample}_{analysis_class}_short_intron40to65_20kb_"
    else:
        root = BASE / "grenedalf_neutral" / analysis_class / site_type
        label = f"{sample}_{analysis_class}_fourfold_20kb_"

    matches = sorted(root.glob(f"diversity-{label}*repeatmasker_P98.csv"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one input for {sample} {analysis_class} "
            f"{site_type}; found {len(matches)}: {matches}"
        )
    return matches[0]


def load_mutation_rates():
    with MUTATION_FILE.open(newline="") as handle:
        rates = {
            row["rate_id"]: float(row["mutation_rate_per_site_per_generation"])
            for row in csv.DictReader(handle, delimiter="\t")
        }

    required = {"primary", "lower_mu", "upper_mu"}
    if not required.issubset(rates):
        raise RuntimeError(f"Missing mutation rates: {required - set(rates)}")

    if not (0 < rates["lower_mu"] < rates["primary"] < rates["upper_mu"]):
        raise RuntimeError(f"Unexpected mutation rate order: {rates}")

    return rates


def load_blocks(path):
    blocks = defaultdict(
        lambda: {
            "pi_num": 0.0,
            "w_num": 0.0,
            "valid": 0,
            "polymorphic": 0,
            "windows": 0,
        }
    )

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise RuntimeError(f"Missing header: {path}")

        pi_cols = [
            name for name in reader.fieldnames
            if name.endswith(".theta_pi")
        ]
        w_cols = [
            name for name in reader.fieldnames
            if name.endswith(".theta_watterson")
        ]
        if len(pi_cols) != 1 or len(w_cols) != 1:
            raise RuntimeError(f"Unexpected theta columns: {path}")

        pi_col, w_col = pi_cols[0], w_cols[0]

        for row in reader:
            start = int(row["start"])
            invariant = int(row["total.invariant"])
            passed = int(row["total.passed"])
            valid = invariant + passed

            if start < 1 or invariant < 0 or passed < 0:
                raise RuntimeError(f"Invalid window counts: {path}: {row}")

            block = blocks[
                (row["chrom"], (start - 1) // BLOCK_SIZE_BP)
            ]
            block["windows"] += 1
            block["polymorphic"] += passed

            if valid > 0:
                pi = finite_number(row[pi_col])
                w = finite_number(row[w_col])

                if pi is None or w is None:
                    raise RuntimeError(
                        f"Non-finite theta at callable sites: "
                        f"{path}: {row['chrom']}:{start}"
                    )

                block["valid"] += valid
                block["pi_num"] += pi * valid
                block["w_num"] += w * valid

    if not blocks or sum(block["valid"] for block in blocks.values()) == 0:
        raise RuntimeError(f"No callable neutral sites: {path}")

    return blocks


def weighted_theta(blocks, numerator):
    return (
        sum(block[numerator] for block in blocks)
        / sum(block["valid"] for block in blocks)
    )


def bootstrap_w(blocks, rng):
    by_chromosome = defaultdict(list)

    for (chrom, _), block in blocks.items():
        if block["valid"] > 0:
            by_chromosome[chrom].append(block)

    estimates = []

    for _ in range(BOOTSTRAP_REPLICATES):
        numerator = 0.0
        denominator = 0.0

        for chrom in sorted(by_chromosome):
            sampled = rng.choices(
                by_chromosome[chrom],
                k=len(by_chromosome[chrom]),
            )
            numerator += sum(block["w_num"] for block in sampled)
            denominator += sum(block["valid"] for block in sampled)

        estimates.append(numerator / denominator)

    return (
        percentile(estimates, 0.025),
        percentile(estimates, 0.975),
    )


def summarize(sample, analysis_class, site_type, rates, rng):
    path = input_path(sample, analysis_class, site_type)
    blocks = load_blocks(path)
    values = list(blocks.values())

    theta_pi = weighted_theta(values, "pi_num")
    theta_w = weighted_theta(values, "w_num")
    theta_low, theta_high = bootstrap_w(blocks, rng)

    factor = 4.0 if analysis_class == "autosome" else 3.0
    primary_mu = rates["primary"]

    return {
        "sample": sample,
        "species": "Dsub" if sample == "Dsub_L5_7" else "Dsuz",
        "analysis_class": analysis_class,
        "site_type": site_type,
        "block_size_bp": BLOCK_SIZE_BP,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "genomic_blocks": len(blocks),
        "windows": sum(block["windows"] for block in values),
        "valid_loci": sum(block["valid"] for block in values),
        "polymorphic_loci": sum(block["polymorphic"] for block in values),
        "theta_pi_per_site": theta_pi,
        "theta_watterson_per_site": theta_w,
        "theta_watterson_bootstrap_95CI_low": theta_low,
        "theta_watterson_bootstrap_95CI_high": theta_high,
        "chromosome_factor": factor,
        "primary_mutation_rate": primary_mu,
        "Ne_primary": theta_w / (factor * primary_mu),
        "Ne_primary_bootstrap_95CI_low": (
            theta_low / (factor * primary_mu)
        ),
        "Ne_primary_bootstrap_95CI_high": (
            theta_high / (factor * primary_mu)
        ),
        "lower_mutation_rate": rates["lower_mu"],
        "upper_mutation_rate": rates["upper_mu"],
        "Ne_mutation_rate_sensitivity_min": (
            theta_w / (factor * rates["upper_mu"])
        ),
        "Ne_mutation_rate_sensitivity_max": (
            theta_w / (factor * rates["lower_mu"])
        ),
        "source_file": str(path),
    }


def main():
    if OUTPUT.exists():
        raise RuntimeError(
            f"Output already exists; will not overwrite: {OUTPUT}"
        )

    rates = load_mutation_rates()
    rng = random.Random(RANDOM_SEED)
    rows = []

    # Match the old script's input order for reproducible bootstrap draws.
    for analysis_class in CLASSES:
        for site_type in SITE_TYPES:
            for sample in SAMPLES:
                rows.append(
                    summarize(
                        sample, analysis_class, site_type, rates, rng
                    )
                )

    if len(rows) != 20:
        raise RuntimeError(f"Expected 20 results; found {len(rows)}")

    rows.sort(
        key=lambda row: (
            row["species"],
            row["sample"],
            row["analysis_class"],
            row["site_type"],
        )
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"input_files\t{len(rows)}")
    print(f"output_rows\t{len(rows)}")
    print(f"bootstrap_replicates\t{BOOTSTRAP_REPLICATES}")
    print(f"block_size_bp\t{BLOCK_SIZE_BP}")
    print(f"random_seed\t{RANDOM_SEED}")
    print(f"mutation_file\t{MUTATION_FILE}")
    print(f"output\t{OUTPUT}")
    print("FINAL: PASS")


if __name__ == "__main__":
    main()
