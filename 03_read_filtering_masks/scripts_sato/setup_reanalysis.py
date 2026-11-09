#!/usr/bin/env python3
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -J prep_RMmask
#SBATCH -o /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20261012_repeatmasker_p98/logs/prep_RMmask.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20261012_repeatmasker_p98/logs/prep_RMmask.%j.err

"""Prepare a separate RepeatMasker + indel, P5/P98 grenedalf rerun on node184."""

from __future__ import annotations

import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

ROOT = Path('/ru-auth/local/home/asato/scratch/population_polymorphism')
OLD = ROOT / 'reanalysis_20260821'
RUN = ROOT / 'reanalysis_20261012_repeatmasker_p98'
SYNTENY = Path('/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824')
SAMPLES = {'Dsub_L5_7': 'Dsub', 'JP-Kan': 'Dsuz', 'JP-Sap': 'Dsuz',
           'JP-Sap2': 'Dsuz', 'JP-Tok': 'Dsuz'}
CLASSES = ('autosome', 'X')


def require(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f'Missing or empty: {path}')
    return path


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def find_repeatmasker_output(species: str) -> Path:
    accession = ('GCF_014743375.2' if species == 'Dsub'
                 else 'GCF_043229965.1')
    return require(SYNTENY / 'repeatmasker' / species /
                   f'{species}.{accession}.fna.out')


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    result = []
    for line_number, line in enumerate(require(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith('#'):
            continue
        columns = line.split('\t')
        if len(columns) < 3:
            raise ValueError(f'Invalid BED: {path}:{line_number}')
        chrom, start, end = columns[:3]
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError(f'Invalid BED coordinates: {path}:{line_number}')
        result.append((chrom, start, end))
    return result


def merge(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    result = []
    for chrom, start, end in sorted(intervals):
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result


def save_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as stream:
        for chrom, start, end in intervals:
            stream.write(f'{chrom}\t{start}\t{end}\n')


def parse_repeatmasker(path: Path):
    """RepeatMasker .out query start/end are 1-based inclusive."""
    intervals = []
    classes = Counter()
    for line_number, line in enumerate(require(path).open(), 1):
        fields = line.split()
        if len(fields) < 11 or not fields[0].isdigit():
            continue  # RepeatMasker header lines
        chrom = fields[4]
        start, end = int(fields[5]) - 1, int(fields[6])
        if start < 0 or end <= start:
            raise ValueError(f'Invalid RepeatMasker position: {path}:{line_number}')
        intervals.append((chrom, start, end))
        classes[fields[10]] += 1
    if not intervals:
        raise ValueError(f'No RepeatMasker hits in {path}')
    return intervals, classes


def clip(intervals, regions):
    by_chrom = {}
    for chrom, start, end in regions:
        by_chrom.setdefault(chrom, []).append((start, end))
    clipped = []
    for chrom, start, end in intervals:
        for left, right in by_chrom.get(chrom, ()):
            a, b = max(start, left), min(end, right)
            if a < b:
                clipped.append((chrom, a, b))
    return merge(clipped)


def percentile_covered(histogram: Path, percentile: float) -> int:
    counts = {}
    with require(histogram).open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            depth, number = int(row['depth']), int(row['site_count'])
            if depth > 0:
                counts[depth] = counts.get(depth, 0) + number
    n = sum(counts.values())
    if not n:
        raise ValueError(f'No covered sites: {histogram}')
    target = math.ceil(percentile * n)
    cumulative = 0
    for depth, number in sorted(counts.items()):
        cumulative += number
        if cumulative >= target:
            return depth
    raise AssertionError('Unreachable percentile calculation')


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f'Expected one instance of {old!r}; found {count}')
    return text.replace(old, new)


def write_job(original: str, output: str, replacements: list[tuple[str, str]]) -> None:
    text = require(OLD / 'scripts' / original).read_text()
    for old, new in replacements:
        text = replace_once(text, old, new)
    path = RUN / 'scripts' / output
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def main() -> None:
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / 'logs').mkdir(exist_ok=True)
    old_thresholds = {}
    with require(OLD / 'config/depth_thresholds.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            old_thresholds[(row['sample'], row['analysis_class'])] = int(row['minimum_depth'])

    repeat_qc = []
    for species in ('Dsub', 'Dsuz'):
        reference = (Path('/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna.fai')
                     if species == 'Dsub' else
                     Path('/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna.fai'))
        lengths = {parts[0]: int(parts[1]) for line in require(reference).read_text().splitlines()
                   if len(parts := line.split('\t')) >= 2}
        source = find_repeatmasker_output(species)
        accession = ('GCF_014743375.2' if species == 'Dsub'
                     else 'GCF_043229965.1')
        input_fasta = require(SYNTENY / 'input_genomes' /
                              f'{species}.{accession}.fna')
        reference_fasta = require(Path(str(reference).removesuffix('.fai')))
        if checksum(input_fasta) != checksum(reference_fasta):
            raise ValueError(f'{species}: RepeatMasker FASTA differs from BAM reference FASTA')
        hits, categories = parse_repeatmasker(source)
        repeat_qc.append((species, 'source', str(source), checksum(source),
                          len(hits), 0, 0))
        invalid = {chrom for chrom, _, _ in hits if chrom not in lengths}
        overflow = [(chrom, end) for chrom, _, end in hits
                    if chrom in lengths and end > lengths[chrom]]
        if invalid or overflow:
            raise ValueError(f'{species}: reference mismatch; unknown={sorted(invalid)[:5]}, overflow={overflow[:5]}')
        for chrom_class in CLASSES:
            regions = read_bed(OLD / 'config/regions' / f'{species}.{chrom_class}.bed')
            repeats = clip(hits, regions)
            if not repeats:
                raise ValueError(f'{species} {chrom_class}: zero RepeatMasker overlaps')
            repeat_bed = RUN / 'masks/repeats' / f'{species}.{chrom_class}.repeatmasker.bed'
            save_bed(repeat_bed, repeats)
            total = sum(end - start for _, start, end in regions)
            masked = sum(end - start for _, start, end in repeats)
            repeat_qc.append((species, chrom_class, str(repeat_bed), '',
                              len(repeats), masked, round(100 * masked / total, 3)))
            for sample, sample_species in SAMPLES.items():
                if sample_species != species:
                    continue
                indel = OLD / 'indel_masks' / f'{sample}.{chrom_class}.indel_mask.support3.fraction0p05.flank5.bed'
                combined = merge(repeats + clip(read_bed(indel), regions))
                save_bed(RUN / 'masks/diversity' / f'{sample}.{chrom_class}.indel_repeatmasker.bed', combined)

        report = RUN / 'checks' / f'{species}.repeatmasker_categories.tsv'
        report.parent.mkdir(exist_ok=True)
        with report.open('w') as stream:
            stream.write('class_or_family\thit_count_before_merging\n')
            for category, count in sorted(categories.items()):
                stream.write(f'{category}\t{count}\n')

    thresholds = RUN / 'config/depth_thresholds.tsv'
    thresholds.parent.mkdir(exist_ok=True)
    with thresholds.open('w') as stream:
        stream.write('sample\tanalysis_class\tminimum_depth\tmaximum_depth\trule\n')
        for sample in SAMPLES:
            for chrom_class in CLASSES:
                p98 = percentile_covered(OLD / 'checks/depth_qc' /
                                         f'{sample}.{chrom_class}.depth_histogram.tsv', 0.98)
                minimum = old_thresholds[(sample, chrom_class)]
                if p98 <= minimum:
                    raise ValueError(f'{sample} {chrom_class}: P98 <= lower depth')
                stream.write(f'{sample}\t{chrom_class}\t{minimum}\t{p98}\tmax10_P05covered_to_P98covered\n')
    with (RUN / 'checks/repeatmasker_coverage.tsv').open('w') as stream:
        stream.write('species\tclass\tfile\tsha256\tintervals\tmasked_bp\tmasked_percent\n')
        for row in repeat_qc:
            stream.write('\t'.join(map(str, row)) + '\n')

    prefix = str(RUN)
    write_job('08_run_grenedalf_diversity_final_array.sh', '01_run_diversity_array.sh', [
        (f'#SBATCH -o {OLD}/logs/poly_grenedalf.%A_%a.out',
         f'#SBATCH -o {RUN}/logs/poly_grenedalf.%A_%a.out'),
        (f'#SBATCH -e {OLD}/logs/poly_grenedalf.%A_%a.err',
         f'#SBATCH -e {RUN}/logs/poly_grenedalf.%A_%a.err'),
        ('DEPTHTAB=${BASE}/config/depth_thresholds.tsv', f'DEPTHTAB={prefix}/config/depth_thresholds.tsv'),
        ('MASKDIR=${BASE}/indel_masks', f'MASKDIR={prefix}/masks/diversity'),
        ('MASK_BED=${MASKDIR}/${SAMPLE}.${CLASS}.indel_mask.support3.fraction0p05.flank5.bed',
         'MASK_BED=${MASKDIR}/${SAMPLE}.${CLASS}.indel_repeatmasker.bed'),
        ('OUTDIR=${BASE}/grenedalf_results/${CLASS}', f'OUTDIR={prefix}/grenedalf_results/${{CLASS}}'),
        ('_indelS3F0p05flank5\n', '_indelS3F0p05flank5_repeatmasker_P98\n'),
    ])
    write_job('10_build_Dsuz_Fst_sample_masks_array.sh', '02_build_fst_masks_array.sh', [
        ('DEPTHTAB=${BASE}/config/depth_thresholds.tsv', f'DEPTHTAB={prefix}/config/depth_thresholds.tsv'),
        ('DEPTH_MASK_DIR=${BASE}/fst_masks/depth_invalid', f'DEPTH_MASK_DIR={prefix}/masks/fst_depth_invalid'),
        ('FINAL_MASK_DIR=${BASE}/fst_masks/combined', f'FINAL_MASK_DIR={prefix}/masks/fst'),
        ('CHECKDIR=${BASE}/checks/fst_masks', f'CHECKDIR={prefix}/checks/fst_masks'),
        ('LOGDIR=${BASE}/logs', f'LOGDIR={prefix}/logs'),
        ('    cat "${INDEL_BED}"\n} |',
         f'    cat "${{INDEL_BED}}"\n    cat "{prefix}/masks/repeats/Dsuz.${{CLASS}}.repeatmasker.bed"\n}} |'),
    ])
    write_job('11_run_grenedalf_Dsuz_Fst_array.sh', '03_run_fst_array.sh', [
        ('MASKDIR=${BASE}/fst_masks/combined', f'MASKDIR={prefix}/masks/fst'),
        ('CONFIGDIR=${BASE}/config/fst', f'CONFIGDIR={prefix}/config/fst'),
        ('CHECKDIR=${BASE}/checks/fst', f'CHECKDIR={prefix}/checks/fst'),
        ('LOGDIR=${BASE}/logs', f'LOGDIR={prefix}/logs'),
        ('OUTDIR=${BASE}/fst_results/${CLASS}', f'OUTDIR={prefix}/fst_results/${{CLASS}}'),
        ('_unbiasedNei"', '_unbiasedNei_repeatmasker_P98"'),
    ])
    print(f'Prepared: {RUN}')
    print('Check RepeatMasker categories, masked fraction, and depth thresholds before sbatch.')
    print('FINAL: PASS')


if __name__ == '__main__':
    main()
