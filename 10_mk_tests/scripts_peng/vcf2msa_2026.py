#!/usr/bin/env python3
"""
vcf2msa.py

Convert VCF SNPs to per-CDS MSAs by applying SNP genotypes onto CDS sequences
whose FASTA headers encode genomic coordinates (possibly multi-exon), e.g.:

>NC_050608.1-:64988-65192;NC_050608.1-:65249-65385;...
ATG...

Behavior:
- Uses *real* coordinates directly (no coordinate liftover/mapping).
- Extracts **biallelic SNPs only** (skips indels and multi-ALT).
- Generates 2 * Nsamples haplotypes (diploid) per CDS, plus a reference row.
- Writes an index CSV aligned to input FASTA order:
    index,header,SNP
  where SNP is 1 if CDS has ≥1 overlapping SNP else 0.
- To save time/space, writes MSAs **only** for CDS with SNPs.
  Output filenames keep the original FASTA index: msa_<index>.fasta

Usage:
  python vcf2msa.py -vcf in.vcf.gz -fasta cds.fa -out outdir -species myspp

Notes:
- VCF POS is 1-based (standard). CDS header coordinates are assumed 1-based too.
- Missing genotype (./.) -> reference allele for both haplotypes.
- Allele index > 1 (unexpected/multiallelic encoding) -> reference allele (conservative).
"""

import os
import sys
import gzip
import argparse
import time
from multiprocessing import Pool
from collections import Counter, defaultdict
import bisect

WCBP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

WORKER_RECORDS = None
WORKER_PARSED_HEADERS = None
WORKER_SNP_IN_RECORD = None
WORKER_NPOPULATIONS = None
WORKER_OUTDIR = None
WORKER_SPECIES = None
WORKER_RESUME = None


def format_seconds(seconds):
    seconds = int(round(seconds))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m {secs}s"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def readlines(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_cds_fasta_records(fasta_path):
    """
    Minimal FASTA reader returning ordered [(header_key, sequence), ...].
    header_key is the header (without '>') up to first whitespace.
    Duplicate coordinate headers are intentionally preserved because different
    CDS records can share the same coordinates but have different protein IDs.
    """
    records = []
    key = None
    chunks = []
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if key is not None:
                    records.append((key, "".join(chunks)))
                key = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if key is not None:
            records.append((key, "".join(chunks)))
    return records


def parse_cds_header(header):
    """
    Parse a CDS header key into:
      - contig (e.g., "NC_050608.1")
      - strand ('+' or '-')
      - coor_range: list of [start,end] exon segments, in header order

    Header examples:
      NC_050608.1-:64988-65192;NC_050608.1-:65249-65385;...
      chr2L+:100-200
    """
    h = header[1:] if header.startswith(">") else header
    parts = h.split(";")
    coor_range = []
    strand = None
    contig = None

    for seg in parts:
        chrom_part, coor_part = seg.split(":")
        start_s, end_s = coor_part.split("-")

        if chrom_part.endswith(("+", "-")):
            seg_strand = chrom_part[-1]
            chrom = chrom_part[:-1]
        else:
            seg_strand = "+"
            chrom = chrom_part

        if strand is None:
            strand = seg_strand
        if contig is None:
            contig = chrom

        start = int(start_s)
        end = int(end_s)
        coor_range.append([start, end])

    if strand == "-":
        coor_range = coor_range[::-1]

    return contig, strand, coor_range


class SNPRecord:
    __slots__ = ("chrom", "pos", "variants")
    def __init__(self, chrom, pos, variants):
        self.chrom = chrom
        self.pos = pos
        self.variants = variants  # list length = 2*nsamples


def extract_haplotype_alleles(sample_field, ref, alt, sep_hint=None):
    """
    sample_field: e.g. "0/1:..." or "0|1:..."
    Return two alleles (a0,a1) as chars.
    Missing or unexpected -> ref allele.
    """
    gt = sample_field.split(":")[0]
    sep = sep_hint if sep_hint is not None else ("|" if "|" in gt else "/")

    if sep not in gt:
        return ref, ref

    a, b = gt.split(sep)[:2]
    if a == "." or b == ".":
        return ref, ref

    # only handle biallelic indexes 0 and 1; else fallback to ref
    a0 = alt if a == "1" else ref
    a1 = alt if b == "1" else ref
    return a0, a1


def build_interval_index(records, parsed_headers):
    """
    Build per-contig interval lists for overlap queries.
    Returns:
      intervals[contig] = list of (start, end, record_idx)
      starts[contig]    = sorted start positions aligned to intervals list
      max_ends[contig]  = prefix max end aligned to intervals list
    """
    intervals = defaultdict(list)

    for record_idx, (contig, strand, coor_range) in enumerate(parsed_headers):
        for (s, e) in coor_range:
            start, end = (s, e) if s <= e else (e, s)
            intervals[contig].append((start, end, record_idx))

    starts = {}
    max_ends = {}
    for contig in intervals:
        intervals[contig].sort(key=lambda x: x[0])
        starts[contig] = [t[0] for t in intervals[contig]]
        running_max = 0
        max_ends[contig] = []
        for _start, end, _record_idx in intervals[contig]:
            running_max = max(running_max, end)
            max_ends[contig].append(running_max)

    return intervals, starts, max_ends


def find_overlapping_cds(contig, pos, intervals, starts, max_ends):
    """
    Return list of record indices whose interval on contig overlaps pos.
    """
    if contig not in intervals:
        return []

    ints = intervals[contig]
    st = starts[contig]
    mx = max_ends[contig]

    idx = bisect.bisect_right(st, pos) - 1
    if idx < 0:
        return []

    hits = []
    j = idx
    while j >= 0 and mx[j] >= pos:
        s, e, record_idx = ints[j]
        if pos <= e:
            hits.append(record_idx)
        j -= 1
    return hits


def snp_cds_index(pos, strand, coor_range):
    index = 0
    for (s, e) in coor_range:
        lo, hi = (s, e) if s <= e else (e, s)
        seg_len = hi - lo + 1
        if lo <= pos <= hi:
            if strand == "+":
                return index + (pos - lo)
            return index + (hi - pos)
        index += seg_len
    return None


def write_msa_for_record(
    record_idx,
    records,
    parsed_headers,
    snp_in_record,
    npopulations,
    outdir,
    species,
    resume=False,
):
    snps = snp_in_record.get(record_idx, [])
    if not snps:
        return "skipped_no_snp"

    fasta_path = os.path.join(outdir, f"msa_{record_idx}.fasta")
    if resume and os.path.exists(fasta_path):
        return "skipped_existing"

    cds_key, original_cds = records[record_idx]
    original_cds = original_cds.upper()
    contig, strand, coor_range = parsed_headers[record_idx]

    snp_index = []
    snp_vars = []
    for pos, rec in sorted(snps, key=lambda x: x[0]):
        idx = snp_cds_index(pos, strand, coor_range)
        if idx is None:
            raise ValueError(
                f"SNP {rec.chrom}:{pos} was assigned to record {record_idx}, "
                f"but is not inside CDS header {cds_key}"
            )
        if idx >= len(original_cds):
            raise IndexError(
                f"SNP index {idx} out of range for record {record_idx} "
                f"(CDS length {len(original_cds)}): {cds_key}"
            )
        if len(rec.variants) < npopulations:
            raise ValueError(
                f"Variant vector too short for record {record_idx}: "
                f"{len(rec.variants)} alleles for {npopulations} haplotypes"
            )
        snp_index.append(idx)
        snp_vars.append(rec.variants)

    os.makedirs(outdir, exist_ok=True)
    tmp_path = f"{fasta_path}.tmp.{os.getpid()}"
    with open(tmp_path, "w") as f:
        f.write(f">{species}_0\n{original_cds}\n")
        for i in range(npopulations):
            alt_cds = list(original_cds)
            for idx, varlist in zip(snp_index, snp_vars):
                var = varlist[i].upper()
                if strand == "-":
                    var = WCBP.get(var, "N")
                alt_cds[idx] = var
            f.write(f">{species}_{i+1}\n{''.join(alt_cds)}\n")
    os.replace(tmp_path, fasta_path)
    return "written"


def init_msa_worker(records, parsed_headers, snp_in_record, npopulations, outdir, species, resume):
    global WORKER_RECORDS
    global WORKER_PARSED_HEADERS
    global WORKER_SNP_IN_RECORD
    global WORKER_NPOPULATIONS
    global WORKER_OUTDIR
    global WORKER_SPECIES
    global WORKER_RESUME

    WORKER_RECORDS = records
    WORKER_PARSED_HEADERS = parsed_headers
    WORKER_SNP_IN_RECORD = snp_in_record
    WORKER_NPOPULATIONS = npopulations
    WORKER_OUTDIR = outdir
    WORKER_SPECIES = species
    WORKER_RESUME = resume


def msa_worker(record_idx):
    return write_msa_for_record(
        record_idx,
        WORKER_RECORDS,
        WORKER_PARSED_HEADERS,
        WORKER_SNP_IN_RECORD,
        WORKER_NPOPULATIONS,
        WORKER_OUTDIR,
        WORKER_SPECIES,
        resume=WORKER_RESUME,
    )


class GenomeReader:
    def __init__(self):
        self.records = []
        self.cds_order = []
        self.parsed_headers = []
        self.snp_in_record = {}   # record_idx -> list[(pos, SNPRecord)]
        self.npopulations = 0     # haplotypes = 2 * nsamples

    def read_genome(self, cds_fasta):
        t0 = time.perf_counter()
        print(f"Opening {cds_fasta} to read")
        self.records = read_cds_fasta_records(cds_fasta)
        t_read = time.perf_counter()
        self.cds_order = [header for header, _seq in self.records]
        self.parsed_headers = [parse_cds_header(header) for header in self.cds_order]
        t_parse = time.perf_counter()

        counts = Counter(self.cds_order)
        n_duplicate_records = sum(count - 1 for count in counts.values() if count > 1)
        print(f"Finished reading {cds_fasta}; CDS records: {len(self.records)}")
        if n_duplicate_records:
            print(
                f"Preserved {n_duplicate_records} duplicate coordinate-header records "
                f"across {sum(1 for count in counts.values() if count > 1)} duplicate groups"
            )
        print(
            "Timing read_genome: "
            f"read_fasta={format_seconds(t_read - t0)}, "
            f"parse_headers={format_seconds(t_parse - t_read)}, "
            f"total={format_seconds(t_parse - t0)}"
        )

    def get_snp(self, vcf_path, delimiter_hint=None):
        """
        Read VCF and assign SNPs to CDS by interval overlap.
        Always sets npopulations based on VCF header, even if 0 overlaps.
        """
        t0 = time.perf_counter()
        intervals, starts, max_ends = build_interval_index(self.records, self.parsed_headers)
        t_index = time.perf_counter()
        snp_in_record = defaultdict(list)
        nsamples = None

        print(f"Opening {vcf_path} to read")
        fh = readlines(vcf_path)

        n_read = 0
        n_used = 0
        for line in fh:
            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                cols = line.strip().split()
                nsamples = max(0, len(cols) - 9)
                continue

            if line.startswith("#"):
                continue

            n_read += 1
            info = line.split()
            if len(info) < 10:
                continue

            chrom = info[0]
            pos = int(info[1])
            ref = info[3].upper()
            alt = info[4].upper()

            # biallelic SNPs only
            if len(ref) != 1 or len(alt) != 1:
                continue
            if "," in alt:
                continue

            samples = info[9:]
            if nsamples is None:
                nsamples = len(samples)

            variants = []
            for s in samples:
                a0, a1 = extract_haplotype_alleles(s, ref, alt, sep_hint=delimiter_hint)
                variants.extend([a0.upper(), a1.upper()])

            hits = find_overlapping_cds(chrom, pos, intervals, starts, max_ends)
            if not hits:
                continue

            rec = SNPRecord(chrom=chrom, pos=pos, variants=variants)
            for record_idx in hits:
                snp_in_record[record_idx].append((pos, rec))
            n_used += 1

            if n_read % 200000 == 0:
                print(f"  processed {n_read:,} variants; matched {n_used:,} SNPs to CDS intervals")

        fh.close()

        if nsamples is None or nsamples == 0:
            print("No samples found in VCF (#CHROM header missing or no sample columns). Exiting.")
            sys.exit(1)

        self.snp_in_record = dict(snp_in_record)
        self.npopulations = 2 * nsamples
        t_scan = time.perf_counter()
        print(f"Finished VCF. Samples: {nsamples}, haplotypes: {self.npopulations}")
        print(f"CDS with ≥1 SNP: {len(self.snp_in_record)}")
        print(
            "Timing get_snp: "
            f"build_interval_index={format_seconds(t_index - t0)}, "
            f"scan_vcf={format_seconds(t_scan - t_index)}, "
            f"total={format_seconds(t_scan - t0)}"
        )

    def write_cds_snp_index_csv(self, outdir, filename="cds_snp_index.csv"):
        """
        Write CSV aligned to input FASTA order:
          index,header,SNP
        where SNP=1 if CDS has ≥1 overlapping SNP else 0.
        """
        t0 = time.perf_counter()
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, filename)
        with open(path, "w") as f:
            f.write("index,header,SNP\n")
            for i, cds_id in enumerate(self.cds_order):
                has = 1 if i in self.snp_in_record else 0
                # Note: header assumed to have no commas. If it can, switch to csv.writer.
                f.write(f"{i},{cds_id},{has}\n")
        print(f"Timing write_cds_snp_index_csv: {format_seconds(time.perf_counter() - t0)}")
        return path

    def mutate_one(self, record_idx, outdir, species, resume=False):
        """
        Write MSA fasta for one CDS (which must have ≥1 SNP).
        Output file: msa_<record_idx>.fasta  (record_idx is the original FASTA index)
        """
        return write_msa_for_record(
            record_idx,
            self.records,
            self.parsed_headers,
            self.snp_in_record,
            self.npopulations,
            outdir,
            species,
            resume=resume,
        )

    def mutate(self, outdir, ncpu=1, species="species", resume=False):
        """
        Write MSAs only for CDS with SNPs, but keep the original FASTA index in filenames.
        """
        t0 = time.perf_counter()
        os.makedirs(outdir, exist_ok=True)

        tasks = []
        for idx in range(len(self.records)):
            if idx in self.snp_in_record:
                tasks.append(idx)

        print(f"MSA tasks to write: {len(tasks)}")
        if resume:
            print("Resume mode: existing msa_<index>.fasta files will be kept")

        counts = Counter()
        if ncpu <= 1:
            for n_done, record_idx in enumerate(tasks, 1):
                counts[self.mutate_one(record_idx, outdir, species, resume=resume)] += 1
                if n_done % 1000 == 0:
                    elapsed = time.perf_counter() - t0
                    rate = n_done / elapsed if elapsed else 0
                    print(
                        f"  wrote/checked {n_done:,}/{len(tasks):,} MSA tasks "
                        f"({rate:.1f} tasks/s, elapsed {format_seconds(elapsed)})"
                    )
        else:
            chunksize = max(1, min(100, len(tasks) // (ncpu * 8) if tasks else 1))
            with Pool(
                ncpu,
                initializer=init_msa_worker,
                initargs=(
                    self.records,
                    self.parsed_headers,
                    self.snp_in_record,
                    self.npopulations,
                    outdir,
                    species,
                    resume,
                ),
            ) as pool:
                for n_done, status in enumerate(pool.imap_unordered(msa_worker, tasks, chunksize), 1):
                    counts[status] += 1
                    if n_done % 1000 == 0:
                        elapsed = time.perf_counter() - t0
                        rate = n_done / elapsed if elapsed else 0
                        print(
                            f"  wrote/checked {n_done:,}/{len(tasks):,} MSA tasks "
                            f"({rate:.1f} tasks/s, elapsed {format_seconds(elapsed)})"
                        )

        print(
            "MSA output summary: "
            + ", ".join(f"{key}={counts[key]}" for key in sorted(counts))
        )
        elapsed = time.perf_counter() - t0
        rate = len(tasks) / elapsed if elapsed else 0
        print(f"Timing mutate: total={format_seconds(elapsed)}, rate={rate:.1f} tasks/s")


def parse_args():
    parser = argparse.ArgumentParser(description="Convert VCF SNPs to CDS MSAs")
    parser.add_argument("-vcf", required=True, help="Input VCF file (can be .gz)")
    parser.add_argument("-fasta", required=True, help="CDS FASTA with coordinate headers")
    parser.add_argument("-out", required=True, help="Output directory")
    parser.add_argument("-species", required=True, help="Species name prefix for MSA headers")
    parser.add_argument("-ncpu", type=int, default=1, help="Number of processes (default: 1)")
    parser.add_argument(
        "-gtsep",
        choices=["/", "|", "auto"],
        default="auto",
        help="Genotype separator: '/' unphased, '|' phased, or auto-detect (default: auto)",
    )
    parser.add_argument(
        "-index_csv",
        default="cds_snp_index.csv",
        help="Index CSV filename to write under -out (default: cds_snp_index.csv)",
    )
    parser.add_argument(
        "-resume",
        action="store_true",
        help="Keep existing msa_<index>.fasta files and write only missing MSA files",
    )
    return parser.parse_args()


def main():
    t0 = time.perf_counter()
    args = parse_args()
    delimiter_hint = None if args.gtsep == "auto" else args.gtsep

    genome = GenomeReader()
    genome.read_genome(args.fasta)
    genome.get_snp(vcf_path=args.vcf, delimiter_hint=delimiter_hint)

    csv_path = genome.write_cds_snp_index_csv(args.out, filename=args.index_csv)
    print(f"Wrote index CSV: {csv_path}")

    genome.mutate(args.out, ncpu=args.ncpu, species=args.species, resume=args.resume)
    print(f"Timing total: {format_seconds(time.perf_counter() - t0)}")
    print("Done.")


if __name__ == "__main__":
    main()
