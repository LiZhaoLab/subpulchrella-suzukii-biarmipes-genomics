#!/usr/bin/env python3
"""Per-site allele counts by strand from a pooled BAM, for one region.

Runs samtools mpileup with the read and base filters of the manuscript and writes,
for each reference position in the region:
  counts[i, 0:8] = A+, A-, C+, C-, G+, G-, T+, T-   (uint16; base quality >= 20)
and a list of candidate-indel intervals to mask (indel allele supported by >= 3 reads
and >= 5% of usable depth; interval plus 5 bp on each side).

Usage: 01_pileup_counts.py SAMPLE:CHROM:START-END   (1-based, inclusive)
"""
import os
import re
import subprocess
import sys
import time

import numpy as np

BAMDIR = '/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup'
REF = {
    'Dsub': '/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna',
    'Dsuz': '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna',
}
OUTDIR = os.environ.get('QC_COUNTS_DIR', '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/counts')
MIN_MAPQ, MIN_BQ = 20, 20
INDEL_MIN_READS, INDEL_MIN_FRAC, INDEL_FLANK = 3, 0.05, 5

BASE_COL = {'A': 0, 'C': 2, 'G': 4, 'T': 6}
RE_START = re.compile(r'\^.')
RE_INDEL = re.compile(r'([+-])(\d+)')


def parse_bases(ref, s):
    """Return (8 counts, {indel allele: n reads}) for one pileup column."""
    if '^' in s:
        s = RE_START.sub('', s)
    indels = {}
    if '+' in s or '-' in s:
        parts = []
        i = 0
        for m in RE_INDEL.finditer(s):
            if m.start() < i:
                continue
            n = int(m.group(2))
            seq = s[m.end():m.end() + n]
            key = (m.group(1), seq.upper())
            indels[key] = indels.get(key, 0) + 1
            parts.append(s[i:m.start()])
            i = m.end() + n
        parts.append(s[i:])
        s = ''.join(parts)
    c = [0] * 8
    r = BASE_COL.get(ref)
    if r is not None:
        c[r] = s.count('.')
        c[r + 1] = s.count(',')
    for b, k in BASE_COL.items():
        if b != ref:
            c[k] = s.count(b)
            c[k + 1] = s.count(b.lower())
    return c, indels


def main():
    sample, chrom, span = sys.argv[1].split(':')
    start, end = (int(x) for x in span.split('-'))
    assert end >= start >= 1
    species = 'Dsub' if sample.startswith('Dsub') else 'Dsuz'
    bam = f'{BAMDIR}/{sample}.rmdup.bam'
    assert os.path.exists(bam), bam
    out = f'{OUTDIR}/{sample}/{chrom}_{start}_{end}.npz'
    os.makedirs(os.path.dirname(out), exist_ok=True)

    cmd = ['samtools', 'mpileup', '-B', '-d', '100000',
           '-q', str(MIN_MAPQ), '-Q', str(MIN_BQ),
           '--rf', '2', '--ff', '3844',
           '-f', REF[species], '-r', f'{chrom}:{start}-{end}', bam]
    L = end - start + 1
    counts = np.zeros((L, 8), dtype=np.uint16)
    masks = []
    t0 = time.time()
    n = 0
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20) as p:
        for line in p.stdout:
            f = line.split('\t')
            pos = int(f[1])
            ref = f[2].upper()
            c, ind = parse_bases(ref, f[4])
            counts[pos - start] = np.minimum(c, 65535)
            if ind:
                depth = sum(c)
                for (kind, seq), k in ind.items():
                    if k >= INDEL_MIN_READS and depth > 0 and k / depth >= INDEL_MIN_FRAC:
                        # deletion removes pos+1..pos+len; insertion sits between pos and pos+1
                        last = pos + len(seq) if kind == '-' else pos
                        masks.append((pos + 1 - INDEL_FLANK, last + INDEL_FLANK))
            n += 1
    assert p.returncode == 0, f'samtools exit {p.returncode}'
    masks = np.array(masks, dtype=np.int64).reshape(-1, 2)
    np.savez_compressed(out + '.tmp.npz', counts=counts, indel_mask=masks,
                        chrom=chrom, start=start, end=end, sample=sample,
                        cmd=' '.join(cmd))
    os.replace(out + '.tmp.npz', out)
    print(f'{sample} {chrom}:{start}-{end} columns={n} masks={len(masks)} '
          f'seconds={time.time() - t0:.0f}', flush=True)


if __name__ == '__main__':
    main()
