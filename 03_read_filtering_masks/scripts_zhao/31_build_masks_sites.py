#!/usr/bin/env python3
"""Build per-chromosome boolean arrays: RepeatMasker mask, short-intron sites, four-fold sites.

Independent of Airi's BED files: parses the RepeatMasker .out and RefSeq GFF directly.
  short_intron: positions 8-30 (from the 5' splice site, strand-aware) of introns 30-65 bp
                in any mRNA; positions overlapping any annotated exon are removed.
  fourfold:     third positions of CCN, ACN, GTN, GCN, GGN codons in the longest complete
                CDS (length divisible by 3) of each protein-coding gene.
Output: QC_DIR/masks_sites/{species}.npz with arrays '{chrom}|rm', '{chrom}|si', '{chrom}|ff'.
Usage: 31_build_masks_sites.py SPECIES   (Dsub or Dsuz)
Environment variable SI_MIN (default 30) sets the lower intron-length bound. With SI_MIN=40
(the definition of Airi Sato's final Ne run, Ne_reanalysis_20261023_repeatmasker_p98) the output goes to
QC_DIR/masks_sites_si40/{species}.npz; the default run and its output are unchanged.
"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

A = '/ru-auth/local/home/asato/scratch'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
CFG = {
    'Dsub': dict(acc='GCF_014743375.2',
                 fna=f'{A}/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna',
                 gff=f'{A}/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff'),
    'Dsuz': dict(acc='GCF_043229965.1',
                 fna=f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna',
                 gff=f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'),
}
REGIONS = f'{A}/population_polymorphism/reanalysis_20260821/config/regions'
FF_PREFIX = {'CC', 'AC', 'GT', 'GC', 'GG'}
SI_MIN, SI_MAX = int(os.environ.get('SI_MIN', 30)), 65
OUTDIR = 'masks_sites' if SI_MIN == 30 else f'masks_sites_si{SI_MIN}'
COMP = str.maketrans('ACGTN', 'TGCAN')


def read_chroms(sp):
    chroms = {}
    for cls in ('autosome', 'X'):
        for line in open(f'{REGIONS}/{sp}.{cls}.bed'):
            c, s, e = line.split()[:3]
            chroms[c] = int(e)
    return chroms


def read_genome(fna, keep):
    seqs, name, buf = {}, None, []
    with open(fna) as f:
        for line in f:
            if line.startswith('>'):
                if name in keep:
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in keep:
                buf.append(line.strip())
    if name in keep:
        seqs[name] = ''.join(buf).upper()
    return seqs


def main():
    sp = sys.argv[1]
    cfg = CFG[sp]
    chroms = read_chroms(sp)
    arr = {c: {k: np.zeros(L + 1, dtype=bool) for k in ('rm', 'si', 'ff', 'exon')}
           for c, L in chroms.items()}           # index = 1-based position

    # RepeatMasker .out: columns 5-7 = query, begin, end (1-based, inclusive)
    rm = f'{A}/synteny_4species/reanalysis_20260824/repeatmasker/{sp}/{sp}.{cfg["acc"]}.fna.out'
    with open(rm) as f:
        for line in f:
            t = line.split()
            if len(t) < 11 or not t[0].isdigit():
                continue
            c, s, e = t[4], int(t[5]), int(t[6])
            if c in arr:
                arr[c]['rm'][s:e + 1] = True

    # GFF: genes, mRNAs, exons, CDS
    parent_of, gene_biotype = {}, {}
    exons = defaultdict(list)       # mRNA id -> [(chrom, s, e, strand)]
    cds = defaultdict(list)         # mRNA id -> [(chrom, s, e, strand, phase)]
    with open(cfg['gff']) as f:
        for line in f:
            if line.startswith('#'):
                continue
            t = line.rstrip('\n').split('\t')
            if len(t) < 9:
                continue
            c, typ, s, e, strand, phase, attr = t[0], t[2], int(t[3]), int(t[4]), t[6], t[7], t[8]
            a = dict(kv.split('=', 1) for kv in attr.split(';') if '=' in kv)
            if typ == 'gene':
                gene_biotype[a['ID']] = a.get('gene_biotype', '')
            elif typ == 'mRNA':
                parent_of[a['ID']] = a.get('Parent', '')
            elif typ == 'exon':
                if c in arr:
                    arr[c]['exon'][s:e + 1] = True
                exons[a.get('Parent', '')].append((c, s, e, strand))
            elif typ == 'CDS':
                cds[a.get('Parent', '')].append((c, s, e, strand, phase))

    # short introns
    n_introns = 0
    for tid, ex in exons.items():
        if tid not in parent_of:
            continue
        ex = sorted(ex, key=lambda x: x[1])
        for (c, s1, e1, st), (c2, s2, e2, st2) in zip(ex, ex[1:]):
            if c not in arr or c != c2:
                continue
            ilen = s2 - e1 - 1
            if SI_MIN <= ilen <= SI_MAX:
                n_introns += 1
                if st == '+':
                    a0, b0 = e1 + 8, e1 + 30              # intron positions 8..30 from donor
                else:
                    a0, b0 = s2 - 30, s2 - 8
                arr[c]['si'][a0:b0 + 1] = True
    for c in arr:
        arr[c]['si'] &= ~arr[c]['exon']

    # four-fold sites from the longest complete CDS per gene
    genome = read_genome(cfg['fna'], set(chroms))
    best = {}
    for tid, parts in cds.items():
        g = parent_of.get(tid)
        if not g or gene_biotype.get(g) != 'protein_coding' or parts[0][0] not in arr:
            continue
        L = sum(e - s + 1 for _, s, e, _, _ in parts)
        if L % 3 == 0 and (g not in best or L > best[g][0] or (L == best[g][0] and tid < best[g][1])):
            best[g] = (L, tid)
    n_ff_genes = 0
    for g, (L, tid) in best.items():
        parts = cds[tid]
        c, strand = parts[0][0], parts[0][3]
        if any(p[0] != c or p[3] != strand for p in parts):
            continue
        pos = []
        for _, s, e, _, _ in sorted(parts, key=lambda x: x[1]):
            pos.extend(range(s, e + 1))
        if strand == '-':
            pos = pos[::-1]
        seq = ''.join(genome[c][p - 1] for p in pos)
        if strand == '-':
            seq = seq.translate(COMP)
        n_ff_genes += 1
        for i in range(0, len(seq) - 2, 3):
            if seq[i:i + 2] in FF_PREFIX and seq[i + 2] in 'ACGT':
                arr[c]['ff'][pos[i + 2]] = True

    out = {}
    for c in arr:
        for k in ('rm', 'si', 'ff'):
            out[f'{c}|{k}'] = arr[c][k]
    os.makedirs(f'{QC}/{OUTDIR}', exist_ok=True)
    np.savez_compressed(f'{QC}/{OUTDIR}/{sp}.npz', **out)
    tot = sum(chroms.values())
    print(f'{sp}: introns {SI_MIN}-{SI_MAX} bp -> {OUTDIR}; chromosomes={len(chroms)} length={tot} '
          f'RM={sum(arr[c]["rm"].sum() for c in arr)} ({100 * sum(arr[c]["rm"].sum() for c in arr) / tot:.2f}%) '
          f'short_introns={n_introns} SI_sites={sum(arr[c]["si"].sum() for c in arr)} '
          f'FF_genes={n_ff_genes} FF_sites={sum(arr[c]["ff"].sum() for c in arr)}')


if __name__ == '__main__':
    main()
