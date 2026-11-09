#!/usr/bin/env python3
"""Script 122: numbers for the rewritten Results paragraph on synteny.

A. Lineage of the D. suzukii-D. subpulchrella breakpoints (same rule as script 91 section C, main parameters,
   outgroup D. biarmipes, tolerance 100 kb), and which D. subpulchrella-lineage breakpoints lie within 100 kb of an
   end of the three inversions still polymorphic in the D. subpulchrella pool (Results_INV document).
B. Centre of D. suzukii chromosome 3 (NC_092082.1, 27-62 Mb; flanks 0-27 Mb = 3L-like, 62 Mb-end = 3R-like):
   PROmer coverage by blocks with any D. subpulchrella sequence (unplaced scaffolds included).
C. Unmasked minimap2 alignment of NC_092082.1 (script 121) to Dsuz_RU_1.0 and to RU_Dsub_v1.1: merged query
   coverage per region for all primary alignments, MAPQ >= 20, and alignment blocks >= 5 kb; target sequences hit.
Output: results/synteny_paragraph_numbers.log (stdout). Python: anaconda3/envs/jupyter-ml.
"""
import os

import numpy as np
import pandas as pd

RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
SYN = '/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824/results'
MM = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/synteny_chr3'
TOL = 100_000
INV = {'NC_050613.1': ('X', 9_926_301, 16_512_084), 'NC_050612.1': ('3L', 9_965_340, 27_321_995),
       'NC_050611.1': ('2R', 16_050_306, 16_287_188)}
CHR3, CHR3_LEN = 'NC_092082.1', 93_255_888
REGIONS = {'3L-like 0-27 Mb': (0, 27_000_000), 'centre 27-62 Mb': (27_000_000, 62_000_000),
           '3R-like 62 Mb-end': (62_000_000, CHR3_LEN)}


def merged_bp(iv, a, b):
    """Merged bp of intervals iv (list of (s, e)) clipped to [a, b)."""
    iv = sorted((max(s, a), min(e, b)) for s, e in iv if e > a and s < b)
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce:
            tot += 0 if cs is None else ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (0 if cs is None else ce - cs)


# ---------------------------------------------------------------- A
SEG = pd.read_csv(f'{RES}/synteny_segments.tsv', sep='\t')
a = SEG[(SEG.params == 'main') & (SEG.comp == 'Dsuz_vs_Dsub')]
o = SEG[(SEG.params == 'main') & (SEG.comp == 'Dsub_vs_Dbia')]
rows = []
for q, g in a.groupby('query'):
    g = g.sort_values('qs')
    bps = [(x1.qe + x2.qs) / 2 for x1, x2 in zip(g.iloc[:-1].itertuples(), g.iloc[1:].itertuples())]
    og = o[o.ref == q].sort_values('rs')
    obps = np.array([(y1.re + y2.rs) / 2 for y1, y2 in zip(og.iloc[:-1].itertuples(), og.iloc[1:].itertuples())])
    for bp in bps:
        if len(obps) and np.min(np.abs(obps - bp)) <= TOL:
            lin = 'Dsub'
        elif ((og.rs <= bp - TOL) & (og.re >= bp + TOL)).any():
            lin = 'Dsuz'
        else:
            lin = 'unassigned'
        end = ''
        if q in INV:
            arm, s, e = INV[q]
            for p in (s, e):
                if abs(bp - p) <= TOL:
                    end = f'{arm} inversion end {p}'
        rows.append(dict(query=q, pos=int(bp), lineage=lin, polymorphic_inv_end=end))
R = pd.DataFrame(rows)
print('== A. D. suzukii-D. subpulchrella breakpoints (main parameters)')
print(R.lineage.value_counts().to_dict())
print(R[R.polymorphic_inv_end != ''].to_string(index=False))
print('D. subpulchrella-lineage breakpoints at a polymorphic-inversion end:',
      int(((R.lineage == 'Dsub') & (R.polymorphic_inv_end != '')).sum()))
print('2L breakpoints:', R[R['query'] == 'NC_050610.1'][['pos', 'lineage']].to_dict('records'))

# ---------------------------------------------------------------- B
COLS = ['S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', 'IDY', 'SIM', 'STP', 'LENR', 'LENQ', 'COVR', 'COVQ',
        'FRM1', 'FRM2', 'R', 'Q']
t = pd.read_csv(f'{SYN}/primary_one_to_one/Dsuz_vs_Dsub.promer.one_to_one.L500.I50.coords.tsv', sep='\t',
                header=None, names=COLS)
t = t[t.R == CHR3]
iv = list(zip(t[['S1', 'E1']].min(1) - 1, t[['S1', 'E1']].max(1)))
print('\n== B. PROmer coverage of D. suzukii chromosome 3 by any D. subpulchrella sequence (L500, I50)')
for name, (s, e) in REGIONS.items():
    sub = t[(t[['S1', 'E1']].max(1) > s) & (t[['S1', 'E1']].min(1) <= e)]
    print(f'{name}: {merged_bp(iv, s, e) / (e - s):.2%} covered; blocks {len(sub)}; '
          f'to unplaced D. subpulchrella sequences {(~sub.Q.str.startswith("NC_") & (sub.Q != "NW_023665577.1")).sum()}')

# ---------------------------------------------------------------- C
print('\n== C. unmasked minimap2 asm20 alignment of NC_092082.1 (primary alignments)')
P = ['qn', 'ql', 'qs', 'qe', 'strand', 'tn', 'tl', 'ts', 'te', 'match', 'blk', 'mapq']
for tgt in ('RU1', 'Dsub'):
    f = f'{MM}/chr3_vs_{tgt}.paf'
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        print(f'{tgt}: PAF not available')
        continue
    p = pd.read_csv(f, sep='\t', header=None, usecols=range(12), names=P)
    assert (p.qn == CHR3).all()
    for label, sel in [('all', p), ('MAPQ>=20', p[p.mapq >= 20]), ('block>=5kb', p[p.blk >= 5000])]:
        ivs = list(zip(sel.qs, sel.qe))
        print(f'{tgt} {label}: ' + '; '.join(f'{n} {merged_bp(ivs, s, e) / (e - s):.1%}'
                                              for n, (s, e) in REGIONS.items()))
    c = p[(p.qe > 27_000_000) & (p.qs < 62_000_000) & (p.mapq >= 20)]
    tops = c.groupby('tn').apply(lambda x: merged_bp(list(zip(x.qs, x.qe)), 27_000_000, 62_000_000)) \
        .sort_values(ascending=False).head(6)
    print(f'{tgt} centre, MAPQ>=20, query bp by target sequence (top 6, Mb):', (tops / 1e6).round(2).to_dict())

# ---------------------------------------------------------------- D (added the same day, after C)
# D1. Size of the unplaced part of each target assembly. D2. Conservative minimap2 measure: residue matches
# (PAF column 10; without -c this is the seed-based count, a lower bound) of primary alignments with MAPQ >= 20,
# divided by region length. D3. Genes of the centre: CDS bp (RefSeq GCF_043229965.1-RS_2025_01, isoforms merged)
# inside PROmer blocks (masked genomes) to D. subpulchrella chromosomes vs unplaced scaffolds, per region.
import re
print('\n== D1. unplaced sequence in the target assemblies')
for name, fai in [('RU1', '/ru-auth/local/home/asato/scratch/genome/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.fna.fai'),
                  ('Dsub', '/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824/repeatmasker/Dsub/Dsub.GCF_014743375.2.fna.masked.fai')]:
    f = pd.read_csv(fai, sep='\t', header=None, usecols=[0, 1], names=['seq', 'len'])
    unpl = ~f.seq.str.startswith('NC_') & (f.seq != 'NW_023665577.1')
    print(f'{name}: total {f.len.sum() / 1e6:.1f} Mb; unplaced {f[unpl].len.sum() / 1e6:.1f} Mb in {unpl.sum()} sequences')
print('\n== D2. residue matches of primary MAPQ >= 20 alignments / region length (lower bound)')
for tgt in ('RU1', 'Dsub'):
    p = pd.read_csv(f'{MM}/chr3_vs_{tgt}.paf', sep='\t', header=None, usecols=range(12), names=P)
    p = p[p.mapq >= 20]
    mid = (p.qs + p.qe) / 2
    unpl = ~p.tn.str.startswith('NC_') & (p.tn != 'NW_023665577.1')
    out = []
    for n, (s, e) in REGIONS.items():
        m = (mid >= s) & (mid < e)
        out.append(f'{n} {p[m].match.sum() / (e - s):.1%} (to unplaced {p[m & unpl].match.sum() / (e - s):.1%}; '
                   f'median match/block {np.median(p[m].match / p[m].blk):.2f})')
    print(f'{tgt}: ' + '; '.join(out))
print('\n== D3. CDS of D. suzukii chromosome 3 inside PROmer blocks to D. subpulchrella')
gff = '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'
cds = []
with open(gff) as fh:
    for line in fh:
        if line.startswith(CHR3 + '\t'):
            f = line.split('\t')
            if f[2] == 'CDS':
                g = re.search(r'gene=([^;]+)', f[8])
                cds.append((int(f[3]) - 1, int(f[4]), g.group(1) if g else 'NA'))
cds = pd.DataFrame(cds, columns=['s', 'e', 'gene'])
unplq = ~t.Q.str.startswith('NC_') & (t.Q != 'NW_023665577.1')
iv_chr = list(zip(t[~unplq][['S1', 'E1']].min(1) - 1, t[~unplq][['S1', 'E1']].max(1)))
iv_unp = list(zip(t[unplq][['S1', 'E1']].min(1) - 1, t[unplq][['S1', 'E1']].max(1)))


def cov_of(cds_iv, blocks):
    """bp of merged CDS intervals that lie inside merged blocks."""
    tot = 0
    for s, e in cds_iv:
        tot += merged_bp(blocks, s, e)
    return tot


for n, (s, e) in REGIONS.items():
    c = cds[(cds.s >= s) & (cds.e <= e)]
    m = []
    for a, b in sorted(zip(c.s, c.e)):
        if m and a <= m[-1][1]:
            m[-1][1] = max(m[-1][1], b)
        else:
            m.append([a, b])
    tot = sum(b - a for a, b in m)
    print(f'{n}: {c.gene.nunique()} genes, CDS {tot / 1e6:.2f} Mb; inside blocks to D. subpulchrella chromosomes '
          f'{cov_of(m, iv_chr) / tot:.1%}, to unplaced scaffolds {cov_of(m, iv_unp) / tot:.1%}')
