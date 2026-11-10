#!/usr/bin/env python3
"""Analysis F of the depth plan: SweepFinder2 input for one pool, genome-wide.

Sites: the P98_RM configuration of scripts/32 (depth max(10, P05)..P98, outside indel masks and
the RepeatMasker mask), major chromosomes only.
Alleles (shared with the DFE analysis, scripts/102): an allele passes if it has >= 3 reads,
>= 1 read on each strand and >= 1% of the site depth; the two most frequent passing alleles are
kept and all other reads are dropped. Then N_SUB = 20 reads are drawn without replacement
(hypergeometric; seed 20261030 + chromosome index) from the kept reads; sites with fewer than 20
kept reads are dropped. The count of the more frequent pool allele in the subsample is x
(folded spectrum; no outgroup for non-coding sites).
Outputs in QC/sweeps/{POOL}/:
  {chrom}.sf2.txt     SweepFinder2 input (position, x, n, folded = 1), polymorphic sites only
  {cls}.all.sf2.txt   all polymorphic sites of the class (autosome or X), for the background SFS
  summary.npz         per chromosome: callable sites (subsampled), folded SFS of all, short-intron
                      and four-fold sites (index 0..10; index 0 = monomorphic)
Usage: 98_sweep_input.py POOL
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
A = '/ru-auth/local/home/asato/scratch'
XCHROM = {'NC_050613.1', 'NC_092084.1'}
N_SUB, MIN_READS, MIN_FRAC, SEED = 20, 3, 0.01, 20261030


def two_alleles(c8):
    """Allele filter and biallelic reduction. c8: (L, 8) counts A+ A- C+ C- G+ G- T+ T-.
    Returns (k1, k2): read counts of the major and second passing allele (k2 = 0 if none)."""
    tot = c8[:, 0::2].astype(np.int64) + c8[:, 1::2]
    depth = tot.sum(1)
    ok = (tot >= MIN_READS) & (c8[:, 0::2] > 0) & (c8[:, 1::2] > 0) & (tot >= MIN_FRAC * depth[:, None])
    t = np.where(ok, tot, 0)
    order = np.argsort(-t, axis=1, kind='stable')
    k1 = np.take_along_axis(t, order[:, :1], 1)[:, 0]
    k2 = np.take_along_axis(t, order[:, 1:2], 1)[:, 0]
    return k1, k2, order[:, 0], order[:, 1]


def subsample(k1, k2, n, rng):
    """Hypergeometric draw of n reads from k1 + k2; returns the count of allele 1 (-1 if < n reads)."""
    m = k1 + k2
    x = np.full(len(m), -1)
    ok = m >= n
    x[ok] = rng.hypergeometric(k1[ok], k2[ok], n)
    return x


def main():
    pool = sys.argv[1]
    sp = 'Dsub' if pool.startswith('Dsub') else 'Dsuz'
    out = f'{QC}/sweeps/{pool}' + ('_test' if os.environ.get('TEST_CHROM') else '')
    os.makedirs(out, exist_ok=True)
    regions = f'{A}/population_polymorphism/reanalysis_20260821/config/regions'
    chroms = {cls: [l.split()[0] for l in open(f'{regions}/{sp}.{cls}.bed')] for cls in ('autosome', 'X')}
    thr = pd.read_csv(f'{RES}/depth_thresholds.tsv', sep='\t')
    thr = thr[thr['sample'] == pool].set_index('cls')
    ms = np.load(f'{QC}/masks_sites/{sp}.npz')
    summ = {}
    for cls, chs in chroms.items():
        dmin, dmax = int(thr.loc[cls, 'min_depth']), int(thr.loc[cls, 'p98'])
        allf = open(f'{out}/{cls}.all.sf2.txt', 'w')
        allf.write('position\tx\tn\tfolded\n')
        for ci, chrom in enumerate(chs):
            files = sorted(glob.glob(f'{QC}/counts/{pool}/{chrom}_*.npz'), key=lambda f: int(f.split('_')[-2]))
            if not files or (os.environ.get('TEST_CHROM') and chrom != os.environ['TEST_CHROM']):
                continue
            rng = np.random.default_rng(SEED + 100 * (cls == 'X') + ci)
            xs, vs = [], []
            for f in files:                     # one 10-Mb count file at a time (memory)
                d = np.load(f)
                s0, e0 = int(d['start']), int(d['end'])
                m = np.zeros(e0 - s0 + 2, dtype=np.int32)
                for s, e in d['indel_mask']:
                    s, e = max(s, s0), min(e, e0)
                    if s <= e:
                        m[s - s0] += 1
                        m[e - s0 + 1] -= 1
                c8 = d['counts']
                depth = c8.sum(1, dtype=np.int64)
                vs.append((depth >= dmin) & (depth <= dmax) & ~(np.cumsum(m)[:-1] > 0)
                          & ~ms[f'{chrom}|rm'][s0:e0 + 1])
                k1, k2, _, _ = two_alleles(c8)
                xs.append(subsample(k1, k2, N_SUB, rng).astype(np.int8))
            x, valid = np.concatenate(xs), np.concatenate(vs)
            L = len(x)
            si, ff = (ms[f'{chrom}|{k}'][1:L + 1] for k in ('si', 'ff'))
            use = valid & (x >= 0)
            poly = use & (x > 0) & (x < N_SUB)
            pos = np.flatnonzero(poly) + 1
            x = x.astype(np.int64)
            fold = np.minimum(x, N_SUB - x)
            with open(f'{out}/{chrom}.sf2.txt', 'w') as fo:
                fo.write('position\tx\tn\tfolded\n')
                for p, xi in zip(pos, x[poly]):
                    fo.write(f'{p}\t{xi}\t{N_SUB}\t1\n')
                    allf.write(f'{p}\t{xi}\t{N_SUB}\t1\n')
            for name, msk in (('all', use), ('si', use & si), ('ff', use & ff)):
                summ[f'{chrom}|{name}'] = np.bincount(fold[msk], minlength=N_SUB // 2 + 1)
            summ[f'{chrom}|length'] = np.array([L])
            print(f'{pool} {chrom}: length {L} callable {use.sum()} polymorphic {poly.sum()} '
                  f'pi_all {np.mean(2 * x[use] * (N_SUB - x[use]) / (N_SUB * (N_SUB - 1))):.5f}',
                  flush=True)
        allf.close()
    np.savez(f'{out}/summary.npz', **summ)


if __name__ == '__main__':
    main()
