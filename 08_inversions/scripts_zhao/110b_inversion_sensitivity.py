#!/usr/bin/env python3
"""Inversion analysis, approach 3: sensitivity of the frequency scan (scripts/110) in the
D. subpulchrella pool.

Synthetic inversion blocks are added to the real window data: in 1-Mb blocks (10 windows) at 50 random
autosomal positions per scenario, each window receives Poisson(lambda) extra SNPs carried by the
inverted arrangement at pool frequency f. Each extra SNP gets a depth drawn from the pool's autosomal
depth distribution (callable range), reads of the arrangement-specific allele ~ Binomial(depth, f), and
passes the filter of scripts/110 (>= 3 reads of each allele, MAF >= 0.05). Window scores are recomputed,
and the block counts as detected when a run or a cluster of scripts/110 overlaps it at the expected
band (|peak MAF - min(f, 1 - f)| <= 0.05).
lambda (extra SNPs per 100 kb) is converted to net divergence between arrangements:
d = lambda / median callable sites per window.
Output: results/INV_sensitivity_Dsub_L5_7.tsv and .log. Seed 20261030.
"""
import importlib.util

import numpy as np
import pandas as pd

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('scan', f'{P}/scripts/110_inversion_freq_scan.py')
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)
POOL, REPS, SEED = 'Dsub_L5_7', 50, 20261030
FS = [0.1, 0.2, 0.3, 0.5]
LAMBDAS = [20, 40, 80, 160, 320]


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(f'{sc.RES}/INV_freq_windows_{POOL}.tsv', sep='\t')
    au = df[df.cls == 'autosome'].reset_index(drop=True)
    bcols = [f'b{k}' for k in range(18)]
    h0 = au[bcols].values.astype(float)
    q = h0.sum(0) / h0.sum()
    thr = pd.read_csv(f'{sc.RES}/depth_thresholds.tsv', sep='\t')
    t = thr[(thr['sample'] == POOL) & (thr.cls == 'autosome')].iloc[0]
    hist = np.load(f'{sc.QC}/depth_hist_{POOL}_autosome.npy').astype(float)
    depths = np.arange(len(hist))
    p = np.where((depths >= t.min_depth) & (depths <= t.p98), hist, 0)
    p /= p.sum()
    med_call = float(np.median(au.callable[au.snps >= 50]))
    starts = [(c, w) for c, g in au.groupby('chrom') for w in g.w.values if (g.w.values == w + 9).any()]
    band_f = sc.EDGES[1:-1]
    rows = []
    for f in FS:
        target = min(f, 1 - f)
        for lam in LAMBDAS:
            hits = 0
            for _ in range(REPS):
                chrom, w0 = starts[rng.integers(len(starts))]
                h = h0.copy()
                idx = au.index[(au.chrom == chrom) & (au.w >= w0) & (au.w <= w0 + 9)].values
                for i in idx:
                    n = rng.poisson(lam)
                    d = rng.choice(depths, size=n, p=p)
                    k = rng.binomial(d, f)
                    keep = (k >= sc.MIN_READS) & (d - k >= sc.MIN_READS)
                    maf = np.minimum(k, d - k)[keep] / d[keep]
                    maf = maf[maf >= sc.MIN_MAF]
                    b = np.clip(np.searchsorted(sc.EDGES, maf, side='right') - 1, 0, 17)
                    np.add.at(h[i], b, 1)
                # rescore the affected chromosome only (q fixed at the genome-wide value)
                cidx = au.index[au.chrom == chrom].values
                E, z = sc.window_scores(h[cidx], q)
                g = au.loc[cidx, ['chrom', 'w', 'start']].copy()
                has = np.isfinite(z).any(1)
                zz = np.where(np.isfinite(z), z, -np.inf)
                g['peak_z'] = np.where(has, zz.max(1), np.nan)
                g['peak_f'] = np.where(has, band_f[zz.argmax(1)], np.nan)
                g = g[g.peak_z.notna()]
                found = False
                for run in sc.runs_of(g):
                    ws = [r.w for r in run]
                    if max(ws) >= w0 and min(ws) <= w0 + 9 and abs(np.mean([r.peak_f for r in run]) - target) <= 0.05:
                        found = True
                for c in sc.clusters_of(g):
                    if c[3] >= w0 and c[2] <= w0 + 9 and abs(c[1] - target) <= 0.05:
                        found = True
                hits += found
            rows.append(dict(f=f, lam_per_100kb=lam, net_divergence=lam / med_call, detected=hits / REPS))
            print(f'f {f} lambda {lam} (d {lam / med_call:.5f}): detected {hits}/{REPS}', flush=True)
    r = pd.DataFrame(rows)
    r.to_csv(f'{sc.RES}/INV_sensitivity_{POOL}.tsv', sep='\t', index=False)
    with open(f'{sc.RES}/INV_sensitivity_{POOL}.log', 'w') as fo:
        fo.write(f'median callable sites per 100-kb window {med_call:.0f}; median SNPs (MAF >= 0.05) per window '
                 f'{np.median(au.snps[au.snps >= 50]):.0f}\n')
        fo.write(r.pivot(index='lam_per_100kb', columns='f', values='detected').to_string() + '\n')
    print(open(f'{sc.RES}/INV_sensitivity_{POOL}.log').read())


if __name__ == '__main__':
    main()
