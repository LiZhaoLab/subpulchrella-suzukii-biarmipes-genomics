#!/usr/bin/env python3
"""Analysis A3 of the depth plan: piN/piS per species and its relation to neutral
diversity along the genome.

piN/piS = pi at 0-fold sites / pi at 4-fold sites (per-site pi in 20 reads without singletons,
scripts/depth_common.py), on orthologous genes present in all five Step-0 tables.
1. Per arm and for the autosomes: piN/piS per pool with a 1-Mb block bootstrap (2,000 replicates;
   seed 20261030), and the ratio D. subpulchrella / D. suzukii pool.
2. Across 1-Mb and 2-Mb windows in D. suzukii coordinates (orthologous genes): regression of
   y = log(piN/piS)_Dsub - log(piN/piS)_Dsuz on x = log(pi_neutral)_Dsub - log(pi_neutral)_Dsuz.
   To avoid a slope created by shared noise (piS in y and in x), x never uses the sites in y:
   (a) x from short introns, y from 0-fold / all 4-fold sites;
   (b) x from 4-fold sites of odd codons, y from 0-fold / 4-fold sites of even codons.
   Under a gamma DFE with shape beta, piN/piS scales approximately as Ne^-beta (Welch et al. 2008;
   Galtier 2016), so the expected slope is about -beta (beta from analysis A1, model G). Noise in
   x biases the OLS slope toward zero (regression dilution); windows are resampled for the 95%
   interval (2,000 replicates).
Outputs: results/A3_pnps_arms.tsv, results/A3_pnps_windows.tsv,
         results/A3_pnps.log
"""
import glob
import time

import numpy as np
import pandas as pd

import depth_common as dc

ARMS = ['2L', '2R', '3L', '3R', '3R_scaffold', 'X']
N_BOOT = int(__import__('os').environ.get('N_BOOT', 2000))


def pnps(s):
    return (s['f0_pi'] / s['f0_n']) / (s['f4_pi'] / s['f4_n'])


def main():
    t0 = time.time()
    logf = open(f'{dc.RES}/A3_pnps.log', 'w')

    def log(x):
        print(x, flush=True)
        logf.write(x + '\n')

    log(f'A3 piN/piS; seed {dc.SEED}; bootstraps {N_BOOT}; started {time.ctime()}')
    d, zs = dc.gene_tables()
    sums = {}
    for p in dc.POOLS:
        rng = np.random.default_rng(dc.SEED)
        sums[p] = dc.gene_sums(zs[p], rng).iloc[d[p].values].reset_index(drop=True)
    rng = np.random.default_rng(dc.SEED + 2)
    rows, boot = [], {}
    groups = {a: (d.arm == a).values for a in ARMS}
    groups['autosomes'] = d.arm.isin(['2L', '2R', '3L', '3R', '3R_scaffold']).values
    for grp, sel in groups.items():
        for p in dc.POOLS:
            g = sums[p]
            blocks = (g.chrom.astype(str) + ':' + (g.pos // 1_000_000).astype(int).astype(str)).values[sel]
            gi = np.flatnonzero(sel)
            ub = np.unique(blocks)
            bsum = np.array([g.iloc[gi[blocks == b]][['f0_pi', 'f0_n', 'f4_pi', 'f4_n']].sum().values for b in ub])
            est = pnps(dict(zip(['f0_pi', 'f0_n', 'f4_pi', 'f4_n'], bsum.sum(0))))
            bs = []
            for _ in range(N_BOOT):
                t = bsum[rng.integers(0, len(ub), len(ub))].sum(0)
                bs.append((t[0] / t[1]) / (t[2] / t[3]))
            boot[(grp, p)] = np.array(bs)
            rows.append(dict(group=grp, pool=p, genes=int(sel.sum()), pnps=est,
                             lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
    ar = pd.DataFrame(rows)
    for grp in groups:
        s = ar[(ar.group == grp) & (ar.pool == 'Dsub_L5_7')].pnps.iloc[0]
        for p in dc.DSUZ:
            r = boot[(grp, 'Dsub_L5_7')] / boot[(grp, p)]
            i = (ar.group == grp) & (ar.pool == p)
            ar.loc[i, 'ratio_dsub_over_pool'] = s / ar.loc[i, 'pnps'].iloc[0]
            ar.loc[i, 'ratio_lo'], ar.loc[i, 'ratio_hi'] = np.percentile(r, [2.5, 97.5])
    ar.to_csv(f'{dc.RES}/A3_pnps_arms.tsv', sep='\t', index=False)
    log('\n== piN/piS (0-fold / 4-fold pi) per arm [95% block bootstrap]; ratio Dsub / pool')
    for r in ar.itertuples():
        extra = '' if r.pool == 'Dsub_L5_7' else f' ratio {r.ratio_dsub_over_pool:.3f} [{r.ratio_lo:.3f}-{r.ratio_hi:.3f}]'
        log(f'{r.group:12s} {r.pool:10s} {r.pnps:.4f} [{r.lo:.4f}-{r.hi:.4f}]{extra}')
    # beta from A1 (model G, autosomes), if available
    beta = {}
    for f in glob.glob(f'{dc.RES}/A1_dfe_estimates_*.tsv'):
        e = pd.read_csv(f, sep='\t')
        e = e[(e.stratum == 'autosomes') & (e.model == 'G') & (e.setting == 'primary')]
        for r in e.itertuples():
            beta[r.pool] = r.b
    log(f'gamma shape beta from A1 model G (autosomes): {beta if beta else "not available yet"}')
    # windows
    kan = sums['JP-Kan']
    wrows = []
    auto = d.arm.isin(['2L', '2R', '3L', '3R', '3R_scaffold']).values
    for size in (1_000_000, 2_000_000):
        win = (kan.chrom.astype(str) + ':' + (kan.pos // size).astype(int).astype(str)).values
        for w in np.unique(win[auto]):
            m = (win == w) & auto
            r = dict(size_mb=size // 1_000_000, window=w, genes=int(m.sum()))
            for p in dc.POOLS:
                s = sums[p][m].sum(numeric_only=True)
                for c in ('f0', 'f4', 'f4even', 'f4odd', 'si'):
                    r[f'{p}_{c}_n'] = s[f'{c}_n']
                    r[f'{p}_{c}_pi'] = s[f'{c}_pi'] / max(s[f'{c}_n'], 1)
            wrows.append(r)
    wd = pd.DataFrame(wrows)
    wd.to_csv(f'{dc.RES}/A3_pnps_windows.tsv', sep='\t', index=False)
    log('\n== across windows: slope of y = log ratio (Dsub / pool) of piN/piS on x = log ratio of neutral pi')
    srng = np.random.default_rng(dc.SEED + 4)
    for size in (1, 2):
        w = wd[(wd.size_mb == size) & (wd['Dsub_L5_7_si_n'] >= 1500) & (wd['Dsub_L5_7_f4_n'] >= 3000)
               & (wd['JP-Kan_f4_n'] >= 3000)]
        for p in dc.DSUZ:
            for design in ('a_introns', 'b_split'):
                if design == 'a_introns':
                    y = np.log(w['Dsub_L5_7_f0_pi'] / w['Dsub_L5_7_f4_pi']) - np.log(w[f'{p}_f0_pi'] / w[f'{p}_f4_pi'])
                    x = np.log(w['Dsub_L5_7_si_pi']) - np.log(w[f'{p}_si_pi'])
                else:
                    y = np.log(w['Dsub_L5_7_f0_pi'] / w['Dsub_L5_7_f4even_pi']) - np.log(w[f'{p}_f0_pi'] / w[f'{p}_f4even_pi'])
                    x = np.log(w['Dsub_L5_7_f4odd_pi']) - np.log(w[f'{p}_f4odd_pi'])
                ok = np.isfinite(x) & np.isfinite(y)
                x, y = x[ok].values, y[ok].values
                slope = np.polyfit(x, y, 1)[0]
                bs = []
                for _ in range(N_BOOT):
                    i = srng.integers(0, len(x), len(x))
                    bs.append(np.polyfit(x[i], y[i], 1)[0])
                r_ = np.corrcoef(x, y)[0, 1]
                log(f'{size} Mb windows (n {len(x)}), Dsub vs {p}, design {design}: slope {slope:.3f} '
                    f'[{np.percentile(bs, 2.5):.3f} to {np.percentile(bs, 97.5):.3f}], Pearson r {r_:.3f}')
    log(f'finished {time.ctime()} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
