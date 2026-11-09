#!/usr/bin/env python3
"""Inversion analysis, approach 3: genome scan of one pool for blocks of SNPs that share an
intermediate frequency, the signature of an inversion polymorphism in pooled data.

Rationale: SNPs that differ between two arrangements of an inversion at pool frequency f are carried by
the same sampled chromosomes, so they share the frequency f (folded: min(f, 1 - f)) up to read-sampling
noise, across the inverted region and most strongly near the breakpoints. Pool data have no haplotypes,
so this frequency signal is the main SNP-based evidence available.

Sites: P98_RM configuration of scripts/32 (depth max(10, P05)..P98, indel and RepeatMasker masks).
Alleles: >= 3 reads, >= 1 read per strand; minor allele frequency (MAF) = second / (first + second
passing allele) at full depth; SNPs with MAF >= 0.05 are used.
Windows: 100 kb. Per window, counts in 18 MAF bins of width 0.025 (0.05-0.50). Expected counts: the
genome-wide spectrum of the chromosome class, with a per-window log-linear tilt in MAF (Poisson GLM:
log E_b = log(S_w q_b) + a + c m_b), so that broad changes of shape are absorbed and only narrow excess
remains. Band score: z = (O - E) / sqrt(E) for bands of two adjacent bins (width 0.05, 17 bands); peak
band = the band with the largest z.
Runs: >= 3 consecutive windows with peak z >= Z_MIN and peak frequencies within 0.05 of each other.
Clusters (added because differentiation inside old inversions is eroded by gene flux and
the signal need not be contiguous): spans of 10 windows (1 Mb) with >= CL_MIN windows of peak z >= Z_MIN
whose peak frequencies lie within 0.025 of a common band; overlapping spans with the same band merged.
Null: windows shuffled within chromosome (1,000 permutations; seed 20261030) -> expected numbers of runs
and clusters.
Outputs: results/INV_freq_windows_{POOL}.tsv, results/INV_freq_runs_{POOL}.tsv,
         results/INV_freq_scan_{POOL}.log; SNP table QC/inv/{POOL}.snps.npz (for approach 2 tests)
Usage: 110_inversion_freq_scan.py POOL
"""
import glob
import os
import sys
import time

import numpy as np
import pandas as pd
import statsmodels.api as sm

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
A = '/ru-auth/local/home/asato/scratch'
XCHROM = {'NC_050613.1', 'NC_092084.1'}
DOT = {'NC_050608.1', 'NC_092083.1'}
WIN, MIN_MAF, MIN_READS = 100_000, 0.05, 3
EDGES = np.round(np.arange(0.05, 0.5001, 0.025), 4)          # 19 edges -> 18 bins
MID = (EDGES[:-1] + EDGES[1:]) / 2
Z_MIN, RUN_MIN, F_TOL, N_PERM, SEED = 4.0, 3, 0.05, 1000, 20261030
SPAN, CL_MIN = 10, 4


def load_chrom(pool, chrom, dmin, dmax, ms):
    files = sorted(glob.glob(f'{QC}/counts/{pool}/{chrom}_*.npz'), key=lambda f: int(f.split('_')[-2]))
    pos_all, maf_all, n_call = [], [], []
    for f in files:
        d = np.load(f)
        s0, e0 = int(d['start']), int(d['end'])
        m = np.zeros(e0 - s0 + 2, dtype=np.int32)
        for s, e in d['indel_mask']:
            s, e = max(s, s0), min(e, e0)
            if s <= e:
                m[s - s0] += 1
                m[e - s0 + 1] -= 1
        c8 = d['counts']
        tot = c8[:, 0::2].astype(np.int64) + c8[:, 1::2]
        depth = tot.sum(1)
        valid = (depth >= dmin) & (depth <= dmax) & ~(np.cumsum(m)[:-1] > 0) & ~ms[f'{chrom}|rm'][s0:e0 + 1]
        ok = (tot >= MIN_READS) & (c8[:, 0::2] > 0) & (c8[:, 1::2] > 0)
        t = np.sort(np.where(ok, tot, 0), axis=1)
        k1, k2 = t[:, 3], t[:, 2]
        maf = np.where(k1 + k2 > 0, k2 / np.maximum(k1 + k2, 1), 0)
        pos = np.arange(s0, e0 + 1)
        n_call.append(np.bincount((pos[valid] - 1) // WIN, minlength=(e0 - 1) // WIN + 1))
        sel = valid & (maf >= MIN_MAF)
        pos_all.append(pos[sel])
        maf_all.append(maf[sel])
    nwin = max(len(x) for x in n_call)
    ncall = np.zeros(nwin, dtype=np.int64)
    for x in n_call:
        ncall[:len(x)] += x
    return np.concatenate(pos_all), np.concatenate(maf_all), ncall


def window_scores(hist, q):
    """hist (W, 18) observed counts; q (18,) class proportions. Returns E (W, 18), z (W, 17)."""
    W = hist.shape[0]
    E = np.zeros_like(hist, dtype=float)
    X = np.column_stack([np.ones(18), MID])
    for w in range(W):
        S = hist[w].sum()
        if S < 50:
            E[w] = np.nan
            continue
        off = np.log(S * q)
        try:
            f = sm.GLM(hist[w], X, family=sm.families.Poisson(), offset=off).fit()
            E[w] = f.mu
        except Exception:
            E[w] = S * q
    O2 = hist[:, :-1] + hist[:, 1:]
    E2 = E[:, :-1] + E[:, 1:]
    z = (O2 - E2) / np.sqrt(E2)
    return E, z


def runs_of(df):
    """Runs of >= RUN_MIN consecutive windows (per chromosome) with peak z >= Z_MIN and peak frequency
    within F_TOL of the run's first window."""
    out = []
    for chrom, g in df.groupby('chrom', sort=False):
        g = g.sort_values('w')
        cur = []
        for r in g.itertuples():
            hit = r.peak_z >= Z_MIN
            if hit and cur and r.w == cur[-1].w + 1 and abs(r.peak_f - cur[0].peak_f) <= F_TOL:
                cur.append(r)
                continue
            if len(cur) >= RUN_MIN:
                out.append(cur)
            cur = [r] if hit else []
        if len(cur) >= RUN_MIN:
            out.append(cur)
    return out


def clusters_of(df):
    """Merged 1-Mb spans with >= CL_MIN strong windows at a common band (+/- 0.025)."""
    out = []
    bands = np.unique(df.peak_f.dropna())
    for chrom, g in df.groupby('chrom', sort=False):
        strong = g[g.peak_z >= Z_MIN]
        if len(strong) < CL_MIN:
            continue
        wmax = int(g.w.max())
        for f in bands:
            ws = np.sort(strong.w[np.abs(strong.peak_f - f) <= 0.025 + 1e-9].values)
            if len(ws) < CL_MIN:
                continue
            spans = [(i, i + SPAN - 1) for i in range(0, wmax + 1)
                     if ((ws >= i) & (ws <= i + SPAN - 1)).sum() >= CL_MIN]
            merged = []
            for a, b in spans:
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = b
                else:
                    merged.append([a, b])
            for a, b in merged:
                inside = ws[(ws >= a) & (ws <= b)]
                out.append((chrom, float(f), int(inside.min()), int(inside.max()), len(inside)))
    # the same windows can satisfy neighbouring bands: keep the band with most windows per region
    out.sort(key=lambda x: (x[0], x[2], -x[4]))
    keep = []
    for c in out:
        if keep and keep[-1][0] == c[0] and c[2] <= keep[-1][3]:
            continue
        keep.append(c)
    return keep


def main():
    pool = sys.argv[1]
    t0 = time.time()
    sp = 'Dsub' if pool.startswith('Dsub') else 'Dsuz'
    regions = f'{A}/population_polymorphism/reanalysis_20260821/config/regions'
    chroms = {cls: [l.split()[0] for l in open(f'{regions}/{sp}.{cls}.bed')] for cls in ('autosome', 'X')}
    thr = pd.read_csv(f'{RES}/depth_thresholds.tsv', sep='\t')
    thr = thr[thr['sample'] == pool].set_index('cls')
    ms = np.load(f'{QC}/masks_sites/{sp}.npz')
    rows, snp_store = [], {}
    for cls, chs in chroms.items():
        dmin, dmax = int(thr.loc[cls, 'min_depth']), int(thr.loc[cls, 'p98'])
        for chrom in chs:
            if chrom in DOT or not glob.glob(f'{QC}/counts/{pool}/{chrom}_*.npz'):
                continue
            pos, maf, ncall = load_chrom(pool, chrom, dmin, dmax, ms)
            snp_store[f'{chrom}|pos'], snp_store[f'{chrom}|maf'] = pos.astype(np.int32), maf.astype(np.float32)
            w = (pos - 1) // WIN
            b = np.clip(np.searchsorted(EDGES, maf, side='right') - 1, 0, 17)
            hist = np.zeros((len(ncall), 18), dtype=np.int64)
            np.add.at(hist, (w, b), 1)
            for wi in range(len(ncall)):
                rows.append(dict(chrom=chrom, cls=cls, w=wi, start=wi * WIN + 1, callable=int(ncall[wi]),
                                 snps=int(hist[wi].sum()), **{f'b{k}': int(hist[wi, k]) for k in range(18)}))
            print(f'{pool} {chrom}: SNPs (MAF >= {MIN_MAF}) {len(pos)}', flush=True)
    df = pd.DataFrame(rows)
    bcols = [f'b{k}' for k in range(18)]
    zs = []
    for cls, g in df.groupby('cls'):
        h = g[bcols].values
        q = h.sum(0) / h.sum()
        E, z = window_scores(h, q)
        zz = pd.DataFrame(z, index=g.index, columns=[f'z{k}' for k in range(17)])
        zs.append(zz)
    df = df.join(pd.concat(zs))
    zc = [f'z{k}' for k in range(17)]
    band_f = (EDGES[1:-1])                       # band k spans bins k, k+1: centre = EDGES[k+1]
    zv = df[zc].values
    has = np.isfinite(zv).any(1)
    df['peak_z'] = np.where(has, np.nanmax(np.where(np.isfinite(zv), zv, -np.inf), 1), np.nan)
    df['peak_f'] = np.where(has, band_f[np.argmax(np.where(np.isfinite(zv), zv, -np.inf), 1)], np.nan)
    df['snp_density'] = df.snps / df.callable.clip(lower=1)
    os.makedirs(f'{QC}/inv', exist_ok=True)
    np.savez_compressed(f'{QC}/inv/{pool}.snps.npz', **snp_store)
    df.to_csv(f'{RES}/INV_freq_windows_{pool}.tsv', sep='\t', index=False)
    ok = df[df.peak_z.notna()].copy()
    obs = runs_of(ok)
    rng = np.random.default_rng(SEED)
    obs_cl = clusters_of(ok)
    perm, perm_cl = [], []
    for _ in range(N_PERM):
        sh = ok.copy()
        for chrom, g in sh.groupby('chrom'):
            idx = g.index.values
            sh.loc[idx, ['peak_z', 'peak_f']] = g[['peak_z', 'peak_f']].values[rng.permutation(len(idx))]
        perm.append(len(runs_of(sh)))
        perm_cl.append(len(clusters_of(sh)))
    perm, perm_cl = np.array(perm), np.array(perm_cl)
    pd.DataFrame(obs_cl, columns=['chrom', 'band_maf', 'first_w', 'last_w', 'strong_windows']).assign(
        start=lambda x: x.first_w * WIN + 1, end=lambda x: (x.last_w + 1) * WIN).to_csv(
        f'{RES}/INV_freq_clusters_{pool}.tsv', sep='\t', index=False)
    rr = []
    for run in obs:
        rr.append(dict(chrom=run[0].chrom, start=run[0].start, end=run[-1].start + WIN - 1, windows=len(run),
                       mean_peak_f=np.mean([r.peak_f for r in run]), min_z=min(r.peak_z for r in run),
                       max_z=max(r.peak_z for r in run)))
    rr = pd.DataFrame(rr)
    rr.to_csv(f'{RES}/INV_freq_runs_{pool}.tsv', sep='\t', index=False)
    with open(f'{RES}/INV_freq_scan_{pool}.log', 'w') as f:
        f.write(f'{pool}: windows with >= 50 SNPs {len(ok)}; windows with peak z >= {Z_MIN}: '
                f'{(ok.peak_z >= Z_MIN).sum()} ({(ok.peak_z >= Z_MIN).mean():.3f})\n')
        f.write(f'runs of >= {RUN_MIN} consecutive windows at a consistent frequency: observed {len(obs)}; '
                f'expected after shuffling windows within chromosome {perm.mean():.2f} (95% {np.percentile(perm, 2.5):.0f}-'
                f'{np.percentile(perm, 97.5):.0f}); P(perm >= obs) {(perm >= len(obs)).mean():.3f}\n')
        for r in rr.itertuples():
            f.write(f'  run {r.chrom}:{r.start}-{r.end} windows {r.windows} MAF {r.mean_peak_f:.3f} z {r.min_z:.1f}-{r.max_z:.1f}\n')
        f.write(f'clusters (>= {CL_MIN} strong windows at a common band within {SPAN} windows): observed {len(obs_cl)}; '
                f'expected {perm_cl.mean():.2f} (95% {np.percentile(perm_cl, 2.5):.0f}-{np.percentile(perm_cl, 97.5):.0f}); '
                f'P(perm >= obs) {(perm_cl >= len(obs_cl)).mean():.3f}\n')
        for c in obs_cl:
            f.write(f'  cluster {c[0]}:{c[2] * WIN + 1}-{(c[3] + 1) * WIN} band MAF {c[1]:.3f} strong windows {c[4]}\n')
        f.write(f'finished {time.ctime()} ({time.time() - t0:.0f} s)\n')
    print(open(f'{RES}/INV_freq_scan_{pool}.log').read())


if __name__ == '__main__':
    main()
