#!/usr/bin/env python3
"""Purity of the four D. suzukii pools (analysis 4.2). Mirror of 65_summarize_contamination.py.

Method 1, nuclear species-specific derived sites (77_dsuz_pool_purity_sites.py): at CDS sites where
the D. subpulchrella base differs from D. suzukii while D. biarmipes matches D. suzukii, the rate of
the D. subpulchrella base among pool reads minus the error rate to one specific base (half of the
reads with either remaining base); sites with > 20% D. subpulchrella base are counted separately.
95% CI from 2,000 bootstraps over 1-Mb blocks. Zero-class estimator (as in script 65): uniform
contamination at rate c gives Poisson(c * depth) species-base reads per site, so the fraction of
sites with zero such reads determines lambda, an upper bound for the uniform component; lambda of
the D. biarmipes- and D. melanogaster-specific control sites is subtracted.
Method 3, mitochondrial diagnostic sites (78_mito_purity.sh): positions of the D. suzukii mtDNA
(NC_060762.1) where the D. subpulchrella pool carries a different base at >= 95% (depth >= 20);
frequency of that base in each D. suzukii pool (depth >= 100), minus the error rate to one base.
Fly numbers per pool (haploid autosomal sample / 2, script 32): Kan 212, Sap 50, Sap2 151, Tok 50.
Seed 20261023. Writes results/purity_dsuz_*.tsv
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
MITO = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/mito'
POPS = ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')
FLIES = {'JP-Kan': 212, 'JP-Sap': 50, 'JP-Sap2': 151, 'JP-Tok': 50}
XCHROM = 'NC_092084.1'
rng = np.random.default_rng(20261023)
pd.set_option('display.width', 250)


def est(d):
    e = d.n_other_two_bases.sum() / 2 / d.depth.sum()
    return d.n_species_allele.sum() / d.depth.sum() - e, e


rows, zc, pc = [], [], []
for pop in POPS:
    c = pd.read_csv(f'{P}/results/purity_dsuz/{pop}.sites.tsv.gz', sep='\t')
    c['freq'] = c.n_species_allele / c.depth
    c['cls'] = np.where(c.chrom == XCHROM, 'X', 'autosome')
    c['block'] = c.chrom + ':' + (c.pos // 1_000_000).astype(str)
    for (sp, cls), d0 in c.groupby(['species', 'cls']):
        hi = d0.freq > 0.2
        d = d0[~hi]
        cc, e = est(d)
        blocks = d.groupby('block')[['n_species_allele', 'n_other_two_bases', 'depth']].sum()
        bs = []
        for _ in range(2000):
            b = blocks.iloc[rng.integers(0, len(blocks), len(blocks))]
            bs.append(b.n_species_allele.sum() / b.depth.sum() - b.n_other_two_bases.sum() / 2 / b.depth.sum())
        dep = d0.depth.values
        n0 = int((d0.n_species_allele == 0).sum())
        lam = brentq(lambda l: np.exp(-l * dep).sum() - n0, 1e-9, 1.0)
        rows.append(dict(pool=pop, species=sp, cls=cls, sites=len(d0), sites_freq_gt_0p2=int(hi.sum()),
                         mean_depth=round(d.depth.mean(), 1),
                         species_base_rate_pct=round(100 * d.n_species_allele.sum() / d.depth.sum(), 3),
                         error_rate_one_base_pct=round(100 * e, 3),
                         species_base_excess_pct=round(100 * cc, 3),
                         ci_lo_pct=round(100 * np.percentile(bs, 2.5), 3),
                         ci_hi_pct=round(100 * np.percentile(bs, 97.5), 3),
                         excess_fly_equivalents=round(cc * FLIES[pop], 2),
                         zero_class_obs=round(n0 / len(d0), 4), lambda_pct=round(100 * lam, 4),
                         zero_class_if_one_fly=round(np.exp(-(lam + 1 / FLIES[pop]) * dep).mean(), 4)))
    d = c[(c.species == 'dsub') & (c.freq <= 0.2)]
    for chrom, g in d.groupby('chrom'):
        cc, _ = est(g)
        pc.append(dict(pool=pop, chrom=chrom, sites=len(g), species_base_excess_pct=round(100 * cc, 3)))
    hist = np.histogram(c[c.species == 'dsub'].freq, bins=[0, 1e-9, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 1.0001])[0]
    print(f'{pop} D. subpulchrella-base frequency at dsub sites [0, (0-1%], (1-2%], (2-3%], (3-5%], (5-10%], '
          f'(10-20%], (20-100%]]: {hist.tolist()}')
t = pd.DataFrame(rows)
for (pop, cls), g in t.groupby(['pool', 'cls']):
    ctrl = g[g.species != 'dsub'].lambda_pct.max()
    t.loc[g.index[g.species == 'dsub'], 'lambda_excess_over_control_pct'] = g[g.species == 'dsub'].lambda_pct - ctrl
t.to_csv(f'{P}/results/purity_dsuz_nuclear.tsv', sep='\t', index=False)
print('\n== Method 1: species-specific derived bases in the D. suzukii pools')
print(t.to_string(index=False))
pcs = pd.DataFrame(pc)
pcs.to_csv(f'{P}/results/purity_dsuz_by_chromosome.tsv', sep='\t', index=False)
print('\n== D. subpulchrella-base excess (%) per chromosome')
print(pcs.pivot(index='chrom', columns='pool', values='species_base_excess_pct').to_string())


# ------------------------------------------------------------------ Method 3: mitochondrial sites
def pileup_counts(path):
    out = {}
    for line in open(path):
        chrom, pos, ref, depth, bases = line.rstrip('\n').split('\t')[:5]
        ref = ref.upper()
        s, i, n = bases, 0, {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        while i < len(s):
            ch = s[i]
            if ch == '^':
                i += 2
                continue
            if ch in '+-':
                j = i + 1
                while j < len(s) and s[j].isdigit():
                    j += 1
                i = j + int(s[i + 1:j])
                continue
            if ch in '.,':
                n[ref] += 1
            elif ch.upper() in n:
                n[ch.upper()] += 1
            i += 1
        out[int(pos)] = (ref, n)
    return out


sub = pileup_counts(f'{MITO}/Dsub_L5_7.mito.pileup')
diag = {}
for pos, (ref, n) in sub.items():
    d = sum(n.values())
    if d < 20 or ref not in 'ACGT':
        continue
    b = max(n, key=n.get)
    if b != ref and n[b] / d >= 0.95:
        diag[pos] = (ref, b, d)
print(f'\n== Method 3: mitochondrial diagnostic sites: {len(diag)} (D. subpulchrella pool depth median '
      f'{np.median([v[2] for v in diag.values()]) if diag else 0:.0f})')
mrows, msites = [], []
for pop in POPS:
    pc_ = pileup_counts(f'{MITO}/{pop}.mito.pileup')
    for pos, (ref, b, _) in diag.items():
        if pos not in pc_:
            continue
        n = pc_[pos][1]
        d = sum(n.values())
        if d < 100:
            continue
        other = d - n[ref] - n[b]
        msites.append(dict(pool=pop, pos=pos, depth=d, n_dsuz=n[ref], n_dsub=n[b], n_other=other,
                           dsub_freq=n[b] / d, err_one_base=other / 2 / d))
m = pd.DataFrame(msites)
m.to_csv(f'{P}/results/purity_dsuz_mito_sites.tsv', sep='\t', index=False)
for pop, g in m.groupby('pool'):
    ex = g.dsub_freq - g.err_one_base
    bs = [ex.iloc[rng.integers(0, len(ex), len(ex))].median() for _ in range(2000)]
    mrows.append(dict(pool=pop, sites=len(g), median_depth=int(g.depth.median()),
                      dsub_base_median_pct=round(100 * g.dsub_freq.median(), 3),
                      dsub_base_pooled_pct=round(100 * g.n_dsub.sum() / g.depth.sum(), 3),
                      error_one_base_pct=round(100 * g.err_one_base.mean(), 3),
                      excess_median_pct=round(100 * ex.median(), 3),
                      ci_lo_pct=round(100 * np.percentile(bs, 2.5), 3), ci_hi_pct=round(100 * np.percentile(bs, 97.5), 3),
                      sites_dsub_freq_gt_10pct=int((g.dsub_freq > 0.10).sum()),
                      iqr_pct=f'{100 * g.dsub_freq.quantile(0.25):.2f}-{100 * g.dsub_freq.quantile(0.75):.2f}'))
mt = pd.DataFrame(mrows)
# Sites where the D. suzukii reference mtDNA carries a base that is rare in the Japanese pools are not
# diagnostic (the D. subpulchrella base is then common in every pool), so the median above is not an
# estimate of contamination. Two robust statistics:
# floor: uniform contamination adds the same fraction to every diagnostic site, so the 10th
#   percentile and the minimum over sites are lower bounds for it (sampling SD at depth 6,000 and 1%:
#   0.13 percentage points);
# leave-one-out: sites where the D. subpulchrella base is below 0.5% in at least two of the other
#   three pools; median frequency in the focal pool minus the error rate to one base.
w = m.pivot(index='pos', columns='pool', values='dsub_freq')
for i, pop in enumerate(mt.pool):
    g = m[m.pool == pop].set_index('pos')
    mt.loc[i, 'floor_min_pct'] = round(100 * g.dsub_freq.min(), 3)
    mt.loc[i, 'floor_p10_pct'] = round(100 * g.dsub_freq.quantile(0.10), 3)
    others = [p for p in POPS if p != pop]
    keep = (w[others] < 0.005).sum(axis=1) >= 2
    gk = g.loc[keep[keep].index]
    exk = gk.dsub_freq - gk.err_one_base
    bs = [exk.iloc[rng.integers(0, len(exk), len(exk))].median() for _ in range(2000)]
    mt.loc[i, 'loo_sites'] = len(gk)
    mt.loc[i, 'loo_excess_median_pct'] = round(100 * exk.median(), 3)
    mt.loc[i, 'loo_ci'] = f'{100 * np.percentile(bs, 2.5):.3f}-{100 * np.percentile(bs, 97.5):.3f}'
mt.to_csv(f'{P}/results/purity_dsuz_mito.tsv', sep='\t', index=False)
print(mt.to_string(index=False))
