#!/usr/bin/env python3
"""Contamination of the D. subpulchrella pool (Dsub_L5_7).

Method 1, allele frequencies at species-specific derived sites (from 60_polarized_mk.py):
  at callable CDS positions where the D. suzukii base differs from D. subpulchrella while
  D. biarmipes matches D. subpulchrella (D. suzukii-specific), and analogous sites for
  D. biarmipes and D. melanogaster. Rate of the species-specific base among pool reads,
  minus the sequencing-error rate to one specific base (half of the reads carrying either
  of the two remaining bases). Sites where the species base exceeds 20% of reads are counted
  separately (shared polymorphism or paralog mismapping, not uniform contamination). The pooled
  excess (species_base_excess_pct) mixes shared polymorphism with contamination; the zero-class
  estimator below gives the upper bound for uniform contamination.
  95% CI from 2,000 bootstraps over 1-Mb blocks; per-chromosome estimates test uniformity;
  X versus autosomes indicates the sex of contaminating flies (all pool flies are female).
Method 2, competitive mapping (63_competitive_mapping.sh): fraction of MAPQ >= 20 reads
  assigned to each species, corrected with misassignment rates from simulated reads.
Writes results/L3_contamination_*.tsv.
"""
import numpy as np
import pandas as pd

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
XCHROM = 'NC_050613.1'
FLIES = 47
rng = np.random.default_rng(20261021)
c = pd.read_csv(f'{P}/results/polarized_mk/Dsub_L5_7.contamination_sites.tsv.gz', sep='\t')
c['freq'] = c.n_species_allele / c.depth
c['cls'] = np.where(c.chrom == XCHROM, 'X', 'autosome')
c['block'] = c.chrom + ':' + (c.pos // 1_000_000).astype(str)


def est(d):
    e = d.n_other_two_bases.sum() / 2 / d.depth.sum()
    return d.n_species_allele.sum() / d.depth.sum() - e, e


rows = []
for (sp, cls), d0 in c.groupby(['species', 'cls']):
    hi = d0.freq > 0.2
    d = d0[~hi]
    cc, e = est(d)
    blocks = d.groupby('block')[['n_species_allele', 'n_other_two_bases', 'depth']].sum()
    bs = []
    for _ in range(2000):
        b = blocks.iloc[rng.integers(0, len(blocks), len(blocks))]
        bs.append(b.n_species_allele.sum() / b.depth.sum() - b.n_other_two_bases.sum() / 2 / b.depth.sum())
    # fraction of sites with >= 1 species-base read: observed vs expected without contamination
    p0 = 1 - (1 - e) ** d.depth
    pc = 1 - (1 - e - max(cc, 0)) ** d.depth
    rows.append(dict(species=sp, cls=cls, sites=len(d0), sites_freq_gt_0p2=int(hi.sum()),
                     mean_depth=round(d.depth.mean(), 1),
                     species_base_rate_pct=round(100 * d.n_species_allele.sum() / d.depth.sum(), 3),
                     error_rate_one_base_pct=round(100 * e, 3),
                     species_base_excess_pct=round(100 * cc, 3),
                     ci_lo_pct=round(100 * np.percentile(bs, 2.5), 3),
                     ci_hi_pct=round(100 * np.percentile(bs, 97.5), 3),
                     excess_fly_equivalents=round(cc * FLIES, 2),
                     sites_with_species_read_obs=round((d.n_species_allele > 0).mean(), 3),
                     expected_if_error_only=round(p0.mean(), 3),
                     expected_if_uniform_contamination=round(pc.mean(), 3)))
t = pd.DataFrame(rows)
t.to_csv(f'{P}/results/L3_contamination_allele_freq.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print('== Method 1: species-specific derived bases in the D. subpulchrella pool')
print(t.to_string(index=False))

# Zero-class estimator. Uniform contamination at rate c (plus a uniform error rate) gives each
# site Poisson(lambda * depth) species-base reads, so the number of sites with zero such reads
# determines lambda. Polymorphic sites only remove sites from the zero class, so lambda is an
# upper bound for the uniform component. The D. biarmipes- and D. melanogaster-specific site
# sets are controls (errors and homoplasy, no plausible contamination); lambda(dsuz) minus
# lambda(control) bounds uniform D. suzukii contamination. Also shown: the zero-class fraction
# expected if one of the 47 flies were D. suzukii (c = 1/47).
from scipy.optimize import brentq
zc = []
for (sp, cls), d0 in c.groupby(['species', 'cls']):
    dep = d0.depth.values
    n0 = int((d0.n_species_allele == 0).sum())
    lam = brentq(lambda l: np.exp(-l * dep).sum() - n0, 1e-9, 1.0)
    zc.append(dict(species=sp, cls=cls, sites=len(d0), zero_class_obs=round(n0 / len(d0), 4),
                   lambda_upper_pct=round(100 * lam, 4),
                   zero_class_if_one_fly=round(np.exp(-(lam + 1 / FLIES) * dep).mean(), 4),
                   zero_class_if_half_fly=round(np.exp(-(lam + 0.5 / FLIES) * dep).mean(), 4)))
z = pd.DataFrame(zc)
for cls in ('autosome', 'X'):
    ctrl = z[(z.cls == cls) & (z.species != 'dsuz')].lambda_upper_pct.max()
    z.loc[(z.cls == cls) & (z.species == 'dsuz'), 'excess_over_control_pct'] = \
        z[(z.cls == cls) & (z.species == 'dsuz')].lambda_upper_pct - ctrl
z['excess_fly_equivalents'] = (z.excess_over_control_pct / 100 * FLIES).round(3)
z.to_csv(f'{P}/results/L3_contamination_zero_class.tsv', sep='\t', index=False)
print('\n== Zero-class estimator (upper bound for uniform contamination; controls dbia, dmel)')
print(z.to_string(index=False))

d = c[(c.species == 'dsuz') & (c.freq <= 0.2)]
pc = []
for chrom, g in d.groupby('chrom'):
    cc, e = est(g)
    pc.append(dict(chrom=chrom, sites=len(g), species_base_excess_pct=round(100 * cc, 3)))
blk = d.groupby('block').apply(lambda g: est(g)[0] if len(g) >= 50 else np.nan).dropna()
pcs = pd.DataFrame(pc)
pcs.to_csv(f'{P}/results/L3_contamination_by_chromosome.tsv', sep='\t', index=False)
print('\n== D. suzukii-specific base: per chromosome')
print(pcs.to_string(index=False))
print(f'1-Mb blocks with >= 50 sites: n = {len(blk)}, species-base excess % median {100 * blk.median():.3f}, '
      f'10th-90th percentile {100 * blk.quantile(0.1):.3f}-{100 * blk.quantile(0.9):.3f}')
hist = np.histogram(d.freq, bins=[0, 1e-9, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2])[0]
print('sites by D. suzukii-base frequency [0, (0-1%], (1-2%], (2-3%], (3-5%], (5-10%], (10-20%]]:', hist)

# Method 2: competitive mapping
try:
    m = pd.read_csv(f'{P}/results/compmap_tally.tsv', sep='\t', header=None,
                    names=['sample', 'what', 'species', 'n'], engine='python')
    q = m[m.what == 'MAPQ>=20'].pivot_table(index='sample', columns='species', values='n', fill_value=0)
    frac = q.div(q.sum(axis=1), axis=0)
    print('\n== Method 2: fraction of MAPQ >= 20 primary reads by species')
    print((100 * frac).round(3).to_string())
    r_sub = frac.loc['simDsub', 'Dsuz']           # Dsub reads misassigned to Dsuz
    r_suz = frac.loc['simDsuz', 'Dsub']           # Dsuz reads misassigned to Dsub
    out = []
    for lane in ('L5', 'L7'):
        f = frac.loc[lane, 'Dsuz'] / (frac.loc[lane, 'Dsub'] + frac.loc[lane, 'Dsuz'])
        corr = (f - r_sub) / (1 - r_sub - r_suz)
        # Upper bound only: simulated reads come from the Dsub assembly, so the correction cannot
        # remove real Dsub reads from sequence missing in that assembly, whose closest match is
        # the D. suzukii homolog (see 66_compmap_diagnostics.sh: edit distances, COI haplotypes).
        out.append(dict(lane=lane, dsuz_fraction_of_dsub_plus_dsuz_pct=round(100 * f, 3),
                        upper_bound_pct=round(100 * corr, 3),
                        upper_bound_fly_equivalents=round(corr * FLIES, 2),
                        dbia_pct=round(100 * frac.loc[lane].get('Dbia', 0), 3),
                        dmel_pct=round(100 * frac.loc[lane].get('Dmel', 0), 3)))
    o = pd.DataFrame(out)
    o.to_csv(f'{P}/results/L3_contamination_compmap.tsv', sep='\t', index=False)
    print(f'misassignment from simulations: Dsub->Dsuz {100 * r_sub:.3f}%, Dsuz->Dsub {100 * r_suz:.3f}%')
    print(o.to_string(index=False))
except FileNotFoundError:
    print('competitive mapping results not available yet')
