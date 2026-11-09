#!/usr/bin/env python3
"""How much do D. subpulchrella reads in a D. suzukii pool inflate its neutral diversity? (4.2 follow-up)

At four-fold sites (D. suzukii site mask, autosomes) where D. subpulchrella carries a derived base
(77_dsuz_pool_purity_sites.py, class 'dsub'; CDS of gene trios only), per-site pi and Watterson
theta are computed with the estimator of 32_chunk_stats.py (alleles need >= 2 reads; pool sizes as
in script 32) from the counts of the D. suzukii reference base and the D. subpulchrella base (reads
with the two other bases are ignored). Their sum per valid four-fold site is compared with the pool's
total four-fold pi and theta (P99_noRM blocks of script 32). Contamination adds the difference between
a contaminated pool and clean pools of similar pool size and depth (JP-Sap: n = 100 like JP-Tok).
The site set covers only CDS of genes with alignments, so the estimate is a lower bound for the
contribution at all four-fold sites of these genes; divergent sites outside these CDS (introns,
intergenic sequence) are affected in the same way but are not counted here.
Output: results/purity_dsuz_diversity_effect.tsv
"""
import glob
import importlib.util
import os

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
spec = importlib.util.spec_from_file_location('cs', f'{PROJ}/scripts/32_chunk_stats.py')
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)
POPS = ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')
ms = np.load(f'{QC}/masks_sites/Dsuz.npz')
b = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob(f'{QC}/stats/div/JP-*.blocks.tsv.gz')])
b = b[(b.config == 'P99_noRM') & (b.site_class == 'fourfold') & (b.cls == 'autosome')]
tot = b.groupby('sample')[['n_valid', 'sum_pi', 'sum_tw']].sum()

rows = []
for pop in POPS:
    f = f'{PROJ}/results/purity_dsuz/{pop}.sites.tsv.gz'
    if not os.path.exists(f):
        print('not yet available:', pop)
        continue
    c = pd.read_csv(f, sep='\t')
    c = c[(c.species == 'dsub') & (c.chrom != 'NC_092084.1')]
    ff = np.concatenate([ms[f'{ch}|ff'][g.pos.values] for ch, g in c.groupby('chrom', sort=False)])
    c = pd.concat([g for _, g in c.groupby('chrom', sort=False)])[ff]
    n = cs.POOL[pop][0]
    dmin, dmax = cs.thresholds(pop, 'autosome')
    c = c[(c.depth >= dmin) & (c.depth <= dmax['p99'])]
    cw, cp = cs.corrections(n, dmax['p99'])
    a = np.where(c.n_focal_allele >= cs.B, c.n_focal_allele, 0)
    y = np.where(c.n_species_allele >= cs.B, c.n_species_allele, 0)
    M = a + y
    snp = (a > 0) & (y > 0)
    Ms = np.where(snp, M, 4)
    h = Ms / (Ms - 1) * (1 - (a / np.maximum(M, 1)) ** 2 - (y / np.maximum(M, 1)) ** 2)
    pi = np.where(snp, h / cp[np.minimum(Ms, len(cp) - 1)], 0).sum()
    tw = np.where(snp, 1 / cw[np.minimum(Ms, len(cw) - 1)], 0).sum()
    t = tot.loc[pop]
    rows.append(dict(pool=pop, dsub_fourfold_sites=len(c), snp_fraction_at_these_sites=round(snp.mean(), 4),
                     pi_share_pct=round(100 * pi / t.sum_pi, 2), theta_share_pct=round(100 * tw / t.sum_tw, 2),
                     fourfold_pi=round(t.sum_pi / t.n_valid, 5), fourfold_theta=round(t.sum_tw / t.n_valid, 5)))
r = pd.DataFrame(rows)
for k in ('pi_share_pct', 'theta_share_pct', 'snp_fraction_at_these_sites'):
    clean = r.loc[r.pool == 'JP-Sap', k].iloc[0]
    r[f'{k}_minus_JP-Sap'] = (r[k] - clean).round(3)
r.to_csv(f'{PROJ}/results/purity_dsuz_diversity_effect.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(r.to_string(index=False))
