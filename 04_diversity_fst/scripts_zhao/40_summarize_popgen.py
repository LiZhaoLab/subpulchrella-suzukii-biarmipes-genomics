#!/usr/bin/env python3
"""Summarize my independent diversity, Ne and FST estimates and compare with the manuscript.

Inputs: QC_DIR/stats/div/*.tsv.gz and QC_DIR/stats/fst/*.tsv.gz from 32_chunk_stats.py.
Outputs (results/):
  L2_diversity_means.tsv   unweighted means of finite 20-kb window estimates (pi, thetaW, D)
  L2_Ne.tsv                callable-weighted theta per site, Ne (thetaW- and pi-based),
                           95% CI from 10,000 chromosome-stratified 1-Mb block bootstraps
  L2_fst_means.tsv         mean Hudson FST of finite windows, per pair and class
  L2_highfst_rbh.tsv       top 1% / 5% FST windows: overlapping genes, % with direct Dmel RBH
"""
import glob

import numpy as np
import pandas as pd

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
A = '/ru-auth/local/home/asato/scratch'
MU = 3.3e-9
SEED = 20261021
NBOOT = 10000
pd.set_option('display.width', 250)


def tajima_d(sum_pi, sum_tw, S, n):
    """Tajima's D with the classical variance, n = empirical minimum read depth of the window."""
    n = np.asarray(n, dtype=float)
    out = np.full(len(n), np.nan)
    for i, (p, t, s, nn) in enumerate(zip(sum_pi, sum_tw, S, n)):
        if s < 1 or nn < 2:
            continue
        k = np.arange(1, int(nn))
        a1, a2 = (1 / k).sum(), (1 / k ** 2).sum()
        b1 = (nn + 1) / (3 * (nn - 1))
        b2 = 2 * (nn ** 2 + nn + 3) / (9 * nn * (nn - 1))
        c1 = b1 - 1 / a1
        c2 = b2 - (nn + 2) / (a1 * nn) + a2 / a1 ** 2
        e1, e2 = c1 / a1, c2 / (a1 ** 2 + a2)
        out[i] = (p - t) / np.sqrt(e1 * s + e2 * s * (s - 1))
    return out


# ---------------- diversity windows ----------------
w = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{QC}/stats/div/*.windows.tsv.gz'))
w = w[w.n_valid > 0].copy()
w['pi'] = w.sum_pi / w.n_valid
w['tw'] = w.sum_tw / w.n_valid
w['D'] = tajima_d(w.sum_pi.values, w.sum_tw.values, w.S.values, w.min_depth.values)
div = (w.groupby(['config', 'cls', 'sample'])
       .agg(windows=('pi', 'size'), pi=('pi', 'mean'), tw=('tw', 'mean'),
            D=('D', 'mean'), D_windows=('D', 'count')).reset_index())
div.to_csv(f'{RES}/L2_diversity_means.tsv', sep='\t', index=False)
print('== window means (unweighted, finite windows)')
print(div.round(5).to_string(index=False))

# ---------------- Ne from neutral sites ----------------
b = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{QC}/stats/div/*.blocks.tsv.gz'))
rng = np.random.default_rng(SEED)
rows = []
for (cfg, cls, sample, sc), g in b.groupby(['config', 'cls', 'sample', 'site_class']):
    g = g[g.n_valid > 0]
    nv, stw, spi = g.n_valid.values, g.sum_tw.values, g.sum_pi.values
    tw = stw.sum() / nv.sum()
    pi = spi.sum() / nv.sum()
    factor = 3 if cls == 'X' else 4
    boot_tw, boot_pi = np.empty(NBOOT), np.empty(NBOOT)
    idx_by_chrom = [np.flatnonzero(g.chrom.values == c) for c in np.unique(g.chrom.values)]
    for r in range(NBOOT):
        pick = np.concatenate([rng.choice(ix, size=len(ix), replace=True) for ix in idx_by_chrom])
        boot_tw[r] = stw[pick].sum() / nv[pick].sum()
        boot_pi[r] = spi[pick].sum() / nv[pick].sum()
    rows.append(dict(config=cfg, cls=cls, sample=sample, site_class=sc, blocks=len(g),
                     valid_sites=int(nv.sum()), snps=int(g.S.sum()), pi=pi, tw=tw,
                     Ne_tw=tw / (factor * MU),
                     Ne_tw_lo=np.percentile(boot_tw, 2.5) / (factor * MU),
                     Ne_tw_hi=np.percentile(boot_tw, 97.5) / (factor * MU),
                     Ne_pi=pi / (factor * MU),
                     Ne_pi_lo=np.percentile(boot_pi, 2.5) / (factor * MU),
                     Ne_pi_hi=np.percentile(boot_pi, 97.5) / (factor * MU)))
ne = pd.DataFrame(rows)
ne.to_csv(f'{RES}/L2_Ne.tsv', sep='\t', index=False)
print('\n== Ne (mu = 3.3e-9; autosome 4mu, X 3mu); seed', SEED)
show = ne[ne.site_class != 'all'].copy()
for c in ('Ne_tw', 'Ne_tw_lo', 'Ne_tw_hi', 'Ne_pi', 'Ne_pi_lo', 'Ne_pi_hi'):
    show[c] = (show[c] / 1e6).round(2)
print(show.drop(columns=['blocks']).round(5).to_string(index=False))
print('\n== genome-wide callable-weighted theta (site_class = all)')
print(ne[ne.site_class == 'all'][['config', 'cls', 'sample', 'valid_sites', 'pi', 'tw']]
      .round(5).to_string(index=False))

# ---------------- FST ----------------
f = pd.concat(pd.read_csv(x, sep='\t') for x in glob.glob(f'{QC}/stats/fst/*.windows.tsv.gz'))
f = f[f.sum_pb > 0].copy()
f['fst'] = 1 - (f.sum_pw1 + f.sum_pw2) / 2 / f.sum_pb
fm = f.groupby(['config', 'cls', 'pair']).agg(windows=('fst', 'size'), fst=('fst', 'mean')).reset_index()
fm.to_csv(f'{RES}/L2_fst_means.tsv', sep='\t', index=False)
print('\n== mean Hudson FST of finite windows')
print(fm.pivot_table(index=['cls', 'pair'], columns='config', values='fst').round(5).to_string())

# ---------------- high-FST windows: genes and direct RBH ----------------
genes = []
for line in open(f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'):
    t = line.split('\t')
    if len(t) > 8 and t[2] == 'gene' and 'gene_biotype=protein_coding' in t[8]:
        gid = t[8].split('ID=')[1].split(';')[0]
        genes.append((t[0], int(t[3]), int(t[4]), gid))
genes = pd.DataFrame(genes, columns=['chrom', 'gene_start', 'gene_end', 'gene'])
MK = f'{A}/MKtest/reanalysis_longest_fix_20260817'
gene_of_prot = {}
for line in open(f'{MK}/pep_longest/Dsuz.IsoJpt1.0.longest.pep.fa'):
    if line.startswith('>'):
        p, g = line[1:].split()[:2]
        gene_of_prot[p] = g.replace('gene=', '')
rbh = {gene_of_prot[p] for p in pd.read_csv(f'{MK}/orth_longest/Dsuz_Dmel.direct.RBH.longest.protein.tsv',
                                            sep='\t', header=None)[1]}
rows = []
for (cfg, cls, pair), g in f.groupby(['config', 'cls', 'pair']):
    elig = g[g.n_snp >= np.percentile(g.n_snp, 10)]
    for top in (1, 5):
        k = int(np.ceil(top / 100 * len(elig)))
        hi = elig.nlargest(k, 'fst')
        hit = set()
        for _, r in hi.iterrows():
            gg = genes[(genes.chrom == r.chrom) & (genes.gene_start <= r['end']) & (genes.gene_end >= r['start'])]
            hit |= set(gg.gene)
        rows.append(dict(config=cfg, cls=cls, pair=pair, top=top, windows=k, genes=len(hit),
                         rbh_genes=len(hit & rbh)))
hf = pd.DataFrame(rows)
agg = hf.groupby(['config', 'top', 'pair'])[['windows', 'genes', 'rbh_genes']].sum().reset_index()
agg['rbh_pct'] = (100 * agg.rbh_genes / agg.genes).round(1)
agg.to_csv(f'{RES}/L2_highfst_rbh.tsv', sep='\t', index=False)
print('\n== high-FST windows (autosome + X): protein-coding genes with a direct Dmel RBH '
      f'(genome baseline {100 * len(rbh & set(genes.gene)) / len(set(genes.gene)):.1f}% of '
      f'{len(set(genes.gene))} protein-coding genes)')
print(agg.to_string(index=False))
