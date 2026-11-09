#!/usr/bin/env python3
"""Does genome composition explain why genome-wide pi is similar between species while
neutral-site pi is ~20% lower in D. subpulchrella?

Uses 1-Mb blocks from 32_chunk_stats.py (all sites). Gene density of a block = number of
annotated four-fold sites (valid or not is irrelevant here; we use valid four-fold sites
of the same sample, which tracks gene density). Reports, per sample and config:
  callable sites and callable-weighted pi in gene-poor (< 5,000 four-fold sites per Mb)
  and gene-rich (>= 5,000) blocks, and pi by gene-density quartile (quartiles pooled
  across samples of both species).
Writes results/L2_composition_test.tsv.
"""
import glob

import numpy as np
import pandas as pd

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
b = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{QC}/stats/div/*.blocks.tsv.gz'))
b = b[b.cls == 'autosome']
key = ['config', 'sample', 'chrom', 'block']
allb = b[b.site_class == 'all'].set_index(key)
ffb = b[b.site_class == 'fourfold'].set_index(key)
d = allb[['n_valid', 'sum_pi', 'sum_tw']].join(ffb[['n_valid']].rename(columns={'n_valid': 'ff_valid'}))
d = d.reset_index()
d = d[d.n_valid > 0]
d['class'] = np.where(d.ff_valid >= 5000, 'gene_rich', 'gene_poor')
rows = []
for (cfg, s, c), g in d.groupby(['config', 'sample', 'class']):
    rows.append(dict(config=cfg, sample=s, block_class=c, blocks=len(g),
                     callable_Mb=round(g.n_valid.sum() / 1e6, 1),
                     pi=g.sum_pi.sum() / g.n_valid.sum(), tw=g.sum_tw.sum() / g.n_valid.sum()))
t = pd.DataFrame(rows)
q = np.quantile(d[d.config == 'P99_noRM'].ff_valid, [0.25, 0.5, 0.75])
d['ffq'] = np.digitize(d.ff_valid, q)
for (cfg, s, k), g in d.groupby(['config', 'sample', 'ffq']):
    rows.append(dict(config=cfg, sample=s, block_class=f'ff_quartile_{k + 1}', blocks=len(g),
                     callable_Mb=round(g.n_valid.sum() / 1e6, 1),
                     pi=g.sum_pi.sum() / g.n_valid.sum(), tw=g.sum_tw.sum() / g.n_valid.sum()))
t = pd.DataFrame(rows)
t.to_csv(f'{RES}/L2_composition_test.tsv', sep='\t', index=False)
pd.set_option('display.width', 200)
x = t[t.config == 'P99_noRM']
print('four-fold-site quartile cut points (sites per Mb):', q)
print(x.pivot_table(index='block_class', columns='sample', values='pi').round(5).to_string())
print(x.pivot_table(index='block_class', columns='sample', values='callable_Mb').to_string())

# Window means (the manuscript's statistic) restricted to windows inside gene-rich blocks.
w = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{QC}/stats/div/*.windows.tsv.gz'))
w = w[(w.cls == 'autosome') & (w.n_valid > 0)].copy()
w['block'] = (w.start - 1) // 1_000_000
w = w.merge(d[['config', 'sample', 'chrom', 'block', 'class']], on=['config', 'sample', 'chrom', 'block'])
w['pi'] = w.sum_pi / w.n_valid
w['tw'] = w.sum_tw / w.n_valid
src = open('/ru-auth/local/home/lzhao/vscode_projects/suzukii/scripts/40_summarize_popgen.py').read()
ns = {'np': np}
exec(src[src.index('def tajima_d'):src.index('# ---------------- diversity windows')], ns)
w['D'] = ns['tajima_d'](w.sum_pi.values, w.sum_tw.values, w.S.values, w.min_depth.values)
wm = w.groupby(['config', 'class', 'sample'])[['pi', 'tw', 'D']].mean().reset_index()
wm.to_csv(f'{RES}/L2_composition_windows.tsv', sep='\t', index=False)
print('\nunweighted 20-kb window means by block class (autosomes)')
print(wm[wm.config == 'P99_noRM'].round(5).to_string(index=False))
