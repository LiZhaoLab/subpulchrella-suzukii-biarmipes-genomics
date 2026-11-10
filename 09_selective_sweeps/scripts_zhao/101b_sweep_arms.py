#!/usr/bin/env python3
"""F follow-up: sweep-peak density per chromosome arm in both species, with D. suzukii
chromosome 3 split into 3L and 3R through orthologs, and the strongest peaks with their genes.

Arm of a D. suzukii 1-Mb window = the D. subpulchrella arm (scripts/depth_common.ARM_DSUB) of most of
the orthologous genes in that window (Step-0 gene positions, JP-Kan table); windows without
orthologous genes are 'unassigned' (mainly the gene-poor centre of chromosome 3). Peaks are assigned
by their maximum position; scanned length per arm = 2 kb x grid points of the windows.
Inputs: results/F_sweep_peaks.tsv, the CLR files. Outputs: results/F_sweep_arms.tsv,
results/F_sweep_arms.log
"""
import glob
import os

import numpy as np
import pandas as pd

import depth_common as dc

QC = f'{dc.QC}/sweeps'
DOT = {'NC_050608.1', 'NC_092083.1'}
d, zs = dc.gene_tables()
z = zs['JP-Kan']
g = z['cds_gene']
n = len(z['genes'])
pos = np.bincount(g, weights=z['cds_pos'].astype(float), minlength=n) / np.maximum(np.bincount(g, minlength=n), 1)
kan = pd.DataFrame({'chrom': z['chroms'][z['gene_chrom']], 'pos': pos}).iloc[d['JP-Kan'].values].reset_index(drop=True)
kan['arm'] = d.arm.values
kan['w'] = (kan.pos // 1_000_000).astype(int)
warm = kan.groupby(['chrom', 'w']).arm.agg(lambda a: a.mode().iloc[0]).to_dict()
pk = pd.read_csv(f'{dc.RES}/F_sweep_peaks.tsv', sep='\t')
rows = []
for pool in dc.POOLS:
    dsub = pool.startswith('Dsub')
    length = {}
    for f in glob.glob(f'{QC}/{pool}/clr/*.txt'):
        chrom = os.path.basename(f).rsplit('_', 1)[0]
        if chrom in DOT:
            continue
        c = pd.read_csv(f, sep='\t')
        if dsub:
            arms = np.full(len(c), dc.ARM_DSUB[chrom])
        else:
            arms = [warm.get((chrom, int(x // 1_000_000)), 'unassigned') for x in c.location]
        for a, k in pd.Series(arms).value_counts().items():
            length[a] = length.get(a, 0) + k * 2000 / 1e6
    p = pk[pk.pool == pool]
    if dsub:
        parm = p.chrom.map(dc.ARM_DSUB)
    else:
        parm = [warm.get((c, int(x // 1_000_000)), 'unassigned') for c, x in zip(p.chrom, p.max_pos)]
    cnt = pd.Series(parm).value_counts()
    for a, mb in length.items():
        rows.append(dict(pool=pool, arm=a, scanned_mb=mb, peaks=int(cnt.get(a, 0)), peaks_per_mb=cnt.get(a, 0) / mb))
r = pd.DataFrame(rows)
r.to_csv(f'{dc.RES}/F_sweep_arms.tsv', sep='\t', index=False)
out = ['Peaks per Mb by arm (D. suzukii chromosome 3 split through orthologs); within a pool the threshold is fixed per chromosome class']
piv = r.pivot(index='arm', columns='pool', values='peaks_per_mb').round(2)
out.append(piv.to_string())
# relative to the pool's mean over 2L, 2R, 3L, 3R
rel = piv.loc[['2L', '2R', '3L', '3R']].mean()
out.append('\nRelative to the pool mean over 2L, 2R, 3L, 3R:')
out.append((piv / rel).round(2).to_string())
out.append('\nStrongest 8 peaks per pool (CLR, position, genes with D. melanogaster symbols):')
for pool in dc.POOLS:
    p = pk[pk.pool == pool].sort_values('max_clr', ascending=False).head(8)
    out.append(f'-- {pool}')
    for x in p.itertuples():
        out.append(f'   {x.chrom}:{x.start}-{x.end} CLR {x.max_clr:.1f} {x.cls} genes: {x.symbols if isinstance(x.symbols, str) else "-"}')
open(f'{dc.RES}/F_sweep_arms.log', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
