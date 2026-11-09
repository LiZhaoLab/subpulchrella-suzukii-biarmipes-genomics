#!/usr/bin/env python3
"""Inversion analysis, approach 3: summary of the frequency clusters of scripts/110 across the
four D. suzukii pools, and their overlap with the sweep peaks of analysis F.

- Replication: a cluster of one pool is replicated when another pool has a cluster overlapping it with a
  band frequency within 0.05.
- Flags: inside the gene-poor, TE-rich centre of D. suzukii chromosome 3 (27-62 Mb; QC report section 12);
  band MAF >= 0.40 (frequencies near 0.5 are typical of collapsed duplicated sequence).
- Sweep peaks (results/F_sweep_peaks.tsv) whose maximum lies inside a cluster of the same pool.
Output: results/INV_cluster_summary.tsv and .log
"""
import numpy as np
import pandas as pd

RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
DSUZ = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
ARM = {'NC_092080.1': '2L', 'NC_092081.1': '2R', 'NC_092082.1': '3', 'NC_092084.1': 'X'}
cl = pd.concat([pd.read_csv(f'{RES}/INV_freq_clusters_{p}.tsv', sep='\t').assign(pool=p) for p in DSUZ],
               ignore_index=True)
pk = pd.read_csv(f'{RES}/F_sweep_peaks.tsv', sep='\t')
rep, flags, nsweep = [], [], []
for r in cl.itertuples():
    o = cl[(cl.pool != r.pool) & (cl.chrom == r.chrom) & (cl.start <= r.end) & (cl.end >= r.start)
           & ((cl.band_maf - r.band_maf).abs() <= 0.05 + 1e-9)]
    rep.append(','.join(sorted(set(o.pool))))
    f = []
    if r.chrom == 'NC_092082.1' and r.start <= 62_000_000 and r.end >= 27_000_000:
        f.append('chr3_gene_poor_centre')
    if r.band_maf >= 0.40:
        f.append('MAF_near_0.5')
    flags.append(','.join(f))
    s = pk[(pk.pool == r.pool) & (pk.chrom == r.chrom) & (pk.max_pos >= r.start) & (pk.max_pos <= r.end)]
    nsweep.append(len(s))
cl['replicated_in'] = rep
cl['n_replicating_pools'] = [len(x.split(',')) if x else 0 for x in rep]
cl['flag'] = flags
cl['sweep_peaks_inside'] = nsweep
cl['arm'] = cl.chrom.map(ARM)
cl.to_csv(f'{RES}/INV_cluster_summary.tsv', sep='\t', index=False)
log = [f'D. suzukii frequency clusters: {len(cl)} in 4 pools; replicated in >= 1 other pool: '
       f'{(cl.n_replicating_pools >= 1).sum()}; in all 3 other pools: {(cl.n_replicating_pools == 3).sum()}; '
       f'flagged (gene-poor centre or MAF near 0.5): {(cl.flag != "").sum()}']
log.append('\nunflagged clusters replicated in >= 2 other pools (merged regions):')
good = cl[(cl.flag == '') & (cl.n_replicating_pools >= 2)].sort_values(['chrom', 'start'])
cur = None
for r in good.itertuples():
    if cur and r.chrom == cur['chrom'] and r.start <= cur['end']:
        cur['end'] = max(cur['end'], r.end)
        cur['maf'].append(r.band_maf)
        cur['pools'].add(r.pool)
        cur['sweeps'] += r.sweep_peaks_inside
    else:
        if cur:
            log.append(f"  {ARM[cur['chrom']]} {cur['chrom']}:{cur['start']}-{cur['end']} band MAF "
                       f"{min(cur['maf']):.3f}-{max(cur['maf']):.3f} pools {sorted(cur['pools'])} sweep peaks inside {cur['sweeps']}")
        cur = dict(chrom=r.chrom, start=r.start, end=r.end, maf=[r.band_maf], pools={r.pool}, sweeps=r.sweep_peaks_inside)
if cur:
    log.append(f"  {ARM[cur['chrom']]} {cur['chrom']}:{cur['start']}-{cur['end']} band MAF "
               f"{min(cur['maf']):.3f}-{max(cur['maf']):.3f} pools {sorted(cur['pools'])} sweep peaks inside {cur['sweeps']}")
tot = pk[pk.pool.isin(DSUZ)].groupby('pool').size()
ins = cl.groupby('pool').sweep_peaks_inside.sum()
cov = cl.assign(L=cl.end - cl.start + 1).groupby('pool').L.sum() / 1e6
log.append('\nsweep peaks (F) inside a frequency cluster of the same pool: ' +
           '; '.join(f'{p} {int(ins.get(p, 0))}/{int(tot.get(p, 0))} peaks in {cov.get(p, 0):.1f} Mb of clusters' for p in DSUZ))
open(f'{RES}/INV_cluster_summary.log', 'w').write('\n'.join(log) + '\n')
print('\n'.join(log))
