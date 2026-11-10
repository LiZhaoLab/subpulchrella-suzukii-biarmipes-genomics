#!/usr/bin/env python3
"""Does a category's lineage difference in MK-positive genes come from divergence or from polymorphism?

Why: in script 130, chemosensory genes and GO 'extracellular region' genes were more
frequent among D. suzukii-lineage-only MK-positive genes, and testis-maximal genes among
D. subpulchrella-lineage-only genes. Two explanations make different predictions:
  (1) lineage-specific adaptive substitution: the category's Dn/Ds is higher on one lineage than in
      the genome, with polymorphism (Pn/Ps) changed no more than in the genome;
  (2) a polymorphism effect: D. subpulchrella (smaller Ne) carries more segregating nonsynonymous
      variants in weakly constrained genes, which lowers its MK power; the category's Pn/Ps ratio
      (D. subpulchrella / D. suzukii) is then higher than in the genome, with Dn/Ds unchanged.
Statistics per category, summed over genes tested on both lineages (results/lineage_genes.tsv,
DAF >= 15%, D. suzukii counts = median over pools):
  RD = (Dn/Ds on the D. subpulchrella lineage) / (Dn/Ds on the D. suzukii lineage);
  RP = (Pn/Ps in D. subpulchrella) / (Pn/Ps in D. suzukii);
  alpha on each lineage (1 - Ds Pn / (Dn Ps), ratio of sums).
Each is divided by the genome value (all genes in the table); 95% intervals from 2,000 bootstrap
replicates over the category's genes (genome value held fixed); seed 20261104. Two-sided bootstrap P
for ratio-to-genome = 1 (2 x the smaller tail fraction). Categories: chemosensory and immunity (script
130 definitions), GO:0005576 extracellular region and GO:0016020 membrane (direct annotations, as in
script 70), FlyAtlas2 testis-maximal genes.
Added after the first run: DnR = Dn_sub/Dn_suz and DsR = Ds_sub/Ds_suz (RD = DnR / DsR), so that a change
in RD can be traced to nonsynonymous or synonymous fixed differences (the residual D. subpulchrella-lineage
Ds excess on chromosome 3 lowers RD by raising Ds); and a matched random-set test: 10,000 random gene sets
of the same size drawn without replacement within strata of chromosome arm x decile of informative counts
(Pn+Ps+Dn+Ds on both lineages), matching the category's stratum counts; two-sided P for RD, RP, DnR,
DsR and delta alpha; per-arm values (descriptive).
Output: results/lineage_category_decomposition.tsv and .log
"""
import gzip

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
GO = f'{PROJ}/inputs/GO_20261022'
NB, NRAND, SEED = 2000, 10000, 20261104


def stats(a):
    pn_s, ps_s, dn_s, ds_s, pn_z, ps_z, dn_z, ds_z = a.sum(0)
    rd = (dn_s / ds_s) / (dn_z / ds_z)
    rp = (pn_s / ps_s) / (pn_z / ps_z)
    al_s = 1 - ds_s * pn_s / (dn_s * ps_s)
    al_z = 1 - ds_z * pn_z / (dn_z * ps_z)
    return np.array([rd, rp, al_s, al_z, dn_s / ds_s, dn_z / ds_z, pn_s / ps_s, pn_z / ps_z, dn_s / dn_z,
                     ds_s / ds_z])


def main():
    t = pd.read_csv(f'{RES}/lineage_genes.tsv', sep='\t')
    go = {'extracellular region': 'GO:0005576', 'membrane': 'GO:0016020'}
    ann = {k: set() for k in go}
    with gzip.open(f'{GO}/fb.gaf.gz', 'rt') as f:
        for line in f:
            if line.startswith('!'):
                continue
            c = line.split('\t')
            if c[0] != 'FB' or 'NOT' in c[3].split('|'):
                continue
            for k, term in go.items():
                if c[4] == term:
                    ann[k].add(c[1])
    cats = {'chemosensory': t.eco.fillna('').str.contains('chemosensory'),
            'immunity': t.eco.fillna('').str.contains('immunity'),
            'extracellular region (GO)': t.FBgn.isin(ann['extracellular region']),
            'membrane (GO)': t.FBgn.isin(ann['membrane']),
            'testis-maximal (FlyAtlas2)': t.FlyAtlas2_max.eq('Male_Testis')}
    cols = ['Pn_sub', 'Ps_sub', 'Dn_sub', 'Ds_sub', 'Pn_suz', 'Ps_suz', 'Dn_suz', 'Ds_suz']
    A = t[cols].to_numpy(float)
    g = stats(A)
    rng = np.random.default_rng(SEED)
    names = ['RD', 'RP', 'alpha_sub', 'alpha_suz', 'DnDs_sub', 'DnDs_suz', 'PnPs_sub', 'PnPs_suz', 'DnR', 'DsR']
    info = A.sum(1)
    dec = pd.qcut(pd.Series(info).rank(method='first'), 10, labels=False).to_numpy()
    strat = t.arm.fillna('NA').astype(str).to_numpy() + '_' + dec.astype(str)
    by_st = {k: np.flatnonzero(strat == k) for k in np.unique(strat)}
    rows = [dict(category='genome', genes=len(t), **dict(zip(names, g)))]
    for cat, m in cats.items():
        idx = np.flatnonzero(m.to_numpy())
        s = stats(A[idx])
        bs = np.array([stats(A[rng.choice(idx, len(idx))]) for _ in range(NB)])
        row = dict(category=cat, genes=len(idx), **dict(zip(names, s)))
        for j, n in enumerate(names[:2]):
            v = bs[:, j] / g[j]
            v = v[np.isfinite(v)]
            row[f'{n}_rel'] = s[j] / g[j]
            row[f'{n}_rel_lo'], row[f'{n}_rel_hi'] = np.percentile(v, [2.5, 97.5])
            row[f'{n}_rel_P'] = min(1.0, 2 * min((v <= 1).mean(), (v >= 1).mean()))
        for j, n in ((2, 'alpha_sub'), (3, 'alpha_suz')):
            v = bs[:, j][np.isfinite(bs[:, j])]
            row[f'{n}_lo'], row[f'{n}_hi'] = np.percentile(v, [2.5, 97.5])
        d = bs[:, 3] - bs[:, 2]
        d = d[np.isfinite(d)]
        row['dalpha'] = s[3] - s[2]
        row['dalpha_lo'], row['dalpha_hi'] = np.percentile(d, [2.5, 97.5])
        row['DnR_rel'], row['DsR_rel'] = s[8] / g[8], s[9] / g[9]
        need = pd.Series(strat[idx]).value_counts()
        rnd = np.array([stats(A[np.concatenate([rng.choice(by_st[k], n, replace=False) for k, n in need.items()])])
                        for _ in range(NRAND)])
        obs = dict(RD=s[0], RP=s[1], DnR=s[8], DsR=s[9], dalpha=s[3] - s[2])
        rv = dict(RD=rnd[:, 0], RP=rnd[:, 1], DnR=rnd[:, 8], DsR=rnd[:, 9], dalpha=rnd[:, 3] - rnd[:, 2])
        for k, o in obs.items():
            v = rv[k][np.isfinite(rv[k])]
            row[f'{k}_matched_median'] = float(np.median(v))
            row[f'{k}_matched_P'] = min(1.0, 2 * min(((v <= o).sum() + 1) / (len(v) + 1), ((v >= o).sum() + 1) / (len(v) + 1)))
        arms = {}
        for a in sorted(t.arm.dropna().unique()):
            j = np.intersect1d(idx, np.flatnonzero(t.arm.eq(a).to_numpy()))
            if len(j) >= 10:
                sa, ga = stats(A[j]), stats(A[np.flatnonzero(t.arm.eq(a).to_numpy())])
                arms[a] = f'{len(j)}:RD {sa[0] / ga[0]:.2f},RP {sa[1] / ga[1]:.2f}'
        row['per_arm_rel_to_arm_genome'] = '; '.join(f'{a} {v}' for a, v in arms.items())
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(f'{RES}/lineage_category_decomposition.tsv', sep='\t', index=False)
    with open(f'{RES}/lineage_category_decomposition.log', 'w') as log:
        s = (f'Seed {SEED}, {NB} bootstrap replicates over genes. RD = Dn/Ds Dsub lineage / Dn/Ds Dsuz lineage; '
             f'RP = Pn/Ps Dsub / Pn/Ps Dsuz; *_rel = category / genome.\n')
        s += out.to_string(index=False, float_format=lambda v: f'{v:.3f}')
        print(s)
        log.write(s + '\n')


if __name__ == '__main__':
    main()
