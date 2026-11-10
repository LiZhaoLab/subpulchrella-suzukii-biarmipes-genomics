#!/usr/bin/env python3
"""Summarize lineage-specific (polarized) MK tests from 60_polarized_mk.py.

Per focal pool (Dsub_L5_7; JP-Kan, JP-Sap, JP-Sap2, JP-Tok), polymorphisms are kept at
derived-allele frequency (DAF) >= threshold (5% or 15%; unpolarized polymorphisms use the
minor-allele frequency instead). Per gene: Fisher's exact test on
[[Pn, Ps], [Dn_lineage, Ds_lineage]]; alpha = 1 - Ds*Pn/(Dn*Ps). Genes are kept if no row or
column of the table is empty. MK-positive: P < 0.01 and alpha > 0; also BH FDR < 0.05.
Genome-wide alpha: (i) ratio of sums with 1,000 gene bootstraps; (ii) asymptotic MK
(Messer and Petrov 2013): alpha(x) from DAF bins 0.05-0.90 (width 0.05), fitted with
a + b*exp(-c*x), alpha_asym = fit at x = 1; 500 gene bootstraps. Seed 20261021.
Outputs: results/L3_polarized_mk_*.tsv
PMK_VARIANT=excl: read results/polarized_mk_excl (60_polarized_mk.py with
EXCLUDE_SHARED=1) and write the L3 tables to results/polarized_mk_excl/summary/.
"""
import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import fisher_exact, hypergeom

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
VARIANT = os.environ.get('PMK_VARIANT', '')
R = f'{P}/results/polarized_mk' + (f'_{VARIANT}' if VARIANT else '')
OUTD = f'{R}/summary' if VARIANT else f'{P}/results'
os.makedirs(OUTD, exist_ok=True)
MK = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817'
POPS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
rng = np.random.default_rng(20261021)


def bh(p):
    p = np.asarray(p)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def gene_table(focal, thr):
    g = pd.read_csv(f'{R}/{focal}.genes.tsv', sep='\t').set_index('gene')
    p = pd.read_csv(f'{R}/{focal}.polymorphisms.tsv.gz', sep='\t')
    f = p.daf.where(p.daf.notna(), p.maf)
    p = p[f >= thr]
    cnt = p.groupby(['gene', 'type']).size().unstack(fill_value=0).reindex(g.index, fill_value=0)
    t = pd.DataFrame({'Pn': cnt.get('n', 0), 'Ps': cnt.get('s', 0),
                      'Dn': g.Dn_focal.round(), 'Ds': g.Ds_focal.round(), 'callable': g.callable})
    t = t[((t.Pn + t.Dn) > 0) & ((t.Ps + t.Ds) > 0) & ((t.Pn + t.Ps) > 0) & ((t.Dn + t.Ds) > 0)]
    t['P'] = [fisher_exact([[a, b], [c, d]])[1] for a, b, c, d in zip(t.Pn, t.Ps, t.Dn, t.Ds)]
    with np.errstate(divide='ignore', invalid='ignore'):
        t['alpha'] = 1 - (t.Ds * t.Pn) / (t.Dn * t.Ps)
    t['q'] = bh(t.P.values)
    return t


def alpha_sum(t, idx=None):
    x = t if idx is None else t.iloc[idx]
    return 1 - (x.Ds.sum() * x.Pn.sum()) / (x.Dn.sum() * x.Ps.sum())


def asymptotic(focal, genes, idx=None):
    g = pd.read_csv(f'{R}/{focal}.genes.tsv', sep='\t').set_index('gene').loc[genes]
    p = pd.read_csv(f'{R}/{focal}.polymorphisms.tsv.gz', sep='\t')
    p = p[p.daf.notna() & p.gene.isin(genes)]
    if idx is not None:
        pick = pd.Series(genes).iloc[idx]
        w = pick.value_counts()
        g = g.loc[pick]
        p = p.merge(w.rename('w'), left_on='gene', right_index=True)
    else:
        p = p.assign(w=1)
    Dn, Ds = g.Dn_focal.sum(), g.Ds_focal.sum()
    edges = np.arange(0.05, 0.9001, 0.05)
    xs, al = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = p[(p.daf >= lo) & (p.daf < hi)]
        pn, ps = b.w[b.type == 'n'].sum(), b.w[b.type == 's'].sum()
        if ps > 0:
            xs.append((lo + hi) / 2)
            al.append(1 - (Ds / Dn) * (pn / ps))
    xs, al = np.array(xs), np.array(al)
    try:
        par, _ = curve_fit(lambda x, a, b, c: a + b * np.exp(-c * x), xs, al, p0=(al[-1], -0.5, 5),
                           maxfev=20000)
        return par[0] + par[1] * np.exp(-par[2]), xs, al
    except RuntimeError:
        k = np.polyfit(xs, al, 1)
        return k[0] + k[1], xs, al


def airi_pos(path):
    a = pd.read_csv(path)
    a['g'] = a.Name.str.replace(r'^rna-', '', regex=True).str.replace(r'\.best$', '', regex=True)
    return set(a.g[(a.Pvalue < 0.01) & (a.Alpha > 0)]), set(a.g)


rows, tabs, genomewide, curves = [], {}, [], []
for focal in ['Dsub_L5_7'] + POPS:
    for thr in (0.05, 0.15):
        t = gene_table(focal, thr)
        tabs[(focal, thr)] = t
        pos = set(t.index[(t.P < 0.01) & (t.alpha > 0)])
        fdr = set(t.index[(t.q < 0.05) & (t.alpha > 0)])
        if focal == 'Dsub_L5_7':
            ap, at = airi_pos(f'{MK}/cmktest_results/Dsub_pairwise/Dsub_vs_Dsuz_mktest.filtered.csv')
        else:
            ap, at = airi_pos(f'{MK}/cmktest_results/Dsuz_pairwise/{focal}/Dsuz_vs_Dsub_mktest.filtered.csv')
        both = set(t.index) & at
        rows.append(dict(focal=focal, DAF_threshold=thr, genes_tested=len(t), MK_positive=len(pos),
                         FDR05_positive=len(fdr), airi_unpolarized_positive=len(ap),
                         shared_with_airi=len(pos & ap),
                         airi_positive_tested_here=len(ap & both),
                         pct_airi_recovered=round(100 * len(pos & ap) / max(1, len(ap & both)), 1)))
    t = tabs[(focal, 0.15)]
    a0 = alpha_sum(t)
    bs = [alpha_sum(t, rng.integers(0, len(t), len(t))) for _ in range(1000)]
    aa, xs, al = asymptotic(focal, list(t.index))
    ba = [asymptotic(focal, list(t.index), rng.integers(0, len(t), len(t)))[0] for _ in range(500)]
    genomewide.append(dict(focal=focal, genes=len(t), Pn=int(t.Pn.sum()), Ps=int(t.Ps.sum()),
                           Dn=int(t.Dn.sum()), Ds=int(t.Ds.sum()),
                           alpha_DAF15=round(a0, 3), alpha_lo=round(np.percentile(bs, 2.5), 3),
                           alpha_hi=round(np.percentile(bs, 97.5), 3),
                           alpha_asymptotic=round(aa, 3), asym_lo=round(np.percentile(ba, 2.5), 3),
                           asym_hi=round(np.percentile(ba, 97.5), 3)))
    curves += [dict(focal=focal, daf=x, alpha=y) for x, y in zip(xs, al)]

res = pd.DataFrame(rows)
res.to_csv(f'{OUTD}/L3_polarized_mk_counts.tsv', sep='\t', index=False)
gw = pd.DataFrame(genomewide)
gw.to_csv(f'{OUTD}/L3_polarized_mk_alpha.tsv', sep='\t', index=False)
pd.DataFrame(curves).to_csv(f'{OUTD}/L3_polarized_mk_alpha_by_daf.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print('== lineage-specific MK-positive genes (P < 0.01, alpha > 0)')
print(res.to_string(index=False))
print('\n== genome-wide alpha (DAF >= 15%; 95% CI from gene bootstraps)')
print(gw.to_string(index=False))

# genes positive on both lineages; Dsuz lineage = union over populations
o = pd.read_csv(f'{MK}/orth_longest/Dsuz_Dsub.direct.RBH.longest.mRNA.tsv', sep='\t', header=None)
s2z = dict(zip(o[0].str.replace('rna-', ''), o[1].str.replace('rna-', '')))
out = []
for thr in (0.05, 0.15):
    ts = tabs[('Dsub_L5_7', thr)]
    tz = {p: tabs[(p, thr)] for p in POPS}
    zpos = set.union(*(set(t.index[(t.P < 0.01) & (t.alpha > 0)]) for t in tz.values()))
    ztested = set.union(*(set(t.index) for t in tz.values()))
    univ = [g for g in ts.index if s2z.get(g) in ztested]
    sp = {g for g in univ if ts.loc[g, 'P'] < 0.01 and ts.loc[g, 'alpha'] > 0}
    zp = {g for g in univ if s2z[g] in zpos}
    k = len(sp & zp)
    exp = len(sp) * len(zp) / len(univ)
    p_hyp = hypergeom.sf(k - 1, len(univ), len(sp), len(zp))
    # power-matched expectation: permute Dsub-lineage labels within deciles of total counts
    tot = ts.loc[univ, ['Pn', 'Ps', 'Dn', 'Ds']].sum(axis=1)
    dec = pd.qcut(tot.rank(method='first'), 10, labels=False)
    lab = pd.Series([g in sp for g in univ], index=univ)
    inz = pd.Series([g in zp for g in univ], index=univ)
    perm = []
    for _ in range(2000):
        sh = lab.groupby(dec).transform(lambda v: rng.permutation(v.values))
        perm.append(int((sh & inz).sum()))
    out.append(dict(DAF_threshold=thr, universe=len(univ), Dsub_lineage_positive=len(sp),
                    Dsuz_lineage_positive_union=len(zp), both=k, expected_random=round(exp, 1),
                    fold=round(k / exp, 2), P_hypergeometric=p_hyp,
                    expected_power_matched=round(np.mean(perm), 1),
                    fold_power_matched=round(k / np.mean(perm), 2),
                    P_power_matched=(1 + sum(x >= k for x in perm)) / (1 + len(perm))))
ov = pd.DataFrame(out)
ov.to_csv(f'{OUTD}/L3_polarized_mk_both_lineages.tsv', sep='\t', index=False)
print('\n== genes MK-positive on both lineages (separate fixed differences on each lineage)')
print(ov.to_string(index=False))
for (focal, thr), t in tabs.items():
    t.to_csv(f'{R}/{focal}.mk_tests.DAF{int(thr * 100)}.tsv', sep='\t')
