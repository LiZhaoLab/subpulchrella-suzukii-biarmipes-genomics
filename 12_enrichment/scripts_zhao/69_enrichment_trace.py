#!/usr/bin/env python3
"""Level 1 trace and power check: Zhao & Begun (2017) overlap and FlyAtlas2 tissue enrichment.

Trace: recomputes the manuscript numbers from Airi's cmktest tables
(reanalysis_longest_fix_20260817, *.filtered.csv) and her FBgn conversion tables (david_mapping),
following the Methods:
  - MK-positive = P < 0.01 and Alpha > 0.
  - Zhao & Begun: universe = tested transcripts with an FBgn; one transcript per FBgn (smallest P,
    then larger Alpha, then transcript ID); one-sided Fisher's exact test; BH (Benjamini-Hochberg).
  - FlyAtlas2: FLYA2_TAU_Genome_NS.txt, maximum-tissue TPM > 1 (secondary: also Tau >= 0.8);
    target and background are FBgn sets (duplicates collapsed) restricted to FlyAtlas2 genes;
    one-sided hypergeometric test per tissue (33 categories), OR = ad/bc, BH across tissues.
    Sets: Dsub-focal (Dsub vs Dsuz); Dsuz union (positive in >= 1 population, background eligible
    in >= 1) and common4 (positive in all four, background eligible in all four).
Power check: MK power = informative sites Pn + Ps + Dn + Ds (union: maximum over populations;
common4: minimum; FBgn with several transcripts: maximum). Background genes are split into power
deciles; with the number of MK-positive genes fixed in each decile, the overlap with the reference
set or tissue is a sum of independent hypergeometric variables. Its exact distribution (convolution)
gives the power-matched expectation and one-sided P. No random numbers are used.
Outputs: results/L1_zhao_begun.tsv, results/L1_flyatlas2.tsv; log to stdout.
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

A = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817'
CM = f'{A}/cmktest_results'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
POPS = ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')


def mk_table(path):
    t = pd.read_csv(path)
    t['id'] = t.Name.str.replace(r'^rna-', '', regex=True).str.replace(r'\.best$', '', regex=True)
    assert t.id.is_unique, path
    t['power'] = t.Pn + t.Ps + t.Dn + t.Ds
    t['pos'] = (t.Pvalue < 0.01) & (t.Alpha > 0)
    t['impossible'] = (t.Dn + t.Ds > t.Ncodon * 3) | (t.Pn + t.Ps > t.Ncodon * 3)
    return t.set_index('id')


def fbgn_map(sp):
    m = pd.read_csv(f'{A}/david_mapping/{sp}_to_Dmel_FBgn.tsv', sep='\t')
    assert m.focal_mrna.is_unique and m.Dmel_FBgn.notna().all()
    return dict(zip(m.focal_mrna, m.Dmel_FBgn))


def matched_test(power, pos, member):
    """Exact power-matched expectation and P(overlap >= observed); arrays over background genes."""
    dec = pd.qcut(pd.Series(power).rank(method='first'), 10, labels=False).to_numpy()
    pmf = np.array([1.0])
    for g in range(10):
        i = dec == g
        n, k, m = int(i.sum()), int(pos[i].sum()), int(member[i].sum())
        pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(k, m) + 1), n, m, k))
    obs = int((pos & member).sum())
    exp = float((np.arange(len(pmf)) * pmf).sum())
    return exp, obs / exp if exp else np.nan, float(pmf[obs:].sum())


Q = {'Dsub': mk_table(f'{CM}/Dsub_pairwise/Dsub_vs_Dsuz_mktest.filtered.csv')}
for p in POPS:
    Q[p] = mk_table(f'{CM}/Dsuz_pairwise/{p}/Dsuz_vs_Dsub_mktest.filtered.csv')
FB = {'Dsub': fbgn_map('Dsub'), 'Dsuz': fbgn_map('Dsuz')}
sp = lambda c: 'Dsub' if c == 'Dsub' else 'Dsuz'
for c, t in Q.items():
    print(f'{c}: tested {len(t)}, MK-positive {int(t.pos.sum())}, with FBgn {t.index.isin(FB[sp(c)]).sum()}, '
          f'impossible-count rows {int(t.impossible.sum())}')

# ---------------------------------------------------------------- Zhao & Begun (2017)
REF = {n: set(pd.read_csv(f'{PROJ}/inputs/ZhaoBegun2017_{f}.tsv', sep='\t').FBgn)
       for n, f in (('638', 'S13_DOS_638'), ('66', 'S15_both_66'))}
TEXT_ZB = {('638', 'Dsub'): (601, 93), ('638', 'JP-Kan'): (None, 124), ('638', 'JP-Sap'): (None, 146),
           ('638', 'JP-Sap2'): (None, 133), ('638', 'JP-Tok'): (None, 145), ('66', 'Dsub'): (64, 20)}
rows = []
for c, t in Q.items():
    fb = FB[sp(c)]
    u = t[t.index.isin(fb)].reset_index()
    u['FBgn'] = u.id.map(fb)
    u = u.sort_values(['FBgn', 'Pvalue', 'Alpha', 'id'], ascending=[True, True, False, True])
    ndup = int(u.FBgn.duplicated().sum())
    u = u.drop_duplicates('FBgn')
    U, pos = set(u.FBgn), set(u.FBgn[u.pos])
    for n, ref in REF.items():
        test = ref & U
        k = len(test & pos)
        b, cc = len(test) - k, len(pos) - k
        d = len(U) - len(pos) - b
        exp = len(test) * len(pos) / len(U)
        member = u.FBgn.isin(ref).to_numpy()
        mexp, mfold, mp = matched_test(u.power.to_numpy(), u.pos.to_numpy(), member)
        rows.append(dict(set=n, comparison=c, universe=len(U), dup_fbgn_dropped=ndup, positive=len(pos),
                         testable=len(test), overlap=k, pct=round(100 * k / len(test), 1),
                         expected=round(exp, 1), fold=round(k / exp, 2),
                         P=fisher_exact([[k, b], [cc, d]], alternative='greater')[1],
                         median_power_ref=u.power[member].median(), median_power_other=u.power[~member].median(),
                         matched_expected=round(mexp, 1), matched_fold=round(mfold, 2), matched_P=mp))
zb = pd.DataFrame(rows)
zb['BH_all10'] = multipletests(zb.P, method='fdr_bh')[1]
zb['BH_within_set'] = zb.groupby('set').P.transform(lambda p: multipletests(p, method='fdr_bh')[1])
zb['matched_BH_all10'] = multipletests(zb.matched_P, method='fdr_bh')[1]
zb.to_csv(f'{RES}/L1_zhao_begun.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print('\n## Zhao & Begun overlap (manuscript: 638 set Dsub 601 testable, 93 overlap, 15.5%, 2.66-fold, '
      'BH P 4.30e-18; Dsuz 596-598 testable, 124/146/133/145, 20.8-24.4%, 1.93-2.15-fold, BH P <= 2.72e-13; '
      '66 set: 64 testable, Dsub 20, Dsuz 24-26, 31.3% and 37.5-40.6%, 3.11-5.37-fold, BH P <= 1.54e-7)')
print(zb.drop(columns=['median_power_ref', 'median_power_other']).to_string(index=False, float_format=lambda x: f'{x:.3g}'))
print(zb[['set', 'comparison', 'median_power_ref', 'median_power_other']].to_string(index=False))

# ---------------------------------------------------------------- FlyAtlas2
fa = pd.read_csv(f'{PROJ}/FLYA2_TAU_Genome_NS.txt', sep='\t')
assert fa.Gene_ID.is_unique and fa.Sex_Tissue.nunique() == 33
TISSUES = sorted(fa.Sex_Tissue.unique())
DATA = {'primary_TPM>1': fa[fa.TPM > 1].set_index('Gene_ID').Sex_Tissue,
        'secondary_Tau>=0.8': fa[(fa.TPM > 1) & (fa.Tau >= 0.8)].set_index('Gene_ID').Sex_Tissue}


def fbgn_level(ids_power, ids_pos, fb):
    """ids_power: Series transcript -> power over the eligible background; ids_pos: positive transcripts."""
    s = ids_power[ids_power.index.isin(fb)]
    g = pd.DataFrame({'FBgn': s.index.map(fb), 'power': s.values, 'pos': s.index.isin(ids_pos)})
    return g.groupby('FBgn').agg(power=('power', 'max'), pos=('pos', 'any'))


dsuz = pd.concat({p: Q[p][['power', 'pos']] for p in POPS}, axis=1)
elig_any, elig_all = dsuz.index, dsuz.dropna().index
SETS = {'Dsub_focal': fbgn_level(Q['Dsub'].power, set(Q['Dsub'].index[Q['Dsub'].pos]), FB['Dsub']),
        'Dsuz_union': fbgn_level(dsuz.xs('power', axis=1, level=1).max(axis=1),
                                 set(elig_any[dsuz.xs('pos', axis=1, level=1).fillna(False).any(axis=1)]), FB['Dsuz']),
        'Dsuz_common4': fbgn_level(dsuz.loc[elig_all].xs('power', axis=1, level=1).min(axis=1),
                                   set(elig_all[dsuz.loc[elig_all].xs('pos', axis=1, level=1).astype(bool).all(axis=1)]),
                                   FB['Dsuz'])}
rows = []
for dname, tis in DATA.items():
    for sname, g in SETS.items():
        g = g[g.index.isin(tis.index)]
        tissue = tis.reindex(g.index).to_numpy()
        pos, M, N = g.pos.to_numpy(), len(g), int(g.pos.sum())
        for T in TISSUES:
            member = tissue == T
            a, n = int((pos & member).sum()), int(member.sum())
            b, c = N - a, n - a
            d = M - N - c
            mexp, mfold, mp = matched_test(g.power.to_numpy(), pos, member)
            rows.append(dict(dataset=dname, set=sname, tissue=T, background=M, target=N, a=a, b=b, c=c, d=d,
                             OR=a * d / (b * c) if b * c else np.nan, P=hypergeom.sf(a - 1, M, n, N),
                             median_power_tissue=np.median(g.power[member]) if n else np.nan,
                             median_power_other=np.median(g.power[~member]),
                             matched_expected=mexp, matched_fold=mfold, matched_P=mp))
fl = pd.DataFrame(rows)
fl['FDR'] = fl.groupby(['dataset', 'set']).P.transform(lambda p: multipletests(p, method='fdr_bh')[1])
fl['matched_FDR'] = fl.groupby(['dataset', 'set']).matched_P.transform(lambda p: multipletests(p, method='fdr_bh')[1])
fl.to_csv(f'{RES}/L1_flyatlas2.tsv', sep='\t', index=False)
print('\n## FlyAtlas2 (manuscript: Dsub Male_Testis 229, OR 2.46, P 1.18e-23, FDR 3.91e-22; union Male_Midgut 85, '
      'OR 1.77, P 1.42e-5, FDR 4.69e-4; common4 Male_Crop 29, OR 1.82, P 0.00357, FDR 0.0729; '
      'common4 Male_Testis 177, OR 1.27, P 0.00442, FDR 0.0729)')
cols = ['set', 'tissue', 'background', 'target', 'a', 'OR', 'P', 'FDR', 'median_power_tissue', 'median_power_other',
        'matched_expected', 'matched_fold', 'matched_P', 'matched_FDR']
for dname in DATA:
    x = fl[(fl.dataset == dname) & ((fl.FDR < 0.1) | (fl.matched_FDR < 0.1))].sort_values(['set', 'P'])
    print(f'\n# {dname}: tissues with FDR < 0.1 (either test)')
    print(x[cols].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
