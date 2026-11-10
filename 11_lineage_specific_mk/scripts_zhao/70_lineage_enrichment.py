#!/usr/bin/env python3
"""Lineage-specific functional enrichment (report section 7.3, analysis 1).

Question: the draft contrasts the D. subpulchrella-focal set (testis; cilium and axoneme GO terms)
with the D. suzukii-focal sets (midgut; membrane and olfactory terms). Those sets come from
unpolarized MK tests that share the same fixed differences. Does the contrast hold when MK-positive
genes are defined per lineage?

Gene sets (MK-positive = P < 0.01 and alpha > 0; FBgn via Airi's david_mapping, duplicates collapsed):
  unpolarized (Airi's cmktest tables): Dsub-focal; Dsuz union; Dsuz common4.
  lineage-specific (64_summarize_polarized_mk.py tables; DAF >= 15% primary, >= 5% sensitivity):
  Dsub lineage; Dsuz lineage union (>= 1 of four pools); Dsuz lineage common4 (all four pools).
  Background = genes tested in that analysis (union: in >= 1 pool; common4: in all four).
Tests per set, unadjusted (one-sided hypergeometric) and power-matched (exact; deciles of
Pn + Ps + Dn + Ds, as in 69_enrichment_trace.py):
  FlyAtlas2 maximum tissue (TPM > 1; 33 tissues; BH across tissues); Zhao & Begun 638 and 66;
  GO direct annotations (FlyBase GAF from geneontology.org; NOT qualifiers removed; all evidence
  codes; background restricted to genes with >= 1 annotation in the aspect; terms with >= 2 target
  genes; BH within aspect).
Direct contrast: genes tested on both lineages (D. subpulchrella transcript and its D. suzukii RBH
ortholog), FBgn from the D. subpulchrella mapping. Among genes MK-positive on exactly one lineage,
table [in category, not] x [Dsub-lineage only, Dsuz-lineage only]; OR > 1 = category more frequent
among Dsub-lineage-only genes (conditional logistic regression for matched pairs reduces to this
table). Two-sided Fisher's exact test; Mantel-Haenszel OR and CMH test stratified by deciles of
summed power on the two lineages. Dsuz lineage = union of pools, and each pool separately.
Outputs: results/L3_lineage_{flyatlas2,zhao_begun,go,contrast}.tsv; log to stdout.
"""
import os
import gzip

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

A = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817'
CM = f'{A}/cmktest_results'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
VARIANT = os.environ.get('PMK_VARIANT', '')   # 'excl' = 60_polarized_mk.py EXCLUDE_SHARED=1
PMK = f'{RES}/polarized_mk' + (f'_{VARIANT}' if VARIANT else '')
VRES = f'{PMK}/summary' if VARIANT else RES
GO = f'{PROJ}/inputs/GO_20261022'
POPS = ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')
KEY_TISSUES = ('Male_Testis', 'Male_Midgut', 'Male_Crop', 'Female_Hindgut')
KEY_GO = ('axoneme', 'sperm flagellum', 'cilium movement', 'cilium movement involved in cell motility',
          'inner dynein arm', 'axonemal dynein complex', 'transmembrane transport', 'plasma membrane',
          'membrane', 'detection of chemical stimulus involved in sensory perception of smell',
          'olfactory receptor activity', 'ionotropic olfactory receptor activity', 'midgut development',
          'cell adhesion', 'G protein-coupled receptor signaling pathway', 'proteolysis')
pd.set_option('display.width', 250)
fmt = lambda v: f'{v:.3g}'


def deciles(power):
    return pd.qcut(pd.Series(power).rank(method='first'), 10, labels=False).to_numpy()


def matched_test(dec, pos, member):
    """Exact power-matched expectation and P(overlap >= observed); see 69_enrichment_trace.py."""
    pmf = np.array([1.0])
    for g in range(10):
        i = dec == g
        n, k, m = int(i.sum()), int(pos[i].sum()), int(member[i].sum())
        pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(k, m) + 1), n, m, k))
    obs = int((pos & member).sum())
    exp = float((np.arange(len(pmf)) * pmf).sum())
    return exp, obs / exp if exp else np.nan, float(pmf[obs:].sum())


def enrich(g, dec, member):
    pos = g.pos.to_numpy()
    M, N, n, a = len(g), int(pos.sum()), int(member.sum()), int((pos & member).sum())
    exp = N * n / M
    mexp, mfold, mp = matched_test(dec, pos, member)
    return dict(background=M, target=N, in_category=n, overlap=a, expected=exp,
                fold=a / exp if exp else np.nan, P=hypergeom.sf(a - 1, M, n, N),
                matched_expected=mexp, matched_fold=mfold, matched_P=mp)


def fbgn_map(sp):
    m = pd.read_csv(f'{A}/david_mapping/{sp}_to_Dmel_FBgn.tsv', sep='\t')
    assert m.focal_mrna.is_unique and m.Dmel_FBgn.notna().all()
    return dict(zip(m.focal_mrna, m.Dmel_FBgn))


def fb_level(power, pos_ids, fb):
    s = power[power.index.isin(fb)]
    g = pd.DataFrame({'FBgn': s.index.map(fb), 'power': s.values, 'pos': s.index.isin(pos_ids)})
    return g.groupby('FBgn').agg(power=('power', 'max'), pos=('pos', 'any'))


def cmk(path):
    t = pd.read_csv(path)
    t.index = t.Name.str.replace(r'^rna-', '', regex=True).str.replace(r'\.best$', '', regex=True)
    assert t.index.is_unique, path
    return pd.DataFrame({'power': t.Pn + t.Ps + t.Dn + t.Ds, 'pos': (t.Pvalue < 0.01) & (t.Alpha > 0)})


def pmk(focal, daf):
    t = pd.read_csv(f'{PMK}/{focal}.mk_tests.DAF{daf}.tsv', sep='\t', index_col=0)
    assert t.index.is_unique
    return pd.DataFrame({'power': t.Pn + t.Ps + t.Dn + t.Ds, 'pos': (t.P < 0.01) & (t.alpha > 0)})


def single(t, fb):
    return fb_level(t.power, set(t.index[t.pos]), fb)


def union_common(tabs, fb):
    pw = pd.concat({p: t.power for p, t in tabs.items()}, axis=1)
    ps = pd.concat({p: t.pos for p, t in tabs.items()}, axis=1).fillna(False).astype(bool)
    allidx = pw.dropna().index
    u = fb_level(pw.max(axis=1), set(pw.index[ps.any(axis=1)]), fb)
    c = fb_level(pw.loc[allidx].min(axis=1), set(allidx[ps.loc[allidx].all(axis=1)]), fb)
    return u, c


FB = {'Dsub': fbgn_map('Dsub'), 'Dsuz': fbgn_map('Dsuz')}
SETS = {('unpolarized', 'Dsub_focal'): single(cmk(f'{CM}/Dsub_pairwise/Dsub_vs_Dsuz_mktest.filtered.csv'), FB['Dsub'])}
SETS[('unpolarized', 'Dsuz_union')], SETS[('unpolarized', 'Dsuz_common4')] = union_common(
    {p: cmk(f'{CM}/Dsuz_pairwise/{p}/Dsuz_vs_Dsub_mktest.filtered.csv') for p in POPS}, FB['Dsuz'])
for daf in (15, 5):
    SETS[(f'DAF{daf}', 'Dsub_lineage')] = single(pmk('Dsub_L5_7', daf), FB['Dsub'])
    SETS[(f'DAF{daf}', 'Dsuz_lineage_union')], SETS[(f'DAF{daf}', 'Dsuz_lineage_common4')] = union_common(
        {p: pmk(p, daf) for p in POPS}, FB['Dsuz'])
print('## Gene sets (FBgn level): background, MK-positive')
for k, g in SETS.items():
    print(f'  {k[0]:12s} {k[1]:22s} {len(g):6d} {int(g.pos.sum()):5d}')

# ---------------------------------------------------------------- annotations
fa = pd.read_csv(f'{PROJ}/FLYA2_TAU_Genome_NS.txt', sep='\t')
assert fa.Gene_ID.is_unique and fa.Sex_Tissue.nunique() == 33
TIS = fa[fa.TPM > 1].set_index('Gene_ID').Sex_Tissue
ZB = {n: set(pd.read_csv(f'{PROJ}/inputs/ZhaoBegun2017_{f}.tsv', sep='\t').FBgn)
      for n, f in (('638', 'S13_DOS_638'), ('66', 'S15_both_66'))}

names, space, cur = {}, {}, None
for line in open(f'{GO}/go-basic.obo'):
    line = line.rstrip('\n')
    if line == '[Term]':
        cur = {}
    elif line.startswith('[') or not line:
        if cur and 'id' in cur:
            names[cur['id']], space[cur['id']] = cur.get('name'), cur.get('namespace')
            for alt in cur.get('alt', []):
                names.setdefault(alt, cur.get('name')); space.setdefault(alt, cur.get('namespace'))
        cur = None if line.startswith('[') else cur
    elif cur is not None:
        k, _, v = line.partition(': ')
        if k == 'alt_id':
            cur.setdefault('alt', []).append(v)
        elif k in ('id', 'name', 'namespace'):
            cur[k] = v
ann, nnot, gaf_date = {'P': {}, 'F': {}, 'C': {}}, 0, None
with gzip.open(f'{GO}/fb.gaf.gz', 'rt') as f:
    for line in f:
        if line.startswith('!'):
            if line.startswith('!date-generated') and gaf_date is None:
                gaf_date = line.split(':', 1)[1].strip()
            continue
        c = line.rstrip('\n').split('\t')
        if c[0] != 'FB' or 'taxon:7227' not in c[12]:
            continue
        if 'NOT' in c[3].split('|'):
            nnot += 1
            continue
        ann[c[8]].setdefault(c[4], set()).add(c[1])
missing = {t for a in ann.values() for t in a if t not in names}
print(f'\nGO: GAF date-generated {gaf_date}; terms P/F/C {[len(ann[a]) for a in "PFC"]}; '
      f'NOT annotations removed {nnot}; GAF terms absent from go-basic.obo {len(missing)}')
GENES_ASPECT = {a: set().union(*d.values()) for a, d in ann.items()}
name2id = {}
for a, d in ann.items():
    for t in d:
        name2id.setdefault(names.get(t), (a, t))
for n in KEY_GO:
    if n not in name2id:
        print(f'  key GO term not found by name in direct annotations: {n}')

# ---------------------------------------------------------------- per-set enrichment
fl_rows, zb_rows, go_rows = [], [], []
for (ver, sname), g in SETS.items():
    gf = g[g.index.isin(TIS.index)]
    dec, tis = deciles(gf.power.to_numpy()), TIS.reindex(gf.index).to_numpy()
    for T in sorted(fa.Sex_Tissue.unique()):
        fl_rows.append(dict(version=ver, set=sname, tissue=T, **enrich(gf, dec, tis == T)))
    dec = deciles(g.power.to_numpy())
    for n, ref in ZB.items():
        zb_rows.append(dict(version=ver, set=sname, ref=n, **enrich(g, dec, g.index.isin(ref))))
    for asp, terms in ann.items():
        ga = g[g.index.isin(GENES_ASPECT[asp])]
        dec, posg = deciles(ga.power.to_numpy()), set(ga.index[ga.pos])
        for t, genes in terms.items():
            if len(genes & posg) >= 2:
                go_rows.append(dict(version=ver, set=sname, aspect=asp, term=t, name=names.get(t),
                                    **enrich(ga, dec, ga.index.isin(genes))))
fl, zb, go = pd.DataFrame(fl_rows), pd.DataFrame(zb_rows), pd.DataFrame(go_rows)
bh = lambda p: multipletests(p, method='fdr_bh')[1]
for d, keys in ((fl, ['version', 'set']), (go, ['version', 'set', 'aspect'])):
    d['FDR'] = d.groupby(keys).P.transform(bh)
    d['matched_FDR'] = d.groupby(keys).matched_P.transform(bh)
zb['matched_BH'] = bh(zb.matched_P)
for d, n in ((fl, 'flyatlas2'), (zb, 'zhao_begun'), (go, 'go')):
    d.to_csv(f'{VRES}/L3_lineage_{n}.tsv', sep='\t', index=False)

cols = ['version', 'set', 'tissue', 'overlap', 'fold', 'FDR', 'matched_fold', 'matched_FDR']
print('\n## FlyAtlas2: key tissues in every set')
print(fl[fl.tissue.isin(KEY_TISSUES)].sort_values(['tissue', 'version', 'set'])[cols].to_string(index=False, float_format=fmt))
print('\n## FlyAtlas2: other tissues with FDR < 0.05 in either test')
x = fl[~fl.tissue.isin(KEY_TISSUES) & ((fl.FDR < 0.05) | (fl.matched_FDR < 0.05))]
print(x[cols].to_string(index=False, float_format=fmt) if len(x) else '  none')
print('\n## Zhao & Begun')
print(zb[['version', 'set', 'ref', 'background', 'target', 'overlap', 'fold', 'P', 'matched_fold', 'matched_P',
          'matched_BH']].to_string(index=False, float_format=fmt))
print('\n## GO: number of terms with FDR < 0.05 (unadjusted test / power-matched test), by aspect')
s = go.groupby(['version', 'set', 'aspect']).agg(tested=('P', 'size'), sig=('FDR', lambda v: int((v < 0.05).sum())),
                                                 sig_matched=('matched_FDR', lambda v: int((v < 0.05).sum())))
print(s.unstack('aspect').to_string())
print('\n## GO: key terms named in the draft')
k = go[go.name.isin(KEY_GO)].sort_values(['name', 'version', 'set'])
print(k[['name', 'aspect', 'version', 'set', 'in_category', 'overlap', 'fold', 'FDR', 'matched_fold', 'matched_FDR']]
      .to_string(index=False, float_format=fmt))
print('\n## GO: top 6 terms per set by power-matched P')
top = go.sort_values('matched_P').groupby(['version', 'set']).head(6).sort_values(['version', 'set', 'matched_P'])
print(top[['version', 'set', 'aspect', 'name', 'overlap', 'fold', 'FDR', 'matched_fold', 'matched_FDR']]
      .to_string(index=False, float_format=fmt))

# ---------------------------------------------------------------- direct lineage contrast
o = pd.read_csv(f'{A}/orth_longest/Dsuz_Dsub.direct.RBH.longest.mRNA.tsv', sep='\t', header=None)
s2z = dict(zip(o[0].str.replace('rna-', ''), o[1].str.replace('rna-', '')))
CATS = [('FlyAtlas2', T, lambda u, T=T: u.index.isin(TIS.index), lambda u, T=T: TIS.reindex(u.index).to_numpy() == T)
        for T in KEY_TISSUES]
CATS += [('ZhaoBegun', n, lambda u: np.ones(len(u), bool), lambda u, n=n: u.index.isin(ZB[n])) for n in ZB]
for n in KEY_GO:
    if n in name2id:
        asp, t = name2id[n]
        CATS.append(('GO_' + asp, n, lambda u, asp=asp: u.index.isin(GENES_ASPECT[asp]),
                     lambda u, asp=asp, t=t: u.index.isin(ann[asp][t])))
ct_rows = []
for daf in (15, 5):
    ts = pmk('Dsub_L5_7', daf)
    tz = {p: pmk(p, daf) for p in POPS}
    for zname in ('union',) + POPS:
        zt = list(tz.values()) if zname == 'union' else [tz[zname]]
        zpow = pd.concat([t.power for t in zt], axis=1).max(axis=1)
        zpos = pd.concat([t.pos for t in zt], axis=1).fillna(False).astype(bool).any(axis=1)
        has_z = np.array([s2z.get(x) in zpow.index for x in ts.index])
        u = ts[has_z & ts.index.isin(FB['Dsub'])].copy()
        u['zpos'] = zpos.reindex(u.index.map(s2z)).to_numpy()
        u['zpow'] = zpow.reindex(u.index.map(s2z)).to_numpy()
        u.index = u.index.map(FB['Dsub'])
        ndup = int(u.index.duplicated().sum())
        u = u[~u.index.duplicated()]
        for src, cat, eligible, member in CATS:
            v = u[eligible(u)]
            x = member(v)
            d1, d2 = v.pos.to_numpy() & ~v.zpos.to_numpy(), ~v.pos.to_numpy() & v.zpos.to_numpy()
            tab = [[int((x & d1).sum()), int((x & d2).sum())], [int((~x & d1).sum()), int((~x & d2).sum())]]
            orr, p = fisher_exact(tab)
            dec = deciles((v.power + v.zpow).to_numpy())
            strata = []
            for gq in range(10):
                i = dec == gq
                t2 = np.array([[(x & d1 & i).sum(), (x & d2 & i).sum()], [(~x & d1 & i).sum(), (~x & d2 & i).sum()]])
                if t2.sum(axis=0).min() > 0 and t2.sum(axis=1).min() > 0:
                    strata.append(t2 + 0.0)
            st = StratifiedTable(strata) if strata else None
            ct_rows.append(dict(DAF=daf, Dsuz_lineage=zname, source=src, category=cat, universe=len(v), dup_fbgn=ndup,
                                Dsub_only_in_cat=tab[0][0], Dsuz_only_in_cat=tab[0][1],
                                Dsub_only_other=tab[1][0], Dsuz_only_other=tab[1][1], OR=orr, P=p,
                                MH_OR=st.oddsratio_pooled if st else np.nan,
                                MH_P=st.test_null_odds(correction=True).pvalue if st else np.nan))
ct = pd.DataFrame(ct_rows)
ct['BH'] = ct.groupby(['DAF', 'Dsuz_lineage']).P.transform(bh)
ct['MH_BH'] = ct.groupby(['DAF', 'Dsuz_lineage']).MH_P.transform(lambda p: bh(p.fillna(1)))
ct.to_csv(f'{VRES}/L3_lineage_contrast.tsv', sep='\t', index=False)
print('\n## Direct lineage contrast, Dsuz lineage = union of pools (OR > 1: more frequent among Dsub-lineage-only genes)')
c = ['DAF', 'source', 'category', 'universe', 'Dsub_only_in_cat', 'Dsuz_only_in_cat', 'Dsub_only_other',
     'Dsuz_only_other', 'OR', 'P', 'BH', 'MH_OR', 'MH_P', 'MH_BH']
print(ct[ct.Dsuz_lineage == 'union'][c].to_string(index=False, float_format=fmt))
print('\n## Direct lineage contrast per D. suzukii pool (DAF >= 15%): OR (MH_OR), BH')
w = ct[(ct.DAF == 15) & (ct.Dsuz_lineage != 'union')]
print(w.pivot_table(index='category', columns='Dsuz_lineage', values=['MH_OR', 'MH_BH'], sort=False)
      .to_string(float_format=fmt))
