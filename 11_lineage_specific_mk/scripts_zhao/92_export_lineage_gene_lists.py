#!/usr/bin/env python3
"""Export the gene lists behind the lineage-specific enrichment (script 70; report section 8).

Script 70 builds the MK-positive gene sets in memory and writes only enrichment statistics. This
script rebuilds the same sets with the same rules (MK-positive = P < 0.01 and alpha > 0; FBgn via
Airi's david_mapping; duplicates collapsed to FBgn; union = positive in >= 1 of four pools, background
tested in >= 1 pool; common4 = positive in all four, background tested in all four) and writes:
  results/L3_lineage_gene_sets.tsv
    one row per set x Dmel FBgn in that set's background: MK_positive, power (Pn + Ps + Dn + Ds; max
    over pools for union, min for common4), pools in which a transcript is MK-positive, focal
    transcripts, FlyAtlas2 maximum tissue (TPM > 1).
  results/L3_lineage_category_genes.tsv
    per-transcript MK counts for every tested gene annotated to the categories in the Discussion
    paragraph (script 70 KEY_TISSUES and KEY_GO, plus every GO term whose name matches CHEMO).
    One row per analysis x pool x transcript x category; unpolarized rows are Airi's cmktest
    tables (species = focal species), DAF15/DAF5 rows are the polarized tests (species = lineage).
Checks: background, target, in_category and overlap of every row of the 20261022 FlyAtlas2 and GO
tables are recomputed from the rebuilt sets and must match exactly; MK-positive FBgn in the category
table must equal the rebuilt sets.
"""
import os
import gzip
import re

import pandas as pd

A = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817'
CM = f'{A}/cmktest_results'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
VARIANT = os.environ.get('PMK_VARIANT', '')   # 'excl' = 60_polarized_mk.py EXCLUDE_SHARED=1
PMK = f'{RES}/polarized_mk' + (f'_{VARIANT}' if VARIANT else '')
VRES = f'{PMK}/summary' if VARIANT else RES
GO = f'{PROJ}/inputs/GO_20261022'
POPS = ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')
DSUB = 'Dsub_L5_7'
KEY_TISSUES = ('Male_Testis', 'Male_Midgut', 'Male_Crop', 'Female_Hindgut')
KEY_GO = ('axoneme', 'sperm flagellum', 'cilium movement', 'cilium movement involved in cell motility',
          'inner dynein arm', 'axonemal dynein complex', 'transmembrane transport', 'plasma membrane',
          'membrane', 'detection of chemical stimulus involved in sensory perception of smell',
          'olfactory receptor activity', 'ionotropic olfactory receptor activity', 'midgut development',
          'cell adhesion', 'G protein-coupled receptor signaling pathway', 'proteolysis')
CHEMO = re.compile(r'olfact|smell|odor|taste|gustat|chemosens|pheromone')
SET_NAMES = {'unpolarized': ('Dsub_focal', 'Dsuz_union', 'Dsuz_common4'),
             'DAF15': ('Dsub_lineage', 'Dsuz_lineage_union', 'Dsuz_lineage_common4'),
             'DAF5': ('Dsub_lineage', 'Dsuz_lineage_union', 'Dsuz_lineage_common4')}


def cmk_raw(path):
    t = pd.read_csv(path)
    t.index = t.Name.str.replace(r'^rna-', '', regex=True).str.replace(r'\.best$', '', regex=True)
    assert t.index.is_unique, path
    return pd.DataFrame({'Pn': t.Pn, 'Ps': t.Ps, 'Dn': t.Dn, 'Ds': t.Ds, 'P': t.Pvalue, 'alpha': t.Alpha})


def pmk_raw(focal, daf):
    t = pd.read_csv(f'{PMK}/{focal}.mk_tests.DAF{daf}.tsv', sep='\t', index_col=0)
    assert t.index.is_unique
    return t[['Pn', 'Ps', 'Dn', 'Ds', 'P', 'alpha']]


def flags(t):
    return pd.DataFrame({'power': t.Pn + t.Ps + t.Dn + t.Ds, 'pos': (t.P < 0.01) & (t.alpha > 0)})


def fbgn_map(sp):
    m = pd.read_csv(f'{A}/david_mapping/{sp}_to_Dmel_FBgn.tsv', sep='\t')
    assert m.focal_mrna.is_unique and m.Dmel_FBgn.notna().all()
    return dict(zip(m.focal_mrna, m.Dmel_FBgn)), dict(zip(m.Dmel_FBgn, m.Dmel_symbol_official))


def fb_level(power, pos_ids, fb):
    s = power[power.index.isin(fb)]
    g = pd.DataFrame({'FBgn': s.index.map(fb), 'power': s.values, 'pos': s.index.isin(pos_ids)})
    return g.groupby('FBgn').agg(power=('power', 'max'), pos=('pos', 'any'))


def union_common(tabs, fb):
    pw = pd.concat({p: t.power for p, t in tabs.items()}, axis=1)
    ps = pd.concat({p: t.pos for p, t in tabs.items()}, axis=1).fillna(False).astype(bool)
    allidx = pw.dropna().index
    u = fb_level(pw.max(axis=1), set(pw.index[ps.any(axis=1)]), fb)
    c = fb_level(pw.loc[allidx].min(axis=1), set(allidx[ps.loc[allidx].all(axis=1)]), fb)
    return u, c


(FB_SUB, SYM), (FB_SUZ, SYM2) = fbgn_map('Dsub'), fbgn_map('Dsuz')
SYM.update(SYM2)
FB = {'Dsub': FB_SUB, 'Dsuz': FB_SUZ}
RAW = {('unpolarized', 'Dsub', DSUB): cmk_raw(f'{CM}/Dsub_pairwise/Dsub_vs_Dsuz_mktest.filtered.csv')}
RAW.update({('unpolarized', 'Dsuz', p): cmk_raw(f'{CM}/Dsuz_pairwise/{p}/Dsuz_vs_Dsub_mktest.filtered.csv')
            for p in POPS})
for daf in (15, 5):
    RAW[(f'DAF{daf}', 'Dsub', DSUB)] = pmk_raw(DSUB, daf)
    RAW.update({(f'DAF{daf}', 'Dsuz', p): pmk_raw(p, daf) for p in POPS})
SETS = {}
for ver, (s_sub, s_un, s_c4) in SET_NAMES.items():
    SETS[(ver, s_sub)] = fb_level(*(lambda f: (f.power, set(f.index[f.pos])))(flags(RAW[(ver, 'Dsub', DSUB)])),
                                  FB['Dsub'])
    SETS[(ver, s_un)], SETS[(ver, s_c4)] = union_common({p: flags(RAW[(ver, 'Dsuz', p)]) for p in POPS}, FB['Dsuz'])

# ---------------------------------------------------------------- annotations (as in script 70)
fa = pd.read_csv(f'{PROJ}/FLYA2_TAU_Genome_NS.txt', sep='\t')
assert fa.Gene_ID.is_unique and fa.Sex_Tissue.nunique() == 33
TIS = fa[fa.TPM > 1].set_index('Gene_ID').Sex_Tissue
names, cur = {}, None
for line in open(f'{GO}/go-basic.obo'):
    line = line.rstrip('\n')
    if line == '[Term]':
        cur = {}
    elif line.startswith('[') or not line:
        if cur and 'id' in cur:
            names[cur['id']] = cur.get('name')
            for alt in cur.get('alt', []):
                names.setdefault(alt, cur.get('name'))
        cur = None if line.startswith('[') else cur
    elif cur is not None:
        k, _, v = line.partition(': ')
        if k == 'alt_id':
            cur.setdefault('alt', []).append(v)
        elif k in ('id', 'name'):
            cur[k] = v
ann = {'P': {}, 'F': {}, 'C': {}}
with gzip.open(f'{GO}/fb.gaf.gz', 'rt') as f:
    for line in f:
        if line.startswith('!'):
            continue
        c = line.rstrip('\n').split('\t')
        if c[0] != 'FB' or 'taxon:7227' not in c[12] or 'NOT' in c[3].split('|'):
            continue
        ann[c[8]].setdefault(c[4], set()).add(c[1])
        SYM.setdefault(c[1], c[2])
GENES_ASPECT = {a: set().union(*d.values()) for a, d in ann.items()}

# ---------------------------------------------------------------- check against the 20261022 tables
for kind, path in (('fl', 'L3_lineage_flyatlas2.tsv'), ('go', 'L3_lineage_go.tsv')):
    old = pd.read_csv(f'{VRES}/{path}', sep='\t')
    sub = {}
    for r in old.itertuples(index=False):
        key = (r.version, r.set, r.aspect if kind == 'go' else '')
        if key not in sub:
            g = SETS[(r.version, r.set)]
            sub[key] = g[g.index.isin(GENES_ASPECT[r.aspect] if kind == 'go' else TIS.index)]
        g = sub[key]
        m = g.index.isin(ann[r.aspect][r.term]) if kind == 'go' else TIS.reindex(g.index).to_numpy() == r.tissue
        new = (len(g), int(g.pos.sum()), int(m.sum()), int((g.pos.to_numpy() & m).sum()))
        assert new == (r.background, r.target, r.in_category, r.overlap), (path, r, new)
    print(f'check: all {len(old)} rows of {path} reproduced (background, target, in_category, overlap)')

# ---------------------------------------------------------------- gene-set table
rows = []
for (ver, sname), g in SETS.items():
    sp = 'Dsub' if sname.startswith('Dsub') else 'Dsuz'
    pools = (DSUB,) if sp == 'Dsub' else POPS
    d = pd.concat([flags(RAW[(ver, sp, p)]).assign(pool=p) for p in pools])
    d = d[d.index.isin(FB[sp])]
    d['FBgn'] = d.index.map(FB[sp])
    tx = d.groupby('FBgn').apply(lambda x: ','.join(sorted(set(x.index))))
    pp = d[d.pos].groupby('FBgn').pool.agg(lambda s: ','.join(p for p in pools if p in set(s)))
    rows.append(pd.DataFrame({'version': ver, 'set': sname, 'FBgn': g.index,
                              'symbol': g.index.map(lambda x: SYM.get(x, '')),
                              'MK_positive': g.pos.values, 'power': g.power.values,
                              'positive_pools': pp.reindex(g.index).fillna('').values,
                              'transcripts': tx.reindex(g.index).values,
                              'FlyAtlas2_max_tissue': TIS.reindex(g.index).fillna('').values}))
gs = pd.concat(rows, ignore_index=True)
gs.to_csv(f'{VRES}/L3_lineage_gene_sets.tsv', sep='\t', index=False)
print('\n## Gene sets written (FBgn level): background, MK-positive')
print(gs.groupby(['version', 'set'], sort=False).MK_positive.agg(['size', 'sum']).to_string())

# ---------------------------------------------------------------- category member table
CATS = [('FlyAtlas2', T, set(TIS.index[TIS == T])) for T in KEY_TISSUES]
for asp in 'PFC':
    for t in sorted(ann[asp]):
        n = names.get(t)
        if n and (n in KEY_GO or CHEMO.search(n)):
            CATS.append((f'GO_{asp}', f'{n} ({t})', ann[asp][t]))
rows = []
for (ver, sp, pool), t in RAW.items():
    t = t[t.index.isin(FB[sp])]
    fbg = pd.Series(t.index.map(FB[sp]), index=t.index)
    for src, cat, members in CATS:
        i = fbg.isin(members).to_numpy()
        if i.any():
            d = t[i]
            rows.append(pd.DataFrame({'analysis': ver, 'species': sp, 'pool': pool, 'category_source': src,
                                      'category': cat, 'FBgn': fbg[i].values,
                                      'symbol': [SYM.get(x, '') for x in fbg[i]], 'transcript': d.index,
                                      **{c: d[c].values for c in d.columns},
                                      'MK_positive': ((d.P < 0.01) & (d.alpha > 0)).values}))
cg = pd.concat(rows, ignore_index=True)
for ver, (s_sub, s_un, _) in SET_NAMES.items():
    for sp, sname in (('Dsub', s_sub), ('Dsuz', s_un)):
        g = SETS[(ver, sname)]
        for src, cat, members in CATS:
            x = cg[(cg.analysis == ver) & (cg.species == sp) & (cg.category == cat) & cg.MK_positive]
            assert set(x.FBgn) == set(g.index[g.pos]) & members, (ver, sname, cat)
cg.to_csv(f'{VRES}/L3_lineage_category_genes.tsv', sep='\t', index=False)
print(f'\ncategory table: {len(cg)} rows, {len(CATS)} categories; MK-positive FBgn match the rebuilt sets')

print('\n## MK-positive genes in chemosensory and midgut categories, DAF >= 15% and unpolarized')
for ver in ('DAF15', 'unpolarized'):
    for sp in ('Dsub', 'Dsuz'):
        x = cg[(cg.analysis == ver) & (cg.species == sp) & cg.MK_positive &
               (cg.category.str.contains(CHEMO) | cg.category.str.contains('midgut', case=False))]
        for cat, y in x.groupby('category'):
            npool = y.groupby('symbol').pool.nunique()
            print(f'  {ver:11s} {sp} {cat}: {len(npool)} genes: '
                  + ', '.join(f'{s}({n})' if sp == 'Dsuz' else s for s, n in npool.items()))
