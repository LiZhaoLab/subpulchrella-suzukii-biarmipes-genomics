#!/usr/bin/env python3
"""Lineage-specific candidate genes and a genome-wide lineage contrast, from the polarized MK test
with sister-pool exclusion (results/polarized_mk_excl; 60_polarized_mk.py EXCLUDE_SHARED=1).

Question (Li Zhao): is it true that the polarized MK test shows no functional pattern
other than testis genes? Which individual genes are the strongest candidates for lineage-specific
adaptive protein evolution on the D. subpulchrella and on the D. suzukii lineage?

Data: <pool>.mk_tests.DAF15.tsv (primary; polymorphisms with derived allele frequency >= 15%) and
DAF5 (sensitivity). MK-positive = Fisher P < 0.01 and alpha > 0 (manuscript definition).
Genes = D. subpulchrella transcripts tested on the D. subpulchrella lineage whose reciprocal-best-hit
D. suzukii ortholog (Airi's orth_longest RBH table) is tested in >= 1 D. suzukii pool; D. melanogaster
FBgn from Airi's david_mapping (D. subpulchrella side).

A. Per-gene table: counts, P, alpha, BH q (over genes, per pool) and DoS on both lineages; D. suzukii
   lineage per pool and the number of pools MK-positive. Fixed differences on the D. suzukii lineage
   are nearly the same in all pools; polymorphism differs, so a D. suzukii-lineage gene claim needs
   >= 3 of 4 pools.
B. Gene-level lineage test: Fisher's exact test (two-sided) of [Dn, Ds] on the D. subpulchrella
   lineage vs [Dn, Ds] on the D. suzukii lineage (median over pools in which the gene is tested,
   rounded). It asks whether the nonsynonymous fraction of fixed differences differs between the two
   lineages. BH over all genes.
   Candidate classes (fixed here, before any list was read):
   Dsub-lineage: q < 0.05 and alpha > 0 on the D. subpulchrella lineage; no D. suzukii pool with
     P < 0.05 and alpha > 0; lineage test P < 0.05 with the higher Dn/Ds on the D. subpulchrella lineage.
   Dsuz-lineage: tested in >= 3 pools; MK-positive in >= 3 pools; median q < 0.05; D. subpulchrella
     lineage not P < 0.05 with alpha > 0; lineage test P < 0.05 with the higher Dn/Ds on the
     D. suzukii lineage.
   Both: q < 0.05 and alpha > 0 on the D. subpulchrella lineage, MK-positive in >= 3 pools and median
     q < 0.05.
C. Genome-wide GO contrast (script 70 tested only 16 preselected terms): universe and discordant
   genes as in script 70 (genes MK-positive on exactly one lineage); every direct GO annotation
   (FlyBase GAF, NOT removed) with >= 5 discordant annotated genes; two-sided Fisher and
   Mantel-Haenszel (strata = deciles of summed Pn+Ps+Dn+Ds on both lineages); BH within aspect,
   D. suzukii lineage = union of pools and each pool.
D. Per-pool power-matched GO enrichment on the D. suzukii lineage (script 70 reports only the union
   and the 4-pool intersection): same exact power-matched test as script 70; BH within aspect.
E. Ecological categories, fixed here before any list was read (D. melanogaster symbols and GO terms
   with is_a/part_of descendants): chemosensory, mechanosensory, cuticle/chitin, detoxification,
   immunity, temperature/circadian/diapause, oviposition/female reproduction, pigmentation, male
   reproduction (GO only). Direct contrast as in C (Fisher, MH), BH over categories.
Sweep peaks: transcripts within 5 kb of a CLR peak (results/F_sweep_peaks.tsv), D. subpulchrella
pool for the D. subpulchrella lineage, number of D. suzukii pools for the D. suzukii lineage.
Outputs: results/lineage_genes.tsv (all genes), lineage_candidates.tsv,
lineage_go_contrast_all.tsv, lineage_go_by_pool.tsv,
lineage_ecological_categories.tsv, lineage_candidates.log.
QC (added after the first run): Ds balance between lineages (two-sided binomial vs the genome-wide share),
Dn excess on the focal lineage (one-sided binomial), callable-codon ratio outside the 1-99% range, divergence
per codon above the 99.5% quantile. robust = candidate class, no QC flag, Dn excess P < 0.05 for one-lineage
classes. No random numbers are used.
"""
import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binom, binom_test, fisher_exact, hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

A = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
PMK = f'{RES}/polarized_mk_excl'
GO = f'{PROJ}/inputs/GO_20261022'
POPS = ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')
ARM = {'NC_050608.1': '4', 'NC_050609.1': '3R', 'NC_050610.1': '2L', 'NC_050611.1': '2R',
       'NC_050612.1': '3L', 'NC_050613.1': 'X', 'NW_023665577.1': '3R_scaffold'}
MECHANO = ['nompC', 'iav', 'nan', 'pain', 'pyx', 'TrpA1', 'trp', 'trpl', 'Trpgamma', 'wtrw', 'Trpm',
           'Trpml', 'Pkd2', 'brv1', 'brv2', 'brv3', 'Piezo', 'Tmc']
# ecological categories: (symbol regex or None, symbol list, GO root names)
ECO = {
    'chemosensory': (r'^(Or\d|Orco$|Gr\d|Ir\d|Obp|ppk|lush$)', [], ['sensory perception of chemical stimulus']),
    'mechanosensory': (None, MECHANO, ['sensory perception of mechanical stimulus']),
    'cuticle_chitin': (r'^(Cpr|Lcp|Edg|Twdl|Ccp)', [], ['structural constituent of cuticle', 'chitin metabolic process']),
    'detoxification': (r'^(Cyp|Gst|Ugt)', [], ['response to xenobiotic stimulus', 'xenobiotic transport']),
    'immunity': (None, [], ['immune response']),
    'temperature_circadian_diapause': (None, [], ['response to temperature stimulus', 'circadian rhythm', '*diapause*']),
    'oviposition_female_reproduction': (None, [], ['egg-laying behavior', 'female mating behavior']),
    'pigmentation': (None, [], ['pigmentation']),
    'male_reproduction_GO': (None, [], ['spermatogenesis', 'sperm competition', 'insemination']),
}
pd.set_option('display.width', 250)
pd.set_option('display.max_colwidth', 60)
LOG = open(f'{RES}/lineage_candidates.log', 'w')


def log(s=''):
    print(s, flush=True)
    LOG.write(s + '\n')


bh = lambda p: multipletests(np.asarray(p, float), method='fdr_bh')[1]


def deciles(v):
    return pd.qcut(pd.Series(v).rank(method='first'), 10, labels=False).to_numpy()


def matched_test(dec, pos, member):
    """Exact power-matched expectation and P(overlap >= observed); identical to script 70."""
    pmf = np.array([1.0])
    for g in range(10):
        i = dec == g
        n, k, m = int(i.sum()), int(pos[i].sum()), int(member[i].sum())
        pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(k, m) + 1), n, m, k))
    obs = int((pos & member).sum())
    exp = float((np.arange(len(pmf)) * pmf).sum())
    return exp, obs / exp if exp else np.nan, float(pmf[obs:].sum())


def contrast(x, d1, d2, dec):
    """x: category membership; d1/d2: Dsub-only / Dsuz-only; returns counts, OR, P, MH OR, MH P."""
    tab = [[int((x & d1).sum()), int((x & d2).sum())], [int((~x & d1).sum()), int((~x & d2).sum())]]
    orr, p = fisher_exact(tab)
    strata = []
    for q in range(10):
        i = dec == q
        t2 = np.array([[(x & d1 & i).sum(), (x & d2 & i).sum()], [(~x & d1 & i).sum(), (~x & d2 & i).sum()]])
        if t2.sum(axis=0).min() > 0 and t2.sum(axis=1).min() > 0:
            strata.append(t2 + 0.0)
    st = StratifiedTable(strata) if strata else None
    return dict(Dsub_only_in_cat=tab[0][0], Dsuz_only_in_cat=tab[0][1], Dsub_only_other=tab[1][0],
                Dsuz_only_other=tab[1][1], OR=orr, P=p, MH_OR=st.oddsratio_pooled if st else np.nan,
                MH_P=st.test_null_odds(correction=True).pvalue if st else np.nan)


# ------------------------------------------------------------------ inputs
def mk(focal, daf):
    t = pd.read_csv(f'{PMK}/{focal}.mk_tests.DAF{daf}.tsv', sep='\t', index_col=0)
    assert t.index.is_unique and {'Pn', 'Ps', 'Dn', 'Ds', 'P', 'alpha', 'q'} <= set(t.columns)
    return t


m = pd.read_csv(f'{A}/david_mapping/Dsub_to_Dmel_FBgn.tsv', sep='\t')
assert m.focal_mrna.is_unique and m.Dmel_FBgn.notna().all()
FB = dict(zip(m.focal_mrna, m.Dmel_FBgn))
o = pd.read_csv(f'{A}/orth_longest/Dsuz_Dsub.direct.RBH.longest.mRNA.tsv', sep='\t', header=None)
S2Z = dict(zip(o[0].str.replace('rna-', ''), o[1].str.replace('rna-', '')))
assert len(S2Z) == len(o)
chrom = pd.read_csv(f'{PMK}/Dsub_L5_7.genes.tsv', sep='\t', usecols=['gene', 'chrom']).set_index('gene').chrom

# GO: names, DAG, annotations
names, parents, cur = {}, defaultdict(set), None
for line in open(f'{GO}/go-basic.obo'):
    line = line.rstrip('\n')
    if line == '[Term]':
        cur = None
    elif line.startswith('id: GO:'):
        cur = line[4:]
    elif cur and line.startswith('name: '):
        names[cur] = line[6:]
    elif cur and line.startswith('is_a: '):
        parents[cur].add(line[6:16])
    elif cur and line.startswith('relationship: part_of '):
        parents[cur].add(line[22:32])
    elif line.startswith('[Typedef]'):
        cur = None
children = defaultdict(set)
for c, ps in parents.items():
    for p in ps:
        children[p].add(c)


def descendants(root):
    out, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


direct = {'P': defaultdict(set), 'F': defaultdict(set), 'C': defaultdict(set)}
term_genes = defaultdict(set)
SYM, FULL = {}, {}
with gzip.open(f'{GO}/fb.gaf.gz', 'rt') as f:
    for line in f:
        if line.startswith('!'):
            continue
        c = line.rstrip('\n').split('\t')
        if c[0] != 'FB' or 'taxon:7227' not in c[12]:
            continue
        SYM.setdefault(c[1], c[2])
        FULL.setdefault(c[1], c[9])
        if 'NOT' in c[3].split('|'):
            continue
        direct[c[8]][c[4]].add(c[1])
        term_genes[c[4]].add(c[1])
fa = pd.read_csv(f'{PROJ}/FLYA2_TAU_Genome_NS.txt', sep='\t')
assert fa.Gene_ID.is_unique
TIS = fa[fa.TPM > 1].set_index('Gene_ID').Sex_Tissue
TAU = fa.set_index('Gene_ID').Tau
ZB638 = set(pd.read_csv(f'{PROJ}/inputs/ZhaoBegun2017_S13_DOS_638.tsv', sep='\t').FBgn)
gs = pd.read_csv(f'{RES}/L3_lineage_gene_sets.tsv', sep='\t')
for fbg, s in zip(gs.FBgn, gs.symbol.astype(str)):
    SYM.setdefault(fbg, s)

# ecological category membership by FBgn (fixed before reading any MK result)
name2id = defaultdict(list)
for t, n in names.items():
    name2id[n].append(t)
ECO_FB, ECO_RULE = {}, {}
all_fb = set(FB.values())
for cat, (pat, lst, roots) in ECO.items():
    terms = set()
    for r in roots:
        if r.startswith('*'):
            ids = [t for t, n in names.items() if r.strip('*') in n]
        else:
            ids = name2id.get(r, [])
        assert ids, f'GO term not found: {r}'
        for t in ids:
            terms |= descendants(t)
    genes = set().union(*(term_genes[t] for t in terms)) if terms else set()
    if pat:
        genes |= {g for g in all_fb if re.search(pat, SYM.get(g, ''))}
    genes |= {g for g in all_fb if SYM.get(g, '') in lst}
    ECO_FB[cat] = genes & all_fb
    ECO_RULE[cat] = f'regex {pat}; list {len(lst)}; GO roots {roots} ({len(terms)} terms incl. descendants)'
log(f'GO: go-basic.obo terms {len(names)}; GAF FlyBase taxon 7227 genes {len(SYM)}')
log('Ecological categories (D. melanogaster genes with a D. subpulchrella transcript in the FBgn map):')
for cat in ECO:
    log(f'  {cat:34s} {len(ECO_FB[cat]):5d}  {ECO_RULE[cat]}')

# sweep peaks
pk = pd.read_csv(f'{RES}/F_sweep_peaks.tsv', sep='\t')
PEAK = {p: set(x for t in g.transcripts.fillna('').astype(str) for x in t.split(',') if x)
        for p, g in pk.groupby('pool')}
log(f'Sweep-peak transcripts per pool: ' + ', '.join(f'{p} {len(v)}' for p, v in PEAK.items()))


# ------------------------------------------------------------------ A, B: per-gene table
def gene_table(daf):
    ts = mk('Dsub_L5_7', daf)
    tz = {p: mk(p, daf) for p in POPS}
    rows = []
    for g in ts.index:
        z = S2Z.get(g)
        if z is None or g not in FB:
            continue
        r = dict(dsub_tx=g, dsuz_tx=z, FBgn=FB[g])
        s = ts.loc[g]
        r.update(Pn_sub=s.Pn, Ps_sub=s.Ps, Dn_sub=s.Dn, Ds_sub=s.Ds, P_sub=s.P, alpha_sub=s.alpha, q_sub=s.q)
        npos = ntest = 0
        zc = []
        for p in POPS:
            if z in tz[p].index:
                x = tz[p].loc[z]
                ntest += 1
                pos = (x.P < 0.01) and (x.alpha > 0)
                npos += pos
                zc.append((x.Pn, x.Ps, x.Dn, x.Ds, x.P, x.alpha, x.q, x.callable))
                r.update({f'P_{p}': x.P, f'alpha_{p}': x.alpha, f'q_{p}': x.q})
            else:
                r.update({f'P_{p}': np.nan, f'alpha_{p}': np.nan, f'q_{p}': np.nan})
        if ntest == 0:
            continue
        zc = np.array(zc, float)
        med = np.median(zc, axis=0)
        r.update(n_pools_tested=ntest, n_pools_pos=npos,
                 n_pools_P05=int(((zc[:, 4] < 0.05) & (zc[:, 5] > 0)).sum()),
                 Pn_suz=med[0], Ps_suz=med[1], Dn_suz=np.round(med[2]), Ds_suz=np.round(med[3]),
                 P_suz_median=med[4], q_suz_median=med[6], alpha_suz_median=med[5],
                 callable_sub=s.callable, callable_suz=med[7])
        rows.append(r)
    t = pd.DataFrame(rows)
    assert t.dsub_tx.is_unique
    with np.errstate(divide='ignore', invalid='ignore'):
        for sfx in ('sub', 'suz'):
            t[f'DoS_{sfx}'] = t[f'Dn_{sfx}'] / (t[f'Dn_{sfx}'] + t[f'Ds_{sfx}']) - t[f'Pn_{sfx}'] / (t[f'Pn_{sfx}'] + t[f'Ps_{sfx}'])
        lp, lo = [], []
        for a, b, c, d in zip(t.Dn_sub, t.Ds_sub, t.Dn_suz, t.Ds_suz):
            lp.append(fisher_exact([[int(a), int(b)], [int(c), int(d)]])[1] if a + b > 0 and c + d > 0 else np.nan)
            lo.append(((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5)))
    t['lineage_P'] = lp
    t['lineage_OR'] = lo   # (Dn/Ds Dsub lineage) / (Dn/Ds Dsuz lineage), 0.5 added to each count
    ok = t.lineage_P.notna()
    t['lineage_BH'] = np.nan
    t.loc[ok, 'lineage_BH'] = bh(t.loc[ok, 'lineage_P'])
    pos_sub = (t.P_sub < 0.01) & (t.alpha_sub > 0)
    sub_strong = (t.q_sub < 0.05) & (t.alpha_sub > 0)
    sub_weak = (t.P_sub < 0.05) & (t.alpha_sub > 0)
    suz_strong = (t.n_pools_tested >= 3) & (t.n_pools_pos >= 3) & (t.q_suz_median < 0.05)
    t['MK_pos_sub'] = pos_sub
    t['MK_pos_suz_union'] = t.n_pools_pos >= 1
    cls = np.full(len(t), '', object)
    cls[(sub_strong & (t.n_pools_P05 == 0) & (t.lineage_P < 0.05) & (t.lineage_OR > 1)).to_numpy()] = 'Dsub_lineage'
    cls[(suz_strong & ~sub_weak & (t.lineage_P < 0.05) & (t.lineage_OR < 1)).to_numpy()] = 'Dsuz_lineage'
    cls[(sub_strong & suz_strong).to_numpy()] = 'both'
    t['candidate_class'] = cls
    # QC (added after CG30427 showed 266 D. suzukii-lineage fixed differences in 471 codons).
    # Synonymous fixed differences should split between sister lineages as in the genome (p0_s); a gene
    # whose MK result rests on a chance (or artefactual) Ds deficit is less credible than one with an
    # excess of Dn at balanced Ds. Binomial tests against the genome-wide shares of this gene set.
    p0_s = t.Ds_sub.sum() / (t.Ds_sub.sum() + t.Ds_suz.sum())
    p0_n = t.Dn_sub.sum() / (t.Dn_sub.sum() + t.Dn_suz.sum())
    ns, nn = (t.Ds_sub + t.Ds_suz).astype(int), (t.Dn_sub + t.Dn_suz).astype(int)
    t['Ds_balance_P'] = [binom_test(int(a), int(n), p0_s) if n > 0 else np.nan for a, n in zip(t.Ds_sub, ns)]
    t['Dn_excess_sub_P'] = [binom.sf(int(a) - 1, int(n), p0_n) if n > 0 else np.nan for a, n in zip(t.Dn_sub, nn)]
    t['Dn_excess_suz_P'] = [binom.cdf(int(a), int(n), p0_n) if n > 0 else np.nan for a, n in zip(t.Dn_sub, nn)]
    t['codon_ratio'] = t.callable_suz / t.callable_sub
    lo_r, hi_r = t.codon_ratio.quantile([0.01, 0.99])
    for sfx in ('sub', 'suz'):
        t[f'div_per_codon_{sfx}'] = (t[f'Dn_{sfx}'] + t[f'Ds_{sfx}']) / t[f'callable_{sfx}']
    hi_div = max(t.div_per_codon_sub.quantile(0.995), t.div_per_codon_suz.quantile(0.995))
    flag = []
    for r in t.itertuples():
        f = []
        if not (lo_r <= r.codon_ratio <= hi_r):
            f.append('callable_mismatch')
        if max(r.div_per_codon_sub, r.div_per_codon_suz) > hi_div:
            f.append('extreme_divergence')
        if r.Ds_balance_P < 0.05:
            f.append('Ds_imbalance')
        flag.append(','.join(f))
    t['qc_flags'] = flag
    dn_ok = np.where(t.candidate_class == 'Dsub_lineage', t.Dn_excess_sub_P < 0.05,
                     np.where(t.candidate_class == 'Dsuz_lineage', t.Dn_excess_suz_P < 0.05, True))
    t['robust'] = (t.candidate_class != '') & (t.qc_flags == '') & dn_ok
    t.attrs['qc'] = (f'p0_s {p0_s:.4f}, p0_n {p0_n:.4f}, codon ratio 1-99% {lo_r:.3f}-{hi_r:.3f}, '
                     f'divergence per codon 99.5% {hi_div:.4f}')
    t['arm'] = t.dsub_tx.map(chrom).map(ARM)
    t['symbol'] = t.FBgn.map(SYM)
    t['name'] = t.FBgn.map(FULL)
    t['FlyAtlas2_max'] = t.FBgn.map(TIS)
    t['tau'] = t.FBgn.map(TAU)
    t['ZhaoBegun638'] = t.FBgn.isin(ZB638)
    t['eco'] = [','.join(c for c in ECO if f in ECO_FB[c]) for f in t.FBgn]
    t['sweep_Dsub'] = t.dsub_tx.isin(PEAK.get('Dsub_L5_7', set()))
    t['sweep_Dsuz_pools'] = [sum(z in PEAK.get(p, set()) for p in POPS) for z in t.dsuz_tx]
    return t


T = {daf: gene_table(daf) for daf in (15, 5)}
t = T[15]
t5 = T[5].set_index('dsub_tx')
t['class_DAF5'] = t.dsub_tx.map(t5.candidate_class).fillna('untested')
t.to_csv(f'{RES}/lineage_genes.tsv', sep='\t', index=False)
log(f'\n## A/B. Genes tested on both lineages (DAF >= 15%): {len(t)} transcripts, {t.FBgn.nunique()} FBgn; '
    f'tested in 4 pools {int((t.n_pools_tested == 4).sum())}')
log(f'MK-positive Dsub lineage {int(t.MK_pos_sub.sum())}; Dsuz lineage >=1 pool {int((t.n_pools_pos >= 1).sum())}, '
    f'>=3 pools {int((t.n_pools_pos >= 3).sum())}, 4 pools {int((t.n_pools_pos == 4).sum())}')
log(f'Lineage test (Dn/Ds differs between lineages): tested {int(t.lineage_P.notna().sum())}; P < 0.05 '
    f'{int((t.lineage_P < 0.05).sum())} (Dsub higher {int(((t.lineage_P < 0.05) & (t.lineage_OR > 1)).sum())}, '
    f'Dsuz higher {int(((t.lineage_P < 0.05) & (t.lineage_OR < 1)).sum())}); BH < 0.05 {int((t.lineage_BH < 0.05).sum())}; '
    f'BH < 0.1 {int((t.lineage_BH < 0.1).sum())}')
exp05 = 0.05 * t.lineage_P.notna().sum()
log(f'  expected P < 0.05 under the null (uniform P): {exp05:.0f}; Fisher P is conservative for small counts')
for c in ('Dsub_lineage', 'Dsuz_lineage', 'both'):
    x = t[t.candidate_class == c]
    log(f'Candidate class {c}: {len(x)} genes; same class at DAF >= 5%: {int((x.class_DAF5 == c).sum())}; '
        f'lineage BH < 0.1: {int((x.lineage_BH < 0.1).sum())}')
log(f'QC thresholds: {t.attrs.get("qc")}')
for c in ('Dsub_lineage', 'Dsuz_lineage', 'both'):
    x = t[t.candidate_class == c]
    log(f'  {c}: flags ' + ', '.join(f'{k} {int(x.qc_flags.str.contains(k).sum())}' for k in
        ('callable_mismatch', 'extreme_divergence', 'Ds_imbalance')) + f'; robust {int(x.robust.sum())}')
cand = t[t.candidate_class != ''].copy()
cols = ['candidate_class', 'symbol', 'name', 'arm', 'Pn_sub', 'Ps_sub', 'Dn_sub', 'Ds_sub', 'P_sub', 'q_sub',
        'Pn_suz', 'Ps_suz', 'Dn_suz', 'Ds_suz', 'n_pools_pos', 'q_suz_median', 'lineage_OR', 'lineage_P',
        'lineage_BH', 'FlyAtlas2_max', 'tau', 'ZhaoBegun638', 'eco', 'sweep_Dsub', 'sweep_Dsuz_pools',
        'class_DAF5', 'Ds_balance_P', 'Dn_excess_sub_P', 'Dn_excess_suz_P', 'codon_ratio', 'div_per_codon_sub',
        'div_per_codon_suz', 'qc_flags', 'robust', 'FBgn', 'dsub_tx', 'dsuz_tx']
cand = cand.sort_values(['candidate_class', 'robust', 'lineage_P'], ascending=[True, False, True])[cols]
cand.to_csv(f'{RES}/lineage_candidates.tsv', sep='\t', index=False)
short = ['symbol', 'arm', 'Pn_sub', 'Ps_sub', 'Dn_sub', 'Ds_sub', 'q_sub', 'Pn_suz', 'Ps_suz', 'Dn_suz', 'Ds_suz',
         'n_pools_pos', 'q_suz_median', 'lineage_OR', 'lineage_P', 'FlyAtlas2_max', 'eco', 'sweep_Dsub',
         'sweep_Dsuz_pools', 'class_DAF5', 'Ds_balance_P', 'qc_flags', 'robust']
for c in ('Dsub_lineage', 'Dsuz_lineage', 'both'):
    log(f'\n### {c} candidates')
    log(cand[cand.candidate_class == c][short].to_string(index=False, float_format=lambda v: f'{v:.3g}'))

# tissue and category composition of the candidate classes (descriptive)
log('\n## FlyAtlas2 maximum tissue of candidate classes (genes with TPM > 1 in some tissue)')
log(pd.crosstab(t.FlyAtlas2_max.fillna('none'), t.candidate_class.replace('', 'other')).loc[
    lambda d: d.drop(columns='other').sum(axis=1) > 0].to_string())

# ------------------------------------------------------------------ C, E: genome-wide and ecological contrasts
GENES_ASPECT = {a: set().union(*d.values()) for a, d in direct.items()}


def contrast_universe(daf, zname):
    tt = T[daf]
    pcols = [f'P_{p}' for p in POPS] if zname == 'union' else [f'P_{zname}']
    acols = [f'alpha_{p}' for p in POPS] if zname == 'union' else [f'alpha_{zname}']
    tested = tt[pcols].notna().any(axis=1)
    zpos = ((tt[pcols].to_numpy() < 0.01) & (tt[acols].to_numpy() > 0)).any(axis=1)
    u = tt[tested].copy()
    u['zpos'] = zpos[tested.to_numpy()]
    u['pow'] = u.Pn_sub + u.Ps_sub + u.Dn_sub + u.Ds_sub + u.Pn_suz + u.Ps_suz + u.Dn_suz + u.Ds_suz
    u = u.drop_duplicates('FBgn').set_index('FBgn')
    return u


go_rows, eco_rows = [], []
for daf in (15, 5):
    for zname in ('union',) + POPS:
        u = contrast_universe(daf, zname)
        d1 = (u.MK_pos_sub & ~u.zpos).to_numpy()
        d2 = (~u.MK_pos_sub & u.zpos).to_numpy()
        dec = deciles(u['pow'].to_numpy())
        disc = d1 | d2
        for asp, terms in direct.items():
            el = u.index.isin(GENES_ASPECT[asp])
            for term, genes in terms.items():
                x = u.index.isin(genes)
                if (x & disc).sum() < 5:
                    continue
                r = contrast(x[el], d1[el], d2[el], dec[el])
                go_rows.append(dict(DAF=daf, Dsuz_lineage=zname, aspect=asp, term=term, name=names.get(term), **r))
        for cat in ECO:
            x = u.index.isin(ECO_FB[cat])
            r = contrast(x, d1, d2, dec)
            eco_rows.append(dict(DAF=daf, Dsuz_lineage=zname, category=cat, genes_in_universe=int(x.sum()), **r))
gc = pd.DataFrame(go_rows)
gc['BH'] = gc.groupby(['DAF', 'Dsuz_lineage', 'aspect']).P.transform(bh)
gc['MH_BH'] = gc.groupby(['DAF', 'Dsuz_lineage', 'aspect']).MH_P.transform(lambda p: bh(p.fillna(1)))
gc.to_csv(f'{RES}/lineage_go_contrast_all.tsv', sep='\t', index=False)
ec = pd.DataFrame(eco_rows)
ec['BH'] = ec.groupby(['DAF', 'Dsuz_lineage']).P.transform(bh)
ec['MH_BH'] = ec.groupby(['DAF', 'Dsuz_lineage']).MH_P.transform(lambda p: bh(p.fillna(1)))
ec.to_csv(f'{RES}/lineage_ecological_categories.tsv', sep='\t', index=False)
log('\n## C. Genome-wide GO contrast (direct annotations, >= 5 discordant genes): terms tested and BH < 0.05')
log(gc.groupby(['DAF', 'Dsuz_lineage', 'aspect']).agg(tested=('P', 'size'), BH05=('BH', lambda v: int((v < 0.05).sum())),
                                                       MH_BH05=('MH_BH', lambda v: int((v < 0.05).sum())),
                                                       minBH=('BH', 'min')).to_string())
log('\nTop 15 terms, DAF >= 15%, Dsuz lineage = union (OR > 1: more frequent among Dsub-lineage-only genes)')
c2 = ['aspect', 'name', 'Dsub_only_in_cat', 'Dsuz_only_in_cat', 'OR', 'P', 'BH', 'MH_OR', 'MH_BH']
log(gc[(gc.DAF == 15) & (gc.Dsuz_lineage == 'union')].sort_values('P').head(15)[c2].to_string(
    index=False, float_format=lambda v: f'{v:.3g}'))
log('\n## E. Ecological categories, direct lineage contrast (OR > 1: more frequent among Dsub-lineage-only genes)')
c3 = ['DAF', 'Dsuz_lineage', 'category', 'genes_in_universe', 'Dsub_only_in_cat', 'Dsuz_only_in_cat', 'Dsub_only_other',
      'Dsuz_only_other', 'OR', 'P', 'BH', 'MH_OR', 'MH_BH']
log(ec[ec.Dsuz_lineage == 'union'][c3].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
log('\nPer pool, DAF >= 15%: MH OR (MH BH)')
w = ec[(ec.DAF == 15) & (ec.Dsuz_lineage != 'union')]
log(w.pivot_table(index='category', columns='Dsuz_lineage', values=['MH_OR', 'MH_BH'], sort=False).to_string(
    float_format=lambda v: f'{v:.3g}'))

# ------------------------------------------------------------------ D: per-pool GO enrichment, Dsuz lineage
m2 = pd.read_csv(f'{A}/david_mapping/Dsuz_to_Dmel_FBgn.tsv', sep='\t')
FBZ = dict(zip(m2.focal_mrna, m2.Dmel_FBgn))
by_pool = []
for p in POPS:
    x = mk(p, 15)
    x = x[x.index.isin(FBZ)]
    g = pd.DataFrame({'FBgn': x.index.map(FBZ), 'power': (x.Pn + x.Ps + x.Dn + x.Ds).values,
                      'pos': ((x.P < 0.01) & (x.alpha > 0)).values})
    g = g.groupby('FBgn').agg(power=('power', 'max'), pos=('pos', 'any'))
    for asp, terms in direct.items():
        ga = g[g.index.isin(GENES_ASPECT[asp])]
        dec, pos = deciles(ga.power.to_numpy()), ga.pos.to_numpy()
        posg = set(ga.index[pos])
        for term, genes in terms.items():
            if len(genes & posg) < 2:
                continue
            mem = ga.index.isin(genes)
            mexp, mfold, mp = matched_test(dec, pos, mem)
            by_pool.append(dict(pool=p, aspect=asp, term=term, name=names.get(term), in_category=int(mem.sum()),
                                overlap=int((mem & pos).sum()), matched_expected=mexp, matched_fold=mfold, matched_P=mp))
bp = pd.DataFrame(by_pool)
bp['matched_FDR'] = bp.groupby(['pool', 'aspect']).matched_P.transform(bh)
bp.to_csv(f'{RES}/lineage_go_by_pool.tsv', sep='\t', index=False)
u70 = pd.read_csv(f'{PMK}/summary/L3_lineage_go.tsv', sep='\t')
sig = u70[(u70.version == 'DAF15') & (u70.set.isin(['Dsuz_lineage_union', 'Dsub_lineage'])) & (u70.matched_FDR < 0.05)]
log('\n## D. GO terms with power-matched FDR < 0.05 in script 70 (DAF >= 15%), and the same term per D. suzukii pool')
for _, r in sig.sort_values(['set', 'matched_FDR']).iterrows():
    q = bp[bp.term == r.term].set_index('pool')
    cells = ', '.join(f'{p} {q.loc[p, "overlap"]:.0f}/{q.loc[p, "in_category"]:.0f} fold {q.loc[p, "matched_fold"]:.2f} '
                      f'FDR {q.loc[p, "matched_FDR"]:.3g}' if p in q.index else f'{p} n.t.' for p in POPS)
    sub = u70[(u70.version == 'DAF15') & (u70.set == 'Dsub_lineage') & (u70.term == r.term)]
    sc = f'Dsub fold {sub.matched_fold.iloc[0]:.2f} FDR {sub.matched_FDR.iloc[0]:.3g}' if len(sub) else 'Dsub < 2 genes'
    log(f'  [{r.set}] {r["name"]} ({r.overlap}/{r.in_category}, fold {r.matched_fold:.2f}, FDR {r.matched_FDR:.3g}); '
        f'{sc}; {cells}')
log('\nPer-pool GO terms with matched FDR < 0.05, Dsuz lineage, by pool:')
for p in POPS:
    s = bp[(bp.pool == p) & (bp.matched_FDR < 0.05)].sort_values('matched_FDR')
    log(f'  {p}: ' + '; '.join(f'{n} ({o}/{i}, {fo:.2f}, {fd:.2g})' for n, o, i, fo, fd in
                              zip(s.name, s.overlap, s.in_category, s.matched_fold, s.matched_FDR)))

# ------------------------------------------------------------------ genes of the ecological categories
log('\n## Candidate genes in ecological categories')
for cat in ECO:
    x = cand[cand.eco.str.contains(cat)]
    log(f'  {cat}: ' + ('; '.join(f'{s} [{c}{", robust" if rb else ""}]' for s, c, rb in zip(x.symbol, x.candidate_class, x.robust))
                        if len(x) else 'none'))
LOG.close()
