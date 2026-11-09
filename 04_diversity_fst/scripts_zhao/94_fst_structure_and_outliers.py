#!/usr/bin/env python3
"""Which population or split drives D. suzukii FST, and what lies in the outlier windows.

Input: Airi's pool-seq Hudson FST in 20-kb windows (grenedalf, run 20261012, filtered sync, RepeatMasker +
indel + depth masks, MAF 0.01; the run behind the FST paragraph), her per-pair passed-site minimum (10th
percentile, MAF01_highFST_thresholds.tsv), her per-sample diversity windows of the same run, the
D. suzukii RefSeq GFF and her Dsuz -> Dmel direct-RBH mapping (david_mapping).

Windows: finite FST in all six pairs and passed sites >= the pair's 10th-percentile minimum in all six.
FST values are used directly as distances (all window means < 0.03; -ln(1 - FST) differs by < 2%).

Per window:
  branch_P = mean over the three pairs {Q, R} of the other populations of (F_PQ + F_PR - F_QR) / 2
             (population-specific branch length, PBS-like; one value per population)
  internal branch of split AB|CD = mean(four cross pairs) - (F_AB + F_CD) / 2  (four-point method)
  splits: study = {Kan, Sap2 | Sap, Tok} (Feng et al. 2024, 2019, NovaSeq 151 bp | Olazcuaga et al. 2020,
          2014/2016, HiSeq 125 bp); {Kan, Sap | Sap2, Tok}; {Kan, Tok | Sap, Sap2}.
Tails: top 1% of windows for each branch and for the study split; their position, callable sites, gene
content, direct-RBH fraction, clustering and per-population pi relative to the genome-wide median.

Writes results/fst_structure.log, results/fst_windows_branches.tsv.gz,
results/fst_outlier_genes.tsv.
"""
import os
import re

import numpy as np
import pandas as pd

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(PROJ, 'results')
AP = '/lustre/fs4/zhao_lab/scratch/asato/population_polymorphism/reanalysis_20261012_repeatmasker_p98'
FSTD = f'{AP}/fst_results_sync_filtered_hudson/MAF01'
THR = f'{AP}/results/high_fst_genes_hudson_MAF01/MAF01_highFST_thresholds.tsv'
DIV = f'{AP}/grenedalf_results'
GFF = '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'
MAP = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/david_mapping/Dsuz_to_Dmel_FBgn.tsv'
CHR = {'NC_092080.1': '2L', 'NC_092081.1': '2R', 'NC_092082.1': '3', 'NC_092083.1': '4', 'NC_092084.1': 'X'}
POPS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
SHORT = {'JP-Kan': 'Kan', 'JP-Sap': 'Sap', 'JP-Sap2': 'Sap2', 'JP-Tok': 'Tok'}
PAIRS = ['JP-Kan:JP-Sap', 'JP-Kan:JP-Sap2', 'JP-Kan:JP-Tok', 'JP-Sap:JP-Sap2', 'JP-Sap:JP-Tok', 'JP-Sap2:JP-Tok']
SPLITS = {'study {Kan,Sap2|Sap,Tok}': (('JP-Kan', 'JP-Sap2'), ('JP-Sap', 'JP-Tok')),
          '{Kan,Sap|Sap2,Tok}': (('JP-Kan', 'JP-Sap'), ('JP-Sap2', 'JP-Tok')),
          '{Kan,Tok|Sap,Sap2}': (('JP-Kan', 'JP-Tok'), ('JP-Sap', 'JP-Sap2'))}
END_MB = 2.0           # "near a chromosome end"
C3_CENTRE = (27e6, 62e6)  # TE-rich centre of D. suzukii chromosome 3 (script 91)
TOP = 0.01

out = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    out.append(s)


def pkey(a, b):
    k = f'{a}:{b}'
    return k if k in PAIRS else f'{b}:{a}'


# ---------------- FST windows
thr = pd.read_csv(THR, sep='\t')
frames = []
for cls in ('autosome', 'X'):
    f = [x for x in os.listdir(f'{FSTD}/{cls}') if x.startswith('fst-') and x.endswith('.csv')]
    assert len(f) == 1, f
    d = pd.read_csv(f'{FSTD}/{cls}/{f[0]}')
    d['cls'] = cls
    frames.append(d)
w = pd.concat(frames, ignore_index=True)
w['chr'] = w.chrom.map(CHR)
assert w.chr.notna().all(), w.chrom.unique()
for p in PAIRS:
    w[p] = pd.to_numeric(w[f'{p}.fst'], errors='coerce')

# per-chromosome means over each pair's finite windows (as in the manuscript text)
say('== Mean FST per chromosome, finite windows of each pair (manuscript convention)')
tab = w.groupby('chr')[PAIRS].mean()
tab.columns = [c.replace('JP-', '') for c in tab.columns]
say(tab.round(4).to_string())
sapcols =[c for c in tab.columns if 'Sap' in c.replace('Sap2', '')]
oth = [c for c in tab.columns if c not in sapcols]
for ch, r in tab.iterrows():
    say(f'  {ch}: min over Sap pairs {r[sapcols].min():.4f} vs max over other pairs {r[oth].max():.4f} '
        f'-> all Sap pairs higher: {r[sapcols].min() > r[oth].max()}')

# joint eligibility
ok = np.ones(len(w), bool)
for p in PAIRS:
    for cls in ('autosome', 'X'):
        m = thr[(thr['class'] == cls) & (thr.pair == p)]
        assert len(m) == 1
        sel = (w.cls == cls).values
        ok[sel] &= (w.loc[sel, f'{p}.passed'] >= m.passed_P10_minimum.iloc[0]).values
    ok &= np.isfinite(w[p].values)
w = w[ok & (w.chr != '4')].copy()
say(f'\n== Joint windows (finite in all 6 pairs, passed >= P10 in all 6; dot chromosome excluded): {len(w)}')
say(w.groupby('chr').size().to_string())

# branches
for P in POPS:
    others = [q for q in POPS if q != P]
    vals = []
    for i in range(3):
        for j in range(i + 1, 3):
            Q, R = others[i], others[j]
            vals.append((w[pkey(P, Q)] + w[pkey(P, R)] - w[pkey(Q, R)]) / 2)
    w[f'br_{SHORT[P]}'] = sum(vals) / 3
for name, ((A, B), (C, D)) in SPLITS.items():
    cross = (w[pkey(A, C)] + w[pkey(A, D)] + w[pkey(B, C)] + w[pkey(B, D)]) / 4
    w[f'int_{name}'] = cross - (w[pkey(A, B)] + w[pkey(C, D)]) / 2
BR = [f'br_{SHORT[P]}' for P in POPS]
INT = [f'int_{n}' for n in SPLITS]

say('\n== Branch lengths (FST units): mean / median over windows / fraction of windows > 0')
for cls_lab, sel in (('autosomes', w.cls == 'autosome'), ('X', w.cls == 'X')):
    say(f'  {cls_lab} ({sel.sum()} windows)')
    for c in BR + INT:
        v = w.loc[sel, c]
        say(f'    {c:32s} mean {v.mean():.4f}  median {v.median():.4f}  >0 {np.mean(v > 0):.3f}')

say('\n== Sap branch by chromosome: mean, median, fraction > 0, and Sap minus mean of other branches (median)')
oth_br = [b for b in BR if b != 'br_Sap']
for ch, g in w.groupby('chr'):
    diff = g.br_Sap - g[oth_br].mean(axis=1)
    say(f'  {ch}: mean {g.br_Sap.mean():.4f} median {g.br_Sap.median():.4f} >0 {np.mean(g.br_Sap > 0):.3f}; '
        f'median(Sap - others) {diff.median():.4f}, fraction of windows Sap longest {np.mean(g[BR].idxmax(axis=1) == "br_Sap"):.3f}')

# ---------------- window annotation
chrlen = w.groupby('chrom').end.max()
w['mid'] = (w.start + w.end) / 2
w['dist_end_Mb'] = np.minimum(w.mid, w.chrom.map(chrlen) - w.mid) / 1e6
w['near_end'] = w.dist_end_Mb < END_MB
w['c3_centre'] = (w.chr == '3') & (w.mid >= C3_CENTRE[0]) & (w.mid < C3_CENTRE[1])
w['passed_mean'] = w[[f'{p}.passed' for p in PAIRS]].mean(axis=1)

# per-population pi of the same run
for P in POPS:
    parts = []
    for cls in ('autosome', 'X'):
        f = [x for x in os.listdir(f'{DIV}/{cls}') if x.startswith(f'diversity-{P}_') and x.endswith('.csv')]
        assert len(f) == 1, (P, cls, f)
        d = pd.read_csv(f'{DIV}/{cls}/{f[0]}')
        col = [c for c in d.columns if c.endswith('.theta_pi')]
        assert len(col) == 1
        parts.append(d[['chrom', 'start', col[0]]].rename(columns={col[0]: f'pi_{SHORT[P]}'}))
    w = w.merge(pd.concat(parts), on=['chrom', 'start'], how='left')
PI = [f'pi_{SHORT[P]}' for P in POPS]
assert w[PI].notna().mean().min() > 0.95, w[PI].notna().mean()
for c in PI:  # relative to the chromosome-class median
    w[c + '_rel'] = w[c] / w.groupby('cls')[c].transform('median')

# genes
genes, mrna2gene = [], {}
with open(GFF) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        t = line.rstrip('\n').split('\t')
        if len(t) < 9 or t[0] not in CHR:
            continue
        if t[2] == 'gene' and 'gene_biotype=protein_coding' in t[8]:
            at = dict(x.split('=', 1) for x in t[8].split(';') if '=' in x)
            genes.append((t[0], int(t[3]), int(t[4]), at['ID'], at.get('Name', ''), at.get('description', '')))
        elif t[2] == 'mRNA':
            m1, m2 = re.search(r'ID=rna-([^;]+)', t[8]), re.search(r'Parent=([^;]+)', t[8])
            if m1 and m2:
                mrna2gene[m1.group(1)] = m2.group(1)
G = pd.DataFrame(genes, columns=['chrom', 'gstart', 'gend', 'gene_id', 'gene_name', 'description'])
G['description'] = G.description.str.replace('%2C', ',').str.replace('%3B', ';')
mp = pd.read_csv(MAP, sep='\t')
mp['gene_id'] = mp.focal_mrna.map(mrna2gene)
assert mp.gene_id.notna().mean() > 0.95, mp.gene_id.notna().mean()  # rest are mRNAs on unplaced scaffolds
g2fb = mp.dropna(subset=['gene_id']).groupby('gene_id').agg(Dmel_FBgn=('Dmel_FBgn', 'first'),
                                                            Dmel_symbol=('Dmel_symbol_official', 'first'))
G = G.merge(g2fb, left_on='gene_id', right_index=True, how='left')
say(f'\n== Protein-coding genes on major chromosomes: {len(G)}; with a direct RBH to Dmel: '
    f'{G.Dmel_FBgn.notna().sum()} ({G.Dmel_FBgn.notna().mean():.3f})')

# gene count per window
w['n_genes'] = 0
for ch, gg in G.groupby('chrom'):
    sel = w.chrom == ch
    s, e = w.loc[sel, 'start'].values, w.loc[sel, 'end'].values
    gs, ge = np.sort(gg.gstart.values), np.sort(gg.gend.values)
    # genes overlapping [s, e]: starts <= e minus ends < s
    w.loc[sel, 'n_genes'] = np.searchsorted(gs, e, side='right') - np.searchsorted(ge, s, side='left')


def describe(mask, label):
    t = w[mask]
    runs = 0
    for ch, g in t.groupby('chrom'):
        st = np.sort(g.start.values)
        runs += 1 + int(np.sum(np.diff(st) > 20000))
    genes_in = overlapping_genes(t)
    rbh = genes_in.Dmel_FBgn.notna().mean() if len(genes_in) else np.nan
    pis = ' '.join(f'{SHORT[P]} {t[f"pi_{SHORT[P]}_rel"].median():.2f}' for P in POPS)
    say(f'  {label}: {len(t)} windows ({"; ".join(f"{k} {v}" for k, v in t.chr.value_counts().items())}); '
        f'{runs} runs of adjacent windows; near end {t.near_end.mean():.2f} (all {w.near_end.mean():.2f}); '
        f'chr3 centre {t.c3_centre.mean():.2f} (all {w.c3_centre.mean():.2f})')
    say(f'      passed sites median {t.passed_mean.median():.0f} (all {w.passed_mean.median():.0f}); genes/window '
        f'{t.n_genes.mean():.2f} (all {w.n_genes.mean():.2f}); genes {len(genes_in)}, direct RBH {rbh:.2f} '
        f'(genome {G.Dmel_FBgn.notna().mean():.2f}); median pi relative to genome median: {pis}')
    return genes_in


def overlapping_genes(t):
    rows = []
    for _, r in t.iterrows():
        gg = G[(G.chrom == r.chrom) & (G.gstart <= r.end) & (G.gend >= r.start)]
        if len(gg):
            rows.append(gg.assign(window=f'{r.chrom}:{int(r.start)}-{int(r.end)}'))
    if not rows:
        return G.iloc[0:0].assign(window='')
    return pd.concat(rows).drop_duplicates('gene_id')


say(f'\n== Top {TOP:.0%} windows of each statistic (autosomes and X separately, then pooled)')
gene_rows = []
for c in ['br_Sap', 'br_Tok', 'br_Kan', 'br_Sap2', 'int_study {Kan,Sap2|Sap,Tok}']:
    m = np.zeros(len(w), bool)
    for cls in ('autosome', 'X'):
        sel = (w.cls == cls).values
        m[sel] = (w.loc[sel, c] >= w.loc[sel, c].quantile(1 - TOP)).values
    gi = describe(m, c)
    t = w[m].set_index(w[m].chrom + ':' + w[m].start.astype(int).astype(str) + '-' + w[m].end.astype(int).astype(str))
    gi = gi.assign(statistic=c, window_value=gi.window.map(t[c]),
                   window_FST_mean=gi.window.map(t[PAIRS].mean(axis=1)),
                   **{f'{p}_rel_window': gi.window.map(t[f'{p}_rel']) for p in PI})
    gene_rows.append(gi)

# overlap between the two strongest tail sets
say('\n== Pairwise FST >= 0.05 (joint windows): within-study pairs vs cross-study pairs')
for p in PAIRS:
    kind = 'within' if p in ('JP-Kan:JP-Sap2', 'JP-Sap:JP-Tok') else 'cross'
    say(f'  {p:16s} {kind:6s} {int((w[p] >= 0.05).sum()):4d} windows')

gr = pd.concat(gene_rows, ignore_index=True)
gr['chr'] = gr.chrom.map(CHR)
cols = ['statistic', 'chr', 'chrom', 'gstart', 'gend', 'gene_name', 'description', 'Dmel_FBgn', 'Dmel_symbol',
        'window', 'window_value', 'window_FST_mean'] + [f'{p}_rel_window' for p in PI]
gr[cols].sort_values(['statistic', 'window_value'], ascending=[True, False]).to_csv(
    os.path.join(RES, 'fst_outlier_genes.tsv'), sep='\t', index=False)
keep = ['chrom', 'chr', 'cls', 'start', 'end'] + PAIRS + BR + INT + PI + ['passed_mean', 'n_genes', 'near_end',
                                                                           'c3_centre']
w[keep].to_csv(os.path.join(RES, 'fst_windows_branches.tsv.gz'), sep='\t', index=False)
with open(os.path.join(RES, 'fst_structure.log'), 'w') as fh:
    fh.write('\n'.join(out) + '\n')

# ---------------- added: X/autosome ratio of the Sap excess, dxy in tail windows, para window
# Drift alone gives X/A FST ratio = Ne_A / Ne_X, between 8/9 and 16/9 (= 1.78) for any breeding sex ratio.
# Extra sampling noise from a pool with fewer independent chromosomes than stated (relatives, fewer flies)
# scales with 1/n: with Sap's stated n (autosomes 100, X 54) the ratio is 100/54 = 1.85, and higher if
# the pooled males are brothers (brothers share the maternal X with probability 1/2, autosomes 1/4).
rng = np.random.default_rng(20261029)
w['block'] = w.chrom + ':' + (w.start // 1_000_000).astype(str)
w['sap_excess'] = w.br_Sap - w[oth_br].mean(axis=1)
say('\n== Sap excess branch (Sap minus mean of the other three branches): X / autosome ratio, 1-Mb block bootstrap')
blk = {cls: w[w.cls == cls].groupby('block').agg(s=('sap_excess', 'sum'), n=('sap_excess', 'size'))
       for cls in ('autosome', 'X')}
est = {cls: b.s.sum() / b.n.sum() for cls, b in blk.items()}
boot = []
for _ in range(2000):
    r = {}
    for cls, b in blk.items():
        i = rng.integers(0, len(b), len(b))
        r[cls] = b.s.values[i].sum() / b.n.values[i].sum()
    boot.append(r['X'] / r['autosome'])
boot = np.array(boot)
say(f'  autosome excess {est["autosome"]:.4f}, X excess {est["X"]:.4f}, ratio {est["X"] / est["autosome"]:.2f} '
    f'(95% {np.percentile(boot, 2.5):.2f}-{np.percentile(boot, 97.5):.2f}; {len(blk["autosome"])} / {len(blk["X"])} blocks); '
    f'P(ratio <= 1.78) = {np.mean(boot <= 16 / 9):.4f}')
for b in BR:
    a, x = w.loc[w.cls == 'autosome', b].mean(), w.loc[w.cls == 'X', b].mean()
    say(f'  {b}: autosome {a:.4f}, X {x:.4f}, X/A {x / a:.2f}')
nA, nX = 100, 54
for lab, ex, n in (('autosome', est['autosome'], nA), ('X', est['X'], nX)):
    neff = 1 / (2 * ex + 1 / n)  # excess ~ (1/n_eff - 1/n) / 2
    say(f'  if the excess were pool sampling only: {lab} n_eff = {neff:.0f} of stated {n} ({neff / n:.2f})')

# dxy (pi-between) and pi-within in tail windows, relative to the genome-wide median of each pair
pb, pw = [], []
for cls in ('autosome', 'X'):
    for kind, store in (('pi-between', pb), ('pi-within', pw)):
        f = [x for x in os.listdir(f'{FSTD}/{cls}') if x.startswith(f'{kind}-') and x.endswith('.csv')]
        assert len(f) == 1
        d = pd.read_csv(f'{FSTD}/{cls}/{f[0]}')
        store.append(d[['chrom', 'start'] + [f'{p}.{kind}' for p in PAIRS]])
w = w.merge(pd.concat(pb), on=['chrom', 'start'], how='left').merge(pd.concat(pw), on=['chrom', 'start'], how='left')
say('\n== Tail windows: median dxy (pi-between) and pi-within relative to the genome median (same class), '
    'mean over the pairs listed')
tails = {'br_Sap': [p for p in PAIRS if 'JP-Sap:' in p or p.endswith(':JP-Sap')],
         'int_study {Kan,Sap2|Sap,Tok}': ['JP-Kan:JP-Sap', 'JP-Kan:JP-Tok', 'JP-Sap:JP-Sap2', 'JP-Sap2:JP-Tok'],
         'br_Tok': [p for p in PAIRS if 'Tok' in p]}
for c, prs in tails.items():
    m = np.zeros(len(w), bool)
    for cls in ('autosome', 'X'):
        sel = (w.cls == cls).values
        m[sel] = (w.loc[sel, c] >= w.loc[sel, c].quantile(1 - TOP)).values
    rel_b = np.mean([(w[f'{p}.pi-between'] / w.groupby('cls')[f'{p}.pi-between'].transform('median'))[m].median()
                     for p in prs])
    rel_w = np.mean([(w[f'{p}.pi-within'] / w.groupby('cls')[f'{p}.pi-within'].transform('median'))[m].median()
                     for p in prs])
    say(f'  {c}: dxy {rel_b:.2f}, pi-within {rel_w:.2f} ({m.sum()} windows; pairs {", ".join(x.replace("JP-", "") for x in prs)})')

say('\n== para window(s)')
para = G[G.gene_name == 'para']
for _, r in para.iterrows():
    t = w[(w.chrom == r.chrom) & (w.end >= r.gstart) & (w.start <= r.gend)]
    q = {c: (w.loc[w.cls == 'X', c] < t[c].max()).mean() for c in ['int_study {Kan,Sap2|Sap,Tok}']}
    say(f'  para {r.chrom}:{r.gstart}-{r.gend}; {len(t)} joint windows; highest study-split window at X quantile '
        f'{q["int_study {Kan,Sap2|Sap,Tok}"]:.3f}')
    say(t[['start', 'end'] + PAIRS + ['int_study {Kan,Sap2|Sap,Tok}'] + PI].round(4).to_string(index=False))

with open(os.path.join(RES, 'fst_structure.log'), 'w') as fh:
    fh.write('\n'.join(out) + '\n')
