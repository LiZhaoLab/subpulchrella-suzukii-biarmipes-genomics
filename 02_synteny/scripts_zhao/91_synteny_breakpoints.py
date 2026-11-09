#!/usr/bin/env python3
"""Script 91: quantitative summary of the four-species synteny (Results section
"Chromosome-scale genome synteny in the suzukii species subgroup"), from Airi Sato's PROmer run
synteny_4species/reanalysis_20260824 (one-to-one blocks, >= 500 bp, >= 50% amino-acid identity;
the blocks of the draft). No new alignment.

A. Conservation of chromosome arms (Muller elements): share of aligned bp between major chromosomes
   that lies on the homologous pair, per comparison.
B. Collinear segments and breakpoints for each homologous pair. Blocks are sorted along the reference;
   a block joins the current segment when it has the same orientation, continues the query in the
   expected direction (back-step <= 20 kb) and both gaps are <= GAP. Segments shorter than MIN_SEG
   (reference span) or with fewer than 3 blocks are removed (local noise, paralogs) and the remaining
   segments are chained again with the same rule until nothing changes. Breakpoints = adjacent final
   segments along the reference. Two parameter sets show the sensitivity.
C. Lineage of the D. suzukii-D. subpulchrella breakpoints (outgroup D. biarmipes), in D. subpulchrella
   coordinates: a breakpoint that is also a D. subpulchrella-D. biarmipes breakpoint (within TOL) is a
   change on the D. subpulchrella lineage; one inside a D. subpulchrella-D. biarmipes collinear segment
   (>= TOL on both sides) is a change on the D. suzukii lineage; otherwise unassigned. Parsimony with
   one outgroup; reuse of breakpoints and errors in the assemblies are not modelled.
D. The central region of D. suzukii chromosome 3 (NC_092082.1): gene and TE density, synteny coverage
   (Airi's tables), and the share of the gene-poor 1-Mb blocks of the diversity analysis (R1 of the
   Additional Methods document) that lies in it.
Output: results/synteny_summary.log (stdout) and results/synteny_segments.tsv.
Python: anaconda3/envs/jupyter-ml.
"""
import glob

import numpy as np
import pandas as pd

SYN = '/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824/results'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
OUT = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
COLS = ['S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', 'IDY', 'SIM', 'STP', 'LENR', 'LENQ', 'COVR', 'COVQ',
        'FRM1', 'FRM2', 'R', 'Q']
ARM = {  # Muller element labels (D. melanogaster arm names)
    'NC_092080.1': '2L', 'NC_092081.1': '2R', 'NC_092082.1': '3', 'NC_092083.1': '4', 'NC_092084.1': 'X',
    'NC_050610.1': '2L', 'NC_050611.1': '2R', 'NC_050612.1': '3L', 'NC_050609.1': '3R',
    'NW_023665577.1': '3R(NW_023665577.1)', 'NC_050608.1': '4', 'NC_050613.1': 'X',
    'NC_066612.1': '2L', 'NC_066615.1': '2R', 'NC_066613.1': '3L', 'NC_066616.1': '3R', 'NC_066614.1': '4',
    'NC_066611.1': 'X',
    'NT_033779.5': '2L', 'NT_033778.4': '2R', 'NT_037436.4': '3L', 'NT_033777.3': '3R', 'NC_004353.4': '4',
    'NC_004354.4': 'X'}
HOMOLOG = lambda a, b: a == b or (a == '3' and b.startswith('3')) or (b == '3' and a.startswith('3')) or \
    (a.startswith('3R') and b.startswith('3R'))
COMPS = ['Dsuz_vs_Dsub', 'Dsub_vs_Dbia', 'Dbia_vs_Dmel']


def load(comp):
    t = pd.read_csv(f'{SYN}/primary_one_to_one/{comp}.promer.one_to_one.L500.I50.coords.tsv', sep='\t',
                    header=None, names=COLS)
    t['rs'], t['re'] = t[['S1', 'E1']].min(1), t[['S1', 'E1']].max(1)
    t['qs'], t['qe'] = t[['S2', 'E2']].min(1), t[['S2', 'E2']].max(1)
    t['ori'] = np.sign((t.E1 - t.S1) * (t.E2 - t.S2)).astype(int)
    t['ra'], t['qa'] = t.R.map(ARM), t.Q.map(ARM)
    return t


def chain(seg, gap, back=20_000):
    """seg: DataFrame with rs, re, qs, qe, ori, n (sorted by rs). Returns merged segments."""
    out = []
    cur = None
    for s in seg.itertuples(index=False):
        if cur is not None and s.ori == cur['ori'] and s.rs - cur['re'] <= gap:
            if s.ori > 0:
                ok = s.qs >= cur['qe'] - back and s.qs - cur['qe'] <= gap
            else:
                ok = s.qe <= cur['qs'] + back and cur['qs'] - s.qe <= gap
            if ok:
                cur['re'] = max(cur['re'], s.re)
                cur['qs'], cur['qe'] = min(cur['qs'], s.qs), max(cur['qe'], s.qe)
                cur['n'] += s.n
                continue
        if cur is not None:
            out.append(cur)
        cur = dict(rs=s.rs, re=s.re, qs=s.qs, qe=s.qe, ori=s.ori, n=s.n)
    if cur is not None:
        out.append(cur)
    return pd.DataFrame(out)


def segments(blocks, gap, min_seg):
    s = blocks[['rs', 're', 'qs', 'qe', 'ori']].assign(n=1).sort_values('rs')
    s = chain(s, gap)
    while True:
        keep = s[((s.re - s.rs) >= min_seg) & (s.n >= 3)].sort_values('rs')
        new = chain(keep, gap * 5)
        if len(new) == len(keep):
            return new
        s = new


# ------------------------------------------------------------------ A. chromosome-arm conservation
print('== A. aligned bp between major chromosomes lying on the homologous pair')
data = {c: load(c) for c in COMPS}
for comp, t in data.items():
    m = t[t.ra.notna() & t.qa.notna()]
    hom = m[[HOMOLOG(a, b) for a, b in zip(m.ra, m.qa)]]
    print(f'{comp}: {hom.LEN1.sum() / m.LEN1.sum() * 100:.2f}% of {m.LEN1.sum() / 1e6:.1f} Mb '
          f'({len(hom)} of {len(m)} blocks); unplaced-sequence blocks excluded')

# ------------------------------------------------------------------ B. collinear segments, breakpoints
PARAMS = {'main': (100_000, 50_000), 'strict': (100_000, 200_000)}
seg_rows = []
summary = []
for pname, (gap, min_seg) in PARAMS.items():
    for comp, t in data.items():
        m = t[t.ra.notna() & t.qa.notna()]
        m = m[[HOMOLOG(a, b) for a, b in zip(m.ra, m.qa)]]
        for (r, q), b in m.groupby(['R', 'Q']):
            if b.LEN1.sum() < 1_000_000:
                continue
            s = segments(b, gap, min_seg)
            span = (s.re - s.rs).sum()
            summary.append(dict(params=pname, comp=comp, ref=r, query=q, ref_arm=ARM[r], query_arm=ARM[q],
                                blocks=len(b), aligned_Mb=round(b.LEN1.sum() / 1e6, 2), segments=len(s),
                                breakpoints=len(s) - 1, seg_span_Mb=round(span / 1e6, 1),
                                seg_N50_Mb=round(float(np.sort(s.re - s.rs)[::-1][np.searchsorted(
                                    np.cumsum(np.sort(s.re - s.rs)[::-1]), span / 2)] / 1e6), 2),
                                orientation_changes=int((s.ori.diff().fillna(0) != 0).sum())))
            for x in s.itertuples(index=False):
                seg_rows.append(dict(params=pname, comp=comp, ref=r, query=q, **x._asdict()))
S = pd.DataFrame(summary)
pd.set_option('display.width', 250)
print('\n== B. collinear segments (main: gap 100 kb, segments >= 50 kb; strict: segments >= 200 kb)')
print(S.to_string(index=False))
print('\nTotals per comparison (breakpoints; breakpoints per 10 Mb of collinear span):')
for (p, c), g in S.groupby(['params', 'comp'], sort=False):
    print(f'  {p:6s} {c}: {g.breakpoints.sum()} breakpoints over {g.seg_span_Mb.sum():.0f} Mb '
          f'({10 * g.breakpoints.sum() / g.seg_span_Mb.sum():.1f} per 10 Mb); by reference arm: ' +
          ', '.join(f'{a} {x.breakpoints.sum()} ({10 * x.breakpoints.sum() / x.seg_span_Mb.sum():.1f}/10 Mb)'
                    for a, x in g.groupby('query_arm' if c == 'Dsuz_vs_Dsub' else 'ref_arm')))
SEG = pd.DataFrame(seg_rows)
SEG.to_csv(f'{OUT}/synteny_segments.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ C. polarize Dsuz-Dsub breakpoints
TOL = 100_000
print(f'\n== C. lineage of D. suzukii-D. subpulchrella breakpoints (outgroup D. biarmipes; tolerance {TOL // 1000} kb)')
for pname in PARAMS:
    a = SEG[(SEG.params == pname) & (SEG.comp == 'Dsuz_vs_Dsub')]
    o = SEG[(SEG.params == pname) & (SEG.comp == 'Dsub_vs_Dbia')]
    res = []
    for q, g in a.groupby('query'):                     # breakpoints along each D. subpulchrella sequence
        g = g.sort_values('qs')
        bps = [((x1.qe + x2.qs) / 2) for x1, x2 in zip(g.iloc[:-1].itertuples(), g.iloc[1:].itertuples())]
        og = o[o.ref == q].sort_values('rs')
        obps = np.array([(y1.re + y2.rs) / 2 for y1, y2 in zip(og.iloc[:-1].itertuples(), og.iloc[1:].itertuples())])
        for bp in bps:
            if len(obps) and np.min(np.abs(obps - bp)) <= TOL:
                lin = 'D. subpulchrella'
            elif ((og.rs <= bp - TOL) & (og.re >= bp + TOL)).any():
                lin = 'D. suzukii'
            else:
                lin = 'unassigned'
            res.append(dict(query=q, arm=ARM[q], pos=bp, lineage=lin))
    r = pd.DataFrame(res)
    print(f'{pname}: {len(r)} breakpoints;', r.lineage.value_counts().to_dict())
    print(r.groupby(['arm', 'lineage']).size().unstack(fill_value=0).to_string())

# ------------------------------------------------------------------ D. central region of D. suzukii chromosome 3
print('\n== D. D. suzukii chromosome 3 (NC_092082.1)')
cov = pd.read_csv(f'{SYN}/NC092082_fusion/NC_092082.1.arm_synteny_coverage.tsv', sep='\t')
gd = pd.read_csv(f'{SYN}/NC092082_fusion/NC_092082.1.region_gene_density.tsv', sep='\t')
print(gd.to_string(index=False))
print(cov[['region', 'arm', 'covered_percent']].to_string(index=False))
te = pd.read_csv(f'{SYN}/TE_gene_tracks_200kb/Dsuz_Dsub.TE_gene_density.200kb.tsv', sep='\t')
z = te[(te.species == 'Dsuz') & (te.chromosome == 'NC_092082.1')].copy()
z['region'] = np.where(z.window_end <= 27e6, '3L_like', np.where(z.window_start_0based >= 62e6, '3R_like', 'central'))
print('TE coverage (Earl Grey, 200-kb windows):',
      {k: f'{100 * x.TE_covered_bp.sum() / x.window_bp.sum():.1f}% ({x.TE_covered_bp.sum() / 1e6:.1f} Mb)'
       for k, x in z.groupby('region')})
tot = te.groupby('species').TE_covered_bp.sum() / 1e6
print('TE-covered Mb on the major chromosomes:', tot.round(1).to_dict())
b = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{QC}/stats/div/JP-*.blocks.tsv.gz'))
b = b[(b.cls == 'autosome') & (b.config == 'P99_noRM')]
key = ['sample', 'chrom', 'block']
d = b[b.site_class == 'all'].set_index(key)[['n_valid']].join(
    b[b.site_class == 'fourfold'].set_index(key)[['n_valid']].rename(columns={'n_valid': 'ff'})).reset_index()
d = d[d.n_valid > 0]
d['poor'] = d.ff < 5000
d['central'] = (d.chrom == 'NC_092082.1') & (d.block * 1e6 >= 27e6) & (d.block * 1e6 < 62e6)
for s_, g in d.groupby('sample'):
    p = g[g.poor]
    print(f'{s_}: gene-poor callable {p.n_valid.sum() / 1e6:.1f} Mb, of which in the chr-3 central region '
          f'{p[p.central].n_valid.sum() / 1e6:.1f} Mb ({100 * p[p.central].n_valid.sum() / p.n_valid.sum():.0f}%)')

# ------------------------------------------------------------------ E. checks (added the same day)
# E1. Do breakpoints lie near assembly gaps (runs of >= 10 N) of the assembly in which they are scored?
# E2. Binomial test of the lineage counts (unassigned excluded). E3. Segments of the 2L comparison.
import gzip
import re
from scipy.stats import binom_test

FNA = {'Dsuz': '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna',
       'Dsub': '/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna'}


def gaps(path, keep):
    out, name, seq = {}, None, []
    op = gzip.open if path.endswith('.gz') else open
    for line in op(path, 'rt'):
        if line.startswith('>'):
            if name in keep:
                out[name] = [m.start() for m in re.finditer(r'[Nn]{10,}', ''.join(seq))]
            name, seq = line[1:].split()[0], []
        elif name in keep:
            seq.append(line.strip())
    if name in keep:
        out[name] = [m.start() for m in re.finditer(r'[Nn]{10,}', ''.join(seq))]
    return out


a = SEG[(SEG.params == 'main') & (SEG.comp == 'Dsuz_vs_Dsub')]
o = SEG[(SEG.params == 'main') & (SEG.comp == 'Dsub_vs_Dbia')]
G = {'Dsuz': gaps(FNA['Dsuz'], set(a.ref)), 'Dsub': gaps(FNA['Dsub'], set(a['query']))}
print('\n== E1. assembly gaps (>= 10 N) on the compared sequences:',
      {sp: sum(len(v) for v in g.values()) for sp, g in G.items()})
rows = []
for q, g in a.groupby('query'):
    g = g.sort_values('qs')
    og = o[o.ref == q].sort_values('rs')
    obps = np.array([(y1.re + y2.rs) / 2 for y1, y2 in zip(og.iloc[:-1].itertuples(), og.iloc[1:].itertuples())])
    for x1, x2 in zip(g.iloc[:-1].itertuples(), g.iloc[1:].itertuples()):
        bp_q = (x1.qe + x2.qs) / 2
        # the same junction in D. suzukii coordinates: the facing ends of the two segments
        ends1 = (x1.rs, x1.re)
        ends2 = (x2.rs, x2.re)
        bp_r = min(((e1 + e2) / 2, abs(e1 - e2)) for e1 in ends1 for e2 in ends2)[0]
        if len(obps) and np.min(np.abs(obps - bp_q)) <= TOL:
            lin = 'Dsub'
        elif ((og.rs <= bp_q - TOL) & (og.re >= bp_q + TOL)).any():
            lin = 'Dsuz'
        else:
            lin = 'unassigned'
        near_sub = any(abs(p - bp_q) <= TOL for p in G['Dsub'].get(q, []))
        near_suz = any(abs(p - bp_r) <= TOL for p in G['Dsuz'].get(x1.ref, []))
        rows.append(dict(arm=ARM[q], lineage=lin, near_gap_Dsub=near_sub, near_gap_Dsuz=near_suz))
E = pd.DataFrame(rows)
print(E.groupby('lineage')[['near_gap_Dsub', 'near_gap_Dsuz']].agg(['sum', 'count']).to_string())
for pname in PARAMS:
    a2 = SEG[(SEG.params == pname) & (SEG.comp == 'Dsuz_vs_Dsub')]
for pname, (ns, nz) in (('main', (33, 19)), ('strict', (25, 14))):
    print(f'E2. {pname}: D. subpulchrella {ns} vs D. suzukii {nz}: two-sided binomial P = {binom_test(ns, ns + nz):.3f}')
print('\n== E3. collinear segments, D. suzukii 2L vs D. subpulchrella 2L (main parameters)')
print(SEG[(SEG.params == 'main') & (SEG.ref == 'NC_092080.1')][['rs', 're', 'qs', 'qe', 'ori', 'n']].to_string(index=False))
