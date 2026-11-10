#!/usr/bin/env python3
"""F follow-up: nucleotide diversity around the sweep loci named in the Results text,
as a regional alternative to the genome-wide CLR panels of Supplementary Fig. 4.

pi per 5-kb window with the estimator and site filter of scripts/32 (P98_RM: depth max(10, P05)..P98,
indel mask, RepeatMasker mask; Kofler et al. 2011 correction, minimum allele count 2, all reads),
divided by the mean pi of the same sample and chromosome class (autosomes or X; dot chromosome
excluded) from the genome-wide 20-kb windows of scripts/32. Windows with < 10% callable sites are
not plotted. Callable fraction = callable sites / window length.
The CLR was computed from SNP allele frequencies only (background spectrum has zero weight at
frequency 0), so the number of SNPs per callable site, which pi includes, was not used by the scan.
Peak-versus-flank summary: pi of the windows that overlap the peak interval, divided by the median
pi of the windows more than 100 kb from it (same sample, same panel).
Outputs: results/F_regional_pi.tsv, results/F_regional_pi.log,
         results/Fig_F_sweep_loci.{png,svg,pdf} (all panels, first draft),
         results/Fig_F_sweep_main.{png,svg,pdf} (old panels a, b, h; for a main figure) and
         results/Fig_F_sweep_loci_supp.{png,svg,pdf} (old panels c, d, f, g; supplementary figure).
         Old panel e (Kan, chromosome 3 37.7-39.3 Mb; peaks in a gap with almost no callable sites) is in
         neither split figure (choice of Li Zhao).
"""
import glob
import importlib.util
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
A = '/ru-auth/local/home/asato/scratch'
GFF = {'Dsub': f'{A}/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff',
       'Dsuz': f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'}
spec = importlib.util.spec_from_file_location('cs', f'{PROJ}/scripts/32_chunk_stats.py')
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)

WIN, MIN_CALL, FLANK = 5000, 0.10, 100_000
DOT = {'NC_050608.1', 'NC_092083.1'}
DSUZ = ['JP-Kan', 'JP-Sap2', 'JP-Sap', 'JP-Tok']
SHORT = {'JP-Kan': 'Kan', 'JP-Sap2': 'Sap2', 'JP-Sap': 'Sap', 'JP-Tok': 'Tok', 'Dsub_L5_7': 'D. subpulchrella'}
# loci named in Manuscript_Text_Inversions_Sweeps.md, in the order of the text:
# (key, species, chrom, chrom label, start, end, pools, genes to label)
LOCI = [
    ('X5.5', 'Dsuz', 'NC_092084.1', 'X', 5_100_000, 5_900_000, DSUZ, ['B-H1', 'B-H2', 'ENGase', 'Marf', 'swa']),
    ('3_70.7', 'Dsuz', 'NC_092082.1', '3', 70_320_000, 71_120_000, DSUZ, ['Fancd2', 'Syp', 'Takl1']),
    ('3_83.5', 'Dsuz', 'NC_092082.1', '3', 83_140_000, 83_940_000, DSUZ, ['EndoU', 'glob1']),
    ('3_38.0', 'Dsuz', 'NC_092082.1', '3', 37_700_000, 39_300_000, DSUZ, []),
    ('X28.3', 'Dsub', 'NC_050613.1', 'X', 27_900_000, 28_700_000, ['Dsub_L5_7'], ['CoRest']),
    ('2L7.9', 'Dsub', 'NC_050610.1', '2L', 7_520_000, 8_320_000, ['Dsub_L5_7'], ['Tep2', 'Myo28B1']),
    ('2R11.5', 'Dsub', 'NC_050611.1', '2R', 11_150_000, 11_950_000, ['Dsub_L5_7'], ['Cpr49Ah', 'Dera']),
]
# figure layouts: (output name, cells in panel order, columns, row height in inches, wspace, column width
# ratios, compact). Compact (narrow panels): no interval label above the panel (the x axis shows it), smaller
# max-CLR text.
LAYOUTS = [
    ('Fig_F_sweep_loci', ['density'] + [x[0] for x in LOCI], 2, 2.4, 0.28, None, False),
    ('Fig_F_sweep_main', ['density', 'X5.5', '2R11.5'], 3, 2.5, 0.38, [0.85, 1.075, 1.075], True),
    ('Fig_F_sweep_loci_supp', ['3_70.7', '3_83.5', 'X28.3', '2L7.9'], 2, 2.4, 0.28, None, False),
]

FONTS = '/ru-auth/local/home/lzhao/.local/share/fonts'
for f in ('arial', 'arialbd', 'ariali', 'arialbi'):
    font_manager.fontManager.addfont(f'{FONTS}/{f}.ttf')
assert 'arial' in font_manager.findfont('Arial', fallback_to_default=False).lower()
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e',
                     'ytick.color': '#52514e', 'axes.linewidth': 0.6})
SUB, SUZ, GREY, INK, MUTED = '#2a78d6', '#eb6834', '#9a9893', '#0b0b0b', '#52514e'


def class_means():
    """Mean pi per (sample, class) over all genome-wide 20-kb windows of scripts/32 (P98_RM, no dot)."""
    tot = {}
    for f in glob.glob(f'{QC}/stats/div/*.windows.tsv.gz'):
        w = pd.read_csv(f, sep='\t')
        w = w[(w.config == 'P98_RM') & ~w.chrom.isin(DOT)]
        for (s, c), g in w.groupby(['sample', 'cls']):
            a = tot.setdefault((s, c), [0.0, 0])
            a[0] += g.sum_pi.sum()
            a[1] += g.n_valid.sum()
    return {k: v[0] / v[1] for k, v in tot.items()}


def chunk_for(sample, chrom, start, end):
    for f in glob.glob(f'{QC}/counts/{sample}/{chrom}_*.npz'):
        s, e = (int(x) for x in os.path.basename(f)[:-4].rsplit('_', 2)[1:])
        if s <= start and end <= e:
            return s, e
    raise SystemExit(f'no single chunk holds {sample} {chrom}:{start}-{end}')


def regional_pi(sample, chrom, start, end):
    """Windows of WIN bp over [start, end) with n_valid, S, sum_pi (scripts/32 P98_RM)."""
    sp = 'Dsub' if sample.startswith('Dsub') else 'Dsuz'
    cls = 'X' if chrom in cs.XCHROM else 'autosome'
    n = cs.POOL[sample][1 if cls == 'X' else 0]
    dmin, dmax = cs.thresholds(sample, cls)
    cw, cp = cs.corrections(n, dmax['p98'])
    c0, c1 = chunk_for(sample, chrom, start, end)
    c4, indel = cs.load_chunk(sample, chrom, c0, c1)
    rm = np.load(f'{QC}/masks_sites/{sp}.npz')[f'{chrom}|rm'][c0:c1 + 1]
    i0, i1 = start - c0, end - c0          # 0-based offsets of positions start..end-1
    c4, indel, rm = c4[i0:i1], indel[i0:i1], rm[i0:i1]
    depth = c4.sum(axis=1)
    valid = (depth >= dmin) & (depth <= dmax['p98']) & ~indel & ~rm
    snp, pi, tw = cs.site_stats(c4, n, cw, cp)
    widx = np.arange(i1 - i0) // WIN
    s = cs.sums_by(widx, widx[-1] + 1, valid, snp, pi, tw, depth)
    ws = start + np.arange(widx[-1] + 1) * WIN
    return pd.DataFrame(dict(sample=sample, cls=cls, chrom=chrom, start=ws, end=np.minimum(ws + WIN, end) - 1,
                             n_valid=s['n_valid'], S=s['S'], sum_pi=s['sum_pi']))


def clr(sample, chrom, start, end):
    fs = [f for f in glob.glob(f'{QC}/sweeps/{sample}/clr/{chrom}_*.txt')
          if os.path.basename(f).rsplit('_', 1)[0] == chrom]
    d = pd.concat([pd.read_csv(f, sep='\t') for f in fs]).sort_values('location')
    return d[(d.location >= start) & (d.location < end)]


_MRNA = {}


def genes(sp, chrom, start, end, sym):
    if sp not in _MRNA:
        rows = []
        with open(GFF[sp]) as f:
            for line in f:
                t = line.split('\t')
                if len(t) > 8 and t[2] == 'mRNA':
                    a = dict(kv.split('=', 1) for kv in t[8].strip().split(';') if '=' in kv)
                    rows.append((t[0], a.get('gene', ''), int(t[3]), int(t[4]), a.get('ID', '').replace('rna-', '')))
        _MRNA[sp] = pd.DataFrame(rows, columns=['chrom', 'gene', 'start', 'end', 'tx'])
    m = _MRNA[sp]
    g = m[(m.chrom == chrom) & (m.end >= start) & (m.start < end)].copy()
    g['symbol'] = g.tx.map(sym).fillna('')
    return g.groupby('gene').agg(start=('start', 'min'), end=('end', 'max'),
                                 symbol=('symbol', lambda x: next((v for v in x if v), ''))).reset_index()


def main():
    pk = pd.read_csv(f'{RES}/F_sweep_peaks.tsv', sep='\t')
    gs = pd.read_csv(f'{RES}/L3_lineage_gene_sets.tsv', sep='\t')
    sym = {}
    for t, s in zip(gs.transcripts.astype(str), gs.symbol.astype(str)):
        for tx in t.split(','):
            sym[tx.strip()] = s
    cm = class_means()
    log = ['Regional pi around the sweep loci of the Results text (scripts/126). Window 5 kb; '
           f'pi relative to the sample mean of its chromosome class; class means: '
           + ', '.join(f'{s} {c} {v:.5f}' for (s, c), v in sorted(cm.items()))]
    rows, panels = [], []
    for key, sp, chrom, clab, start, end, pools, label in LOCI:
        tracks = {}
        for p in pools:
            w = regional_pi(p, chrom, start, end)
            assert len(w) and (w.n_valid <= WIN).all()
            w['locus'] = key
            w['callable'] = w.n_valid / (w.end - w.start + 1)
            w['pi'] = np.where(w.callable >= MIN_CALL, w.sum_pi / w.n_valid.clip(lower=1), np.nan)
            w['pi_rel'] = w.pi / cm[(p, w.cls.iloc[0])]
            q = pk[(pk.pool == p) & (pk.chrom == chrom) & (pk.end >= start) & (pk.start < end)]
            top = q.sort_values('max_clr', ascending=False).head(1)
            if len(top):
                ps, pe = int(top.start.iloc[0]), int(top.end.iloc[0])
                inpk = (w.end >= ps - 1000) & (w.start <= pe + 1000)
                far = (w.end < ps - FLANK) | (w.start > pe + FLANK)
                r_pk = w.pi[inpk].mean() / np.nanmedian(w.pi[far])
                c_pk = w.callable[inpk].mean() / np.nanmedian(w.callable[far])
                log.append(f'{key} {p}: strongest peak {chrom}:{ps}-{pe} CLR {top.max_clr.iloc[0]:.1f}; '
                           f'pi in peak / median flank pi = {r_pk:.2f}; callable in peak / flank = {c_pk:.2f}; '
                           f'flank pi relative to class mean = {np.nanmedian(w.pi_rel[far]):.2f}')
            else:
                log.append(f'{key} {p}: no peak above threshold in the plotted interval')
            c = clr(p, chrom, start, end)
            t = pk[(pk.pool == p) & (pk.cls == w.cls.iloc[0])].threshold.iloc[0]
            tracks[p] = (w, c, t, q)
            rows.append(w)
        g = genes(sp, chrom, start, end, sym)
        panels.append((key, sp, chrom, clab, start, end, pools, label, tracks, g))
    pd.concat(rows).to_csv(f'{RES}/F_regional_pi.tsv', sep='\t', index=False)
    open(f'{RES}/F_regional_pi.log', 'w').write('\n'.join(log) + '\n')
    print('\n'.join(log))
    for layout in LAYOUTS:
        figure(panels, *layout)


def figure(panels, name, cells, ncols, row_h, wspace, width_ratios, compact):
    inv = pd.read_csv(f'{RES}/INV_cluster_summary.tsv', sep='\t')
    inv = inv[inv.flag.isna() & (inv.n_replicating_pools >= 2)]
    bykey = {p[0]: p for p in panels}
    assert all(c == 'density' or c in bykey for c in cells), cells
    nrows = -(-len(cells) // ncols)
    fig = plt.figure(figsize=(7.2, row_h * nrows))
    outer = fig.add_gridspec(nrows, ncols, hspace=0.42, wspace=wspace, width_ratios=width_ratios)
    letters = 'abcdefgh'
    for k, cell in enumerate(cells):
        r, c = divmod(k, ncols)
        if cell == 'density':
            ax = fig.add_subplot(outer[r, c])
            density_panel(ax)
            ax.set_title(letters[k], loc='left', fontweight='bold', fontsize=8)
            continue
        key, sp, chrom, clab, start, end, pools, label, tracks, g = bykey[cell]
        sub = outer[r, c].subgridspec(4, 1, height_ratios=[1, 1.15, 0.6, 0.45], hspace=0.22)
        axs = [fig.add_subplot(sub[i]) for i in range(4)]
        col = SUB if sp == 'Dsub' else SUZ
        x0, x1 = start / 1e6, end / 1e6
        for a in axs:
            a.set_xlim(x0, x1)
        # candidate inversion regions (D. suzukii clusters replicated in >= 3 populations)
        if sp == 'Dsuz':
            for z in inv[inv.chrom == chrom].itertuples():
                if z.end >= start and z.start < end:
                    for a in axs[:3]:
                        a.axvspan(max(z.start, start) / 1e6, min(z.end, end) / 1e6, color='#ecebe8', lw=0, zorder=0)
        txt = []
        for p in pools:
            w, cl, t, q = tracks[p]
            axs[0].plot(cl.location / 1e6, cl.LR, lw=0.5, color=col, alpha=0.85)
            axs[0].axhline(t, color=col, lw=0.4, ls='--', alpha=0.6)
            mid = (w.start + w.end) / 2e6
            axs[1].step(mid, w.pi_rel, where='mid', lw=0.6, color=col, alpha=0.85)
            axs[2].step(mid, w.callable, where='mid', lw=0.5, color=col, alpha=0.85)
            if sp == 'Dsuz':
                txt.append(f'{SHORT[p]} {cl.LR.max():.0f}')
        axs[0].set_ylim(0, 1.45 * max(tracks[p][1].LR.max() for p in pools))
        axs[0].set_ylabel('CLR', fontsize=6.5)
        if txt:
            axs[0].text(0.01, 0.97, 'max CLR: ' + ', '.join(txt), transform=axs[0].transAxes, ha='left',
                        va='top', fontsize=5.0 if compact else 5.5, color=MUTED)
        axs[1].axhline(1, color=GREY, lw=0.5, ls=':')
        axs[1].set_ylim(0, min(2.5, max(1.6, np.nanmax([tracks[p][0].pi_rel.max() for p in pools]) * 1.05)))
        axs[1].set_ylabel('Relative\nπ', fontsize=6.5)
        axs[2].set_ylim(0, 1.05)
        axs[2].set_yticks([0, 1])
        axs[2].set_ylabel('Callable', fontsize=6)
        # genes: alternate two rows to limit overlap; label only the genes named in the text
        gg = g.sort_values('start').reset_index(drop=True)
        for i, z in gg.iterrows():
            y = i % 2
            axs[3].plot([z.start / 1e6, z.end / 1e6], [y, y], lw=2.2, color=MUTED, solid_capstyle='butt')
        # one label per group of named genes whose centres lie within 20% of the plotted span
        named = gg[gg.symbol.isin(label)].drop_duplicates('symbol')
        named = named.assign(mid=(named.start + named.end) / 2e6).sort_values('mid')
        groups = []
        for z in named.itertuples():
            if groups and z.mid - groups[-1][-1][1] < 0.20 * (x1 - x0):
                groups[-1].append((z.symbol, z.mid))
            else:
                groups.append([(z.symbol, z.mid)])
        for grp in groups:
            axs[3].text(np.mean([m for _, m in grp]), 1.45, ', '.join(s for s, _ in grp), fontsize=5.5,
                        fontstyle='italic', ha='center', va='bottom', color=INK)
        axs[3].set_ylim(-0.6, 3.0)
        axs[3].set_yticks([])
        axs[3].spines['left'].set_visible(False)
        for a in axs[:3]:
            a.tick_params(labelbottom=False, labelsize=6)
        axs[3].tick_params(labelsize=6)
        fig.align_ylabels(axs[:3])
        axs[3].set_xlabel(f'Chromosome {clab} (Mb)' if sp == 'Dsuz' else f'{clab} (Mb)', fontsize=6.5)
        species = 'D. suzukii' if sp == 'Dsuz' else 'D. subpulchrella'
        axs[0].set_title(letters[k], loc='left', fontweight='bold', fontsize=8)
        axs[0].text(0.08 if compact else 0.06, 1.12, species, transform=axs[0].transAxes, fontstyle='italic',
                    fontsize=7, va='bottom')
        if not compact:
            axs[0].text(1.0, 1.12, f'{clab}: {x0:.2f}–{x1:.2f} Mb', transform=axs[0].transAxes, fontsize=6.5,
                        ha='right', va='bottom', color=MUTED)
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(f'{RES}/{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {name}')


def density_panel(ax):
    """Panel c of the earlier figure (scripts/108 fig_f), unchanged: unrounded per-Mb densities."""
    ar = pd.read_csv(f'{RES}/F_sweep_arms.tsv', sep='\t')
    arms = ['2L', '2R', '3L', '3R', '3R_scaffold']
    for i, pool in enumerate(['Dsub_L5_7'] + DSUZ):
        q = ar[ar.pool == pool].set_index('arm')
        ref = q.loc[['2L', '2R', '3L', '3R'], 'peaks_per_mb'].mean()
        ax.plot(np.arange(5) - 0.2 + 0.1 * i, q.loc[arms, 'peaks_per_mb'] / ref, 'o', ms=3.5,
                color=SUB if i == 0 else SUZ,
                label='D. subpulchrella' if i == 0 else ('D. suzukii (4 populations)' if i == 1 else None))
    ax.axhline(1, color=GREY, lw=0.6, ls=':')
    ax.set_ylim(0.2, 2.85)
    ax.set_xticks(np.arange(5), ['2L', '2R', '3L', '3R', '3R\nscaffold'], fontsize=6.5)
    ax.tick_params(labelsize=6)
    ax.set_ylabel('Peaks per Mb, relative to the\nsample mean of 2L, 2R, 3L, 3R', fontsize=6.5)
    ax.legend(frameon=False, prop={'style': 'italic', 'size': 6}, loc='upper left')


if __name__ == '__main__':
    main()
