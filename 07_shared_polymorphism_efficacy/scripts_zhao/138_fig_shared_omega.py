#!/usr/bin/env python3
"""Script 138: figure for the shared-polymorphism paragraph of the manuscript.

Question shown: lineage-specific fixed differences whose derived base also segregates in the other
species ("shared" sites, script 113) - is their 0-fold / four-fold ratio per site like that of
polymorphism (piN/piS) or like that of substitutions (omega of all fixed differences)?

Quantities (autosomes = Muller 2L, 2R, 3; primary settings of A1 and script 113; sharing threshold 5%):
  pi          piN/piS: pi at 0-fold sites / pi at four-fold sites of the focal pool (A1 spectra, 20 reads
              drawn without replacement, seed 20261030; scripts/102)
  shared      (shared 0-fold fixed differences / 0-fold sites) / (shared four-fold fixed differences /
              four-fold sites)  (= script 113 omega_shared)
  all         (0-fold fixed differences / 0-fold sites) / (four-fold fixed differences / four-fold sites)
              (= A1 omega = script 113 omega_obs)
  not_shared  omega after removal of the shared sites (= script 113 omega_removed; tabulated, not plotted)
  shared_over_pi, shared_over_all  ratios of the above (tabulated, not plotted)
Lineages: D. subpulchrella (Dsub_L5_7; shared = derived base present in JP-Kan, JP-Sap or JP-Sap2) and
D. suzukii (JP-Kan, JP-Sap, JP-Sap2; shared = derived base present in Dsub_L5_7). JP-Tok is excluded
because it contains D. subpulchrella DNA (report section 10.2).
95% intervals: percentile intervals from 100 bootstrap replicates over genes, with the same gene weights
as A1 and script 113 (seed 20261030 + 7). All quantities of one pool use the same replicates, so the
ratios are paired within a pool.
Every point estimate is checked by assertion against results/alpha_shared_polymorphism.tsv and
the A1 estimates (results/A1_dfe_estimates_*.tsv).

Added: panels b-d for the efficacy paragraph of the manuscript ("Despite similar
alpha estimates, purifying selection was less effective in D. subpulchrella ..."), per pool (Dsub_L5_7 and
the four D. suzukii pools, JP-Tok included as in that paragraph):
  b  pN/pS = (Pn / Ln) / (Ps / Ls), polymorphisms with derived allele frequency (DAF) >= 15%
  d  dN/dS = (Dn / Ln) / (Ds / Ls), lineage-specific fixed differences (no multiple-hit correction)
     b and d: polarized MK after exclusion of fixed differences shared with the other species
     (results/polarized_mk_excl, scripts 60, 114), genes tested in Table 7, Nei-Gojobori site counts scaled
     to callable codons (scripts/96, PMK_VARIANT=excl; its functions are reused, not copied). 95% intervals:
     1,000 bootstrap replicates over genes (seed 20261030, one generator per pool).
  c  piN/piS = pi at 0-fold / pi at four-fold sites, autosomes, orthologous genes present in all five
     tables, 95% interval from the 1-Mb block bootstrap (2,000 replicates) of analysis A3
     (results/A3_pnps_arms.tsv, scripts/104). Not recomputed.
  The values are checked by assertion against results/polarized_mk_excl/summary/mk_per_site_rates.tsv
  and against the numbers printed in the manuscript paragraph.
Outputs: results/shared_omega_ci.tsv, results/shared_omega_ci.log,
         results/efficacy_per_site_ci.tsv, results/efficacy_per_site_ci.log,
         results/Fig_shared_omega.{pdf,svg,png}          (panel a alone, horizontal grid lines)
         results/Fig_shared_omega_nogrid.{pdf,svg,png}   (panel a alone, no grid lines)
         results/Fig_shared_omega_efficacy.{pdf,svg,png} (panels a-d in a 2 x 2 layout, 180 mm wide, no grid lines)
         results/Fig_shared_omega_efficacy_merged.{pdf,svg,png} (panel a + one panel b holding pN/pS,
           piN/piS and dN/dS on the 0-0.25 axis of panel a; above each group, the range over the four D. suzukii
           pools of 100 x (D. subpulchrella / D. suzukii pool - 1), one decimal, computed from unrounded values)
Usage:   138_fig_shared_omega.py compute      (loads the four site tables; run on a compute node)
         138_fig_shared_omega.py compute_mk   (reads the reference genomes; run on a compute node)
         138_fig_shared_omega.py plot
Python:  envs/popgen.
"""
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{P}/results'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
TSV = f'{RES}/shared_omega_ci.tsv'
FIG = f'{RES}/Fig_shared_omega'
SUB = 'Dsub_L5_7'
POOLS = ['JP-Kan', 'JP-Sap', 'JP-Sap2']
THRESH = 0.05
DSUZ4 = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
TSV_MK = f'{RES}/efficacy_per_site_ci.tsv'
SEED, N_BOOT_MK = 20261030, 1000
AUTO = ['2L', '2R', '3']
N_SUB = 20


def omega(L0, D0, L4, D4):
    return (D0 / L0) / (D4 / L4)


def pi_ratio(sp):
    """sp (..., 2, N_SUB + 1): class 0 = 0-fold, 1 = four-fold; same estimator as scripts/102 pi_from."""
    j = np.arange(N_SUB + 1)
    h = 2 * j * (N_SUB - j) / (N_SUB * (N_SUB - 1))
    pi = (sp * h).sum(-1) / sp.sum(-1)
    return pi[..., 0] / pi[..., 1]


def ci(b):
    return np.percentile(b, 2.5), np.percentile(b, 97.5)


def compute():
    spec = importlib.util.spec_from_file_location('s113', f'{P}/scripts/113_alpha_shared_polymorphism.py')
    s113 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s113)
    assert s113.SUB == SUB and s113.POOLS == POOLS and s113.N_SUB == N_SUB and s113.MIN_FRAC == 0.01
    log = open(f'{RES}/shared_omega_ci.log', 'w')

    def say(s):
        print(s, flush=True)
        log.write(s + '\n')

    say(f'started {time.ctime()}; numpy {np.__version__}; pandas {pd.__version__}; seed {s113.SEED} '
        f'(+7 for gene weights); bootstraps {s113.N_BOOT}; threshold {THRESH}; strata autosomes {AUTO}')
    t113 = pd.read_csv(f'{RES}/alpha_shared_polymorphism.tsv', sep='\t')
    t113 = t113[(t113.thresh == THRESH) & (t113.stratum == 'autosomes')].set_index(['pair', 'focal'])
    est = pd.concat([pd.read_csv(f'{RES}/A1_dfe_estimates_{p}.tsv', sep='\t') for p in [SUB] + POOLS])
    est = est[(est.setting == 'primary') & (est.model == 'D') & (est.stratum == 'autosomes')].set_index('pool')

    tabs, W, pis = {}, {}, {}
    for pool in [SUB] + POOLS:
        z = s113.load(pool)
        tabs[pool] = z
        idx = np.flatnonzero(z['gene_df'].muller.isin(AUTO).values)
        assert len(idx) == est.loc[pool, 'genes'], (pool, len(idx), est.loc[pool, 'genes'])
        w = s113.weights(len(idx))
        W[pool] = (idx, w)
        # per-gene site-frequency spectra exactly as scripts/102 gene_spectra (same j, same seed)
        j, fold = z['j'], z['cds_fold']
        counts = np.zeros((len(z['genes']), 2, N_SUB + 1))
        for ci_, f in ((0, 0), (1, 4)):
            m = (fold == f) & (j >= 0)
            np.add.at(counts, (z['cds_gene'][m], ci_, j[m]), 1)
        sp = counts[idx]
        p0 = pi_ratio(sp.sum(0))
        assert abs(p0 - est.loc[pool, 'pi0'] / est.loc[pool, 'pi4']) < 1e-12
        pb = pi_ratio(np.einsum('rg,gck->rck', w, sp))
        pis[pool] = (p0, pb)
        say(f'{pool}: autosomal genes {len(idx)}; piN/piS {p0:.4f} ({ci(pb)[0]:.4f}-{ci(pb)[1]:.4f}); '
            f'equals A1 pi0/pi4')

    m = np.load(f'{QC}/site_table/{SUB}.dsuz_map.npz')
    rows = []
    for pool in POOLS:
        zs, zp = tabs[SUB], tabs[pool]
        say(f'== pair {SUB} / {pool}')
        i_s, i_p = s113.join(zs, m, zp, say)
        for focal, zf, zo, i_f, i_o in ((SUB, zs, zp, i_s, i_p), (pool, zp, zs, i_p, i_s)):
            d, _ = s113.direction(zf, zo, i_f, i_o, THRESH)
            gs = s113.gene_sums(d, len(zf['genes']))
            idx, w = W[focal]

            def q(t):
                L0, D0, S0 = t[..., 0, 0], t[..., 0, 1], t[..., 0, 3]
                L4, D4, S4 = t[..., 1, 0], t[..., 1, 1], t[..., 1, 3]
                return {'all': omega(L0, D0, L4, D4), 'shared': (S0 / L0) / (S4 / L4),
                        'not_shared': omega(L0 - S0, D0 - S0, L4 - S4, D4 - S4)}

            t = gs[idx].sum(0)
            pt, bt = q(t), q(np.einsum('rg,gck->rck', w, gs[idx]))
            pt['pi'], bt['pi'] = pis[focal]
            pt['shared_over_pi'], bt['shared_over_pi'] = pt['shared'] / pt['pi'], bt['shared'] / bt['pi']
            pt['shared_over_all'], bt['shared_over_all'] = pt['shared'] / pt['all'], bt['shared'] / bt['all']
            r = t113.loc[(pool, focal)]
            for k, col in (('all', 'omega_obs'), ('shared', 'omega_shared'), ('not_shared', 'omega_removed'),
                           ('pi', 'pi0_over_pi4_focal')):
                assert abs(pt[k] - r[col]) < 1e-9, (pool, focal, k, pt[k], r[col])
            (L0, D0, _, S0, _), (L4, D4, _, S4, _) = t[0], t[1]
            assert (int(S0), int(S4), int(D0), int(D4)) == (r.S0, r.S4, r.D0, r.D4)
            for k in ('pi', 'shared', 'all', 'not_shared', 'shared_over_pi', 'shared_over_all'):
                lo, hi = ci(bt[k])
                rows.append(dict(lineage='D. subpulchrella' if focal == SUB else 'D. suzukii', focal=focal,
                                 other=pool if focal == SUB else SUB, quantity=k, value=pt[k], lo=lo, hi=hi,
                                 genes=len(idx), L0=int(L0), L4=int(L4), D0=int(D0), D4=int(D4), S0=int(S0),
                                 S4=int(S4)))
            say(f'  focal {focal}: ' + '; '.join(
                f'{k} {pt[k]:.4f} ({ci(bt[k])[0]:.4f}-{ci(bt[k])[1]:.4f})'
                for k in ('pi', 'shared', 'all', 'not_shared', 'shared_over_pi', 'shared_over_all'))
                + '; point estimates equal script 113')
    pd.DataFrame(rows).to_csv(TSV, sep='\t', index=False)
    say(f'wrote {TSV}; finished {time.ctime()}')


def compute_mk():
    os.environ['PMK_VARIANT'] = 'excl'
    src = open(f'{P}/scripts/96_mk_per_site_rates.py').read()
    head = src[:src.index('\nrows = []\n')]          # constants and functions of script 96, no outputs
    s96 = {'__name__': 's96'}
    exec(compile(head, f'{P}/scripts/96_mk_per_site_rates.py', 'exec'), s96)
    assert s96['R'].endswith('polarized_mk_excl') and s96['POPS'] == DSUZ4
    log = open(f'{RES}/efficacy_per_site_ci.log', 'w')

    def say(x):
        print(x, flush=True)
        log.write(x + '\n')

    say(f'started {time.ctime()}; numpy {np.__version__}; pandas {pd.__version__}; source {s96["R"]}; '
        f'seed {SEED}; bootstraps {N_BOOT_MK} over genes')
    ref = pd.read_csv(f'{s96["OUTD"]}/mk_per_site_rates.tsv', sep='\t')
    ref = ref[ref.set == 'tested (Table 7)'].set_index('focal')
    sites, rows = {}, []
    for focal in [SUB] + DSUZ4:
        sp = 'Dsub' if focal.startswith('Dsub') else 'Dsuz'
        g = pd.read_csv(f'{s96["R"]}/{focal}.genes.tsv', sep='\t').set_index('gene')
        g = g[g.callable > 0]
        miss = [x for x in g.index if sp not in sites or x not in sites[sp].index]
        if miss:
            sites[sp] = pd.concat([sites.get(sp), s96['gene_sites'](sp, miss)])
        st = sites[sp].reindex(g.index)
        assert st.notna().all().all(), focal
        frac = g.callable / g.codons
        g['Ls'], g['Ln'] = st.Ls_ref * frac, st.Ln_ref * frac
        t = pd.read_csv(f'{s96["R"]}/{focal}.mk_tests.DAF15.tsv', sep='\t').set_index('gene')
        M = np.column_stack([t.Pn, t.Ps, t.Dn, t.Ds, g.loc[t.index, 'Ln'], g.loc[t.index, 'Ls']]).astype(float)

        def q(m):
            Pn, Ps, Dn, Ds, Ln, Ls = (m[..., k] for k in range(6))
            return {'pNpS': (Pn / Ln) / (Ps / Ls), 'dNdS': (Dn / Ln) / (Ds / Ls)}

        pt = q(M.sum(0))
        assert round(pt['pNpS'], 4) == ref.loc[focal, 'pNpS'] and round(pt['dNdS'], 4) == ref.loc[focal, 'KaKs']
        assert len(t) == ref.loc[focal, 'genes']
        rng = np.random.default_rng(SEED)
        n = len(t)
        W = np.array([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(N_BOOT_MK)])
        bt = q(W @ M)
        for k in ('pNpS', 'dNdS'):
            lo, hi = np.percentile(bt[k], [2.5, 97.5])
            rows.append(dict(pool=focal, quantity=k, value=pt[k], lo=lo, hi=hi, genes=n,
                             interval='gene bootstrap, 1,000 replicates'))
        say(f'{focal}: tested genes {n}; pN/pS {pt["pNpS"]:.4f} ({np.percentile(bt["pNpS"], 2.5):.4f}-'
            f'{np.percentile(bt["pNpS"], 97.5):.4f}); dN/dS {pt["dNdS"]:.4f} ({np.percentile(bt["dNdS"], 2.5):.4f}-'
            f'{np.percentile(bt["dNdS"], 97.5):.4f}); equal script 96 (excl)')
    a3 = pd.read_csv(f'{RES}/A3_pnps_arms.tsv', sep='\t')
    a3 = a3[a3.group == 'autosomes'].set_index('pool')
    for p in [SUB] + DSUZ4:
        rows.append(dict(pool=p, quantity='piNpiS', value=a3.loc[p, 'pnps'], lo=a3.loc[p, 'lo'],
                         hi=a3.loc[p, 'hi'], genes=int(a3.loc[p, 'genes']),
                         interval='A3 1-Mb block bootstrap, 2,000 replicates'))
    d = pd.DataFrame(rows)
    v = d.pivot(index='pool', columns='quantity', values='value')
    # numbers printed in the manuscript paragraph
    assert round(v.loc[SUB, 'pNpS'], 3) == 0.066 and round(v.loc[SUB, 'dNdS'], 3) == 0.229
    assert (round(v.loc[DSUZ4, 'pNpS'].min(), 3), round(v.loc[DSUZ4, 'pNpS'].max(), 3)) == (0.057, 0.058)
    assert (round(v.loc[DSUZ4, 'dNdS'].min(), 3), round(v.loc[DSUZ4, 'dNdS'].max(), 3)) == (0.193, 0.202)
    r = d[(d.pool == SUB) & (d.quantity == 'piNpiS')].iloc[0]
    assert (round(r.value, 4), round(r.lo, 4), round(r.hi, 4)) == (0.0785, 0.0754, 0.0818)
    assert (round(v.loc[DSUZ4, 'piNpiS'].min(), 4), round(v.loc[DSUZ4, 'piNpiS'].max(), 4)) == (0.0655, 0.068)
    # The manuscript percentages were computed (script 134) from values rounded to 4 decimals; reproduce that,
    # and log the unrounded range. dN/dS JP-Sap: unrounded 18.50% (rounds to 19%), from 4-decimal values 18.497%.
    for k, rng_txt in (('pNpS', (14, 16)), ('dNdS', (13, 18)), ('piNpiS', (15, 20))):
        rr = v.loc[SUB, k] / v.loc[DSUZ4, k]
        r4 = v.loc[SUB, k].round(4) / v.loc[DSUZ4, k].round(4)
        assert (round(100 * (r4.min() - 1)), round(100 * (r4.max() - 1))) == rng_txt, (k, r4)
        say(f'{k}: D. subpulchrella / D. suzukii, unrounded {100 * (rr.min() - 1):.2f}-{100 * (rr.max() - 1):.2f}% '
            f'higher; from 4-decimal values {100 * (r4.min() - 1):.2f}-{100 * (r4.max() - 1):.2f}% '
            f'(paragraph: {rng_txt[0]}-{rng_txt[1]}%)')
    d.to_csv(TSV_MK, sep='\t', index=False)
    say(f'all values equal the manuscript paragraph; wrote {TSV_MK}; finished {time.ctime()}')


def setup_mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fonts = '/ru-auth/local/home/lzhao/.local/share/fonts'
    for f in ('arial', 'arialbd', 'ariali', 'arialbi'):
        font_manager.fontManager.addfont(f'{fonts}/{f}.ttf')
    assert 'arial' in font_manager.findfont('Arial', fallback_to_default=False).lower()
    # style and colours as in scripts/118 (Fig. 4)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e',
                         'ytick.color': '#52514e', 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                         'ytick.major.width': 0.6, 'legend.frameon': False,
                         'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic'})
    return plt


COL = {'D. subpulchrella': '#2a78d6', 'D. suzukii': '#eb6834'}
DOT = dict(fmt='o', ms=3.2, mew=0, elinewidth=0.8, capsize=1.6, capthick=0.8, zorder=3)
PI_TXT = r'$\mathrm{\pi}_\mathrm{N}/\mathrm{\pi}_\mathrm{S}$'


def draw_shared(ax, grid):
    from matplotlib.lines import Line2D
    d = pd.read_csv(TSV, sep='\t')
    cats = [('pi', 'Polymorphism\n(' + PI_TXT + ')'), ('shared', 'Shared fixed\ndifferences'),
            ('all', 'All fixed\ndifferences')]
    for ci_, (q, _) in enumerate(cats):
        for li, lin in enumerate(COL):
            s = d[(d.lineage == lin) & (d.quantity == q)]
            if lin == 'D. subpulchrella':
                s = s.set_index('other').loc[POOLS]
                if q in ('pi', 'all'):      # one D. subpulchrella pool: identical across the three pairs
                    assert s.value.nunique() == 1 and s.lo.nunique() == 1
                    s = s.iloc[:1]
            else:
                s = s.set_index('focal').loc[POOLS]
            x = ci_ + (-0.17 if li == 0 else 0.17) + (np.arange(len(s)) - (len(s) - 1) / 2) * 0.085
            v, lo, hi = s.value.to_numpy(), s.lo.to_numpy(), s.hi.to_numpy()
            ax.errorbar(x, v, yerr=[v - lo, hi - v], color=COL[lin], **DOT)
    ax.set_xticks(range(len(cats)), [lab for _, lab in cats])
    ax.tick_params(axis='x', length=0, pad=4)
    ax.set_xlim(-0.55, len(cats) - 0.45)
    ax.set_ylim(0, 0.25)
    ax.set_yticks(np.arange(0, 0.26, 0.05))
    if grid:
        ax.yaxis.grid(True, color='#e4e3df', lw=0.5)
        ax.set_axisbelow(True)
    ax.set_ylabel('Nonsynonymous / synonymous ratio\n(0-fold / four-fold, per site)')
    handles = [Line2D([], [], marker='o', ls='', ms=3.2, color=COL[k], mew=0) for k in COL]
    ax.legend(handles, [r'$\it{D.\ subpulchrella}$ lineage', r'$\it{D.\ suzukii}$ lineage'],
              loc='upper left', fontsize=6.5, handletextpad=0.2, borderaxespad=0.2)


def draw_pool(ax, d, q, ylab, ymax, step):
    s = d[d.quantity == q].set_index('pool').loc[[SUB] + DSUZ4]
    short = ['Dsub', 'Kan', 'Sap', 'Sap2', 'Tok']
    for i, (p, r) in enumerate(s.iterrows()):
        c = COL['D. subpulchrella'] if p == SUB else COL['D. suzukii']
        ax.errorbar([i], [r.value], yerr=[[r.value - r.lo], [r.hi - r.value]], color=c, **DOT)
    ax.set_xticks(range(5), short)
    for lab, p in zip(ax.get_xticklabels(), s.index):
        lab.set_color(COL['D. subpulchrella'] if p == SUB else COL['D. suzukii'])
    ax.tick_params(axis='x', length=0, pad=3)
    ax.set_xlim(-0.6, 4.6)
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax + 1e-9, step))
    ax.set_ylabel(ylab)


def draw_merged(ax, d):
    cats = [('pNpS', 'pN/pS\n(DAF ≥ 15%)'), ('piNpiS', PI_TXT), ('dNdS', 'dN/dS\n(lineage-specific)')]
    for ci_, (q, _) in enumerate(cats):
        s = d[d.quantity == q].set_index('pool')
        r = s.loc[SUB]
        ax.errorbar([ci_ - 0.2], [r.value], yerr=[[r.value - r.lo], [r.hi - r.value]], color=COL['D. subpulchrella'],
                    **DOT)
        z = s.loc[DSUZ4]
        x = ci_ + 0.17 + (np.arange(4) - 1.5) * 0.085
        v, lo, hi = z.value.to_numpy(), z.lo.to_numpy(), z.hi.to_numpy()
        ax.errorbar(x, v, yerr=[v - lo, hi - v], color=COL['D. suzukii'], **DOT)
        pct = 100 * (r.value / v - 1)
        top = max(r.hi, hi.max())
        ax.text(ci_, top + 0.012, f'+{pct.min():.1f}–{pct.max():.1f}%', ha='center', va='bottom', fontsize=6.5,
                color='#0b0b0b')
    ax.set_xticks(range(len(cats)), [lab for _, lab in cats])
    ax.tick_params(axis='x', length=0, pad=4)
    ax.set_xlim(-0.55, len(cats) - 0.45)
    ax.set_ylim(0, 0.25)
    ax.set_yticks(np.arange(0, 0.26, 0.05))
    ax.set_ylabel('Nonsynonymous / synonymous ratio\n(per site)')


def plot():
    plt = setup_mpl()
    outs = []
    for grid, name in ((True, FIG), (False, f'{FIG}_nogrid')):
        fig, ax = plt.subplots(figsize=(3.46, 2.45))
        draw_shared(ax, grid)
        fig.tight_layout(pad=0.4)
        for ext in ('pdf', 'svg', 'png'):
            fig.savefig(f'{name}.{ext}', dpi=300)
        plt.close(fig)
        outs.append(name)
    d = pd.read_csv(TSV_MK, sep='\t')
    fig = plt.figure(figsize=(7.087, 4.6))
    gs = fig.add_gridspec(2, 2, left=0.085, right=0.985, bottom=0.085, top=0.955, wspace=0.3, hspace=0.42)
    axs = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(2)]
    draw_shared(axs[0], grid=False)
    draw_pool(axs[1], d, 'pNpS', 'pN/pS (polymorphism, DAF ≥ 15%)', 0.10, 0.02)
    draw_pool(axs[2], d, 'piNpiS', PI_TXT + ' (nucleotide diversity)', 0.10, 0.02)
    draw_pool(axs[3], d, 'dNdS', 'dN/dS (lineage-specific divergence)', 0.25, 0.05)
    for ax, lab in zip(axs, 'abcd'):
        ax.text(-0.2, 1.07, lab, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
    name = f'{FIG}_efficacy'
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(f'{name}.{ext}', dpi=300)
    plt.close(fig)
    outs.append(name)
    fig, axs = plt.subplots(1, 2, figsize=(7.087, 2.55))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.19, top=0.93, wspace=0.3)
    draw_shared(axs[0], grid=False)
    draw_merged(axs[1], d)
    for ax, lab in zip(axs, 'ab'):
        ax.text(-0.2, 1.07, lab, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
    name = f'{FIG}_efficacy_merged'
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(f'{name}.{ext}', dpi=300)
    outs.append(name)
    print('wrote ' + ', '.join(f'{o}.{{pdf,svg,png}}' for o in outs))


if __name__ == '__main__':
    {'compute': compute, 'compute_mk': compute_mk, 'plot': plot}[sys.argv[1]]()
