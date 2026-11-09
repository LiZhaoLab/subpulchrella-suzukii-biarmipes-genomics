#!/usr/bin/env python3
"""Analysis A1 of the depth plan: distribution of fitness effects (DFE) of new
nonsynonymous mutations, alpha and omega_a, per lineage, from the Step-0 site tables.

Sites: 0-fold (selected) and 4-fold (neutral) positions in the ancestral codon (site table of
scripts/97). Alleles: the filter of scripts/98 (>= 3 reads, >= 1 read per strand, >= MIN_FRAC of
depth); ancestral base (D. biarmipes parsimony) and the most frequent passing derived base are
kept; sites with two passing derived bases and no ancestral base are dropped. N_SUB reads are
drawn without replacement (hypergeometric); j = derived reads in the subsample. Unfolded
spectrum j = 0..N_SUB per gene; j = N_SUB counts as a lineage-specific fixed difference.
Model: fastDFE 1.5 BaseInference (polyDFE-type likelihood with nuisance parameters r_i taken
from the neutral spectrum), fitted to POLYMORPHISM ONLY: a deleterious-only DFE cannot fit the
adaptive part of the divergence, so divergence enters only the alpha step below.
  D (primary): discrete DFE, fractions of new 0-fold mutations with S = 4Ne*s in (-1e5, -100],
     (-100, -10], (-10, -1], (-1, 0] (fastDFE DiscreteFractionalParametrization; S uniform within
     a bin); ancestral misidentification eps estimated.
  G: gamma DFE of deleterious mutations (S_d = mean, b = shape), eps estimated, fastDFE default
     bound S_d >= -1e5 (the likelihood rises toward the bound; S_d is not identifiable).
  GB: gamma + exponential beneficial component (p_b, S_b), eps = 0.
alpha and omega_a (Eyre-Walker and Keightley 2009): omega = (D_sel / L_sel) / (D_neut / L_neut)
from the lineage-specific fixed differences (j = N_SUB); omega_na = expected relative fixation rate
of new deleterious mutations under the fitted DFE, E[S / (1 - exp(-S))] (for GB times 1 - p_b);
alpha = 1 - omega_na / omega; omega_a = omega - omega_na.
Strata: autosomes (Muller 2L, 2R, 3), X, 2L, 2R, 3, autosomal genes with highest FlyAtlas2
expression in the testis vs other autosomal genes. Dot chromosome excluded.
Uncertainty: 100 bootstrap replicates over genes (resampled within stratum; seed 20261030).
Sensitivity (autosomes, model D): N_SUB = 10; MIN_FRAC = 0.02.
Validation: 20 Poisson-noise replicates simulated from the fitted D. subpulchrella autosomal
DFE (fastDFE Simulation) and refitted.
Outputs: results/A1_dfe_estimates.tsv, results/A1_dfe_spectra.tsv,
         results/A1_dfe.log
Usage: 102_dfe.py [POOL ...]
"""
import logging
import sys
import time

import fastdfe as fd
import numpy as np
import pandas as pd
from scipy import integrate, stats

logging.getLogger('fastdfe').setLevel(logging.ERROR)
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
POOLS = ['Dsub_L5_7', 'JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
MULLER = {'NC_050609.1': '3', 'NC_050610.1': '2L', 'NC_050611.1': '2R', 'NC_050612.1': '3',
          'NC_050613.1': 'X', 'NW_023665577.1': '3', 'NC_050608.1': '4',
          'NC_092080.1': '2L', 'NC_092081.1': '2R', 'NC_092082.1': '3', 'NC_092083.1': '4',
          'NC_092084.1': 'X'}
MIN_READS, SEED, N_BOOT = 3, 20261030, 100
INTERVALS = np.array([-np.inf, -100, -10, -1, 0])
EDGES = np.array([-1e5, -100, -10, -1, 0])
FIXED = {'D': {'all': {'h': 0.5}},
         'G': {'all': {'h': 0.5, 'p_b': 0, 'S_b': 1}},
         'GB': {'all': {'eps': 0, 'h': 0.5}}}


def _fix(S):
    return 1.0 if abs(S) < 1e-12 else S / (-np.expm1(-S))


def _mean_fix(a, b):
    return integrate.quad(_fix, a, b, limit=400)[0] / (b - a)


MEAN_FIX_BIN = np.array([_mean_fix(EDGES[i], EDGES[i + 1]) for i in range(4)])


def omega_na_gamma(S_d, b):
    sc = abs(S_d) / b
    f = lambda x: stats.gamma.pdf(x, b, scale=sc) * (x / np.expm1(x) if x > 1e-12 else 1.0)
    return sum(integrate.quad(f, lo, hi, limit=500)[0]
               for lo, hi in ((0, 1e-6), (1e-6, 1e-3), (1e-3, 1), (1, 10), (10, 100), (100, np.inf)))


def load_table(pool):
    return np.load(f'{QC}/site_table/{pool}.npz')


def derived_counts(z, pre, min_frac, n_sub, rng, sel=None):
    """(gene or intron index, j) for sites usable at depth n_sub; j = -1 if unusable."""
    tot = z[f'{pre}_tot'].astype(np.int64)
    sb = z[f'{pre}_sb']
    anc = z[f'{pre}_anc'].astype(np.int64)
    if sel is not None:
        tot, sb, anc = tot[sel], sb[sel], anc[sel]
    depth = tot.sum(1)
    strand = ((sb[:, None] >> np.arange(4)) & 1).astype(bool)
    ok = (tot >= MIN_READS) & strand & (tot >= min_frac * depth[:, None])
    rows = np.arange(len(anc))
    anc_ok = ok[rows, anc]
    k_anc = np.where(anc_ok, tot[rows, anc], 0)
    t_der = np.where(ok, tot, 0)
    t_der[rows, anc] = 0
    n_der_alleles = (t_der > 0).sum(1)
    der = np.argmax(t_der, 1)
    k_der = t_der[rows, der]
    usable = (k_anc + k_der >= n_sub) & ~(~anc_ok & (n_der_alleles > 1))
    j = np.full(len(anc), -1)
    j[usable] = rng.hypergeometric(k_der[usable], k_anc[usable], n_sub)
    return j


def gene_spectra(pool, min_frac, n_sub):
    """counts (G, 2, n_sub + 1): class 0 = 0-fold (selected), 1 = 4-fold (neutral); gene table."""
    z = load_table(pool)
    rng = np.random.default_rng(SEED)
    fold = z['cds_fold']
    j = derived_counts(z, 'cds', min_frac, n_sub, rng)
    gene = z['cds_gene']
    G = len(z['genes'])
    counts = np.zeros((G, 2, n_sub + 1), dtype=np.int64)
    for ci, f in ((0, 0), (1, 4)):
        m = (fold == f) & (j >= 0)
        np.add.at(counts, (gene[m], ci, j[m]), 1)
    chroms = z['chroms']
    genes = pd.DataFrame({'tx': z['genes'], 'chrom': chroms[z['gene_chrom']]})
    genes['muller'] = genes.chrom.map(MULLER)
    assert genes.muller.notna().all()
    return counts, genes


def testis_map():
    g = pd.read_csv(f'{RES}/L3_lineage_gene_sets.tsv', sep='\t')
    g = g[g.version == 'DAF15']
    m = {}
    for t, tis in zip(g.transcripts.astype(str), g.FlyAtlas2_max_tissue.astype(str)):
        for tx in t.split(','):
            m[tx.strip()] = tis
    return m


def fit(sel, neut, model='D'):
    kw = dict(model=fd.DiscreteFractionalParametrization(intervals=EDGES)) if model == 'D' else {}
    inf = fd.BaseInference(sfs_neut=fd.Spectrum(neut.tolist()), sfs_sel=fd.Spectrum(sel.tolist()),
                           fixed_params=FIXED[model], do_bootstrap=False, n_runs=10 if model != 'GB' else 20,
                           parallelize=False, seed=SEED, **kw)
    inf.run()
    p = inf.params_mle
    disc = inf.get_discretized(intervals=INTERVALS, confidence_intervals=False)[0]
    omega = (sel[-1] / sel.sum()) / (neut[-1] / neut.sum())
    if model == 'D':
        w_na = float(np.dot(disc, MEAN_FIX_BIN))
    else:
        w_na = omega_na_gamma(p['S_d'], p['b']) * (1 - float(p.get('p_b', 0)))
    return dict(S_d=float(p.get('S_d', np.nan)), b=float(p.get('b', np.nan)), p_b=float(p.get('p_b', 0)),
                S_b=float(p.get('S_b', 0)), eps=float(p.get('eps', 0)), lnl=float(inf.likelihood),
                omega=float(omega), omega_na=w_na, alpha=1 - w_na / omega, omega_a=omega - w_na,
                **{f'frac_{lab}': float(v) for lab, v in zip(['lt-100', '-100_-10', '-10_-1', '-1_0'], disc)})


def pi_from(spec):
    n = len(spec) - 1
    j = np.arange(n + 1)
    return float((spec * 2 * j * (n - j) / (n * (n - 1))).sum() / spec.sum())


def main():
    pools = sys.argv[1:] or POOLS
    tm = testis_map()
    rows, spec_rows = [], []
    log = open(f'{RES}/A1_dfe{"_" + "_".join(pools) if sys.argv[1:] else ""}.log', 'w')

    def say(s):
        print(s, flush=True)
        log.write(s + '\n')

    say(f'fastdfe {fd.__version__}; seed {SEED}; bootstraps {N_BOOT}; started {time.ctime()}')
    for pool in pools:
        for n_sub, min_frac, tag in ((20, 0.01, 'primary'), (10, 0.01, 'n10'), (20, 0.02, 'frac02')):
            counts, genes = gene_spectra(pool, min_frac, n_sub)
            genes['testis'] = genes.tx.map(tm).eq('Male_Testis')
            strata = {'autosomes': genes.muller.isin(['2L', '2R', '3'])}
            if tag == 'primary':
                strata.update({'X': genes.muller.eq('X'), '2L': genes.muller.eq('2L'),
                               '2R': genes.muller.eq('2R'), '3': genes.muller.eq('3'),
                               'auto_testis': strata['autosomes'] & genes.testis,
                               'auto_other': strata['autosomes'] & ~genes.testis})
            for st, mask in strata.items():
                idx = np.flatnonzero(mask.values)
                tot = counts[idx].sum(0)
                sel, neut = tot[0], tot[1]
                spec_rows.append(dict(pool=pool, stratum=st, setting=tag, cls='0fold', genes=len(idx),
                                      **{f'j{k}': int(v) for k, v in enumerate(sel)}))
                spec_rows.append(dict(pool=pool, stratum=st, setting=tag, cls='4fold', genes=len(idx),
                                      **{f'j{k}': int(v) for k, v in enumerate(neut)}))
                models = ['D', 'G', 'GB'] if (tag == 'primary' and st == 'autosomes') else ['D']
                for model in models:
                    est = fit(sel, neut, model)
                    boot = []
                    if (model == 'D' and tag == 'primary') or (model == 'G' and st == 'autosomes'):
                        rng = np.random.default_rng(SEED + 7)
                        for _ in range(N_BOOT):
                            w = np.bincount(rng.integers(0, len(idx), len(idx)), minlength=len(idx))
                            bt = np.tensordot(w, counts[idx], axes=(0, 0))
                            boot.append(fit(bt[0], bt[1], model))
                    row = dict(pool=pool, stratum=st, setting=tag, model=model, n_sub=n_sub,
                               min_frac=min_frac, genes=len(idx), L_sel=int(sel.sum()), L_neut=int(neut.sum()),
                               P_sel=int(sel[1:-1].sum()), P_neut=int(neut[1:-1].sum()),
                               D_sel=int(sel[-1]), D_neut=int(neut[-1]),
                               pi0=pi_from(sel), pi4=pi_from(neut), **est)
                    if boot:
                        bd = pd.DataFrame(boot)
                        for k in ('S_d', 'b', 'eps', 'alpha', 'omega_a', 'omega_na', 'frac_lt-100',
                                  'frac_-100_-10', 'frac_-10_-1', 'frac_-1_0'):
                            row[f'{k}_lo'], row[f'{k}_hi'] = np.percentile(bd[k], [2.5, 97.5])
                        bd.assign(pool=pool, stratum=st, model=model).to_csv(
                            f'{QC}/site_table/A1_boot_{pool}_{st}_{model}.tsv', sep='\t', index=False)
                    rows.append(row)
                    say(f"{pool} {st} {tag} {model}: genes {len(idx)} pi0/pi4 {row['pi0'] / row['pi4']:.4f} "
                        f"S_d {est['S_d']:.0f} b {est['b']:.3f} alpha {est['alpha']:.3f} omega {est['omega']:.4f} "
                        f"omega_a {est['omega_a']:.4f} frac(-1,0] {est['frac_-1_0']:.3f} eps {est['eps']:.4f} "
                        f"p_b {est['p_b']:.4f} lnL {est['lnl']:.1f}")
    out = pd.DataFrame(rows)
    suffix = ('_' + '_'.join(pools)) if sys.argv[1:] else ''
    out.to_csv(f'{RES}/A1_dfe_estimates{suffix}.tsv', sep='\t', index=False)
    pd.DataFrame(spec_rows).to_csv(f'{RES}/A1_dfe_spectra{suffix}.tsv', sep='\t', index=False)
    # validation: simulate from the fitted D. subpulchrella autosomal model-D DFE and refit
    if 'Dsub_L5_7' in pools:
        sp = pd.DataFrame(spec_rows)
        pick = (sp.pool == 'Dsub_L5_7') & (sp.stratum == 'autosomes') & (sp.setting == 'primary')
        neut = sp[pick & (sp.cls == '4fold')].filter(regex=r'^j\d+$').values[0].astype(float)
        sel_obs = sp[pick & (sp.cls == '0fold')].filter(regex=r'^j\d+$').values[0].astype(float)
        inf = fd.BaseInference(sfs_neut=fd.Spectrum(neut.tolist()), sfs_sel=fd.Spectrum(sel_obs.tolist()),
                               model=fd.DiscreteFractionalParametrization(intervals=EDGES),
                               fixed_params={'all': {'h': 0.5, 'eps': 0}}, do_bootstrap=False, n_runs=10,
                               parallelize=False, seed=SEED)
        inf.run()
        true = inf.get_discretized(intervals=INTERVALS, confidence_intervals=False)[0]
        scaled = neut * sel_obs.sum() / neut.sum()     # selected spectrum takes its site count from here
        scaled[-1] = 0
        sim = fd.Simulation(sfs_neut=fd.Spectrum(scaled.tolist()), params=inf.params_mle,
                            model=fd.DiscreteFractionalParametrization(intervals=EDGES))
        expected = np.array(sim.run().data, dtype=float)
        rng = np.random.default_rng(SEED + 11)
        rec = []
        for _ in range(20):
            s_ = rng.poisson(expected)
            s_[0] = int(sel_obs.sum() - s_[1:].sum())
            n_ = rng.poisson(neut)
            n_[0] = int(neut.sum() - n_[1:].sum())
            r2 = fd.BaseInference(sfs_neut=fd.Spectrum(n_.tolist()), sfs_sel=fd.Spectrum(s_.tolist()),
                                  model=fd.DiscreteFractionalParametrization(intervals=EDGES),
                                  fixed_params={'all': {'h': 0.5}}, do_bootstrap=False, n_runs=10,
                                  parallelize=False, seed=SEED)
            r2.run()
            rec.append(list(r2.get_discretized(intervals=INTERVALS, confidence_intervals=False)[0])
                       + [r2.params_mle['eps']])
        rec = np.array(rec)
        say('validation (simulated from the fitted Dsub autosomal DFE with eps = 0, Poisson noise, '
            '20 replicates, refitted with eps free): true fractions ' + ' '.join(f'{v:.4f}' for v in true)
            + '; recovered median ' + ' '.join(f'{v:.4f}' for v in np.median(rec[:, :4], 0))
            + '; 2.5-97.5%: ' + ' '.join(f'{a:.4f}-{b:.4f}' for a, b in zip(np.percentile(rec[:, :4], 2.5, 0),
                                                                          np.percentile(rec[:, :4], 97.5, 0)))
            + f'; eps median {np.median(rec[:, 4]):.4f}')
    say(f'finished {time.ctime()}')


if __name__ == '__main__':
    main()
