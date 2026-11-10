#!/usr/bin/env python3
"""Analysis F of the depth plan: neutral null distribution of the SweepFinder2 CLR for
one pool and chromosome class.

1. Demographic fit. The folded read spectrum of all scanned sites (QC/sweeps/{POOL}/summary.npz,
   'all' sites of the class; 20 reads per site) is fitted with a constant-size model and a
   two-epoch model (moments 1.6, Demographics1D; parameters nu = present/ancestral size, T in
   2*N_anc generations). The expected read spectrum is the moments spectrum for the pool's
   haploid sample size C, transformed by a Monte Carlo matrix that reproduces the data pipeline:
   site depth drawn from the pool's depth histogram (callable range of scripts/98, >= 20 reads),
   allele reads binomial from the pool frequency, strands binomial(1/2), allele filter of
   scripts/98 (>= 3 reads, >= 1 read per strand, >= 1% of depth), hypergeometric subsample of
   20 reads. Poisson likelihood over folded classes 1..10; theta per site from the total.
2. Simulation. N_REP segments of 1.2 Mb under the better model (msprime 1.4, Hudson coalescent,
   binary mutation model, N_anc = 10,000 diploids, mutation and recombination rates scaled to
   keep theta per site; recombination / mutation rate ratio R_RATIO = 2, chosen low because less
   recombination gives larger neutral CLR values, so the null is conservative). Each base is callable with the
   pool's callable fraction; each variant goes through the same pool sampling and filter.
3. SweepFinder2 (-lu, grid every 2 kb over the central 1 Mb) on each segment, with the folded
   background spectrum estimated (-f) from all simulated segments.
Outputs: QC/sweeps/{POOL}/null_{CLS}/ (fit.json, expected/observed spectra, clr_{rep}.txt).
Usage: 100_sweep_null.py POOL CLS     Seed 20261030.
"""
import json
import os
import subprocess
import sys

import moments
import msprime
import numpy as np
import pandas as pd
from scipy.optimize import minimize

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
SF2 = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/envs/popgen/bin/SweepFinder2'
POOL = {'Dsub_L5_7': (94, 94), 'JP-Kan': (424, 241), 'JP-Sap': (100, 54),
        'JP-Sap2': (302, 175), 'JP-Tok': (100, 90)}
XCHROM = {'NC_050613.1', 'NC_092084.1'}
N_SUB, MIN_READS, MIN_FRAC = 20, 3, 0.01
N_MC, N_REP, SEG, CORE0, CORE1, GRID = 4000, int(os.environ.get('N_REP', 40)), 1_200_000, 100_000, 1_100_000, 2000
NA_SIM, R_RATIO, SEED = 10_000, 2.0, 20261030
NU_MAX = float(os.environ.get('NU_MAX', 0))


def pipeline_matrix(C, depth_p, depths, rng):
    """T[k, x]: probability that a site with k of C pool chromosomes carrying allele 1 gives folded
    read count x (0..10) after depth sampling, allele filter and subsampling."""
    T = np.zeros((C + 1, N_SUB // 2 + 1))
    T[0, 0] = T[C, 0] = 1
    for k in range(1, C):
        d = rng.choice(depths, size=N_MC, p=depth_p)
        r = rng.binomial(d, k / C)
        rs = rng.binomial(r, 0.5)
        o = d - r
        os_ = rng.binomial(o, 0.5)
        p1 = (r >= MIN_READS) & (rs > 0) & (rs < r) & (r >= MIN_FRAC * d)
        p2 = (o >= MIN_READS) & (os_ > 0) & (os_ < o) & (o >= MIN_FRAC * d)
        kept = np.where(p1 & p2, d, np.where(p1, r, np.where(p2, o, 0)))
        ok = kept >= N_SUB
        x = np.zeros(N_MC, dtype=int)
        both = ok & p1 & p2
        x[both] = rng.hypergeometric(r[both], o[both], N_SUB)
        f = np.minimum(x, N_SUB - x)
        f[~both] = 0
        f = f[ok]
        T[k] = np.bincount(f, minlength=N_SUB // 2 + 1) / max(len(f), 1)
    return T


def expected(params, C, T):
    nu, t = params
    fs = moments.Demographics1D.two_epoch((nu, t), [C]) if nu is not None else \
        moments.Demographics1D.snm([C])
    p = np.array(fs.data)          # per unit theta, entries 0..C (0 and C masked -> 0 here)
    p[0] = p[C] = 0
    return p @ T                   # expected folded read classes per unit theta


def fit(obs, C, T):
    """Poisson fit; obs = folded classes 0..10 (0 = monomorphic). Returns dict."""
    S, L = obs[1:].sum(), obs.sum()

    def nll(e):
        th = S / (L * e[1:].sum())
        mu = L * th * e[1:]
        return float((mu - obs[1:] * np.log(mu)).sum()), th

    e0 = expected((None, None), C, T)
    ll0, th0 = nll(e0)
    best = None
    for start in ([0.0, -1.0], [1.0, -2.0], [-1.0, -1.0], [2.0, -0.5]):
        # NU_MAX (used for JP-Kan X only): bound nu and T at NU_MAX, because the unbounded
        # fit for that class did not finish within 45 min (moments is slow at extreme nu and T)
        bnd = [(-5, np.log(NU_MAX)), (-7, np.log(NU_MAX))] if NU_MAX else None
        r = minimize(lambda z: nll(expected((np.exp(z[0]), np.exp(z[1])), C, T))[0], start,
                     method='Nelder-Mead', bounds=bnd, options=dict(xatol=1e-3, fatol=1e-3, maxiter=400))
        if best is None or r.fun < best.fun:
            best = r
    nu, t = np.exp(best.x)
    e1 = expected((nu, t), C, T)
    _, th1 = nll(e1)
    return dict(S=int(S), L=int(L), snm_nll=ll0, snm_theta=th0, two_epoch_nll=float(best.fun),
                nu=float(nu), T=float(t), theta=float(th1),
                obs=obs.tolist(), exp_snm=(L * th0 * e0).tolist(), exp_two_epoch=(L * th1 * e1).tolist())


def main():
    pool, cls = sys.argv[1], sys.argv[2]
    out = f'{QC}/sweeps/{pool}/null_{cls}'
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED + (cls == 'X'))
    C = POOL[pool][1 if cls == 'X' else 0]
    summ = np.load(f'{QC}/sweeps/{pool}/summary.npz')
    chroms = sorted({k.split('|')[0] for k in summ.files})
    chroms = [c for c in chroms if (c in XCHROM) == (cls == 'X')]
    obs = sum(summ[f'{c}|all'] for c in chroms).astype(float)
    length = sum(int(summ[f'{c}|length'][0]) for c in chroms)
    f_call = obs.sum() / length
    thr = pd.read_csv(f'{RES}/depth_thresholds.tsv', sep='\t')
    thr = thr[(thr['sample'] == pool) & (thr.cls == cls)].iloc[0]
    h = np.load(f'{QC}/depth_hist_{pool}_{cls}.npy').astype(float)
    depths = np.arange(len(h))
    keep = (depths >= max(int(thr.min_depth), N_SUB)) & (depths <= int(thr.p98))
    depth_p = np.where(keep, h, 0)
    depth_p /= depth_p.sum()
    T = pipeline_matrix(C, depth_p, depths, rng)
    res = fit(obs, C, T)
    res.update(pool=pool, cls=cls, C=C, callable_fraction=f_call, chroms=chroms, nu_max=NU_MAX,
               r_over_mu=R_RATIO, n_rep=N_REP, seed=SEED, moments=moments.__version__,
               msprime=msprime.__version__)
    use_two = res['two_epoch_nll'] < res['snm_nll'] - 2
    res['model_used'] = 'two_epoch' if use_two else 'constant'
    json.dump(res, open(f'{out}/fit.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if not k.startswith(('obs', 'exp'))}), flush=True)
    theta = res['theta'] if use_two else res['snm_theta']
    mu = theta / (4 * NA_SIM)
    dem = msprime.Demography()
    if use_two:
        dem.add_population(name='p', initial_size=res['nu'] * NA_SIM)
        dem.add_population_parameters_change(time=res['T'] * 2 * NA_SIM, initial_size=NA_SIM)
    else:
        dem.add_population(name='p', initial_size=NA_SIM)
    allf = open(f'{out}/all.sf2.txt', 'w')
    allf.write('position\tx\tn\tfolded\n')
    for rep in range(N_REP):
        ts = msprime.sim_ancestry(samples=[msprime.SampleSet(C, population='p', ploidy=1)],
                                  demography=dem, ploidy=2, sequence_length=SEG,
                                  recombination_rate=R_RATIO * mu, random_seed=SEED + 1000 * rep + 1)
        ts = msprime.sim_mutations(ts, rate=mu, model=msprime.BinaryMutationModel(),
                                   random_seed=SEED + 1000 * rep + 2)
        gm = ts.genotype_matrix()                         # (sites, C)
        k = gm.sum(1)
        pos = ts.tables.sites.position.astype(np.int64) + 1
        sel = (k > 0) & (k < C) & (rng.random(len(k)) < f_call)
        k, pos = k[sel], pos[sel]
        d = rng.choice(depths, size=len(k), p=depth_p)
        r = rng.binomial(d, k / C)
        o = d - r
        rs, os_ = rng.binomial(r, 0.5), rng.binomial(o, 0.5)
        p1 = (r >= MIN_READS) & (rs > 0) & (rs < r) & (r >= MIN_FRAC * d)
        p2 = (o >= MIN_READS) & (os_ > 0) & (os_ < o) & (o >= MIN_FRAC * d)
        both = p1 & p2
        x = np.zeros(len(k), dtype=int)
        x[both] = rng.hypergeometric(r[both], o[both], N_SUB)
        keep = both & (x > 0) & (x < N_SUB)
        rows = list(zip(pos[keep].tolist(), x[keep].tolist()))
        with open(f'{out}/sim_{rep}.sf2.txt', 'w') as f:
            f.write('position\tx\tn\tfolded\n')
            for p, x in rows:
                f.write(f'{p}\t{x}\t{N_SUB}\t1\n')
                allf.write(f'{p}\t{x}\t{N_SUB}\t1\n')
        print(f'rep {rep}: polymorphic sites {len(rows)}', flush=True)
    allf.close()
    subprocess.run([SF2, '-f', f'{out}/all.sf2.txt', f'{out}/spect.txt'], check=True,
                   stdout=subprocess.DEVNULL)
    with open(f'{out}/grid.txt', 'w') as f:
        for g in range(CORE0, CORE1 + 1, GRID):
            f.write(f'{g}\n')
    for rep in range(N_REP):
        subprocess.run([SF2, '-lu', f'{out}/grid.txt', f'{out}/sim_{rep}.sf2.txt', f'{out}/spect.txt',
                        f'{out}/clr_{rep}.txt'], check=True, stdout=subprocess.DEVNULL)
        print(f'SF2 rep {rep} done', flush=True)


if __name__ == '__main__':
    main()
