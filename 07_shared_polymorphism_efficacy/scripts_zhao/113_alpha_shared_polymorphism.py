#!/usr/bin/env python3
"""Direct test of the alpha counting effect, step 2.

Question: is the lower alpha of the D. subpulchrella lineage (A1, element 3) caused by
lineage-specific fixed differences that are in fact ancestral polymorphisms, sorted in
D. subpulchrella and still segregating in D. suzukii?

Sites: coding sites of the Step-0 site tables (scripts/97), D. subpulchrella sites mapped to
D. suzukii positions by scripts/112, joined with each D. suzukii pool table on the D. suzukii
position. A joined site is kept when the two tables agree on the D. suzukii base and on the
D. subpulchrella base, and the position is unique in both tables. The ancestral bases are not
required to agree: when the other pool's majority allele is the focal derived base, the two tables
polarize the site differently, and such sites are shared polymorphisms that must stay in the test
(a first version that required agreement lost 17% of D. subpulchrella four-fold fixed differences).
All A1-usable focal sites are kept, so 'before' equals A1 exactly; fixed differences at sites not
in the join cannot be tested and stay in (conservative: shared sites are under-removed).
Fixed differences: the A1 primary definition (scripts/102: allele filter >= 3 reads, >= 1 read per
strand, >= 1% of depth; 20 reads drawn without replacement, same random seed, so the calls are
identical to A1; j = 20 is a lineage-specific fixed difference). Derived base = the most frequent
passing non-ancestral base, as in A1.
Shared: the focal derived base passes, in the other species' pool, >= 3 reads, >= 1 read per strand
and >= THRESH of depth (primary 5%; sensitivity 1% and 10%). Both directions are tested:
D. subpulchrella-lineage fixed differences against each D. suzukii pool, and each D. suzukii
pool's lineage fixed differences against the D. subpulchrella pool.
Baseline (independent mutation, sequencing error, mapping error): among focal sites with no derived
read (j = 0), the rate at which the other pool carries base b, by fold class and ancestral base;
expected shared = sum of these rates over the fixed differences.
alpha = 1 - omega_na / omega, with omega_na of the A1 primary DFE fit (model D) of the focal pool
and stratum, and omega from fixed differences on the joined sites, before and after excluding
shared sites. 95% intervals: the 100 A1 bootstrap replicates over genes, reproduced with the same
weights (seed 20261030 + 7) and paired with the A1 bootstrap omega_na.
JP-Tok is excluded: it contains about 1.7% D. subpulchrella DNA (report section 10.2), which would
create false sharing at exactly these sites.
Lineage difference in alpha (D. subpulchrella minus D. suzukii pool): interval from pairing
bootstrap replicate k of the two lineages (independent resampling; conservative).
Outputs: results/alpha_shared_polymorphism.tsv (per lineage),
results/alpha_shared_polymorphism_summary.tsv (per pair), results/alpha_shared_polymorphism.log,
QC/site_table/alpha_shared_boot.npz (bootstrap alpha replicates)
Usage: 113_alpha_shared_polymorphism.py   (run with envs/popgen python)
"""
import importlib.util
import time

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = f'{PROJ}/results'
spec = importlib.util.spec_from_file_location('a1', f'{PROJ}/scripts/102_dfe.py')
a1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a1)
SEED, N_BOOT, MULLER, MIN_READS = a1.SEED, a1.N_BOOT, a1.MULLER, a1.MIN_READS
N_SUB, MIN_FRAC = 20, 0.01
SUB = 'Dsub_L5_7'
POOLS = ['JP-Kan', 'JP-Sap', 'JP-Sap2']
THRESH = [0.05, 0.01, 0.10]
STRATA = ['autosomes', '2L', '2R', '3', 'X']


def load(pool):
    z = dict(np.load(f'{QC}/site_table/{pool}.npz'))
    rng = np.random.default_rng(SEED)
    z['j'] = a1.derived_counts(z, 'cds', MIN_FRAC, N_SUB, rng)
    tot = z['cds_tot'].astype(np.int64)
    anc = z['cds_anc'].astype(np.int64)
    depth = tot.sum(1)
    strand = ((z['cds_sb'][:, None] >> np.arange(4)) & 1).astype(bool)
    ok = (tot >= MIN_READS) & strand & (tot >= MIN_FRAC * depth[:, None])
    t_der = np.where(ok, tot, 0)
    t_der[np.arange(len(anc)), anc] = 0
    z['der'] = np.argmax(t_der, 1)
    z['depth'] = depth
    genes = pd.DataFrame({'tx': z['genes'], 'chrom': z['chroms'][z['gene_chrom']]})
    genes['muller'] = genes.chrom.map(MULLER)
    assert genes.muller.notna().all()
    z['gene_df'] = genes
    return z


def passes(z, rows, base, f):
    """other-pool base `base` at `rows` passes >= 3 reads, >= 1 per strand, >= f of depth."""
    t = z['cds_tot'][rows, base].astype(np.int64)
    s = ((z['cds_sb'][rows] >> base) & 1).astype(bool)
    return (t >= MIN_READS) & s & (t >= f * z['depth'][rows])


def uniq(k):
    valid = np.flatnonzero(k >= 0)
    kk = k[valid]
    u, inv, cnt = np.unique(kk, return_inverse=True, return_counts=True)
    keep = cnt[inv] == 1
    return kk[keep], valid[keep], int((cnt > 1).sum())


def join(zs, m, zp, say):
    names = list(m['dsuz_chroms'])
    ks = np.where(m['dsuz_pos'] >= 0, m['dsuz_chrom'].astype(np.int64) * 10 ** 9 + m['dsuz_pos'], -1)
    pc = np.array([names.index(c) for c in zp['chroms']], dtype=np.int64)
    kp = pc[zp['cds_chrom']] * 10 ** 9 + zp['cds_pos'].astype(np.int64)
    ks_u, rs, ds = uniq(ks)
    kp_u, rp, dp = uniq(kp)
    _, a, b = np.intersect1d(ks_u, kp_u, assume_unique=True, return_indices=True)
    i_s, i_p = rs[a], rp[b]
    same_suz = zs['cds_sis'][i_s] == zp['cds_ref'][i_p]
    same_sub = zp['cds_sis'][i_p] == zs['cds_ref'][i_s]
    same_anc = zs['cds_anc'][i_s] == zp['cds_anc'][i_p]
    same_out = zs['cds_out'][i_s] == zp['cds_out'][i_p]
    keep = same_suz & same_sub
    say(f'  join: duplicated positions dropped Dsub {ds}, pool {dp}; joined {len(i_s)} '
        f'({100 * len(i_s) / len(zs["cds_gene"]):.1f}% of Dsub sites, {100 * len(i_p) / len(zp["cds_gene"]):.1f}% '
        f'of pool sites); agree Dsuz base {same_suz.mean():.5f}, Dsub base {same_sub.mean():.5f}, '
        f'ancestral {same_anc.mean():.5f}, D. biarmipes {same_out.mean():.5f}; kept (both bases agree) '
        f'{int(keep.sum())}')
    return i_s[keep], i_p[keep]


def direction(zf, zo, i_f, i_o, f):
    """per-site arrays over all A1-usable focal sites (so that 'before' equals A1 exactly); the
    sharing test applies to fixed differences at joined sites (testable), the others stay in."""
    j, fold = zf['j'], zf['cds_fold']
    use = (j >= 0) & np.isin(fold, (0, 4))
    other = np.full(len(j), -1, dtype=np.int64)
    other[i_f] = i_o
    anc = zf['cds_anc'].astype(np.int64)
    der = zf['der']
    D = use & (j == N_SUB)
    T = D & (other >= 0)
    S = np.zeros(len(j), bool)
    S[T] = passes(zo, other[T], der[T], f)
    # baseline rates from joined focal sites without a derived read
    rate = np.zeros((2, 4, 4))
    mono = use & (j == 0) & (other >= 0)
    for fi, fv in enumerate((0, 4)):
        for a in range(4):
            sel = mono & (fold == fv) & (anc == a)
            for b in range(4):
                if b != a and sel.any():
                    rate[fi, a, b] = passes(zo, other[sel], np.full(int(sel.sum()), b), f).mean()
    fidx = (fold == 4).astype(int)
    E = np.zeros(len(j))
    E[T] = rate[fidx[T], anc[T], der[T]]
    return dict(gene=zf['cds_gene'][use], cls=fidx[use], D=D[use], T=T[use], S=S[use], E=E[use]), rate


def gene_sums(d, G):
    """(G, 2, 5): per gene and class (0 = 0-fold, 1 = 4-fold): L, D, T (testable D), S (shared), E."""
    out = np.zeros((G, 2, 5))
    for k, v in enumerate((np.ones(len(d['D'])), d['D'], d['T'], d['S'], d['E'])):
        np.add.at(out, (d['gene'], d['cls'], k), v)
    return out


def omega(L0, D0, L4, D4):
    return (D0 / L0) / (D4 / L4)


def full_gene_sums(z):
    """A1 per-gene L and D on all usable sites (to reproduce A1 omega)."""
    use = (z['j'] >= 0) & np.isin(z['cds_fold'], (0, 4))
    G = len(z['genes'])
    out = np.zeros((G, 2, 2))
    g, c = z['cds_gene'][use], (z['cds_fold'][use] == 4).astype(int)
    np.add.at(out, (g, c, 0), 1)
    np.add.at(out, (g, c, 1), (z['j'][use] == N_SUB).astype(float))
    return out


def weights(n):
    rng = np.random.default_rng(SEED + 7)
    return np.array([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(N_BOOT)])


def main():
    log = open(f'{RES}/alpha_shared_polymorphism.log', 'w')

    def say(s):
        print(s, flush=True)
        log.write(s + '\n')

    say(f'started {time.ctime()}; seed {SEED}; bootstraps {N_BOOT}; thresholds {THRESH} (first = primary)')
    est = pd.concat([pd.read_csv(f'{RES}/A1_dfe_estimates_{p}.tsv', sep='\t') for p in [SUB] + POOLS])
    est = est[(est.setting == 'primary') & (est.model == 'D')].set_index(['pool', 'stratum'])
    zs = load(SUB)
    m = np.load(f'{QC}/site_table/{SUB}.dsuz_map.npz')
    tabs = {SUB: zs}
    # reproduction check of A1 omega (point estimate and bootstrap replicates)
    boot_wna, W = {}, {}
    for pool in [SUB] + POOLS:
        z = tabs[pool] if pool in tabs else load(pool)
        tabs[pool] = z
        fs = full_gene_sums(z)
        gdf = z['gene_df']
        for st in STRATA:
            mask = gdf.muller.isin(['2L', '2R', '3']) if st == 'autosomes' else gdf.muller.eq(st)
            idx = np.flatnonzero(mask.values)
            t = fs[idx].sum(0)
            w0 = omega(t[0, 0], t[0, 1], t[1, 0], t[1, 1])
            W[pool, st] = (idx, weights(len(idx)))
            bt = np.einsum('rg,gck->rck', W[pool, st][1], fs[idx])
            wb = omega(bt[:, 0, 0], bt[:, 0, 1], bt[:, 1, 0], bt[:, 1, 1])
            b = pd.read_csv(f'{QC}/site_table/A1_boot_{pool}_{st}_D.tsv', sep='\t')
            assert len(b) == N_BOOT
            boot_wna[pool, st] = b.omega_na.values
            dmax = np.abs(wb - b.omega.values).max()
            say(f'check {pool} {st}: omega {w0:.6f} vs A1 {est.loc[(pool, st), "omega"]:.6f}; '
                f'bootstrap omega max abs difference {dmax:.2e}')
            assert abs(w0 - est.loc[(pool, st), 'omega']) < 1e-9 and dmax < 1e-9
    rows, boots = [], {}
    for pool in POOLS:
        zp = tabs[pool]
        say(f'== pair {SUB} / {pool}')
        i_s, i_p = join(zs, m, zp, say)
        for f in THRESH:
            for focal, zf, zo, i_f, i_o in ((SUB, zs, zp, i_s, i_p), (pool, zp, zs, i_p, i_s)):
                d, rate = direction(zf, zo, i_f, i_o, f)
                gs = gene_sums(d, len(zf['genes']))
                for st in STRATA:
                    idx, w = W[focal, st]
                    t = gs[idx].sum(0)
                    (L0, D0, T0, S0, E0), (L4, D4, T4, S4, E4) = t[0], t[1]
                    w_obs = omega(L0, D0, L4, D4)
                    assert abs(w_obs - est.loc[(focal, st), 'omega']) < 1e-9   # 'before' = A1
                    w_rem = omega(L0 - S0, D0 - S0, L4 - S4, D4 - S4)
                    w_na = est.loc[(focal, st), 'omega_na']
                    bt = np.einsum('rg,gck->rck', w, gs[idx])
                    b_obs = omega(bt[:, 0, 0], bt[:, 0, 1], bt[:, 1, 0], bt[:, 1, 1])
                    b_rem = omega(bt[:, 0, 0] - bt[:, 0, 3], bt[:, 0, 1] - bt[:, 0, 3],
                                  bt[:, 1, 0] - bt[:, 1, 3], bt[:, 1, 1] - bt[:, 1, 3])
                    a_obs, a_rem = 1 - boot_wna[focal, st] / b_obs, 1 - boot_wna[focal, st] / b_rem
                    boots[pool, focal, st, f] = (a_obs, a_rem)
                    sf4 = bt[:, 1, 3] / bt[:, 1, 2]
                    row = dict(pair=pool, focal=focal, other=pool if focal == SUB else SUB, stratum=st,
                               thresh=f, L0=int(L0), D0=int(D0), T0=int(T0), S0=int(S0), E0=E0, L4=int(L4),
                               D4=int(D4), T4=int(T4), S4=int(S4), E4=E4, testable_frac4=T4 / D4,
                               testable_frac0=T0 / D0, shared_frac4=S4 / T4,
                               shared_frac4_lo=np.percentile(sf4, 2.5),
                               shared_frac4_hi=np.percentile(sf4, 97.5), shared_frac0=S0 / T0,
                               obs_over_exp4=S4 / E4 if E4 else np.nan, obs_over_exp0=S0 / E0 if E0 else np.nan,
                               d4_per_site=D4 / L4, d4_per_site_removed=(D4 - S4) / (L4 - S4),
                               omega_obs=w_obs, omega_removed=w_rem,
                               omega_shared=(S0 / L0) / (S4 / L4) if S4 else np.nan,
                               omega_na=w_na, alpha_obs=1 - w_na / w_obs, alpha_removed=1 - w_na / w_rem,
                               alpha_obs_lo=np.percentile(a_obs, 2.5), alpha_obs_hi=np.percentile(a_obs, 97.5),
                               alpha_removed_lo=np.percentile(a_rem, 2.5),
                               alpha_removed_hi=np.percentile(a_rem, 97.5),
                               alpha_A1=est.loc[(focal, st), 'alpha'],
                               pi0_over_pi4_focal=est.loc[(focal, st), 'pi0'] / est.loc[(focal, st), 'pi4'],
                               baseline_rate4=float(rate[1][rate[1] > 0].mean()))
                    rows.append(row)
                    say(f'  f={f} focal {focal} {st}: D4 {int(D4)} testable {100 * T4 / D4:.1f}% shared {int(S4)} '
                        f'({100 * S4 / T4:.2f}%, expected {E4:.0f}); D0 {int(D0)} testable {100 * T0 / D0:.1f}% '
                        f'shared {int(S0)} ({100 * S0 / T0:.2f}%, expected '
                        f'{E0:.0f}); omega shared {row["omega_shared"]:.3f}; alpha {row["alpha_obs"]:.3f} '
                        f'({row["alpha_obs_lo"]:.3f}-{row["alpha_obs_hi"]:.3f}) -> {row["alpha_removed"]:.3f} '
                        f'({row["alpha_removed_lo"]:.3f}-{row["alpha_removed_hi"]:.3f}); A1 {row["alpha_A1"]:.3f}')
    out = pd.DataFrame(rows)
    out.to_csv(f'{RES}/alpha_shared_polymorphism.tsv', sep='\t', index=False)
    np.savez_compressed(f'{QC}/site_table/alpha_shared_boot.npz',
                        **{'|'.join(map(str, k)) + '|' + lab: v[i] for k, v in boots.items()
                           for i, lab in enumerate(('obs', 'removed'))})
    # Lineage difference in alpha. The two lineages' bootstraps resample different gene tables
    # independently; pairing replicate k with replicate k gives a draw of the difference that ignores
    # the positive correlation between orthologous genes, so the interval is conservative (too wide).
    summ = []
    for f in THRESH:
        p = out[out.thresh == f]
        for pool in POOLS:
            for st in STRATA:
                s = p[(p.pair == pool) & (p.focal == SUB) & (p.stratum == st)].iloc[0]
                z = p[(p.pair == pool) & (p.focal == pool) & (p.stratum == st)].iloc[0]
                d_obs = boots[pool, SUB, st, f][0] - boots[pool, pool, st, f][0]
                d_rem = boots[pool, SUB, st, f][1] - boots[pool, pool, st, f][1]
                r = dict(pair=pool, stratum=st, thresh=f,
                         dalpha_obs=s.alpha_obs - z.alpha_obs, dalpha_obs_lo=np.percentile(d_obs, 2.5),
                         dalpha_obs_hi=np.percentile(d_obs, 97.5),
                         dalpha_removed=s.alpha_removed - z.alpha_removed,
                         dalpha_removed_lo=np.percentile(d_rem, 2.5), dalpha_removed_hi=np.percentile(d_rem, 97.5),
                         excess_d4_obs=s.d4_per_site / z.d4_per_site - 1,
                         excess_d4_removed=s.d4_per_site_removed / z.d4_per_site_removed - 1,
                         shared4_per_site_dsub=s.S4 / s.L4, shared4_per_site_dsuz=z.S4 / z.L4,
                         nonshared4_per_site_dsub=(s.D4 - s.S4) / s.L4, nonshared4_per_site_dsuz=(z.D4 - z.S4) / z.L4)
                summ.append(r)
                if f == THRESH[0]:
                    say(f'summary {pool} {st}: alpha Dsub - Dsuz {r["dalpha_obs"]:+.3f} '
                        f'({r["dalpha_obs_lo"]:+.3f} to {r["dalpha_obs_hi"]:+.3f}) -> {r["dalpha_removed"]:+.3f} '
                        f'({r["dalpha_removed_lo"]:+.3f} to {r["dalpha_removed_hi"]:+.3f}) after exclusion; '
                        f'excess 4-fold divergence per site {100 * r["excess_d4_obs"]:+.1f}% -> '
                        f'{100 * r["excess_d4_removed"]:+.1f}%; 4-fold per 1000 sites shared Dsub '
                        f'{1000 * r["shared4_per_site_dsub"]:.2f} / Dsuz {1000 * r["shared4_per_site_dsuz"]:.2f}, '
                        f'not shared Dsub {1000 * r["nonshared4_per_site_dsub"]:.2f} / Dsuz '
                        f'{1000 * r["nonshared4_per_site_dsuz"]:.2f}')
    pd.DataFrame(summ).to_csv(f'{RES}/alpha_shared_polymorphism_summary.tsv', sep='\t', index=False)
    say(f'finished {time.ctime()}')


if __name__ == '__main__':
    main()
