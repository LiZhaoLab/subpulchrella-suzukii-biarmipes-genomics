#!/usr/bin/env python3
"""Inversion analysis: combine breakpoint evidence (approach 2, scripts/109) with the SNP
frequency signal (approach 3, scripts/110) for the D. subpulchrella pool.

1. Delly inversion calls genotyped in the pool and in the reference strain (sv/geno.inv.tsv). Read
   support: split reads (RV of RR + RV) for precise calls, discordant pairs (DV of DR + DV) otherwise;
   pool frequency of the inverted junction f = support / (support + reference reads).
   Candidate polymorphic inversion: major chromosome; length >= 10 kb; pool support >= 4 reads and
   f >= 0.05; reference strain with >= 10 informative reads and at most 1 supporting read (a junction
   seen in the reference strain is an assembly or mapping artefact, because the reference was
   assembled from this strain). The two junctions of one inversion (3to3 and 5to5) are paired when both
   ends lie within 5 kb.
2. Annotation: fraction of +/- 500 bp around each breakpoint in the RepeatMasker mask; distance to the
   nearest D. suzukii-D. subpulchrella synteny breakpoint (results/synteny_segments.tsv).
   Tiers: 'high' = both junctions paired, both breakpoints < 50% repeat, length >= 100 kb, pool frequency
   0.05-0.95; 'single_junction' = only one junction called (typical of TE-mediated artefacts and
   non-reference TE insertions); 'repeat_breakpoint' = a breakpoint >= 50% repeat. Calls supported in the
   reference strain are listed separately as 'reference_strain' (assembly orientation errors, or an
   arrangement carried by the reference strain that is rare in the pool).
3. SNP frequency test per candidate >= 100 kb: fraction of SNPs (MAF >= 0.05, scripts/110) within 0.025
   of the folded junction frequency min(f, 1 - f), inside the inversion and within 200 kb inside each
   breakpoint, against 1,000 random intervals of the same length on the same chromosome (seed 20261030).
Outputs: results/INV_candidates.tsv, results/INV_candidates.log
"""
import numpy as np
import pandas as pd

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
MAJOR = ['NC_050609.1', 'NC_050610.1', 'NC_050611.1', 'NC_050612.1', 'NC_050613.1', 'NW_023665577.1']
SEED = 20261030


def main():
    g = pd.read_csv(f'{QC}/sv/geno.inv.tsv', sep='\t')
    g.columns = ['chrom', 'pos', 'end', 'id', 'filter', 'ct', 'pe', 'sr', 'precise', 'qual',
                 'p_dr', 'p_dv', 'p_rr', 'p_rv', 'r_dr', 'r_dv', 'r_rr', 'r_rv']
    for c in ['pe', 'sr', 'p_dr', 'p_dv', 'p_rr', 'p_rv', 'r_dr', 'r_dv', 'r_rr', 'r_rv']:
        g[c] = pd.to_numeric(g[c], errors='coerce').fillna(0)
    g['precise'] = g.precise.astype(str).isin(['1', 'True'])
    g['length'] = g.end - g.pos
    g['p_sup'] = np.where(g.precise, g.p_rv, g.p_dv)
    g['p_ref'] = np.where(g.precise, g.p_rr, g.p_dr)
    g['r_sup'] = np.where(g.precise, g.r_rv, g.r_dv)
    g['r_ref'] = np.where(g.precise, g.r_rr, g.r_dr)
    g['f_pool'] = g.p_sup / (g.p_sup + g.p_ref).clip(lower=1)
    g['f_ref'] = g.r_sup / (g.r_sup + g.r_ref).clip(lower=1)
    log = [f'Delly inversion junction calls (merged sites genotyped in both BAMs): {len(g)}']
    m = g.chrom.isin(MAJOR) & (g.length >= 10_000)
    log.append(f'on major chromosomes and >= 10 kb: {m.sum()}')
    art = m & (g.r_sup >= 2) & (g.f_ref >= 0.2)
    log.append(f'  supported in the reference strain (>= 2 reads, f >= 0.2; artefacts or fixed differences '
               f'between the strain and the assembly): {art.sum()}')
    refsup = g[art].copy()
    cand = m & (g.p_sup >= 4) & (g.f_pool >= 0.05) & (g.r_sup <= 1) & (g.r_sup + g.r_ref >= 10)
    c = g[cand].copy()
    log.append(f'  candidate polymorphic junctions (pool >= 4 reads, f >= 0.05; reference strain <= 1 read of '
               f'>= 10): {len(c)}; >= 100 kb: {(c.length >= 100_000).sum()}')
    # pair 3to3 with 5to5 junctions
    c['pair'] = ''
    a = c[c.ct == '3to3']
    b = c[c.ct == '5to5']
    print(f'pairing {len(a)} 3to3 with {len(b)} 5to5 junctions', flush=True)
    for i, r in a.iterrows():
        mm = b[(b.chrom == r.chrom) & ((b.pos - r.pos).abs() <= 5000) & ((b.end - r.end).abs() <= 5000)]
        if len(mm):
            c.loc[i, 'pair'] = ','.join(mm.id)
            c.loc[mm.index, 'pair'] = r.id
    log.append(f'  candidates with both junctions (3to3 + 5to5 within 5 kb): {(c.pair != "").sum()} junction calls')
    # repeats around breakpoints
    ms = np.load(f'{QC}/masks_sites/Dsub.npz')
    rm_cache = {ch: ms[f'{ch}|rm'] for ch in c.chrom.unique()}      # decompress each mask once
    rep = []
    for r in c.itertuples():
        rm = rm_cache[r.chrom]
        f1 = rm[max(1, r.pos - 500):r.pos + 501].mean()
        f2 = rm[max(1, r.end - 500):r.end + 501].mean()
        rep.append(max(f1, f2))
    c['repeat_frac_max_bp'] = rep
    c['tier'] = np.where(c.repeat_frac_max_bp >= 0.5, 'repeat_breakpoint',
                np.where(c.pair == '', 'single_junction',
                np.where((c.length >= 100_000) & (c.f_pool <= 0.95), 'high', 'paired_other')))
    log.append('  tiers: ' + ', '.join(f'{k} {v}' for k, v in c.tier.value_counts().items()))
    # synteny breakpoints in D. subpulchrella coordinates
    syn = pd.read_csv(f'{RES}/synteny_segments.tsv', sep='\t')
    syn = syn[(syn.params == 'main') & (syn.comp == 'Dsuz_vs_Dsub')]
    sb = pd.concat([syn[['query', 'qs']].rename(columns={'query': 'chrom', 'qs': 'x'}),
                    syn[['query', 'qe']].rename(columns={'query': 'chrom', 'qe': 'x'})])
    dist = []
    for r in c.itertuples():
        xs = sb[sb.chrom == r.chrom].x.values
        dist.append(min(np.abs(xs - r.pos).min(), np.abs(xs - r.end).min()) if len(xs) else np.nan)
    c['dist_synteny_bp_kb'] = np.array(dist) / 1000
    # SNP frequency test
    snpz = np.load(f'{QC}/inv/Dsub_L5_7.snps.npz')
    snp = {k: (snpz[k].astype(np.int64) if k.endswith('|pos') else snpz[k]) for k in snpz.files}   # decompress once;
    # int64 positions: searchsorted with int64 bounds on an int32 array copies the array at every call
    rng = np.random.default_rng(SEED)
    out = []
    for r in c.itertuples():
        if r.length < 100_000:
            out.append((np.nan,) * 6)
            continue
        pos, maf = snp[f'{r.chrom}|pos'], snp[f'{r.chrom}|maf']
        fb = min(r.f_pool, 1 - r.f_pool)
        cum_near = np.concatenate([[0], np.cumsum(np.abs(maf - fb) <= 0.025)])

        def count(s, e):
            i, j = np.searchsorted(pos, s, 'left'), np.searchsorted(pos, e, 'right')
            return cum_near[j] - cum_near[i], j - i

        def frac(s, e):
            k, n = count(s, e)
            return k / n if n else np.nan

        def frac_bp(s, e):
            w = min(200_000, (e - s) // 2)
            k1, n1 = count(s, s + w)
            k2, n2 = count(e - w, e)
            return (k1 + k2) / (n1 + n2) if n1 + n2 else np.nan

        obs_in, obs_bp = frac(r.pos, r.end), frac_bp(r.pos, r.end)
        L, hi = r.end - r.pos, int(pos.max())
        st = rng.integers(1, max(2, hi - L), 1000)
        rin = np.array([frac(s, s + L) for s in st])
        rbp = np.array([frac_bp(s, s + L) for s in st])
        out.append((fb, obs_in, np.nanmean(rin), (rin >= obs_in).mean(), obs_bp, (rbp >= obs_bp).mean()))
    c[['folded_f', 'snp_frac_inside', 'snp_frac_random', 'P_inside', 'snp_frac_breakpoints', 'P_breakpoints']] = out
    c = c.sort_values(['chrom', 'pos'])
    c.to_csv(f'{RES}/INV_candidates.tsv', sep='\t', index=False)
    big_ref = refsup[(refsup.length >= 100_000)].sort_values('r_sup', ascending=False)
    big_ref.to_csv(f'{RES}/INV_reference_strain_supported.tsv', sep='\t', index=False)
    log.append(f'\nreference-strain-supported junctions >= 100 kb: {len(big_ref)} (results/INV_reference_strain_supported.tsv); largest support:')
    for r in big_ref.head(10).itertuples():
        log.append(f'  {r.chrom}:{r.pos}-{r.end} ({r.length / 1e6:.2f} Mb) {r.ct} ref strain {int(r.r_sup)}/{int(r.r_sup + r.r_ref)} '
                   f'f {r.f_ref:.2f}; pool {int(r.p_sup)}/{int(r.p_sup + r.p_ref)} f {r.f_pool:.2f}')
    log.append('\nhigh-tier candidates (paired junctions, breakpoints outside repeats, >= 100 kb):')
    for r in c[c.tier == 'high'].itertuples():
        log.append(f'{r.chrom}:{r.pos}-{r.end} ({r.length / 1e6:.2f} Mb) {r.ct} {"PRECISE" if r.precise else "IMPRECISE"} '
                   f'pool {int(r.p_sup)}/{int(r.p_sup + r.p_ref)} f {r.f_pool:.3f}; ref strain {int(r.r_sup)}/'
                   f'{int(r.r_sup + r.r_ref)}; paired {bool(r.pair)}; repeats at bp {r.repeat_frac_max_bp:.2f}; '
                   f'synteny bp {r.dist_synteny_bp_kb:.0f} kb; SNPs at f: inside {r.snp_frac_inside:.4f} vs random '
                   f'{r.snp_frac_random:.4f} (P {r.P_inside:.3f}), near breakpoints {r.snp_frac_breakpoints:.4f} '
                   f'(P {r.P_breakpoints:.3f})')
    log.append(f'\ncandidates 10-100 kb: {((c.length >= 10_000) & (c.length < 100_000)).sum()} '
               f'(listed in the table only)')
    open(f'{RES}/INV_candidates.log', 'w').write('\n'.join(log) + '\n')
    print('\n'.join(log))


if __name__ == '__main__':
    main()
