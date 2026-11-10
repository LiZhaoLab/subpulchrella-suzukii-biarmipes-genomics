#!/usr/bin/env python3
"""Analysis F of the depth plan: summary of the SweepFinder2 scans.

Peaks: runs of grid points (2-kb grid) with CLR above a threshold t, merged when closer than
10 kb. Null: neutral simulations of scripts/100 (two-epoch model fitted to each pool and chromosome
class, recombination/mutation rate ratio 2). For thresholds t at null grid-point quantiles
(0.99, 0.999, 0.9999 and the null maximum), the null peak rate per Mb gives the expected number of
false peaks in the scanned genome; estimated false discovery rate (FDR) = expected / observed.
With no null peak above t, the rate is bounded by 3 / (null Mb) (rule of three). The reported
threshold is the lowest with FDR <= 0.10 (else the null maximum).
Genes: RefSeq mRNAs overlapping a peak +/- 5 kb, with D. melanogaster symbols from
results/L3_lineage_gene_sets.tsv. Overlap with lineage-specific MK-positive genes
(DAF >= 15%): Fisher's exact test (genes are not independent; descriptive).
Between species: D. subpulchrella peak genes mapped to D. suzukii orthologs (Step-0 tables);
expected overlap from 1,000 random placements of the peaks (rotation within chromosome).
The dot chromosome (Muller F) is excluded: it barely recombines, and the null does not model it.
Outputs: results/F_sweep_peaks.tsv, results/F_sweep_thresholds.tsv,
         results/F_sweep_summary.log
"""
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
A = '/ru-auth/local/home/asato/scratch'
GFF = {'Dsub': f'{A}/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff',
       'Dsuz': f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'}
POOLS = ['Dsub_L5_7', 'JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
XCHROM = {'NC_050613.1', 'NC_092084.1'}
DOT = {'NC_050608.1', 'NC_092083.1'}      # Muller F: almost no recombination; not calibrated by the null
MERGE, FLANK, SEED = 10_000, 5_000, 20261030


def peaks_from(pos, clr, t):
    above = clr > t
    if not above.any():
        return []
    p, c = pos[above], clr[above]
    out, s, e, m, mp = [], p[0], p[0], c[0], p[0]
    for x, y in zip(p[1:], c[1:]):
        if x - e <= MERGE:
            e = x
            if y > m:
                m, mp = y, x
        else:
            out.append((s, e, m, mp))
            s, e, m, mp = x, x, y, x
    out.append((s, e, m, mp))
    return out


def mrnas(gff):
    rows = []
    with open(gff) as f:
        for line in f:
            t = line.split('\t')
            if len(t) > 8 and t[2] == 'mRNA':
                a = dict(kv.split('=', 1) for kv in t[8].strip().split(';') if '=' in kv)
                rows.append((t[0], int(t[3]), int(t[4]), a.get('ID', '').replace('rna-', ''), a.get('gene', '')))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'tx', 'gene'])


def main():
    t0 = time.time()
    tag = '' if not sys.argv[1:] else '_test'
    logf = open(f'{RES}/F_sweep_summary{tag}.log', 'w')

    def log(s):
        print(s, flush=True)
        logf.write(s + '\n')

    log(f'F sweep summary; started {time.ctime()}')
    gs = pd.read_csv(f'{RES}/L3_lineage_gene_sets.tsv', sep='\t')
    sym = {}
    for t, s in zip(gs.transcripts.astype(str), gs.symbol.astype(str)):
        for tx in t.split(','):
            sym[tx.strip()] = s
    mkpos = {'Dsub': set(), 'Dsuz': set()}
    for setname, sp in (('Dsub_lineage', 'Dsub'), ('Dsuz_lineage_union', 'Dsuz')):
        g = gs[(gs.version == 'DAF15') & (gs.set == setname)]
        for t, flag in zip(g.transcripts.astype(str), g.MK_positive):
            for tx in t.split(','):
                if flag:
                    mkpos[sp].add(tx.strip())
    mk_tested = {sp: set() for sp in mkpos}
    for setname, sp in (('Dsub_lineage', 'Dsub'), ('Dsuz_lineage_union', 'Dsuz')):
        g = gs[(gs.version == 'DAF15') & (gs.set == setname)]
        for t in g.transcripts.astype(str):
            mk_tested[sp].update(x.strip() for x in t.split(','))
    mr = {sp: mrnas(GFF[sp]) for sp in GFF}
    thr_rows, peak_rows = [], []
    pools = sys.argv[1:] or POOLS
    for pool in pools:
        sp = 'Dsub' if pool.startswith('Dsub') else 'Dsuz'
        files = glob.glob(f'{QC}/sweeps/{pool}/clr/*.txt')
        scan = []
        for f in files:
            chrom = os.path.basename(f).rsplit('_', 1)[0]
            d = pd.read_csv(f, sep='\t')
            d['chrom'] = chrom
            scan.append(d)
        scan = pd.concat(scan).sort_values(['chrom', 'location']).reset_index(drop=True)
        scan = scan[~scan.chrom.isin(DOT)]
        scan['cls'] = np.where(scan.chrom.isin(XCHROM), 'X', 'autosome')
        for cls, sc in scan.groupby('cls'):
            nd = f'{QC}/sweeps/{pool}/null_{cls}'
            if len(glob.glob(f'{nd}/clr_*.txt')) < 40:
                log(f'{pool} {cls}: null not complete; skipped')
                continue
            fit = json.load(open(f'{nd}/fit.json'))
            null = [pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{nd}/clr_*.txt'))]
            null_mb = len(null) * 1.0
            allv = np.concatenate([n.LR.values for n in null])
            scanned_mb = len(sc) * 2000 / 1e6
            cands = sorted(set([np.quantile(allv, q) for q in (0.99, 0.999, 0.9999)] + [allv.max()]))
            best = None
            for t in cands:
                n_null = sum(len(peaks_from(n.location.values, n.LR.values, t)) for n in null)
                rate = max(n_null, 3) / null_mb if n_null == 0 else n_null / null_mb
                obs = sum(len(peaks_from(g.location.values, g.LR.values, t)) for _, g in sc.groupby('chrom'))
                exp = rate * scanned_mb
                fdr = min(1.0, exp / obs) if obs else np.nan
                lo_k = chi2.ppf(0.025, 2 * n_null) / 2 if n_null else 0.0
                hi_k = chi2.ppf(0.975, 2 * n_null + 2) / 2
                thr_rows.append(dict(pool=pool, cls=cls, threshold=t, null_peaks=n_null, null_mb=null_mb,
                                     null_rate_per_mb=rate, scanned_mb=scanned_mb, observed_peaks=obs,
                                     expected_false=exp, fdr=fdr,
                                     fdr_lo=min(1.0, lo_k / null_mb * scanned_mb / obs) if obs else np.nan,
                                     fdr_hi=min(1.0, hi_k / null_mb * scanned_mb / obs) if obs else np.nan,
                                     model=fit['model_used'], nu=fit['nu'], T=fit['T']))
                if best is None and obs and fdr <= 0.10:
                    best = (t, fdr)
            if best is None:
                best = (allv.max(), thr_rows[-1]['fdr'])
            t, fdr = best
            tr = [r for r in thr_rows if r['pool'] == pool and r['cls'] == cls and r['threshold'] == t][0]
            log(f'{pool} {cls}: model {fit["model_used"]} (nu {fit["nu"]:.2f}, T {fit["T"]:.3f}); scanned '
                f'{scanned_mb:.0f} Mb; null {null_mb:.0f} Mb; threshold CLR {t:.1f}; peaks {tr["observed_peaks"]}; '
                f'estimated FDR {fdr:.3f} (Poisson 95% {tr["fdr_lo"]:.3f}-{tr["fdr_hi"]:.3f}; null peaks {tr["null_peaks"]})')
            for chrom, g in sc.groupby('chrom'):
                for s, e, m, mp in peaks_from(g.location.values, g.LR.values, t):
                    mm = mr[sp]
                    hit = mm[(mm.chrom == chrom) & (mm.end >= s - FLANK) & (mm.start <= e + FLANK)]
                    txs = sorted(set(hit.tx))
                    peak_rows.append(dict(pool=pool, species=sp, cls=cls, chrom=chrom, start=int(s), end=int(e),
                                          max_clr=m, max_pos=int(mp), threshold=t, fdr=fdr,
                                          transcripts=','.join(txs),
                                          symbols=','.join(sorted({sym.get(x, '') for x in txs} - {''})),
                                          genes=','.join(sorted(set(hit.gene))),
                                          mk_positive=','.join(sorted({sym.get(x, x) for x in txs if x in mkpos[sp]}))))
    th = pd.DataFrame(thr_rows)
    pk = pd.DataFrame(peak_rows)
    # peak density per chromosome (scanned grid length), for comparison with the diversity deficit
    dens = []
    for pool in pools:
        for f in glob.glob(f'{QC}/sweeps/{pool}/clr/*.txt'):
            chrom = os.path.basename(f).rsplit('_', 1)[0]
            if chrom in DOT:
                continue
            dens.append((pool, chrom, len(pd.read_csv(f, sep='\t')) * 2000 / 1e6))
    dens = pd.DataFrame(dens, columns=['pool', 'chrom', 'mb']).groupby(['pool', 'chrom']).mb.sum().reset_index()
    if not pk.empty:
        cnt = pk.groupby(['pool', 'chrom']).agg(peaks=('start', 'size'),
                                                peak_kb=('end', lambda e: 0)).reset_index()
        width = pk.assign(w=(pk.end - pk.start + 2000) / 1000).groupby(['pool', 'chrom']).w.sum().reset_index()
        dens = dens.merge(cnt[['pool', 'chrom', 'peaks']], how='left').merge(width, how='left').fillna(0)
        dens['peaks_per_mb'] = dens.peaks / dens.mb
        dens['peak_fraction'] = dens.w / 1000 / dens.mb
        dens.to_csv(f'{RES}/F_sweep_density{tag}.tsv', sep='\t', index=False)
        log('\n== peaks per Mb by chromosome (dot chromosome excluded)')
        for r in dens.itertuples():
            log(f'{r.pool:10s} {r.chrom:15s} scanned {r.mb:6.1f} Mb peaks {int(r.peaks):4d} per Mb {r.peaks_per_mb:.2f} '
                f'fraction of length in peaks {r.peak_fraction:.4f}')
    th.to_csv(f'{RES}/F_sweep_thresholds{tag}.tsv', sep='\t', index=False)
    pk.to_csv(f'{RES}/F_sweep_peaks{tag}.tsv', sep='\t', index=False)
    if pk.empty or 'Dsub_L5_7' not in set(pk.pool) or len(set(pk.pool)) < 2:
        log(f'finished {time.ctime()} (partial run; overlap analyses skipped)')
        return
    log('\n== peaks per pool (FDR <= 0.10 threshold)')
    for pool, g in pk.groupby('pool'):
        log(f'{pool}: {len(g)} peaks, median width {np.median(g.end - g.start + 2000) / 1000:.0f} kb, '
            f'{g.transcripts.str.len().gt(0).sum()} with genes')
    # MK-positive overlap
    log('\n== overlap of peak genes with lineage-specific MK-positive genes (DAF >= 15%), MK-tested genes only')
    for pool, g in pk.groupby('pool'):
        sp = 'Dsub' if pool.startswith('Dsub') else 'Dsuz'
        inpk = set(x for t in g.transcripts for x in t.split(',') if x) & mk_tested[sp]
        tested = mk_tested[sp]
        a = len(inpk & mkpos[sp])
        b = len(inpk) - a
        c = len(mkpos[sp] & tested) - a
        d = len(tested) - a - b - c
        if inpk:
            orr, p = fisher_exact([[a, b], [c, d]])
            log(f'{pool}: peak genes {len(inpk)}, MK-positive {a} ({a / len(inpk):.3f}) vs genome '
                f'{(a + c) / len(tested):.3f}; odds ratio {orr:.2f}, Fisher P {p:.2g}')
    # replication between D. suzukii pools: fraction of peaks of pool A within 10 kb of a peak of pool B
    log('\n== D. suzukii pools: peaks of A within 10 kb of a peak of B (expected: peaks of A rotated '
        'within chromosome, 200 replicates)')
    rng0 = np.random.default_rng(SEED + 1)
    suz = [p for p in pools if not p.startswith('Dsub')]
    lens_z = {}
    sz = np.load(f'{QC}/sweeps/JP-Kan/summary.npz')
    lens_z = {k.split('|')[0]: int(sz[k][0]) for k in sz.files if k.endswith('|length')}

    def frac_near(ga, gb, shift=None):
        hit = 0
        for chrom, a in ga.groupby('chrom'):
            b = gb[gb.chrom == chrom]
            ca = ((a.start + a.end) // 2).values
            if shift is not None:
                ca = (ca + shift[chrom]) % lens_z[chrom]
            if len(b) == 0:
                continue
            bs, be = b.start.values - MERGE, b.end.values + MERGE
            hit += sum(((x >= bs) & (x <= be)).any() for x in ca)
        return hit / max(len(ga), 1)

    for i, a in enumerate(suz):
        for b in suz[i + 1:]:
            ga, gb = pk[pk.pool == a], pk[pk.pool == b]
            obs = frac_near(ga, gb)
            ex = [frac_near(ga, gb, {c: rng0.integers(0, L) for c, L in lens_z.items()}) for _ in range(200)]
            log(f'{a} peaks near {b} peaks: {obs:.3f} (expected {np.mean(ex):.3f}, 95% {np.percentile(ex, 2.5):.3f}-'
                f'{np.percentile(ex, 97.5):.3f})')
    # between species
    z = np.load(f'{QC}/site_table/Dsub_L5_7.npz')
    sub2suz = dict(zip(z['genes'], z['gene_sister']))
    log('\n== D. subpulchrella peaks with an orthologous gene in a peak of each D. suzukii pool')
    rng = np.random.default_rng(SEED)
    dsub = pk[pk.pool == 'Dsub_L5_7']
    sub_genes = set(sub2suz.get(x) for t in dsub.transcripts for x in t.split(',') if x in sub2suz)
    for pool in POOLS[1:]:
        g = pk[pk.pool == pool]
        suz_genes = set(x for t in g.transcripts for x in t.split(',') if x)
        obs = len(sub_genes & suz_genes)
        # expectation: rotate D. suzukii peaks within chromosomes and count overlapping orthologs
        mm = mr['Dsuz']
        lens = mm.groupby('chrom').end.max()
        arr = {c: (x.start.values, x.end.values, x.tx.values) for c, x in mm.groupby('chrom')}
        exp = []
        for _ in range(1000):
            genes_r = set()
            for chrom, gg in g.groupby('chrom'):
                L = lens.get(chrom, 0)
                if not L:
                    continue
                sh = rng.integers(0, L)
                st, en, tx = arr[chrom]
                for s, e in zip(gg.start, gg.end):
                    s2, e2 = (s + sh) % L, (e + sh) % L
                    if e2 < s2:
                        continue
                    genes_r.update(tx[(en >= s2 - FLANK) & (st <= e2 + FLANK)])
            exp.append(len(sub_genes & genes_r))
        exp = np.array(exp)
        log(f'Dsub vs {pool}: shared orthologous peak genes {obs}; expected {exp.mean():.1f} '
            f'(95% {np.percentile(exp, 2.5):.0f}-{np.percentile(exp, 97.5):.0f}); P(random >= obs) '
            f'{(exp >= obs).mean():.3g}')
    log(f'finished {time.ctime()} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
