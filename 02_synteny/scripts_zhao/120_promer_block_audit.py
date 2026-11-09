#!/usr/bin/env python3
"""Script 120: what one PROmer "alignment block" is, to decide whether the block counts in the
Results sentence (77,307 / 45,151 / 26,249 blocks) are a meaningful quantity to report. Airi Sato's run
synteny_4species/reanalysis_20260824 (MUMmer 4.0.1 PROmer --mumreference on RepeatMasker hard-masked genomes;
delta-filter -1 -l 100; then LEN1 and LEN2 >= 500 bp and amino-acid identity >= 50%). No new alignment.

A. Threshold sensitivity: block counts for minimum length (both sides) x minimum amino-acid identity, from the
   one-to-one L100 coordinate files (the input of Airi's filter).
B. Length and identity of the retained blocks.
C. Overlap between retained blocks: sum of block bp vs merged bp on each side.
D. Sequence classes: blocks between homologous major chromosomes, blocks involving unplaced sequences or mtDNA.
E. Coding content on the reference side (NCBI RefSeq GFF of the reference genome): share of block bp in CDS
   (all isoforms merged), blocks without CDS, genes touched, blocks per gene, number of merged CDS intervals.
F. Fragmentation: blocks adjacent along the reference that continue each other (same query sequence, same
   orientation, query continues in the expected direction; both gaps between -500 bp and 10 kb). Gap length,
   share of gaps containing hard-masked sequence (N) in the reference, share with a reading-frame change.
   Simple chains of such blocks (gap <= 10 kb, then <= 100 kb, no size filter) give the number of runs.
G. Blocks per collinear segment of script 91 (results/synteny_segments.tsv).
Output: results/promer_block_audit.log (stdout) and results/promer_block_thresholds.tsv.
Python: anaconda3/envs/jupyter-ml.
"""
import re
import sys

import numpy as np
import pandas as pd

BASE = '/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824'
GEN = '/ru-auth/local/home/asato/scratch/genome'
OUT = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
COLS = ['S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', 'IDY', 'SIM', 'STP', 'LENR', 'LENQ', 'COVR', 'COVQ',
        'FRM1', 'FRM2', 'R', 'Q']
COMPS = {  # comparison: (reference masked FASTA, reference RefSeq GFF)
    'Dsuz_vs_Dsub': (f'{BASE}/repeatmasker/Dsuz/Dsuz.GCF_043229965.1.fna.masked',
                     f'{GEN}/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'),
    'Dsub_vs_Dbia': (f'{BASE}/repeatmasker/Dsub/Dsub.GCF_014743375.2.fna.masked',
                     f'{GEN}/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff'),
    'Dbia_vs_Dmel': (f'{BASE}/repeatmasker/Dbia/Dbia.GCF_025231255.1.fna.masked',
                     f'{GEN}/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.gff'),
}
MAJOR = {  # major chromosomes (arm label); same table as script 91
    'NC_092080.1': '2L', 'NC_092081.1': '2R', 'NC_092082.1': '3', 'NC_092083.1': '4', 'NC_092084.1': 'X',
    'NC_050610.1': '2L', 'NC_050611.1': '2R', 'NC_050612.1': '3L', 'NC_050609.1': '3R',
    'NW_023665577.1': '3R', 'NC_050608.1': '4', 'NC_050613.1': 'X',
    'NC_066612.1': '2L', 'NC_066615.1': '2R', 'NC_066613.1': '3L', 'NC_066616.1': '3R', 'NC_066614.1': '4',
    'NC_066611.1': 'X',
    'NT_033779.5': '2L', 'NT_033778.4': '2R', 'NT_037436.4': '3L', 'NT_033777.3': '3R', 'NC_004353.4': '4',
    'NC_004354.4': 'X'}
MITO = {'NC_060762.1', 'NC_024511.2'}  # D. suzukii, D. melanogaster (D. subpulchrella, D. biarmipes: none)
HOMOLOG = lambda a, b: a == b or (a == '3' and b.startswith('3')) or (b == '3' and a.startswith('3'))
LMIN = [100, 200, 300, 500, 1000, 2000]
IMIN = [0, 30, 50, 70, 90]


def load(path):
    t = pd.read_csv(path, sep='\t', header=None, names=COLS)
    assert t.shape[1] == 17 and t['IDY'].between(0, 100).all()
    t['rs'], t['re'] = t[['S1', 'E1']].min(1) - 1, t[['S1', 'E1']].max(1)  # 0-based half-open
    t['qs'], t['qe'] = t[['S2', 'E2']].min(1) - 1, t[['S2', 'E2']].max(1)
    t['ori'] = np.sign((t.E1 - t.S1) * (t.E2 - t.S2)).astype(int)
    return t


def merge(iv):
    """iv: DataFrame seq, s, e -> merged intervals (same columns)."""
    out = []
    for seq, g in iv.sort_values(['seq', 's']).groupby('seq', sort=False):
        s, e = g.s.to_numpy(), g.e.to_numpy()
        cur_s, cur_e = s[0], e[0]
        for a, b in zip(s[1:], e[1:]):
            if a <= cur_e:
                cur_e = max(cur_e, b)
            else:
                out.append((seq, cur_s, cur_e))
                cur_s, cur_e = a, b
        out.append((seq, cur_s, cur_e))
    return pd.DataFrame(out, columns=['seq', 's', 'e'])


class Cover:
    """Base pairs of a merged interval set inside [a, b), vectorized per sequence."""

    def __init__(self, merged):
        self.d = {}
        for seq, g in merged.groupby('seq'):
            s, e = g.s.to_numpy(np.int64), g.e.to_numpy(np.int64)
            self.d[seq] = (s, e, np.concatenate([[0], np.cumsum(e - s)]))

    def _upto(self, seq, x):
        s, e, cum = self.d[seq]
        i = np.searchsorted(s, x, side='right')  # intervals starting before or at x
        full = cum[np.maximum(i - 1, 0)]
        part = np.where(i > 0, np.minimum(x, e[np.maximum(i - 1, 0)]) - s[np.maximum(i - 1, 0)], 0)
        return np.where(i > 0, full + np.maximum(part, 0), 0)

    def bp(self, seqs, a, b):
        res = np.zeros(len(a), np.int64)
        seqs = np.asarray(seqs)
        for seq in np.unique(seqs):
            if seq not in self.d:
                continue
            m = seqs == seq
            res[m] = self._upto(seq, b[m]) - self._upto(seq, a[m])
        return res


def n_runs(fasta):
    """Hard-masked (N) runs and sequence lengths of a FASTA."""
    runs, lens, name, buf = [], {}, None, []

    def flush():
        if name is None:
            return
        s = ''.join(buf).upper()
        lens[name] = len(s)
        for m in re.finditer(r'N+', s):
            runs.append((name, m.start(), m.end()))

    with open(fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    flush()
    return pd.DataFrame(runs, columns=['seq', 's', 'e']), lens


def read_cds(gff):
    rows = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            m = re.search(r'(?:^|;)gene=([^;]+)', f[8])
            rows.append((f[0], int(f[3]) - 1, int(f[4]), m.group(1) if m else 'NA'))
    return pd.DataFrame(rows, columns=['seq', 's', 'e', 'gene'])


def genes_per_block(t, cds):
    """Set of genes whose CDS overlaps each block (reference side)."""
    res = [set() for _ in range(len(t))]
    maxlen = int((cds.e - cds.s).max())
    for seq, g in cds.groupby('seq'):
        g = g.sort_values('s')
        s, e, gn = g.s.to_numpy(), g.e.to_numpy(), g.gene.to_numpy()
        idx = np.flatnonzero(t.R.to_numpy() == seq)
        for k in idx:
            a, b = t.rs.iat[k], t.re.iat[k]
            lo, hi = np.searchsorted(s, a - maxlen), np.searchsorted(s, b)
            sel = e[lo:hi] > a
            res[k] = set(gn[lo:hi][sel])
    return res


def neighbours(t, maxgap=10_000, minov=-500):
    """Adjacent blocks along the reference that continue each other in the query."""
    t = t.sort_values(['R', 'rs']).reset_index(drop=True)
    a, b = t.iloc[:-1].reset_index(drop=True), t.iloc[1:].reset_index(drop=True)
    rgap = b.rs - a.re
    qgap = np.where(a.ori > 0, b.qs - a.qe, a.qs - b.qe)
    ok = (a.R == b.R) & (a.Q == b.Q) & (a.ori == b.ori) & rgap.between(minov, maxgap) & \
        (qgap >= minov) & (qgap <= maxgap)
    return t, a[ok], b[ok], rgap[ok].to_numpy(), qgap[ok.to_numpy()]


def chains(t, maxgap):
    t = t.sort_values(['R', 'rs']).reset_index(drop=True)
    a, b = t.iloc[:-1].reset_index(drop=True), t.iloc[1:].reset_index(drop=True)
    rgap = b.rs - a.re
    qgap = np.where(a.ori > 0, b.qs - a.qe, a.qs - b.qe)
    cont = (a.R == b.R) & (a.Q == b.Q) & (a.ori == b.ori) & (rgap <= maxgap) & (qgap >= -20_000) & \
        (qgap <= maxgap)
    return int(len(t) - cont.sum())


def n50(x):
    x = np.sort(np.asarray(x))[::-1]
    return int(x[np.searchsorted(np.cumsum(x), x.sum() / 2)])


def main():
    print(f'python {sys.version.split()[0]}; pandas {pd.__version__}; numpy {np.__version__}')
    thr_rows = []
    seg = pd.read_csv(f'{OUT}/synteny_segments.tsv', sep='\t')
    for comp, (fasta, gff) in COMPS.items():
        print(f'\n################ {comp}')
        raw = load(f'{BASE}/alignments/{comp}/{comp}.promer.one_to_one.L100.coords.tsv')
        t = load(f'{BASE}/results/primary_one_to_one/{comp}.promer.one_to_one.L500.I50.coords.tsv')
        keep = (raw.LEN1 >= 500) & (raw.LEN2 >= 500) & (raw.IDY >= 50)
        assert keep.sum() == len(t), (keep.sum(), len(t))

        # A. thresholds
        print('== A. block counts by minimum length (both sides, bp) x minimum amino-acid identity (%)')
        mlen = raw[['LEN1', 'LEN2']].min(1)
        tab = pd.DataFrame({f'idy>={i}': [int(((mlen >= L) & (raw.IDY >= i)).sum()) for L in LMIN]
                            for i in IMIN}, index=[f'len>={L}' for L in LMIN])
        print(tab.to_string())
        for L in LMIN:
            for i in IMIN:
                thr_rows.append((comp, L, i, int(((mlen >= L) & (raw.IDY >= i)).sum())))
        print(f'input one-to-one blocks (>= 100 bp): {len(raw)}; shortest side < 500 bp: '
              f'{(mlen < 500).mean():.1%}; identity < 50%: {(raw.IDY < 50).mean():.1%}')

        # B. size and identity
        print('== B. retained blocks (>= 500 bp, >= 50% identity)')
        q = t.IDY.quantile([.1, .5, .9]).round(1).tolist()
        print(f'blocks {len(t)}; reference length median {int(t.LEN1.median())} bp, mean {t.LEN1.mean():.0f}, '
              f'N50 {n50(t.LEN1)}, max {t.LEN1.max()}; aa identity 10/50/90% quantiles {q}; '
              f'share >= 90%: {(t.IDY >= 90).mean():.1%}; share of blocks 500-999 bp: '
              f'{(t.LEN1 < 1000).mean():.1%}')

        # C. overlap between retained blocks
        print('== C. overlap of retained blocks')
        for side, sq, s, e in [('reference', 'R', 'rs', 're'), ('query', 'Q', 'qs', 'qe')]:
            iv = t[[sq, s, e]].rename(columns={sq: 'seq', s: 's', e: 'e'})
            mg = merge(iv)
            tot, uni = int((iv.e - iv.s).sum()), int((mg.e - mg.s).sum())
            print(f'{side}: sum of block bp {tot / 1e6:.2f} Mb, merged {uni / 1e6:.2f} Mb '
                  f'(redundancy {tot / uni - 1:.2%}); merged intervals {len(mg)}')

        # D. sequence classes
        print('== D. sequence classes')
        ra, qa = t.R.map(MAJOR), t.Q.map(MAJOR)
        both = ra.notna() & qa.notna()
        homo = both & np.array([HOMOLOG(x, y) for x, y in zip(ra.fillna(''), qa.fillna(''))])
        mito = t.R.isin(MITO) | t.Q.isin(MITO)
        print(f'major-major homologous {homo.sum()} ({homo.mean():.1%}); major-major non-homologous '
              f'{(both & ~homo).sum()}; involving unplaced sequence {(~both & ~mito).sum()} '
              f'({(~both & ~mito).mean():.1%}); involving mtDNA {mito.sum()}')

        # E. coding content (reference side)
        print('== E. coding content, reference side (RefSeq GFF, CDS of all isoforms merged)')
        cds = read_cds(gff)
        cm = merge(cds[['seq', 's', 'e']])
        cov = Cover(cm)
        cbp = cov.bp(t.R.to_numpy(), t.rs.to_numpy(), t.re.to_numpy())
        frac = cbp / (t.re - t.rs).to_numpy()
        blk = merge(t[['R', 'rs', 're']].rename(columns={'R': 'seq', 'rs': 's', 're': 'e'}))
        bcov = Cover(blk)
        cds_in_blocks = bcov.bp(cm.seq.to_numpy(), cm.s.to_numpy(), cm.e.to_numpy()).sum()
        print(f'CDS: {len(cds)} records, {cds.gene.nunique()} genes, {len(cm)} merged intervals, '
              f'{(cm.e - cm.s).sum() / 1e6:.2f} Mb')
        print(f'share of block bp in CDS {cbp.sum() / (t.re - t.rs).sum():.1%}; blocks with no CDS '
              f'{(cbp == 0).mean():.1%}; blocks >= 50% CDS {(frac >= 0.5).mean():.1%}; '
              f'share of all CDS bp inside blocks {cds_in_blocks / (cm.e - cm.s).sum():.1%}')
        gpb = genes_per_block(t, cds)
        ng = np.array([len(g) for g in gpb])
        g2b = pd.Series([g for gs in gpb for g in gs]).value_counts()
        print(f'blocks touching CDS of 0 / 1 / >= 2 genes: {(ng == 0).sum()} / {(ng == 1).sum()} / '
              f'{(ng >= 2).sum()}; genes touched {len(g2b)} of {cds.gene.nunique()}; blocks per touched gene '
              f'median {g2b.median():.0f}, mean {g2b.mean():.2f}, 90% quantile {g2b.quantile(.9):.0f}')

        # F. fragmentation
        print('== F. adjacent blocks that continue each other (gaps -500 bp to 10 kb on both sides)')
        nr, lens = n_runs(fasta)
        nmask = sum(int(e - s) for s, e in zip(nr.s, nr.e))
        print(f'reference hard-masked (N) {nmask / 1e6:.1f} Mb of {sum(lens.values()) / 1e6:.1f} Mb '
              f'({nmask / sum(lens.values()):.1%})')
        ts, a, b, rgap, qgap = neighbours(t)
        ncov = Cover(merge(nr))
        pos = rgap > 0
        nin = np.zeros(len(a), np.int64)
        nin[pos] = ncov.bp(a.R.to_numpy()[pos], a.re.to_numpy()[pos], b.rs.to_numpy()[pos])
        has_nb = np.zeros(len(ts), bool)
        has_nb[a.index] = True
        has_nb[b.index + 1] = True
        print(f'adjacent continuing pairs {len(a)}; blocks with such a neighbour {has_nb.mean():.1%}')
        print(f'reference gap quantiles 10/50/90%: {np.quantile(rgap, [.1, .5, .9]).round(0).tolist()} bp; '
              f'gap <= 0 (touching/overlapping) {(rgap <= 0).mean():.1%}; gap <= 100 bp {(rgap <= 100).mean():.1%}')
        print(f'gaps > 0 containing hard-masked N {(nin[pos] > 0).mean():.1%}; '
              f'gaps > 0 with >= 50% N {(nin[pos] >= 0.5 * rgap[pos]).mean():.1%}; '
              f'reading-frame change between the two blocks {(a.FRM1.to_numpy() != b.FRM1.to_numpy()).mean():.1%}')
        for g in (10_000, 100_000):
            print(f'chains of continuing blocks (gaps <= {g // 1000} kb, no size filter): {chains(t, g)}')

        # G. collinear segments of script 91
        s91 = seg[seg.comp == comp]
        for p, g in s91.groupby('params'):
            print(f'== G. script 91 {p}: {len(g)} collinear segments on homologous major chromosomes, '
                  f'{g.n.sum()} blocks in them (median {g.n.median():.0f} blocks per segment)')
    pd.DataFrame(thr_rows, columns=['comparison', 'min_len_bp', 'min_aa_identity', 'blocks']).to_csv(
        f'{OUT}/promer_block_thresholds.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
