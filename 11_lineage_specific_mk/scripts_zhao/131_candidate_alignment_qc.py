#!/usr/bin/env python3
"""Alignment QC of the lineage candidate genes of script 130 (results/lineage_candidates.tsv).

Why: CG30427 had 266 D. suzukii-lineage fixed differences in 471 codons, and raptor (a conserved TOR
complex gene) 71 nonsynonymous ones; misaligned or misannotated exons produce such counts. Both are
consistent between the D. subpulchrella- and D. suzukii-focal alignments, so the check must look at
the alignment itself.
Method: reference sequences only (no pool data), from Airi's Dsub-focal PRANK codon alignments
(prank_out/Dsub_vs_Dsuz and Dsub_vs_Dbia, joined through the D. subpulchrella reference Dsub_0).
Comparable codon = all three references have no gap and only ACGT. Each differing position is
assigned by parsimony as in script 60 (D. subpulchrella == D. biarmipes -> D. suzukii lineage;
D. suzukii == D. biarmipes -> D. subpulchrella lineage; else skipped); a codon change is
nonsynonymous if the amino acid of the derived codon differs from that of the codon with the
ancestral base(s) restored.
Per gene and lineage: number of nonsynonymous codon changes (reference level); clustering = largest
number of these changes in any 30-codon window, compared with 2,000 uniform placements over the
comparable codons (P_cluster, seed 20261104); indel proximity = fraction of the changes within 3 codons
of a gap in either pairwise alignment, and the same fraction for all comparable codons.
Flags: clustered (P_cluster < 0.01 and >= 5 changes); near_indels (>= 50% of changes near a gap and
>= 3 times the background fraction). Focal lineage: D. subpulchrella for Dsub_lineage, D. suzukii for
Dsuz_lineage, both lineages for 'both' (flag if either is flagged).
Output: results/lineage_candidates_alignment_qc.tsv and a log.
"""
import numpy as np
import pandas as pd

MK = '/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/prank_out'
RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
W, NPERM, SEED = 30, 2000, 20261104
BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}


def read_fas(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def anchored(path, other):
    s = read_fas(path)
    ref, oth = s['Dsub_0'], s[other]
    assert len(ref) == len(oth)
    keep = [i for i, b in enumerate(ref) if b != '-']
    r = ''.join(ref[i] for i in keep)
    o = ''.join(oth[i] for i in keep)
    # gap flag per anchor position: gap in other, or an insertion in other right after this position
    gap = np.zeros(len(keep), bool)
    for j, i in enumerate(keep):
        nxt = keep[j + 1] if j + 1 < len(keep) else len(ref)
        if oth[i] == '-' or nxt - i > 1:
            gap[j] = True
    return r, o, gap


def changes(tx):
    r1, suz, g1 = anchored(f'{MK}/Dsub_vs_Dsuz/rna-{tx}.best.fas', 'dsuz')
    r2, bia, g2 = anchored(f'{MK}/Dsub_vs_Dbia/rna-{tx}.best.fas', 'dbia')
    assert r1 == r2 and len(r1) % 3 == 0, tx
    n = len(r1) // 3
    gapc = (g1 | g2).reshape(n, 3).any(axis=1)
    comp = np.zeros(n, bool)
    ch = {'sub': [], 'suz': []}
    for c in range(n):
        s, z, b = r1[3 * c:3 * c + 3], suz[3 * c:3 * c + 3], bia[3 * c:3 * c + 3]
        if gapc[c] or any(x not in 'ACGT' for x in s + z + b):
            continue
        comp[c] = True
        anc_s, anc_z = list(s), list(z)
        hit = {'sub': False, 'suz': False}
        for k in range(3):
            if s[k] == z[k]:
                continue
            if z[k] == b[k]:
                anc_s[k] = z[k]
                hit['sub'] = True
            elif s[k] == b[k]:
                anc_z[k] = s[k]
                hit['suz'] = True
        if hit['sub'] and CODE[s] != CODE[''.join(anc_s)]:
            ch['sub'].append(c)
        if hit['suz'] and CODE[z] != CODE[''.join(anc_z)]:
            ch['suz'].append(c)
    near = np.zeros(n, bool)
    for c in np.flatnonzero(gapc):
        near[max(0, c - 3):c + 4] = True
    return n, comp, near, ch


def max_window(pos, n):
    if len(pos) == 0:
        return 0
    cnt = np.bincount(pos, minlength=n)
    cs = np.concatenate([[0], np.cumsum(cnt)])
    w = min(W, n)
    return int((cs[w:] - cs[:-w]).max())


def main():
    cand = pd.read_csv(f'{RES}/lineage_candidates.tsv', sep='\t')
    rng = np.random.default_rng(SEED)
    rows = []
    for r in cand.itertuples():
        n, comp, near, ch = changes(r.dsub_tx)
        idx = np.flatnonzero(comp)
        row = dict(dsub_tx=r.dsub_tx, symbol=r.symbol, candidate_class=r.candidate_class, codons=n,
                   comparable=len(idx), bg_near_gap=float(near[idx].mean()) if len(idx) else np.nan)
        for lin in ('sub', 'suz'):
            pos = np.array(ch[lin], int)
            k = len(pos)
            obs = max_window(pos, n)
            if k >= 2 and len(idx) >= k:
                null = np.array([max_window(np.sort(rng.choice(idx, k, replace=False)), n) for _ in range(NPERM)])
                p = float(((null >= obs).sum() + 1) / (NPERM + 1))
            else:
                p = np.nan
            fn = float(near[pos].mean()) if k else np.nan
            row.update({f'ref_Dn_{lin}': k, f'maxwin_{lin}': obs, f'P_cluster_{lin}': p, f'near_gap_{lin}': fn})
            fl = []
            if k >= 5 and p < 0.01:
                fl.append('clustered')
            if k >= 2 and fn >= 0.5 and fn >= 3 * max(row['bg_near_gap'], 1e-9):
                fl.append('near_indels')
            row[f'aln_flags_{lin}'] = ','.join(fl)
        focal = {'Dsub_lineage': ['sub'], 'Dsuz_lineage': ['suz'], 'both': ['sub', 'suz']}[r.candidate_class]
        row['aln_flags'] = ','.join(sorted({f for lin in focal for f in row[f'aln_flags_{lin}'].split(',') if f}))
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(f'{RES}/lineage_candidates_alignment_qc.tsv', sep='\t', index=False)
    m = cand.merge(out[['dsub_tx', 'ref_Dn_sub', 'ref_Dn_suz', 'aln_flags']], on='dsub_tx')
    with open(f'{RES}/lineage_candidates_alignment_qc.log', 'w') as log:
        def p(s=''):
            print(s)
            log.write(s + '\n')
        p(f'Alignment QC, window {W} codons, {NPERM} permutations, seed {SEED}; candidates {len(out)}')
        p(f'Reference-level vs pool-level focal Dn (Spearman): sub {m.ref_Dn_sub.corr(m.Dn_sub, method="spearman"):.3f}, '
          f'suz {m.ref_Dn_suz.corr(m.Dn_suz, method="spearman"):.3f}')
        for c in ('Dsub_lineage', 'Dsuz_lineage', 'both'):
            x = m[m.candidate_class == c]
            p(f'{c}: {len(x)}; alignment-flagged {int((x.aln_flags != "").sum())} '
              f'(clustered {int(x.aln_flags.str.contains("clustered").sum())}, near_indels '
              f'{int(x.aln_flags.str.contains("near_indels").sum())}); robust before {int(x.robust.sum())}, '
              f'robust and not flagged {int((x.robust & (x.aln_flags == "")).sum())}')
        p('\nFlagged genes:')
        f = out[out.aln_flags != ''].merge(cand[['dsub_tx', 'robust']], on='dsub_tx')
        p(f[['symbol', 'candidate_class', 'robust', 'codons', 'ref_Dn_sub', 'maxwin_sub', 'P_cluster_sub', 'near_gap_sub',
             'ref_Dn_suz', 'maxwin_suz', 'P_cluster_suz', 'near_gap_suz', 'bg_near_gap', 'aln_flags']]
          .to_string(index=False, float_format=lambda v: f'{v:.3g}'))


if __name__ == '__main__':
    main()
