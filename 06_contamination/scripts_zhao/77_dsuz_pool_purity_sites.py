#!/usr/bin/env python3
"""Purity of a D. suzukii pool (analysis 4.2): pool read counts at species-specific derived sites.

Mirror of the D. subpulchrella-pool check in 60_polarized_mk.py (contamination block) and
65_summarize_contamination.py, with the roles of the species exchanged. For each D. suzukii pool:
at callable CDS positions (depth 20..P99 of that pool, outside indel masks) of gene trios with a
closed RBH triangle, the D. suzukii reference base r, the D. subpulchrella base y and the
D. biarmipes base o come from Airi's pool-specific PRANK alignments (Dsuz_vs_Dsub, Dsuz_vs_Dbia),
joined through the D. suzukii reference CDS (checked against the genome, as in script 60), and the
D. melanogaster base m from Dsuz_vs_Dmel. Site classes:
  dsub : y != r and o == r  (D. subpulchrella-specific derived base; the contaminant signal)
  dbia : o != r and y == r  (control: contamination by D. biarmipes is not plausible)
  dmel : m != r, y == r, o == r  (control)
Output: results/purity_dsuz/{POOL}.sites.tsv.gz with species, chrom, pos, depth, n_focal_allele
(reads with r), n_species_allele (reads with the other species' base), n_other_two_bases.
Usage: 77_dsuz_pool_purity_sites.py POOL   (JP-Kan | JP-Sap | JP-Sap2 | JP-Tok); TEST_CHROM=... for one
chromosome.
"""
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('pm', f'{PROJ}/scripts/60_polarized_mk.py')
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)
cs = pm.cs


def main():
    focal = sys.argv[1]
    assert focal in ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok'), focal
    base = f'{pm.MK}/prank_out/Dsuz/{focal}'
    anchor = f'{focal}_0'
    regions = f'{pm.A}/population_polymorphism/reanalysis_20260821/config/regions'
    major = {l.split()[0] for cls in ('autosome', 'X') for l in open(f'{regions}/Dsuz.{cls}.bed')}
    tri, _ = pm.orthologs()
    tri = tri['Dsuz']
    models = pm.cds_models(pm.GFF['Dsuz'], major)
    thr = pd.read_csv(f'{PROJ}/results/depth_thresholds.tsv', sep='\t')
    thr = thr[thr['sample'] == focal].set_index('cls')
    qc = defaultdict(int)
    by_chrom = defaultdict(list)
    for tx in tri:
        if tx in models:
            by_chrom[models[tx][0]].append(tx)
    rows = []
    code = {b'A': 0, b'C': 1, b'G': 2, b'T': 3}
    isb = lambda v: np.isin(v, [b'A', b'C', b'G', b'T'])
    test_chrom = os.environ.get('TEST_CHROM')
    for chrom, txs in sorted(by_chrom.items()):
        if test_chrom and chrom != test_chrom:
            continue
        cls = 'X' if chrom in cs.XCHROM else 'autosome'
        dmax = int(thr.loc[cls, 'p99'])
        cnt, indel = pm.load_chrom(focal, chrom)
        gseq = pm.chrom_seq(pm.FNA['Dsuz'], chrom)
        assert len(gseq) == len(cnt), (chrom, len(gseq), len(cnt))
        for tx in txs:
            f_s, f_o, f_m = (f'{base}/Dsuz_vs_{x}/rna-{tx}.best.fas' for x in ('Dsub', 'Dbia', 'Dmel'))
            if not (os.path.exists(f_s) and os.path.exists(f_o)):
                qc['skip_alignment_missing'] += 1
                continue
            ref1, sis = pm.anchored(f_s, anchor, 'dsub')
            ref2, outg = pm.anchored(f_o, anchor, 'dbia')
            dmel = None
            if os.path.exists(f_m):
                ref3, dmel = pm.anchored(f_m, anchor, 'dmel')
                if ref3 != ref1:
                    dmel = None
            _, pos, strand = models[tx]
            if len(pos) == len(ref1) + 3:
                pos = pos[:-3]
            if ref1 != ref2 or len(ref1) != len(pos):
                qc['skip_anchor_mismatch'] += 1
                continue
            gref = ''.join(gseq[p - 1] for p in pos)
            if strand == '-':
                gref = ''.join(pm.COMP.get(b, 'N') for b in gref)
            if gref != ref1.upper():
                qc['skip_anchor_differs_from_genome'] += 1
                continue
            qc['genes_used'] += 1
            c = cnt[pos - 1].astype(np.int64)
            if strand == '-':
                c = c[:, [7, 6, 5, 4, 3, 2, 1, 0]]
            tot = c[:, 0::2] + c[:, 1::2]
            depth = tot.sum(1)
            ok_site = (depth >= pm.MIN_DEPTH) & (depth <= dmax) & ~indel[pos - 1]
            r = np.frombuffer(ref1.upper().encode(), dtype='S1')
            y = np.frombuffer(sis.upper().encode(), dtype='S1')
            o = np.frombuffer(outg.upper().encode(), dtype='S1')
            m = np.frombuffer(dmel.upper().encode(), dtype='S1') if dmel else None
            base_ok = ok_site & isb(r) & isb(y) & isb(o)
            tests = [('dsub', y, base_ok & (y != r) & (o == r)),
                     ('dbia', o, base_ok & (o != r) & (y == r))]
            if m is not None:
                tests.append(('dmel', m, base_ok & isb(m) & (m != r) & (y == r) & (o == r)))
            for label, other, sel in tests:
                for i in np.flatnonzero(sel):
                    ri, oi = code[r[i]], code[other[i]]
                    rows.append((label, chrom, int(pos[i]), int(depth[i]), int(tot[i, ri]),
                                 int(tot[i, oi]), int(depth[i] - tot[i, ri] - tot[i, oi])))
        print(f'{focal} {chrom}: sites so far {len(rows)}', flush=True)
    out = f'{PROJ}/results/purity_dsuz' + ('_test' if test_chrom else '')
    os.makedirs(out, exist_ok=True)
    pd.DataFrame(rows, columns=['species', 'chrom', 'pos', 'depth', 'n_focal_allele', 'n_species_allele',
                                'n_other_two_bases']).to_csv(f'{out}/{focal}.sites.tsv.gz', sep='\t', index=False)
    print(dict(qc))


if __name__ == '__main__':
    main()
