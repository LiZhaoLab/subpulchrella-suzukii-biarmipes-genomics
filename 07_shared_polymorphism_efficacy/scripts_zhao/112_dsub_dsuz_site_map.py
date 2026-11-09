#!/usr/bin/env python3
"""Direct test of the alpha counting effect, step 1: map every coding site of the
D. subpulchrella site table (scripts/97) to its aligned position in the D. suzukii genome.

Question (Results_A_Efficacy_Of_Selection.md, limitations): are the excess
D. subpulchrella-lineage fixed differences ancestral polymorphisms that are still segregating in
D. suzukii? Testing this needs, for each D. subpulchrella site, the D. suzukii position, which
Step 0 did not store.

Mapping: Airi's PRANK codon alignment Dsub_vs_Dsuz of each ortholog pair (the alignment scripts/60
and 97 use for the D. subpulchrella pool). D. subpulchrella CDS index (3 * codon + q) -> D. suzukii
CDS index -> D. suzukii genome position, through the CDS model of the sister transcript.
Checks: per gene, the ungapped D. suzukii sequence of the alignment equals the genome CDS of the
sister model (genes that fail are left unmapped); per site, the D. suzukii base at the mapped
position (transcript orientation) equals cds_sis of the site table.
Output: QC/site_table/Dsub_L5_7.dsuz_map.npz
  dsuz_chroms (names), dsuz_chrom (index, -1 = unmapped), dsuz_pos (1-based, -1 = unmapped)
  and a QC log results/alpha_shared_map.log
Usage: 112_dsub_dsuz_site_map.py
"""
import importlib.util
import os
import time
from collections import defaultdict

import numpy as np

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('s97', f'{PROJ}/scripts/97_site_table.py')
s97 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s97)
p60 = s97.p60
A, MK, QC = s97.A, s97.MK, s97.QC
COMP = p60.COMP
FOCAL = 'Dsub_L5_7'


def main():
    log = open(f'{PROJ}/results/alpha_shared_map.log', 'w')

    def say(s):
        print(s, flush=True)
        log.write(s + '\n')

    say(f'started {time.ctime()}')
    z = np.load(f'{QC}/site_table/{FOCAL}.npz')
    genes, gene_sis = z['genes'], z['gene_sister']
    cg, cc, cq, csis = z['cds_gene'], z['cds_codon'], z['cds_q'], z['cds_sis']
    regions = f'{A}/population_polymorphism/reanalysis_20260821/config/regions'
    major = {l.split()[0] for cls in ('autosome', 'X') for l in open(f'{regions}/Dsuz.{cls}.bed')}
    models = p60.cds_models(p60.GFF['Dsuz'], major)
    say(f'{FOCAL}: {len(genes)} genes, {len(cg)} coding sites; D. suzukii CDS models on major '
        f'chromosomes: {len(models)}')
    gen = s97.read_genome(p60.FNA['Dsuz'])
    chroms = sorted(major)
    cidx = {c: i for i, c in enumerate(chroms)}
    out_chrom = np.full(len(cg), -1, dtype=np.int16)
    out_pos = np.full(len(cg), -1, dtype=np.int64)
    qc = defaultdict(int)
    order = np.argsort(cg, kind='stable')
    bounds = np.searchsorted(cg[order], np.arange(len(genes) + 1))
    sis_dir = f'{MK}/prank_out/Dsub_vs_Dsuz'
    for gi, tx in enumerate(genes):
        rows = order[bounds[gi]:bounds[gi + 1]]
        if len(rows) == 0:
            continue
        stx = str(gene_sis[gi])
        if stx not in models:
            qc['genes_sister_model_missing'] += 1
            continue
        _, sis_u, sis_ix = s97.anchored_index(f'{sis_dir}/rna-{tx}.best.fas', 'Dsub_0', 'dsuz')
        chrom, pos, strand = models[stx]
        g = ''.join(gen[chrom][p - 1] for p in pos)
        if strand == '-':
            g = ''.join(COMP.get(b, 'N') for b in g)
        su = sis_u.upper()
        if len(su) == len(g) + 3:      # alignment includes the stop codon, the model may not
            su = su[:-3]
        if su != g[:len(su)] or len(su) not in (len(g), len(g) - 3):
            qc['genes_sister_seq_differs_from_genome'] += 1
            continue
        i = 3 * cc[rows].astype(np.int64) + cq[rows]
        j = sis_ix[i]
        ok = (j >= 0) & (j < len(pos))
        qc['genes_mapped'] += 1
        qc['sites_unmapped_gap'] += int((~ok).sum())
        r = rows[ok]
        out_chrom[r] = cidx[chrom]
        out_pos[r] = pos[j[ok]]
        base = np.array(['ACGT'.find(su[k]) for k in j[ok]])
        mism = base != csis[r]
        qc['sites_mapped'] += int(ok.sum())
        qc['sites_base_mismatch_cds_sis'] += int(mism.sum())
        out_chrom[r[mism]] = -1         # keep only sites where the mapped base equals cds_sis
        out_pos[r[mism]] = -1
    mapped = int((out_pos >= 0).sum())
    for k in sorted(qc):
        say(f'  {k}: {qc[k]}')
    say(f'  sites kept (mapped, base equals cds_sis): {mapped} of {len(cg)} ({100 * mapped / len(cg):.2f}%)')
    key = out_chrom.astype(np.int64) * 10 ** 9 + out_pos
    u, n = np.unique(key[out_pos >= 0], return_counts=True)
    say(f'  D. suzukii positions hit by more than one D. subpulchrella site: {int((n > 1).sum())}')
    np.savez_compressed(f'{QC}/site_table/{FOCAL}.dsuz_map.npz', dsuz_chroms=np.array(chroms),
                        dsuz_chrom=out_chrom, dsuz_pos=out_pos)
    say(f'finished {time.ctime()}')


if __name__ == '__main__':
    main()
