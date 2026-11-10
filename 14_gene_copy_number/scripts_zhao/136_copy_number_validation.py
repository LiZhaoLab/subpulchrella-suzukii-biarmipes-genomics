#!/usr/bin/env python3
"""Validation of the C2 first-pass copy-number calls (scripts/135) with read depth and miniprot.

Question: which of the annotation-based copy-number differences between D. subpulchrella (Dsub) and
D. suzukii (Dsuz) are real, and which come from assembly duplication, collapsed duplicates or missed
annotation?

Two independent lines of evidence:
 1. Read depth (purge_dups logic). Mean depth over the CDS of every annotated protein-coding gene in
    each pool, divided by the median of strict single-copy genes on the same chromosome class
    (autosome / X). Summed over the copies of one D. melanogaster gene this estimates the copy number
    in the population (dCN). An assembly duplicate (two assembled copies of one locus) gives
    about 0.5 + 0.5; a real duplicate 1 + 1; a duplicate collapsed into one assembled copy about 2.
    BAMs: Airi's rmdup BAMs (bwa mem, MAPQ-0 reads kept): Dsub_L5_7 on RU_Dsub_v1.1; JP-Kan, JP-Sap,
    JP-Sap2, JP-Tok on CBGP_Dsuzu_IsoJpt1.0. Reference-strain reads (SRR12476582) on RU_Dsub_v1.1:
    Data_scratch/zhao/suzukii_QC/sv/refstrain.bam (bwa 0.7.18).
 2. miniprot 0.13-r248 alignment of D. melanogaster longest proteins (one query set for every genome,
    so no genome aligns its own annotation) and of D. biarmipes, D. subpulchrella and D. suzukii
    (IsoJpt1.0) proteins to the four assemblies. Loci = overlapping alignments on one strand; each
    locus is assigned to the D. melanogaster gene of its best-scoring alignment and classified as
    intact / disrupted (frameshift or in-frame stop) / partial.

Modes: 'prep' (login node): CDS BED files and miniprot query FASTA. The Slurm job
scripts/136_copy_number_validation.sbatch runs bedcov and miniprot. 'analyze' is in
scripts/137_copy_number_analysis.py.
"""
import os
import re
import sys
from urllib.parse import unquote

import pandas as pd

OUT = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn'
G = '/ru-auth/local/home/asato/scratch/genome'
GFF = {
    'Dsub': f'{G}/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff',
    'DsuzI': f'{G}/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff',
}


def attrs(s):
    return {k: unquote(v) for k, v in (x.split('=', 1) for x in s.rstrip(';').split(';') if '=' in x)}


def cds_bed(sp):
    genes = pd.read_csv(f'{OUT}/{sp}_genes.tsv.gz', sep='\t', low_memory=False)
    genes = genes[genes.coding]
    want = dict(zip(genes.protein, genes.gene))
    assert len(want) == len(genes), 'protein ids not unique'
    iv = {}
    chrom_class = {}
    with open(GFF[sp]) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            if p[2] == 'region':
                a = attrs(p[8])
                chrom_class[p[0]] = a.get('chromosome', 'Unknown')
                continue
            if p[2] != 'CDS':
                continue
            a = attrs(p[8])
            pid = a.get('protein_id')
            if pid in want:
                iv.setdefault(want[pid], []).append((p[0], int(p[3]) - 1, int(p[4])))
    missing = set(want.values()) - set(iv)
    assert not missing, f'{sp}: {len(missing)} genes without CDS'
    rows = []
    for g, ivs in iv.items():
        seqs = {s for s, _, _ in ivs}
        assert len(seqs) == 1, (g, seqs)
        ivs = sorted(ivs, key=lambda x: x[1])
        merged = []
        for s, a, b in ivs:  # merge overlaps (alternative CDS records of one protein should not overlap)
            if merged and a <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], b)
            else:
                merged.append([s, a, b])
        rows += [(s, a, b, g) for s, a, b in merged]
    bed = pd.DataFrame(rows, columns=['seqid', 'start', 'end', 'gene']).sort_values(['seqid', 'start'])
    bed.to_csv(f'{OUT}/{sp}_cds.bed', sep='\t', header=False, index=False)
    pd.Series(chrom_class, name='chromosome').rename_axis('seqid').to_csv(f'{OUT}/{sp}_chrom_class.tsv', sep='\t')
    print(f'{sp}: {bed.gene.nunique()} genes, {len(bed)} CDS intervals, {(bed.end - bed.start).sum() / 1e6:.1f} Mb; '
          f'chromosome classes: {pd.Series(chrom_class).value_counts().to_dict()}', flush=True)


def queries():
    n = 0
    with open(f'{OUT}/mp_queries.faa', 'w') as fo:
        for sp in ['Dmel', 'Dbia', 'Dsub', 'DsuzI']:
            with open(f'{OUT}/{sp}_longest.faa') as f:
                for line in f:
                    fo.write(line)
                    n += line.startswith('>')
    print(f'miniprot queries: {n} proteins', flush=True)


if __name__ == '__main__':
    assert sys.argv[1] == 'prep'
    for sp in GFF:
        cds_bed(sp)
    queries()
