#!/usr/bin/env python3
"""Read check of disrupted (pseudogene) miniprot loci in candidate gene families.

Question: is a frameshift or premature stop codon found by miniprot in the Dsub or D. suzukii (IsoJpt1.0)
assembly real in the population, or an error in the assembly sequence?
Test: samtools mpileup over the exons of the locus (seed alignment of scripts/137) in the pool(s) mapped to that
assembly (Dsub: Dsub_L5_7 pool and the reference strain; IsoJpt1.0: JP-Kan, JP-Sap, JP-Sap2, JP-Tok). At each
position: fraction of reads with an indel starting there, and fraction with a base other than the reference.
An assembly error at a disrupting site appears as a position where most reads (>= 0.5) disagree with the
reference. If no exon position has a majority non-reference allele, the reference sequence, including its
disruption, is the majority allele in the population. Minor alleles (0.1-0.5) at indels are reported as possible
polymorphic disruptions.
Loci: disrupted loci (no intact alignment) of chemoreceptor, detoxification, mechanosensory and cuticle genes in
the Dsub and IsoJpt1.0 genomes. Output: results/C2_disruption_read_check.tsv
samtools 1.24 (envs/sv); mpileup -B -q 0 -Q 13 -d 100000. No random numbers are used.
Functions moved to module level and the run put under a main guard, so that
scripts/139_solid_copy_number_changes.py can apply the same checks to other loci; outputs unchanged.
"""
import re
import subprocess
from collections import Counter

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
OUT = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn'
S = f'{PROJ}/envs/sv/bin/samtools'
G = '/ru-auth/local/home/asato/scratch/genome'
B = '/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup'
REF = {'Dsub': f'{G}/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna',
       'DsuzI': f'{G}/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna'}
BAMS = {'Dsub': {'Dsub_L5_7': f'{B}/Dsub_L5_7.rmdup.bam',
                 'refstrain': '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sv/refstrain.bam'},
        'DsuzI': {p: f'{B}/{p}.rmdup.bam' for p in ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']}}
FOCUS = re.compile(r'^(Or\d|Gr\d|Ir\d|Obp|ppk|Cyp|Gst|Ugt|alpha-Est|Est-|Cpr|Twdl|Lcp|nompC|iav|nan$|pain|pyx|'
                   r'TrpA1|trp|Trp|wtrw|Pkd2|brv|Piezo|Tmc|nompA)')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CODE = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}


def load_loci():
    LO = pd.read_csv(f'{OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
    dm = pd.read_csv(f'{OUT}/Dmel_genes.tsv.gz', sep='\t', low_memory=False)
    LO['symbol'] = LO.dmel.map(dict(zip(dm.gene, dm.symbol))).astype(str)
    return LO


def load_bed():
    bed = {}
    for gnm in ('Dsub', 'DsuzI'):
        b = pd.read_csv(f'{OUT}/mp/{gnm}_loci_exons.bed', sep='\t', header=None,
                        names=['seqid', 'start', 'end', 'name'])
        bed.update({k: v for k, v in b.groupby('name')})
    return bed


def pileup(ref, bam, regions):
    """Per position: depth, indel-start fraction, non-reference base fraction."""
    rows = []
    for seqid, a, b in regions:
        out = subprocess.run([S, 'mpileup', '-B', '-q', '0', '-Q', '13', '-d', '100000', '-f', ref,
                              '-r', f'{seqid}:{a + 1}-{b}', bam], capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            p = line.split('\t')
            depth, s = int(p[3]), p[4]
            if depth == 0:
                rows.append((int(p[1]), 0, np.nan, np.nan))
                continue
            s = re.sub(r'\^.', '', s).replace('$', '')
            n_indel = len(re.findall(r'[+-]\d+', s))
            s = re.sub(r'[+-](\d+)', lambda m: '\x00' + m.group(1), s)
            # remove inserted/deleted sequence after the length marker
            t, i = [], 0
            while i < len(s):
                if s[i] == '\x00':
                    j = i + 1
                    while j < len(s) and s[j].isdigit():
                        j += 1
                    i = j + int(s[i + 1:j])
                    continue
                t.append(s[i])
                i += 1
            bases = ''.join(t)
            nonref = sum(ch in 'ACGTNacgtn' for ch in bases)
            rows.append((int(p[1]), depth, n_indel / depth, nonref / depth))
    return pd.DataFrame(rows, columns=['pos', 'depth', 'indel', 'nonref'])


def check_loci(sel, bed):
    """Per locus and sample: exon pileup summary (indel and base disagreement with the reference)."""
    res = []
    for r in sel.itertuples():
        ex = bed.get(r.name)
        if ex is None:
            continue
        regions = list(zip(ex.seqid, ex.start, ex.end))
        for smp, bam in BAMS[r.genome].items():
            x = pileup(REF[r.genome], bam, regions)
            res.append(dict(genome=r.genome, locus=r.name, symbol=r.symbol, tname=r.tname, ts=r.ts, te=r.te,
                            seed_q=r.seed_q, seed_cov=r.seed_cov, seed_fs=r.seed_fs, seed_st=r.seed_st, sample=smp,
                            exon_bp=len(x), median_depth=x.depth.median(),
                            max_indel=x.indel.max(), max_nonref=x.nonref.max(),
                            n_majority_indel=int((x.indel >= 0.5).sum()),
                            n_majority_snp=int((x.nonref >= 0.5).sum()),
                            n_minor_indel=int(((x.indel >= 0.1) & (x.indel < 0.5)).sum())))
    return pd.DataFrame(res)


# Assembly-error test. Dsub: the reference-strain reads come from the assembled strain, so any position where
# most of them disagree (indel or base) marks an error in the assembly sequence. IsoJpt1.0: no reads from the
# assembled strain; a position where most reads of all four pools carry an indel marks either an assembly error or
# a minority allele in the assembled strain. Majority SNPs in pools are expected (the reference strain carries the
# minor allele at some polymorphic sites) and are reported but not used.
def summarise(R):
    pool = R[R['sample'] != 'refstrain'].groupby(['genome', 'locus', 'symbol']).agg(
        pool_majority_indel_min=('n_majority_indel', 'min'), pool_majority_indel_max=('n_majority_indel', 'max'),
        pool_majority_snp=('n_majority_snp', 'median'), pool_minor_indel=('n_minor_indel', 'max'),
        depth=('median_depth', 'median'), seed_fs=('seed_fs', 'first'), seed_st=('seed_st', 'first'),
        seed_cov=('seed_cov', 'first')).reset_index()
    ref = R[R['sample'] == 'refstrain'].set_index('locus')
    pool['refstrain_majority_indel'] = pool.locus.map(ref.n_majority_indel)
    pool['refstrain_majority_snp'] = pool.locus.map(ref.n_majority_snp)

    def verdict(r):
        if r.genome == 'Dsub':
            if r.refstrain_majority_indel > 0 or r.refstrain_majority_snp > 0:
                return 'assembly error at >= 1 position (reference-strain reads disagree)'
            return 'assembly sequence confirmed by reference-strain reads' + (
                '; pool majority indel elsewhere' if r.pool_majority_indel_min > 0 else '')
        if r.pool_majority_indel_min > 0:
            return 'majority indel in all four pools (assembly error or minority allele in the assembled strain)'
        return 'no majority indel in any pool' if r.pool_majority_indel_max == 0 else 'majority indel in some pools'
    pool['verdict'] = pool.apply(verdict, axis=1)
    return pool


# Premature stop codons without frameshift: translate the seed exons, locate each internal stop codon, and estimate
# the frequency of the stop allele in each sample as the lowest reference-base frequency over the three codon
# positions (a stop is lost if any of its bases changes; pools carry no haplotype information, so this is an
# approximation).
def ref_seq(ref, seqid, a, b):
    out = subprocess.run([S, 'faidx', ref, f'{seqid}:{a + 1}-{b}'], capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:]).upper()


def ref_freq(ref, bam, seqid, pos):
    out = subprocess.run([S, 'mpileup', '-B', '-q', '0', '-Q', '13', '-d', '100000', '-f', ref, '-r',
                          f'{seqid}:{pos}-{pos}', bam], capture_output=True, text=True, check=True).stdout
    if not out:
        return np.nan, 0
    p = out.split('\t')
    s = re.sub(r'\^.', '', p[4]).replace('$', '')
    s = re.sub(r'[+-](\d+)', lambda m: '#' + m.group(1), s)
    t, i = [], 0
    while i < len(s):
        if s[i] == '#':
            j = i + 1
            while j < len(s) and s[j].isdigit():
                j += 1
            i = j + int(s[i + 1:j])
            continue
        t.append(s[i])
        i += 1
    bases = [ch for ch in t if ch in '.,ACGTNacgtn*']
    n = len(bases)
    return (sum(ch in '.,' for ch in bases) / n if n else np.nan), n


def stop_freqs(sel, bed):
    stops = []
    for r in sel[(sel.seed_fs == 0) & (sel.seed_st > 0)].itertuples():
        ex = bed.get(r.name)
        if ex is None:
            continue
        ex = ex.sort_values('start')
        pos = []  # genomic positions (1-based) of the CDS in transcription order
        seq = ''
        for e in ex.itertuples():
            seq += ref_seq(REF[r.genome], e.seqid, e.start, e.end)
            pos += list(range(e.start + 1, e.end + 1))
        if r.strand == '-':
            seq = seq.translate(COMP)[::-1]
            pos = pos[::-1]
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
        for k, cdn in enumerate(codons[:-1]):  # internal stops only
            if CODE.get(cdn) != '*':
                continue
            cpos = pos[3 * k:3 * k + 3]
            row = dict(genome=r.genome, locus=r.name, symbol=r.symbol, seqid=r.tname, codon_index=k,
                       n_codons=len(codons), stop_codon=cdn, positions=','.join(map(str, cpos)))
            for smp, bam in BAMS[r.genome].items():
                fr = [ref_freq(REF[r.genome], bam, r.tname, p_) for p_ in cpos]
                row[f'stop_freq_{smp}'] = min(f for f, _ in fr)
                row[f'depth_{smp}'] = int(np.median([n for _, n in fr]))
            stops.append(row)
    return pd.DataFrame(stops)


def main():
    LO = load_loci()
    sel = LO[LO.genome.isin(['Dsub', 'DsuzI']) & ~LO.intact & LO.symbol.str.match(FOCUS)].copy()
    bed = load_bed()
    R = check_loci(sel, bed)
    R.to_csv(f'{PROJ}/results/C2_disruption_read_check.tsv', sep='\t', index=False, float_format='%.3f')
    pool = summarise(R)
    pool.to_csv(f'{PROJ}/results/C2_disruption_read_check_summary.tsv', sep='\t', index=False,
                float_format='%.1f')
    print(f'disrupted focus loci checked: {len(pool)} (Dsub {int((pool.genome == "Dsub").sum())}, IsoJpt1.0 '
          f'{int((pool.genome == "DsuzI").sum())})')
    print(dict(Counter(zip(pool.genome, pool.verdict))))
    pd.set_option('display.width', 250)
    print(pool.sort_values(['genome', 'symbol']).to_string(index=False))
    ST = stop_freqs(sel, bed)
    ST.to_csv(f'{PROJ}/results/C2_premature_stop_frequencies.tsv', sep='\t', index=False,
              float_format='%.3f')
    print('\nPremature stop codons (no frameshift) in disrupted focus loci: estimated stop-allele frequency per sample')
    print(ST.drop(columns=['positions', 'seqid']).to_string(index=False))


if __name__ == '__main__':
    main()
