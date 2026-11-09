#!/usr/bin/env python3
"""Script 123: is mitochondrial divergence between D. subpulchrella and D. suzukii unusually low
relative to nuclear divergence? Same design on both sides: the D. subpulchrella pool against the D. suzukii
reference genome, at four-fold degenerate sites.

mtDNA: D. subpulchrella raw-read subsample mapped to the D. suzukii mtDNA NC_060762.1 (script 78; MAPQ >= 20,
base quality >= 20). Four-fold sites = third positions of the 13 protein-coding genes (RefSeq GFF of
GCF_043229965.1) whose first two bases define a four-fold family in the invertebrate mitochondrial code
(NCBI table 5; AGN = Ser is four-fold there). Same measures for the four D. suzukii pools (within-species baseline).
Nuclear: site table of script 97 (Dsub_L5_7.npz): cds_fold == 4 (degeneracy in the ancestral codon), cds_sis = the
aligned D. suzukii reference base (PRANK ortholog alignment).
Measures, at sites with depth >= 20: d = mean over sites of (1 - frequency of the D. suzukii reference base);
fixed = share of sites where the reference base is at <= 5%; pi = mean 1 - sum(p^2) within the pool;
Jukes-Cantor correction of d.
Output: results/mito_vs_nuclear_divergence.log (stdout). Python: anaconda3/envs/jupyter-ml.
"""
import re

import numpy as np

QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
GFF = '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'
MT = 'NC_060762.1'
FOUR = {'CT', 'GT', 'TC', 'CC', 'AC', 'GC', 'CG', 'GG', 'AG'}  # table 5
COMP = str.maketrans('ACGT', 'TGCA')
MIN_DEPTH = 20


def read_fasta(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def pileup_counts(path, L):
    """ACGT counts per position (genomic orientation) from samtools mpileup text."""
    c = np.zeros((L + 1, 4), np.int64)
    for line in open(path):
        f = line.split('\t')
        pos, ref, s = int(f[1]), f[2].upper(), f[4]
        s = re.sub(r'\^.', '', s).replace('$', '')
        while True:  # remove indel sequences
            m = re.search(r'[+-](\d+)', s)
            if not m:
                break
            s = s[:m.start()] + s[m.end() + int(m.group(1)):]
        nref = s.count('.') + s.count(',')
        su = s.upper()
        for i, b in enumerate('ACGT'):
            c[pos, i] = su.count(b) + (nref if b == ref else 0)
    return c


def measures(cnt, refidx):
    dep = cnt.sum(1)
    ok = dep >= MIN_DEPTH
    cnt, refidx, dep = cnt[ok], refidx[ok], dep[ok]
    fref = cnt[np.arange(len(cnt)), refidx] / dep
    p = cnt / dep[:, None]
    d = float(np.mean(1 - fref))
    jc = -0.75 * np.log(1 - 4 * d / 3)
    return dict(sites=int(ok.sum()), d=d, d_JC=float(jc), fixed=float(np.mean(fref <= 0.05)),
                pi=float(np.mean(1 - (p ** 2).sum(1))))


def fmt(m):
    return (f"sites {m['sites']:,}; d {m['d']:.4f} (JC {m['d_JC']:.4f}); fixed differences {m['fixed']:.4f}; "
            f"pi within pool {m['pi']:.4f}")


# ---------------------------------------------------------------- mtDNA four-fold sites
seq = read_fasta(f'{QC}/mito/dsuz_mito.fa')
L = len(seq)
assert L == 16392
four = set()
ncodon = 0
for line in open(GFF):
    if not line.startswith(MT + '\t'):
        continue
    f = line.split('\t')
    if f[2] != 'CDS':
        continue
    s, e, strand = int(f[3]), int(f[4]), f[6]
    assert f[7] == '0'
    cds = seq[s - 1:e] if strand == '+' else seq[s - 1:e][::-1].translate(COMP)
    for k in range(0, len(cds) - 2, 3):
        ncodon += 1
        if cds[k:k + 2] in FOUR:
            q3 = k + 2  # index in transcript orientation -> genomic position
            four.add(s + q3 if strand == '+' else e - q3)
four = np.array(sorted(four))
print(f'mtDNA: {ncodon} codons in 13 CDS; {len(four)} four-fold sites (table 5)')
refidx = np.array(['ACGT'.index(seq[p - 1]) for p in four])
allpos = np.array([p for p in range(1, L + 1) if seq[p - 1] in 'ACGT'])
allref = np.array(['ACGT'.index(seq[p - 1]) for p in allpos])
for pool, path in [('Dsub_L5_7', f'{QC}/mito/Dsub_L5_7.mito.pileup')] + \
        [(p, f'{QC}/mito/{p}.mito.pileup') for p in ('JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok')]:
    c = pileup_counts(path, L)
    print(f'{pool} mtDNA four-fold: ' + fmt(measures(c[four], refidx)))
    if pool == 'Dsub_L5_7':
        m = measures(c[allpos], allref)
        print(f'{pool} mtDNA all positions: ' + fmt(m) +
              f'; fixed-difference positions {round(m["fixed"] * m["sites"])}')

# ---------------------------------------------------------------- nuclear four-fold sites
z = np.load(f'{QC}/site_table/Dsub_L5_7.npz')
chroms = list(z['chroms'])
xi = chroms.index('NC_050613.1')
sel = (z['cds_fold'] == 4) & (z['cds_sis'] >= 0) & (z['cds_sis'] <= 3)
for name, m in [('autosomes', sel & (z['cds_chrom'] != xi)), ('X', sel & (z['cds_chrom'] == xi))]:
    print(f'Dsub_L5_7 nuclear four-fold, {name}: ' +
          fmt(measures(z['cds_tot'][m].astype(np.int64), z['cds_sis'][m].astype(np.int64))))
