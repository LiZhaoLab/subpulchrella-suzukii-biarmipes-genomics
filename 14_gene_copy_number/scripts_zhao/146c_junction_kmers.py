#!/usr/bin/env python3
"""Processed (intronless) copies in D. suzukii pool reads, from exon-exon junction k-mers.

Question: are retrocopies of Sce and betaTub97EF (Chiu et al. 2013, WT3 strain) present in the four Japanese
D. suzukii pools, although neither IsoJpt1.0 nor RU_1.0 contains them?
A processed copy contains the spliced junctions of its parent; the intron-containing parent does not. For each
CDS junction of the parent (IsoJpt1.0 RefSeq models, all transcripts), the k-mers that span the junction with at
least 6 nt on each side are collected (k = 20 and 16), both orientations. Reads (all mapping qualities, primary and
secondary, soft clips included) mapped within +-2 kb of the parent are scanned; a read with a junction k-mer
supports a processed copy. Normalisation: reads containing any exonic k-mer of the same length taken entirely
inside the CDS (number of such reads per junction position = local depth).
Positive controls: the same test on RpL36, CCT1 and VhaAC39-1, whose retrocopies (92.8-97.8% identical to the
parent) are in IsoJpt1.0; their reads are scanned in the retrocopy region with k-mers built from the parent, which
measures sensitivity when the copy has diverged.
samtools 1.24 (envs/sv). No random numbers.
"""
import subprocess
from collections import defaultdict

S = '/lustre/fs8/home/lzhao/vscode_projects/suzukii/envs/sv/bin/samtools'
GI = '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna'
GFF = '/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'
B = '/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup'
POOLS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
COMP = str.maketrans('ACGTN', 'TGCAN')
RETRO = {'RpL36': ('NW_027255939.1', 1918069, 1918413), 'CCT1': ('NW_027255899.1', 1547553, 1549223),
         'VhaAC39-1': ('NW_027255896.1', 24702032, 24703088)}
GENES = ['Sce', 'betaTub97EF', 'RpL36', 'CCT1', 'VhaAC39-1']
PAD = 2000


def rc(s):
    return s.translate(COMP)[::-1]


def fetch(c, a, b):
    out = subprocess.run([S, 'faidx', GI, f'{c}:{a}-{b}'], capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:]).upper()


cds = defaultdict(list)
with open(GFF) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'CDS':
            continue
        at = dict(x.split('=', 1) for x in p[8].split(';') if '=' in x)
        if at.get('gene') in GENES:
            cds[(at['gene'], at['Parent'])].append((p[0], int(p[3]), int(p[4]), p[6]))


def kmers(gene, k):
    """Junction and exonic k-mers over all transcripts of the gene (forward and reverse complement)."""
    junc, exo, region = set(), set(), None
    for (g, tx), ex in cds.items():
        if g != gene:
            continue
        c, strand = ex[0][0], ex[0][3]
        ex = sorted(ex, key=lambda e: e[1], reverse=(strand == '-'))
        segs = [fetch(c, a, b) if strand == '+' else rc(fetch(c, a, b)) for _, a, b, _ in ex]
        seq, bounds = '', []
        for s_ in segs:
            seq += s_
            bounds.append(len(seq))
        for j in bounds[:-1]:
            for i in range(j - k + 6, j - 5):
                junc.add(seq[i:i + k])
        for s_ in segs:
            for i in range(0, len(s_) - k + 1, 7):
                exo.add(s_[i:i + k])
        lo, hi = min(e[1] for e in ex), max(e[2] for e in ex)
        region = (c, lo - PAD, hi + PAD) if region is None else (c, min(region[1], lo - PAD), max(region[2], hi + PAD))
    exo -= junc
    return junc | {rc(x) for x in junc}, exo | {rc(x) for x in exo}, region


def scan(region, junc, exo, k):
    res = {}
    for p in POOLS:
        out = subprocess.run([S, 'view', '-F', '0x4', f'{B}/{p}.rmdup.bam', f'{region[0]}:{region[1]}-{region[2]}'],
                             capture_output=True, text=True, check=True).stdout
        nj = ne = 0
        for line in out.splitlines():
            sq = line.split('\t')[9]
            ks = {sq[i:i + k] for i in range(len(sq) - k + 1)}
            nj += bool(ks & junc)
            ne += bool(ks & exo)
        res[p] = (nj, ne)
    return res


print('gene\tk\tscanned\tjunction_kmers\t' + '\t'.join(f'{p} junction/exonic reads' for p in POOLS))
for gene in GENES:
    for k in (20, 16):
        junc, exo, region = kmers(gene, k)
        r = scan(region, junc, exo, k)
        print(f'{gene}\t{k}\tparent {region[0]}:{region[1]}-{region[2]}\t{len(junc) // 2}\t' +
              '\t'.join(f'{r[p][0]}/{r[p][1]}' for p in POOLS))
        if gene in RETRO:
            c, a, b = RETRO[gene]
            r = scan((c, a - 200, b + 200), junc, exo, k)
            print(f'{gene}\t{k}\tretrocopy {c}:{a}-{b}\t{len(junc) // 2}\t' +
                  '\t'.join(f'{r[p][0]}/{r[p][1]}' for p in POOLS))
