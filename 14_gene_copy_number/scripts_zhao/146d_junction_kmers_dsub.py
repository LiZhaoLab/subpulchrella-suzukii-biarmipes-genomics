#!/usr/bin/env python3
"""Are the D. suzukii retrocopies (RpL36, RpS14, CCT1/T-cp1, VhaAC39-1) also present in D. subpulchrella reads?

Same junction k-mer test as scripts/146c, on the D. subpulchrella pool (Dsub_L5_7, 47 females) and reference-strain
reads (SRR12476582), both mapped to RU_Dsub_v1.1. Junction k-mers come from the D. subpulchrella parent genes
(RefSeq CDS of every transcript whose CDS overlaps the parent locus of scripts/137). Reads of a retrocopy missing
from the assembly would map to the parent with soft clips at exon ends, so reads mapped within +-2 kb of the parent
are scanned (all mapping qualities, unmapped reads excluded).
Positive controls: two D. subpulchrella-lineage retrocopies of the solid table (Cox17, RpLP1), whose reads are
scanned in the retrocopy region with k-mers from the parent. samtools 1.24 (envs/sv). No random numbers.
"""
import subprocess

S = '/lustre/fs8/home/lzhao/vscode_projects/suzukii/envs/sv/bin/samtools'
G = '/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1'
FA = f'{G}/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna'
GFF = f'{G}/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff'
BAMS = {'Dsub_pool': '/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup/Dsub_L5_7.rmdup.bam',
        'Dsub_refstrain': '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sv/refstrain.bam'}
COMP = str.maketrans('ACGTN', 'TGCAN')
PAD = 2000
# label: (parent seqid, start, end, strand), retrocopy region or None (1-based, from scripts/137 loci)
GENES = {'RpL36': (('NC_050613.1', 11819012, 11819446, '-'), None),
         'RpS14': (('NC_050613.1', 10794984, 10796672, '-'), None),
         'CCT1': (('NC_050609.1', 10041994, 10044024, '-'), None),
         'VhaAC39-1': (('NC_050613.1', 15207826, 15209829, '+'), None),
         'control Cox17': (('NC_050613.1', 9630941, 9631250, None), ('NW_023665645.1', 505151, 505387)),
         'control RpLP1': (('NC_050610.1', 25610253, 25610736, None), ('NW_023665596.1', 770246, 770520))}


def rc(s):
    return s.translate(COMP)[::-1]


def fetch(c, a, b):
    out = subprocess.run([S, 'faidx', FA, f'{c}:{a}-{b}'], capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:]).upper()


def transcripts(c, a, b):
    tx = {}
    with open(GFF) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'CDS' or p[0] != c:
                continue
            at = dict(x.split('=', 1) for x in p[8].split(';') if '=' in x)
            tx.setdefault(at['Parent'], []).append((int(p[3]), int(p[4]), p[6]))
    return {t: e for t, e in tx.items() if min(x[0] for x in e) <= b and max(x[1] for x in e) >= a}


def kmers(c, a, b, k):
    junc, exo = set(), set()
    for t, ex in transcripts(c, a, b).items():
        strand = ex[0][2]
        ex = sorted(ex, key=lambda e: e[0], reverse=(strand == '-'))
        segs = [fetch(c, x, y) if strand == '+' else rc(fetch(c, x, y)) for x, y, _ in ex]
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
    exo -= junc
    return junc | {rc(x) for x in junc}, exo | {rc(x) for x in exo}


def scan(c, a, b, junc, exo, k):
    res = {}
    for smp, bam in BAMS.items():
        out = subprocess.run([S, 'view', '-F', '0x4', bam, f'{c}:{a}-{b}'], capture_output=True, text=True,
                             check=True).stdout
        nj = ne = 0
        for line in out.splitlines():
            sq = line.split('\t')[9]
            ks = {sq[i:i + k] for i in range(len(sq) - k + 1)}
            nj += bool(ks & junc)
            ne += bool(ks & exo)
        res[smp] = f'{nj}/{ne}'
    return res


print('gene\tk\tscanned\tn_junction_kmers\t' + '\t'.join(f'{s} junction/exonic reads' for s in BAMS))
for g, ((c, a, b, _), retro) in GENES.items():
    for k in (20, 16):
        junc, exo = kmers(c, a, b, k)
        r = scan(c, a - PAD, b + PAD, junc, exo, k)
        print(f'{g}\t{k}\tparent {c}:{a}-{b}\t{len(junc) // 2}\t' + '\t'.join(r.values()))
        if retro:
            r = scan(retro[0], retro[1] - 200, retro[2] + 200, junc, exo, k)
            print(f'{g}\t{k}\tretrocopy {retro[0]}:{retro[1]}-{retro[2]}\t{len(junc) // 2}\t' + '\t'.join(r.values()))
