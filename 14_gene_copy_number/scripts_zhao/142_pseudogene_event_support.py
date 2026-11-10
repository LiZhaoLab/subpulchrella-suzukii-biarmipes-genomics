#!/usr/bin/env python3
"""Population support for each disrupting mutation in candidate pseudogene loci.

Question: is each frameshift or premature stop codon that miniprot finds in a disrupted locus (a) a real sequence
difference from the intact orthologous copies rather than an alignment artefact, and (b) present at high frequency
in the population samples?
Replaces the locus-level read check of scripts/138 for the solid-change table (scripts/139, section C). That test
asked only whether reads contradicted the assembly anywhere in the exons; it did not locate the disruption, so a
frameshift created by the alignment passed whenever the reads agreed with the assembly, and a disruption carried by
only part of the population passed when the assembly carried it.

Method, per disrupted locus (Dsub on RU_Dsub_v1.1; D. suzukii on IsoJpt1.0):
1. Events. The miniprot CIGAR (cg tag) of each alignment is walked from the alignment start. Op consumption
   (query aa / target nt), checked on 30,000 alignments against the PAF spans: M n/3n; I n/0; D 0/3n; F 0/n;
   G 1/n; N 0/n; U, V 1/n, where for U and V the n nt include the codon split by the intron (1 + 2 nt for U,
   2 + 1 nt for V); the exon boundaries derived this way equal the miniprot CDS features in 60,000 alignments.
   F and G ops are frameshifts; codons of M ops are translated and in-frame stop codons recorded. The numbers of
   frameshifts and stops found equal the miniprot fs and st tags of every seed alignment (assertion in main).
   Events are taken from the seed alignment of the locus (scripts/137).
2. Artefact filters. An event is excluded if
   - it lies within 10 nt of an intron boundary of the seed alignment (a shifted splice junction creates
     frameshifts and stops);
   - it is a frameshift with another frameshift within 30 nt whose combined length change is a multiple of 3 (the
     frame is restored after a local misalignment);
   - it lies in the last 10% of the query protein (little of the protein is lost), or in the first 10% (a
     downstream start codon could restore most of the protein);
   - another alignment of the locus with a different query and >= 0.9 of the seed alignment score places the
     position in an exon without the same event (a frameshift with the same length change modulo 3 within 60 nt;
     a stop codon at the same position), or places it in an intron (a near-equal gene structure that avoids it);
   - DNA test: the event is not confirmed against the intact orthologous copies. The genomic sequence of the
     disrupted locus (90 nt on each side of the event, in the orientation of the gene) is aligned locally
     (Biopython PairwiseAligner: match 2, mismatch -3, gap open -5, gap extend -2) to the genomic region (+-2 kb)
     of every intact locus of the same D. melanogaster gene in D. biarmipes and in the sister species (per genome,
     the locus with the best score). The comparison is limited to the exon that contains the event, up to 30 nt on
     each side (>= 10 nt required). A frameshift is confirmed by a genome when the two sequences differ in length
     over this segment by a number of nt that is not a multiple of 3; a stop codon is confirmed when the aligned
     triplet of the intact copy is not a stop codon. Alignments with < 70% identity over the segment, or not
     covering it, are not used. The event is confirmed when >= 1 genome is usable and every usable genome
     confirms it. (Genomic rather than coding sequence is used because clean alignments made with a sister-species
     protein that ends at the premature stop do not cover the event.)
   An additional locus-level control in scripts/139: the NCBI model at the locus is full length and of normal
   quality.
3. Population support (samtools mpileup -B, MAPQ >= 20, base quality >= 13; >= 10 reads per pool).
   Frameshift (read CIGARs, samtools view): a read carries another allele if it has an indel whose length is not
   a multiple of 3, or a soft clip of >= 5 nt, within 30 nt of the event (reads with the other allele of a long
   indel or of a structural difference are soft-clipped rather than gapped; an mpileup indel count misses them);
   it carries the assembly allele if it spans the event +- 10 nt without either. Frequency = assembly-allele
   reads / all informative reads. At 200 control exon sites per genome (calibrate()), this estimator gave >= 0.8
   at 99.5-100% of sites in every sample and < 0.9 at up to 8% (soft clips unrelated to an event), so the
   frameshift threshold is 0.8. Stop codon: frequency = lowest reference-base frequency over the three codon
   positions (mpileup; pools carry no haplotype information, so this is an approximation); threshold 0.9.
   Dsub: frequency >= threshold in the Dsub pool (47 females), and the reference-strain reads agree with the
   assembly (frequency >= 0.5 when the reference strain has >= 5 reads; a lower value marks an assembly error).
   D. suzukii: frequency >= threshold in >= 3 of the 4 pools with >= 10 reads, and no pool with >= 10 reads below
   0.5. Frequencies between 0.1 and the threshold are reported as segregating and are not counted.
4. Shared events. For an event in a Dsub locus, the same DNA comparison is made with the genomic region of the
   disrupted loci of the same gene in IsoJpt1.0 (and the reverse): 'same' when the sister sequence has no frame
   difference from it over the segment (length difference a multiple of 3) or the same stop codon at the aligned
   position, i.e. the sister carries the same mutation.
A locus is population-supported if >= 1 of its events passes 2 and 3.
Output (when run as a script): the event table for the loci listed in argv[1], written to argv[2].
samtools 1.24 (envs/sv); Biopython 1.88 (envs/popgen). No random numbers are used.
"""
import importlib.util
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
OUT = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn'
spec = importlib.util.spec_from_file_location('rc', f'{PROJ}/scripts/138_disruption_read_check.py')
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)

GENOME_FA = dict(rc.REF, Dbia=f'{rc.G}/Dbia/genome.fa')  # identical to GCF_025231255.1_RU_DBia_V1.1_genomic.fna (cmp), indexed
OWN = {'Dbia': 'Dbia', 'Dsub': 'Dsub', 'DsuzI': 'DsuzI'}  # query set not used on this genome (as scripts/137)
SISTER = {'Dsub': 'DsuzI', 'DsuzI': 'Dsub'}
MIN_COV, SCOV_FULL, MAX_LEN_RATIO = 0.5, 0.5, 1.5  # as scripts/137
INTRON_NT, COMP_NT, REPL_NT, END_FRAC, ALT_AS = 10, 30, 60, 0.9, 0.9
DNA_FLANK, DNA_SIDE, DNA_MIN_ID, Y_PAD = 90, 30, 0.7, 2000
WIN, MIN_DEPTH, MIN_DEPTH_REF, MAPQ, CLIP_MIN = 30, 10, 5, 20, 5
FIX_STOP, FIX_FS, SEG = 0.9, 0.8, 0.1  # FIX_FS from the control-site calibration (calibrate())
POOLS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
STOPS = {'TAA', 'TAG', 'TGA'}
ALIGNER = PairwiseAligner(mode='local', match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)


# ------------------------------------------------------------------ alignments and loci
def query_tables():
    """dmel_of (query key -> D. melanogaster gene) and the set of full-length queries, as in scripts/137."""
    asg = pd.read_csv(f'{OUT}/C2_gene_assignments.tsv.gz', sep='\t', low_memory=False)
    cols = ['q', 's', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']
    scov, lratio = {}, {}
    for sp in ['Dmel', 'Dbia', 'Dsub', 'DsuzI']:
        h = pd.read_csv(f'{OUT}/{sp}_vs_Dmel.tsv', sep='\t', names=cols)
        h = h[h.evalue <= 1e-5].sort_values(['q', 'bitscore'], ascending=[True, False]).drop_duplicates('q')
        scov.update(dict(zip(h.q, h.length / h.slen)))
        lratio.update(dict(zip(h.q, h.qlen / h.slen)))
    asg['key'] = asg.species + '|' + asg.gene
    asg['scov'] = asg.key.map(scov)
    asg['len_ratio'] = asg.key.map(lratio)
    full = set(asg.key[(asg.scov >= SCOV_FULL) & (asg.len_ratio <= MAX_LEN_RATIO)])
    return dict(zip(asg.key, asg.dmel)), full


def read_alignments(gnm, dmel_of, full):
    rows = []
    with open(f'{OUT}/mp/{gnm}.gff') as f:
        for line in f:
            if not line.startswith('##PAF'):
                continue
            p = line.rstrip('\n').split('\t')[1:]
            tags = {t[:2]: t[5:] for t in p[12:]}
            rows.append((p[0], int(p[1]), int(p[2]), int(p[3]), p[4], p[5], int(p[7]), int(p[8]), int(p[9]),
                         int(p[10]), int(tags.get('AS', 0)), int(tags.get('fs', 0)), int(tags.get('st', 0)),
                         tags['cg']))
    d = pd.DataFrame(rows, columns=['q', 'qlen', 'qs', 'qe', 'strand', 'tname', 'ts', 'te', 'nmatch', 'blen', 'AS',
                                    'fs', 'st', 'cg'])
    d['qsp'] = d.q.str.split('|').str[0]
    d['cov'] = (d.qe - d.qs) / d.qlen
    d = d[(d.qsp != OWN[gnm]) & (d['cov'] >= MIN_COV) & ((d.qsp == 'Dmel') | d.q.isin(full))].copy()
    d['dmel'] = np.where(d.qsp == 'Dmel', d.q.str.split('|').str[1], d.q.map(dmel_of))
    return d.reset_index(drop=True)


def locus_alignments(L, A):
    """Alignments of the locus: same sequence and strand, same D. melanogaster gene, >= 50% overlap (of the shorter
    span) with the seed alignment. The seed itself must be found (checks the reconstruction)."""
    x = A[(A.tname == L.tname) & (A.strand == L.strand) & (A.dmel == L.dmel)]
    ov = np.minimum(x.te, L.te) - np.maximum(x.ts, L.ts)
    x = x[(ov > 0) & (ov >= 0.5 * np.minimum(x.te - x.ts, L.te - L.ts))].copy()
    seed = (x.q == L.seed_q) & (x.ts == L.ts) & (x.te == L.te)
    assert seed.sum() >= 1, f'seed alignment of {L["name"]} not found'
    x['is_seed'] = seed & ~seed.duplicated()
    return x.sort_values(['is_seed', 'AS'], ascending=False)


# ------------------------------------------------------------------ events
def transcript_seq(gnm, a):
    s = rc.ref_seq(GENOME_FA[gnm], a.tname, a.ts, a.te)
    return s.translate(rc.COMP)[::-1] if a.strand == '-' else s


def gpos(a, t0, t1):
    """Genomic 1-based closed interval of transcript offsets [t0, t1)."""
    if a.strand == '+':
        return a.ts + t0 + 1, a.ts + t1
    return a.te - t1 + 1, a.te - t0


class Walk:
    """Events (frameshifts, in-frame stops), exons and spliced coding sequence of one alignment."""

    def __init__(self, a, seq):
        t, q = 0, a.qs
        ev, bnd = [], []
        for k, o in re.findall(r'(\d+)([MIDFGNUV])', a.cg):
            k = int(k)
            if o == 'M':
                for i in range(k):
                    if seq[t + 3 * i:t + 3 * i + 3] in STOPS:
                        ev.append(dict(type='stop', t0=t + 3 * i, t1=t + 3 * i + 3, qpos=q + i, shift=0,
                                       tpos=[t + 3 * i, t + 3 * i + 1, t + 3 * i + 2]))
                t, q = t + 3 * k, q + k
            elif o == 'I':
                q += k
            elif o == 'D':
                t += 3 * k
            elif o in 'FG':
                ev.append(dict(type='frameshift', t0=t, t1=t + k, qpos=q, shift=k if o == 'F' else k - 3))
                t, q = t + k, q + (1 if o == 'G' else 0)
            else:  # N; U, V: the k nt include the split codon (1 + 2 nt for U, 2 + 1 nt for V)
                lead = {'N': 0, 'U': 1, 'V': 2}[o]
                bnd += [t + lead, t + k - (3 - lead if lead else 0)]
                if lead:  # codon split by the intron
                    tp = list(range(t, t + lead)) + list(range(t + k - 3 + lead, t + k))
                    if ''.join(seq[i] for i in tp) in STOPS:
                        ev.append(dict(type='stop', t0=t, t1=t + k, qpos=q, shift=0, tpos=tp))
                t, q = t + k, q + (0 if o == 'N' else 1)
        assert t == a.te - a.ts, f'CIGAR walk does not reach the alignment end ({a.q})'
        cut = [0] + bnd + [t]
        self.toff = [(cut[i], cut[i + 1]) for i in range(0, len(cut), 2)]  # exons, transcript offsets
        self.exons = [gpos(a, x0, x1) for x0, x1 in self.toff]
        self.cds = ''.join(seq[x0:x1] for x0, x1 in self.toff)
        ev.sort(key=lambda e: e['t0'])
        for e in ev:
            e['g0'], e['g1'] = gpos(a, e['t0'], e['t1'])
            if e['type'] == 'stop':
                e['gpos'] = sorted(gpos(a, i, i + 1)[0] for i in e['tpos'])
            e['intron_nt'] = min((min(abs(e['t0'] - b), abs(e['t1'] - b)) for b in bnd), default=10 ** 9)
            e['c'] = self.spliced(e['t0'])
            e['exon_t'] = next(((x0, x1) for x0, x1 in self.toff if x0 <= e['t0'] < x1), (e['t0'], e['t0']))
        fsh = [e for e in ev if e['type'] == 'frameshift']
        for e1, e2 in zip(fsh, fsh[1:]):
            if e2['t0'] - e1['t1'] <= COMP_NT and (e1['shift'] + e2['shift']) % 3 == 0:
                e1['compensated'] = e2['compensated'] = True
        self.events, self.a, self.seq = ev, a, seq

    def spliced(self, t):
        n = 0
        for x0, x1 in self.toff:
            if x0 <= t < x1:
                return n + t - x0
            n += x1 - x0
        return None

    def relation(self, g0, g1):
        """Where a genomic interval lies in this alignment: 'exon', 'intron' or 'outside'."""
        if not (self.a.ts < g0 and g1 <= self.a.te):
            return 'outside'
        return 'exon' if any(x0 <= g0 and g1 <= x1 for x0, x1 in self.exons) else 'intron'


def same_event(e, f):
    if e['type'] != f['type']:
        return False
    if e['type'] == 'stop':
        return e['g0'] == f['g0']
    return abs(e['g0'] - f['g0']) <= REPL_NT and (e['shift'] - f['shift']) % 3 == 0


# ------------------------------------------------------------------ DNA comparison with another copy
def dna_compare(xseq, e, ys):
    """Compare the event of the transcript-oriented genomic sequence xseq (seed alignment span) with the genomic
    sequence ys of another copy. Returns (verdict, detail); verdict 'differs' (the other copy lacks the event),
    'same' (it carries it) or 'unusable'. The comparison is limited to the exon that contains the event."""
    t0, t1 = e['t0'], e['t1']
    lo = max(0, t0 - DNA_FLANK)
    w = xseq[lo:t1 + DNA_FLANK]
    x0, x1 = e['exon_t']
    xa, xb = max(x0, t0 - DNA_SIDE), min(x1, t1 + DNA_SIDE)
    if e['type'] == 'stop' and t1 - t0 != 3:
        return 'unusable', 'stop codon split by an intron'
    if t0 - xa < 10 or xb - t1 < 10:
        return 'unusable', 'fewer than 10 exon nt on one side'
    aln = ALIGNER.align(w, ys)[0]
    bx, by = aln.aligned
    xa, xb, cc = xa - lo, xb - lo, t0 - lo

    def mp(x):
        for (u0, u1), (v0, _) in zip(bx, by):
            if u0 <= x < u1:
                return v0 + x - u0
        return None
    if xa < bx[0][0] or xb > bx[-1][1]:
        return 'unusable', 'local alignment does not cover the exon flanks'
    inb = [x for x in range(xa, xb) if mp(x) is not None]
    ident = sum(w[x] == ys[mp(x)] for x in inb) / (xb - xa)
    if ident < DNA_MIN_ID:
        return 'unusable', f'identity {ident:.2f}'
    if e['type'] == 'frameshift':
        xa_, xb_ = min(inb), max(inb)
        delta = (mp(xb_) - mp(xa_)) - (xb_ - xa_)
        return ('differs' if delta % 3 else 'same'), f'length difference {delta} nt, identity {ident:.2f}'
    y0, y2 = mp(cc), mp(cc + 2)
    if y0 is None or y2 is None or y2 - y0 != 2:
        return 'unusable', 'stop codon not aligned without gaps'
    cod = ys[y0:y0 + 3]
    return ('same' if cod in STOPS else 'differs'), f'aligned triplet {cod}, identity {ident:.2f}'


# ------------------------------------------------------------------ reads
def fs_freq(ref, bam, seqid, g0, g1):
    """Frequency of the assembly allele at a frameshift, from read CIGARs (samtools view, MAPQ >= 20, primary
    alignments without duplicates). A read carries another allele if it has an indel whose length is not a
    multiple of 3, or a soft clip of >= 5 nt, at a position within 30 nt of the event (reads with the other allele
    of a long indel are often soft-clipped rather than gapped). A read carries the assembly allele if it spans the
    event +- 10 nt without such an indel or clip."""
    a, b, ia, ib = g0 - WIN, g1 + WIN, g0 - 10, g1 + 10
    out = subprocess.run([rc.S, 'view', '-q', str(MAPQ), '-F', '0xF04', bam, f'{seqid}:{max(1, a)}-{b}'],
                         capture_output=True, text=True, check=True).stdout
    n_ref = n_alt = 0
    for line in out.splitlines():
        p = line.split('\t', 6)
        r = int(p[3])
        start, alt = r, False
        ops = re.findall(r'(\d+)([MIDNSHP=X])', p[5])
        for i, (k, o) in enumerate(ops):
            k = int(k)
            if o == 'S' and k >= CLIP_MIN and a <= r <= b:
                alt = True
            elif o in 'ID' and k % 3 and a <= r <= b:
                alt = True
            if o in 'MDN=X':
                r += k
        if alt:
            n_alt += 1
        elif start <= ia and r - 1 >= ib:
            n_ref += 1
    n = n_ref + n_alt
    return (n_ref / n if n else np.nan), n


def ref_freq_q(ref, bam, seqid, pos):
    out = subprocess.run([rc.S, 'mpileup', '-B', '-q', str(MAPQ), '-Q', '13', '-d', '100000', '-f', ref, '-r',
                          f'{seqid}:{pos}-{pos}', bam], capture_output=True, text=True, check=True).stdout
    if not out:
        return np.nan, 0
    s = re.sub(r'\^.', '', out.split('\t')[4]).replace('$', '')
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
    return (sum(ch in '.,' for ch in bases) / len(bases) if bases else np.nan), len(bases)


def stop_freq(ref, bam, seqid, positions):
    fr = [ref_freq_q(ref, bam, seqid, p) for p in positions]
    f = [x for x, n in fr if n > 0]
    return (min(f) if len(f) == 3 else np.nan), min(n for _, n in fr)


def population(gnm, e, seqid):
    res = {}
    for smp, bam in rc.BAMS[gnm].items():
        res[f'freq_{smp}'], res[f'depth_{smp}'] = (
            fs_freq(rc.REF[gnm], bam, seqid, e['g0'], e['g1']) if e['type'] == 'frameshift'
            else stop_freq(rc.REF[gnm], bam, seqid, e['gpos']))
    fix = FIX_FS if e['type'] == 'frameshift' else FIX_STOP
    if gnm == 'Dsub':
        f, n = res['freq_Dsub_L5_7'], res['depth_Dsub_L5_7']
        fr, nr = res['freq_refstrain'], res['depth_refstrain']
        if nr >= MIN_DEPTH_REF and fr < 0.5:
            v = 'assembly error (reference-strain reads disagree)'
        elif n < MIN_DEPTH:
            v = 'not testable (pool depth < 10)'
        elif f >= fix:
            v = 'supported'
        elif f > SEG:
            v = 'segregating in the pool'
        else:
            v = 'absent from the pool'
    else:
        ok = [res[f'freq_{p}'] for p in POOLS if res[f'depth_{p}'] >= MIN_DEPTH]
        if sum(f >= fix for f in ok) >= 3 and not any(f < 0.5 for f in ok):
            v = 'supported'
        elif len(ok) < 3:
            v = 'not testable (< 3 pools with depth >= 10)'
        elif all(f <= SEG for f in ok):
            v = 'absent from the pools'
        else:
            v = 'segregating or variable among pools'
    res['population'] = v
    return res


# ------------------------------------------------------------------ per locus
class Context:
    def __init__(self, LO):
        self.LO = LO
        self.dmel_of, self.full = query_tables()
        self.A = {g: read_alignments(g, self.dmel_of, self.full) for g in ['Dbia', 'Dsub', 'DsuzI']}
        self._w, self._r = {}, {}

    def seed_walk(self, L):
        if L['name'] not in self._w:
            X = locus_alignments(L, self.A[L.genome])
            walks = [(a.is_seed, Walk(a, transcript_seq(L.genome, a))) for a in X.itertuples()]
            self._w[L['name']] = walks
        return self._w[L['name']]

    def region(self, M):
        if M['name'] not in self._r:
            s_ = rc.ref_seq(GENOME_FA[M.genome], M.tname, max(0, M.ts - Y_PAD), M.te + Y_PAD)
            self._r[M['name']] = s_.translate(rc.COMP)[::-1] if M.strand == '-' else s_
        return self._r[M['name']]

    def best_copy(self, w0, e, loci):
        """Among loci of one genome, the copy whose genomic region (+-2 kb) aligns best to the event window."""
        w = w0.seq[max(0, e['t0'] - DNA_FLANK):e['t1'] + DNA_FLANK]
        best = None
        for _, M in loci.iterrows():
            sc = ALIGNER.score(w, self.region(M))
            if best is None or sc > best[0]:
                best = (sc, M)
        if best is None:
            return None
        v, det = dna_compare(w0.seq, e, self.region(best[1]))
        return best[1]['name'], v, det


def locus_events(L, ctx):
    walks = ctx.seed_walk(L)
    w0 = [w for s, w in walks if s][0]
    a0 = w0.a
    n_fs = sum(e['type'] == 'frameshift' for e in w0.events)
    n_st = sum(e['type'] == 'stop' for e in w0.events)
    assert (n_fs, n_st) == (a0.fs, a0.st), f'{L["name"]}: events found ({n_fs}, {n_st}) != tags ({a0.fs}, {a0.st})'
    same_gene = ctx.LO[(ctx.LO.dmel == L.dmel) & (ctx.LO['name'] != L['name'])]
    rows = []
    for e in w0.events:
        rel = {'exon_same': [], 'exon_without': [], 'intron': []}
        for s, w in walks:
            if s:
                continue
            r = w.relation(e['g0'], e['g1'])
            if r == 'exon':
                rel['exon_same' if any(same_event(e, f) for f in w.events) else 'exon_without'].append(w.a)
            elif r == 'intron':
                rel['intron'].append(w.a)
        rep = sorted({a.q for a in rel['exon_same'] if a.q != a0.q})
        alt = [a for a in rel['exon_without'] + rel['intron'] if a.q != a0.q and a.AS >= ALT_AS * a0.AS]
        dna, dna_txt = [], []
        for g in ['Dbia', SISTER[L.genome]]:
            x = same_gene[(same_gene.genome == g) & same_gene.intact]
            r = ctx.best_copy(w0, e, x) if len(x) else None
            if r:
                dna.append(r[1])
                dna_txt.append(f'{g} {r[1]} ({r[2]})')
        usable = [v for v in dna if v != 'unusable']
        dna_ok = bool(usable) and all(v == 'differs' for v in usable)
        flt = []
        if e['intron_nt'] < INTRON_NT:
            flt.append(f'within {INTRON_NT} nt of an intron boundary')
        if e.get('compensated'):
            flt.append('compensated by a nearby frameshift')
        if e['qpos'] >= END_FRAC * a0.qlen:
            flt.append('in the last 10% of the protein')
        if e['qpos'] < (1 - END_FRAC) * a0.qlen:
            flt.append('in the first 10% of the protein (a downstream start codon could restore most of it)')
        if alt:
            flt.append(f'another query aligns through the position without the event ({len(alt)} alignment(s) '
                       f'with >= {ALT_AS} of the seed score)')
        if not usable:
            flt.append('DNA test not possible (no usable intact copy)')
        elif not dna_ok:
            flt.append('DNA test: an intact copy carries the same sequence state')
        # shared with the sister species: its disrupted copies of the same gene
        xs_ = same_gene[(same_gene.genome == SISTER[L.genome]) & ~same_gene.intact]
        sh = ctx.best_copy(w0, e, xs_) if len(xs_) else None
        row = dict(genome=L.genome, locus=L['name'], dmel=L.dmel, seqid=L.tname, strand=L.strand, type=e['type'],
                   g0=e['g0'], g1=e['g1'], shift=e['shift'], seed_query=a0.q, seed_AS=a0.AS, qpos=e['qpos'] + 1,
                   qlen=a0.qlen, intron_nt=e['intron_nt'] if e['intron_nt'] < 10 ** 9 else np.nan,
                   reproduced_by_other_query=','.join(rep), n_contradicting_alignments=len(alt),
                   dna_test='; '.join(dna_txt), dna_confirmed=dna_ok,
                   sister_disrupted_copy=f'{sh[0]} {sh[1]} ({sh[2]})' if sh else '',
                   shared_with_sister=(sh[1] == 'same') if sh else False,
                   artefact_filter='; '.join(flt))
        row.update(population(L.genome, e, L.tname))
        row['supported'] = (not flt) and row['population'] == 'supported'
        rows.append(row)
    return rows


def calibrate(n_sites=200, seed=20261106):
    """Frameshift-allele estimator (fs_freq) at control sites: random exon positions (>= 60 nt from the exon ends)
    of intact loci of single-locus genes without frameshift or stop, n_sites per genome, numpy default_rng(seed).
    No frameshifting indel is expected there, so the values show the background of the estimator (soft clips and
    indels unrelated to an event)."""
    rng = np.random.default_rng(seed)
    LO = pd.read_csv(f'{OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
    rows = []
    for gnm in ['Dsub', 'DsuzI']:
        L = LO[(LO.genome == gnm) & LO.intact & (LO.seed_fs == 0) & (LO.seed_st == 0)]
        L = L[L.dmel.map(L.groupby('dmel').size()) == 1]
        bed = pd.read_csv(f'{OUT}/mp/{gnm}_loci_exons.bed', sep='\t', header=None, names=['seqid', 'start', 'end', 'name'])
        bed = bed[bed.name.isin(L.name) & (bed.end - bed.start >= 150)]
        for r in bed.iloc[rng.choice(len(bed), n_sites, replace=False)].itertuples():
            pos = int(rng.integers(r.start + 60, r.end - 60))
            row = dict(genome=gnm, seqid=r.seqid, pos=pos)
            for smp, bam in rc.BAMS[gnm].items():
                row[f'freq_{smp}'], row[f'depth_{smp}'] = fs_freq(rc.REF[gnm], bam, r.seqid, pos, pos)
            rows.append(row)
    return pd.DataFrame(rows)


def run(loci_names):
    LO = pd.read_csv(f'{OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
    ctx = Context(LO)
    rows = []
    for _, L in LO[LO.name.isin(loci_names)].iterrows():
        rows += locus_events(L, ctx)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    E = run([x.strip() for x in open(sys.argv[1]) if x.strip()])
    E.to_csv(sys.argv[2], sep='\t', index=False, float_format='%.3f')
    print(E.drop(columns=['seqid', 'strand', 'dna_test', 'sister_disrupted_copy']).to_string(index=False))
