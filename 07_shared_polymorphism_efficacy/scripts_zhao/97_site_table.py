#!/usr/bin/env python3
"""Step 0 of the depth plan: per-site table of coding and short-intron sites for one
focal pool, with the ancestral state from D. biarmipes (Depth_Analysis_Plan.md).

Reuses the ortholog triangles, alignment joins, depth and indel filters and the polarization rule
of scripts/60_polarized_mk.py, but stores sites instead of per-gene sums, and stores read counts
without any frequency threshold, so that later analyses choose their own allele filters and
subsample every pool to a common depth.

Coding sites: every position of a callable codon, defined as in script 60: all three positions
with depth 20..P99 and outside indel masks; sister and D. biarmipes bases in ACGT; focal major
codon, sister codon and ancestral codon not stop codons; the major allele passes the script-60
allele filter (>= 3 reads, >= 1 read per strand, frequency > 0.05); every position polarizable
(focal major == sister, or sister == outgroup -> ancestral = sister; focal major == outgroup ->
ancestral = focal major; otherwise the codon is skipped).
Short-intron sites: positions 8-30 from the donor of introns 30-65 bp of the focal transcript that
lie between two CDS segments, restricted to the short-intron mask of scripts/31 (no exon overlap),
when the sister and the D. biarmipes transcripts have an intron at the same aligned CDS position.
Intron sequences (transcript orientation) are aligned pairwise, focal-sister and focal-outgroup
(Biopython PairwiseAligner, global; match 2, mismatch -3, gap open -5, gap extend -2). A site is
kept when the three bases are ACGT, it passes the depth/indel filter and it is polarizable by the
same rule.

Output: QC/site_table/{FOCAL}.npz
  cds_gene, cds_codon, cds_q, cds_chrom, cds_pos, cds_ref, cds_sis, cds_out, cds_anc  (bases 0-3 = ACGT)
  cds_fold      degeneracy of the site in the ancestral codon (0, 2, 3, 4)
  cds_anc_codon ancestral codon as 16*b1 + 4*b2 + b3
  cds_tot (N, 4) read counts ACGT, both strands (transcript orientation)
  cds_sb  bitmask: bit b set if base b has >= 1 read on each strand
  si_intron, si_chrom, si_pos, si_ref, si_sis, si_out, si_anc, si_tot, si_sb, si_offset (1-based
  position from the donor)
  genes (focal tx), gene_sister, gene_dbia, gene_chrom, chroms, introns (gene index, k)
Usage: 97_site_table.py FOCAL    (TEST_CHROM=<chrom> for a one-chromosome test)
"""
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import Align

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('p60', f'{PROJ}/scripts/60_polarized_mk.py')
p60 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p60)
mk, cs = p60.mk, p60.cs
A, MK, QC = p60.A, p60.MK, p60.QC
FNA = dict(p60.FNA, Dbia=f'{A}/genome/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.fna')
GFF = dict(p60.GFF, Dbia=f'{A}/genome/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.gff')
OUT = f'{QC}/site_table'
BASES = 'ACGT'
BI = {b: i for i, b in enumerate(BASES)}
COMP = p60.COMP
SI_MIN, SI_MAX, SI_FROM, SI_TO = 30, 65, 8, 30

ALN = Align.PairwiseAligner()
ALN.mode = 'global'
ALN.match_score, ALN.mismatch_score = 2, -3
ALN.open_gap_score, ALN.extend_gap_score = -5, -2


def revcomp(s):
    return ''.join(COMP.get(b, 'N') for b in reversed(s))


def read_genome(fna):
    seqs, name, buf = {}, None, []
    with open(fna) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    seqs[name] = ''.join(buf).upper()
    return seqs


def cds_segments(gff, wanted):
    """{tx: (chrom, strand, [(start, end), ...] in transcript order)} for transcripts in `wanted`."""
    parts = defaultdict(list)
    with open(gff) as f:
        for line in f:
            t = line.split('\t')
            if len(t) > 8 and t[2] == 'CDS':
                tx = p60.strip(t[8].split('Parent=')[1].split(';')[0])
                if tx in wanted:
                    parts[tx].append((t[0], int(t[3]), int(t[4]), t[6]))
    out = {}
    for tx, p in parts.items():
        if len({x[0] for x in p}) != 1 or len({x[3] for x in p}) != 1:
            continue
        strand = p[0][3]
        segs = sorted(((s, e) for _, s, e, _ in p), reverse=(strand == '-'))
        out[tx] = (p[0][0], strand, segs)
    return out


def introns_of(model):
    """List of (cds_coord, donor_genomic, length) for introns between CDS segments.
    cds_coord = number of CDS bases before the intron; donor = genomic position of intron base 1."""
    chrom, strand, segs = model
    res, cum = [], 0
    for k in range(len(segs) - 1):
        (s1, e1), (s2, e2) = segs[k], segs[k + 1]
        cum += e1 - s1 + 1
        if strand == '+':
            res.append((cum, e1 + 1, s2 - e1 - 1))
        else:
            res.append((cum, s1 - 1, s1 - e2 - 1))
    return res


def intron_seq(genome, chrom, strand, donor, length):
    if strand == '+':
        return genome[chrom][donor - 1:donor - 1 + length]
    return revcomp(genome[chrom][donor - length:donor])


def anchored_index(path, anchor, other_prefix):
    """(anchor ungapped seq, other ungapped seq, other index per anchor position or -1)."""
    s = mk.read_fasta(path)
    a = s[anchor]
    oth = [v for k, v in s.items() if k.lower().startswith(other_prefix)]
    assert len(oth) == 1, (path, list(s))
    o = oth[0]
    idx, j = [], -1
    for ca, co in zip(a, o):
        if co != '-':
            j += 1
        if ca != '-':
            idx.append(j if co != '-' else -1)
    return a.replace('-', ''), o.replace('-', ''), np.array(idx)


def pair_map(x, y):
    """Global alignment of x and y: array with the y index for each x position, or -1."""
    m = np.full(len(x), -1)
    aln = ALN.align(x, y)[0]
    for (xs, xe), (ys, ye) in zip(*aln.aligned):
        m[xs:xe] = np.arange(ys, ye)
    return m


def fold_of(codon, q):
    aa = mk.CODE[codon]
    syn = sum(mk.CODE[codon[:q] + b + codon[q + 1:]] == aa for b in BASES if b != codon[q])
    return {0: 0, 1: 2, 2: 3, 3: 4}[syn]


def main():
    focal = sys.argv[1]
    sp = 'Dsub' if focal.startswith('Dsub') else 'Dsuz'
    sis_sp = 'Dsuz' if sp == 'Dsub' else 'Dsub'
    if sp == 'Dsub':
        sis_dir, out_dir = f'{MK}/prank_out/Dsub_vs_Dsuz', f'{MK}/prank_out/Dsub_vs_Dbia'
        sis_prefix, anchor = 'dsuz', 'Dsub_0'
    else:
        sis_dir = f'{MK}/prank_out/Dsuz/{focal}/Dsuz_vs_Dsub'
        out_dir = f'{MK}/prank_out/Dsuz/{focal}/Dsuz_vs_Dbia'
        sis_prefix, anchor = 'dsub', f'{focal}_0'
    regions = f'{A}/population_polymorphism/reanalysis_20260821/config/regions'
    major_chroms = {l.split()[0] for cls in ('autosome', 'X') for l in open(f'{regions}/{sp}.{cls}.bed')}
    tri, _ = p60.orthologs()
    tri = tri[sp]
    models = p60.cds_models(p60.GFF[sp], major_chroms)
    fseg = cds_segments(GFF[sp], set(tri))
    sseg = cds_segments(GFF[sis_sp], {v[0] for v in tri.values()})
    bseg = cds_segments(GFF['Dbia'], {v[1] for v in tri.values()})
    print(f'{focal}: loading sister and outgroup genomes', flush=True)
    gen_s, gen_b = read_genome(FNA[sis_sp]), read_genome(FNA['Dbia'])
    masks = np.load(f'{QC}/masks_sites/{sp}.npz')
    thr = pd.read_csv(f'{PROJ}/results/depth_thresholds.tsv', sep='\t')
    thr = thr[thr['sample'] == focal].set_index('cls')
    qc = defaultdict(int)
    by_chrom = defaultdict(list)
    for tx in tri:
        if tx in models:
            by_chrom[models[tx][0]].append(tx)
    test_chrom = os.environ.get('TEST_CHROM')
    C = defaultdict(list)          # coding columns
    S = defaultdict(list)          # short-intron columns
    genes, gene_sis, gene_dbia, gene_chrom, introns, chroms = [], [], [], [], [], []
    for chrom, txs in sorted(by_chrom.items()):
        if test_chrom and chrom != test_chrom:
            continue
        ci = len(chroms)
        chroms.append(chrom)
        cls = 'X' if chrom in cs.XCHROM else 'autosome'
        dmax = int(thr.loc[cls, 'p99'])
        cnt, indel = p60.load_chrom(focal, chrom)
        gseq = p60.chrom_seq(FNA[sp], chrom)
        assert len(gseq) == len(cnt)
        si_mask = masks[f'{chrom}|si']
        for tx in txs:
            f_s, f_o = f'{sis_dir}/rna-{tx}.best.fas', f'{out_dir}/rna-{tx}.best.fas'
            if not (os.path.exists(f_s) and os.path.exists(f_o)):
                qc['skip_alignment_missing'] += 1
                continue
            ref1, sis_u, sis_ix = anchored_index(f_s, anchor, sis_prefix)
            ref2, out_u, out_ix = anchored_index(f_o, anchor, 'dbia')
            sis = ''.join(sis_u[j] if j >= 0 else '-' for j in sis_ix)
            outg = ''.join(out_u[j] if j >= 0 else '-' for j in out_ix)
            _, pos, strand = models[tx]
            if len(pos) == len(ref1) + 3:
                pos = pos[:-3]
            if ref1 != ref2 or len(ref1) != len(pos) or len(ref1) % 3:
                qc['skip_anchor_mismatch'] += 1
                continue
            gref = ''.join(gseq[p - 1] for p in pos)
            if strand == '-':
                gref = ''.join(COMP.get(b, 'N') for b in gref)
            if gref != ref1.upper():
                qc['skip_anchor_differs_from_genome'] += 1
                continue
            gi = len(genes)
            genes.append(tx)
            gene_sis.append(tri[tx][0])
            gene_dbia.append(tri[tx][1])
            gene_chrom.append(ci)
            c = cnt[pos - 1].astype(np.int64)
            if strand == '-':
                c = c[:, [7, 6, 5, 4, 3, 2, 1, 0]]
            tot = c[:, 0::2] + c[:, 1::2]
            sb = ((c[:, 0::2] > 0) & (c[:, 1::2] > 0)).astype(np.uint8) @ (1 << np.arange(4)).astype(np.uint8)
            depth = tot.sum(1)
            ok_site = (depth >= p60.MIN_DEPTH) & (depth <= dmax) & ~indel[pos - 1]
            passed = (tot >= p60.MIN_ALT) & (c[:, 0::2] > 0) & (c[:, 1::2] > 0) & \
                     (tot > p60.MIN_FRAC * np.maximum(depth, 1)[:, None])
            ref_u = ref1.upper()
            for k in range(0, len(pos), 3):
                idx = range(k, k + 3)
                if not ok_site[k:k + 3].all():
                    continue
                if any(sis[i] not in BASES or outg[i] not in BASES for i in idx):
                    continue
                major = ''.join(BASES[int(np.argmax(tot[i]))] for i in idx)
                if mk.CODE[major] == '*' or mk.CODE[sis[k:k + 3]] == '*':
                    continue
                if any(not passed[i, int(np.argmax(tot[i]))] for i in idx):
                    continue
                anc, bad = [], False
                for q, i in enumerate(idx):
                    x, y, o = major[q], sis[i], outg[i]
                    if x == y or y == o:
                        anc.append(y)
                    elif x == o:
                        anc.append(x)
                    else:
                        bad = True
                if bad:
                    qc['codons_unpolarized'] += 1
                    continue
                anc = ''.join(anc)
                if mk.CODE[anc] == '*':
                    continue
                qc['codons_callable'] += 1
                acode = 16 * BI[anc[0]] + 4 * BI[anc[1]] + BI[anc[2]]
                for q, i in enumerate(idx):
                    C['gene'].append(gi)
                    C['codon'].append(k // 3)
                    C['q'].append(q)
                    C['chrom'].append(ci)
                    C['pos'].append(int(pos[i]))
                    C['ref'].append(BI.get(ref_u[i], 4))
                    C['sis'].append(BI[sis[i]])
                    C['out'].append(BI[outg[i]])
                    C['anc'].append(BI[anc[q]])
                    C['fold'].append(fold_of(anc, q))
                    C['anc_codon'].append(acode)
                    C['tot'].append(tot[i])
                    C['sb'].append(sb[i])
            # short introns
            if tx not in fseg or tri[tx][0] not in sseg or tri[tx][1] not in bseg:
                qc['genes_no_segment_model'] += 1
                continue
            fm, smod, bmod = fseg[tx], sseg[tri[tx][0]], bseg[tri[tx][1]]
            s_int = {cc: (d, L) for cc, d, L in introns_of(smod)}
            b_int = {cc: (d, L) for cc, d, L in introns_of(bmod)}
            for kk, (cc, donor, L) in enumerate(introns_of(fm)):
                if not (SI_MIN <= L <= SI_MAX) or cc >= len(ref1):
                    continue
                qc['short_introns_focal'] += 1
                js, jb = sis_ix[cc], out_ix[cc]
                if js < 1 or jb < 1 or sis_ix[cc - 1] != js - 1 or out_ix[cc - 1] != jb - 1:
                    qc['short_introns_flank_gap'] += 1
                    continue
                if js not in s_int or jb not in b_int:
                    qc['short_introns_no_orthologous_intron'] += 1
                    continue
                fseq = intron_seq({chrom: gseq}, chrom, strand, donor, L)
                sseq = intron_seq(gen_s, smod[0], smod[1], *s_int[js])
                bseq = intron_seq(gen_b, bmod[0], bmod[1], *b_int[jb])
                if len(sseq) > 3 * SI_MAX or len(bseq) > 3 * SI_MAX:
                    qc['short_introns_orthologue_too_long'] += 1
                    continue
                qc['short_introns_used'] += 1
                ms, mb = pair_map(fseq, sseq), pair_map(fseq, bseq)
                ii = len(introns)
                introns.append((gi, kk))
                for off in range(SI_FROM, min(SI_TO, L) + 1):
                    gpos = donor + off - 1 if strand == '+' else donor - off + 1
                    if not si_mask[gpos]:
                        continue
                    i0 = off - 1
                    if ms[i0] < 0 or mb[i0] < 0:
                        continue
                    y, o = sseq[ms[i0]], bseq[mb[i0]]
                    cg = cnt[gpos - 1].astype(np.int64)
                    if strand == '-':
                        cg = cg[[7, 6, 5, 4, 3, 2, 1, 0]]
                    tg = cg[0::2] + cg[1::2]
                    d = tg.sum()
                    if not (p60.MIN_DEPTH <= d <= dmax) or indel[gpos - 1]:
                        continue
                    if y not in BASES or o not in BASES or fseq[i0] not in BASES:
                        continue
                    x = BASES[int(np.argmax(tg))]
                    if x == y or y == o:
                        a = y
                    elif x == o:
                        a = x
                    else:
                        continue
                    S['intron'].append(ii)
                    S['offset'].append(off)
                    S['chrom'].append(ci)
                    S['pos'].append(gpos)
                    S['ref'].append(BI[fseq[i0]])
                    S['sis'].append(BI[y])
                    S['out'].append(BI[o])
                    S['anc'].append(BI[a])
                    S['tot'].append(tg)
                    S['sb'].append(int(((cg[0::2] > 0) & (cg[1::2] > 0)) @ (1 << np.arange(4))))
        print(f'{focal} {chrom}: genes {len(genes)} coding sites {len(C["pos"])} '
              f'intron sites {len(S["pos"])}', flush=True)
    os.makedirs(OUT, exist_ok=True)
    arr = {}
    small = {'codon': np.int32, 'q': np.int8, 'chrom': np.int8, 'pos': np.int32, 'ref': np.int8,
             'sis': np.int8, 'out': np.int8, 'anc': np.int8, 'fold': np.int8, 'anc_codon': np.int8,
             'sb': np.uint8, 'gene': np.int32, 'intron': np.int32, 'offset': np.int8}
    for pre, D in (('cds', C), ('si', S)):
        for key, v in D.items():
            arr[f'{pre}_{key}'] = np.array(v, dtype=np.uint16 if key == 'tot' else small[key])
    arr.update(genes=np.array(genes), gene_sister=np.array(gene_sis), gene_dbia=np.array(gene_dbia),
               gene_chrom=np.array(gene_chrom, dtype=np.int8), chroms=np.array(chroms),
               introns=np.array(introns, dtype=np.int32).reshape(-1, 2))
    name = focal + (f'_test_{test_chrom}' if test_chrom else '')
    np.savez_compressed(f'{OUT}/{name}.npz', **arr)
    # validation against script 60: callable codons per gene
    g60 = pd.read_csv(f'{PROJ}/results/polarized_mk/{focal}.genes.tsv', sep='\t').set_index('gene')
    n_cod = pd.Series(arr['cds_codon'][arr['cds_q'] == 0]).groupby(
        arr['cds_gene'][arr['cds_q'] == 0]).size()
    mine = pd.Series(n_cod.values, index=np.array(genes)[n_cod.index])
    comp = g60.loc[g60.index.isin(genes), 'callable'].reindex(genes).fillna(0)
    mine = mine.reindex(genes).fillna(0)
    qc['validation_genes'] = len(genes)
    qc['validation_genes_callable_mismatch'] = int((comp != mine).sum())
    qc['validation_callable_codons_script60'] = int(comp.sum())
    qc['validation_callable_codons_here'] = int(mine.sum())
    pd.Series(qc).to_csv(f'{OUT}/{name}.qc.tsv', sep='\t', header=False)
    print(dict(qc))


if __name__ == '__main__':
    main()
