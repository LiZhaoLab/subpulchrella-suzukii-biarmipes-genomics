#!/usr/bin/env python3
"""Lineage-specific (polarized) McDonald-Kreitman counts with D. biarmipes as outgroup.

For one focal pool (Dsub_L5_7, or one D. suzukii population):
  - sister and outgroup bases at each focal CDS position come from Airi's pairwise PRANK
    codon alignments (focal-vs-sister and focal-vs-D. biarmipes), joined through the focal
    reference sequence (<prefix>_0), which must equal the reference CDS from the GFF;
  - polymorphism and allele frequencies come from my per-site pool allele counts
    (01_pileup_counts.py; MAPQ >= 20, base quality >= 20, duplicates removed);
  - only gene trios with a closed reciprocal-best-hit triangle (focal-sister-D. biarmipes)
    on the major chromosomes are used.
Callable codon: all 3 positions with depth 20..P99, outside indel masks, and sister and
outgroup bases in ACGT. An allele passes if it has >= 3 reads, >= 1 read on each strand and
frequency > 0.05 (the manuscript's MK thresholds). Focal state = major pool allele.
Polarization per position (parsimony): if focal == sister, no change; if sister == outgroup,
the change is on the focal lineage; if focal == outgroup, on the sister lineage; otherwise
the codon is skipped. Polymorphisms are polarized when sister == outgroup and that base is
one of the two pool alleles. Syn/nonsyn by pathway averaging (20_mk_recount.path_counts).
Dsub_L5_7 run also writes contamination evidence: pool read counts at positions where the
D. suzukii (or D. biarmipes, D. melanogaster) base differs from the D. subpulchrella base.

Usage: 60_polarized_mk.py FOCAL   (Dsub_L5_7 | JP-Kan | JP-Sap | JP-Sap2 | JP-Tok)
Outputs in results/polarized_mk/: {FOCAL}.genes.tsv, {FOCAL}.polymorphisms.tsv.gz,
{FOCAL}.qc.tsv, and for Dsub_L5_7 {FOCAL}.contamination_sites.tsv.gz

EXCLUDE_SHARED=1 (Results_Alpha_Shared_Polymorphism.md): a focal-lineage
fixed difference is not counted when the other species' pool carries the focal derived base (same
allele filter: >= 3 reads, >= 1 per strand, > 5%) at the aligned position, i.e. the site is
polymorphic in the other species (the sister state is otherwise one reference genome). Aligned
position: the focal-vs-sister PRANK alignment and the sister CDS model (genes whose aligned sister
sequence differs from the sister genome are not tested). The other pool's position must have depth
20..P99 of that pool and lie outside its indel mask; untestable fixed differences are kept.
Other pools: D. suzukii JP-Kan, JP-Sap, JP-Sap2 for the D. subpulchrella lineage (excluded if the
base passes in any of them; JP-Tok not used, ~1.7% D. subpulchrella DNA), Dsub_L5_7 for the
D. suzukii lineages. Dn_focal / Ds_focal are then the counts after exclusion; Dn_focal_all /
Ds_focal_all the counts without exclusion (equal to the default run); Dn_focal_only_{pool} /
Ds_focal_only_{pool} the counts with exclusion by that pool alone. Outputs go to
results/polarized_mk_excl/; the contamination file is not written.
"""
import glob
import gzip
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

A = '/ru-auth/local/home/asato/scratch'
MK = f'{A}/MKtest/reanalysis_longest_fix_20260817'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
OUT = f'{PROJ}/results/polarized_mk'
spec = importlib.util.spec_from_file_location('mk', f'{PROJ}/scripts/20_mk_recount.py')
mk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mk)
spec = importlib.util.spec_from_file_location('cs', f'{PROJ}/scripts/32_chunk_stats.py')
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)

FNA = {'Dsub': f'{A}/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna',
       'Dsuz': f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna'}
GFF = {'Dsub': f'{A}/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff',
       'Dsuz': f'{A}/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff'}
MIN_DEPTH, MIN_ALT, MIN_FRAC = 20, 3, 0.05
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
BASES = 'ACGT'
EXCL = os.environ.get('EXCLUDE_SHARED') == '1'
OTHER_POOLS = {'Dsub': ['JP-Kan', 'JP-Sap', 'JP-Sap2'], 'Dsuz': ['Dsub_L5_7']}
if EXCL:
    OUT += '_excl'


def strip(x):
    return x.replace('rna-', '')


def orthologs():
    """Return dict focal_species -> {focal tx: (sister tx, dbia tx)} with closed RBH triangles."""
    o = f'{MK}/orth_longest'
    sub_suz = pd.read_csv(f'{o}/Dsuz_Dsub.direct.RBH.longest.mRNA.tsv', sep='\t', header=None)
    sub_bia = pd.read_csv(f'{o}/Dsub_Dbia.RBH.longest.mRNA.tsv', sep='\t', header=None)
    suz_bia = pd.read_csv(f'{o}/Dsuz_Dbia.direct.RBH.longest.mRNA.tsv', sep='\t', header=None)
    s2z = dict(zip(sub_suz[0].map(strip), sub_suz[1].map(strip)))      # Dsub -> Dsuz
    s2b = dict(zip(sub_bia[1].map(strip), sub_bia[0].map(strip)))      # Dsub -> Dbia
    z2b = dict(zip(suz_bia[1].map(strip), suz_bia[0].map(strip)))      # Dsuz -> Dbia
    tri = {s: (z, s2b[s]) for s, z in s2z.items() if s in s2b and z2b.get(z) == s2b[s]}
    return {'Dsub': tri, 'Dsuz': {z: (s, b) for s, (z, b) in tri.items()}}, len(s2z)


def cds_models(gff, chroms):
    parts = defaultdict(list)
    with open(gff) as f:
        for line in f:
            t = line.split('\t')
            if len(t) > 8 and t[2] == 'CDS' and t[0] in chroms:
                par = t[8].split('Parent=')[1].split(';')[0]
                parts[strip(par)].append((t[0], int(t[3]), int(t[4]), t[6]))
    models = {}
    for tx, p in parts.items():
        if len({x[0] for x in p}) != 1 or len({x[3] for x in p}) != 1:
            continue
        pos = []
        for _, s, e, _ in sorted(p, key=lambda x: x[1]):
            pos.extend(range(s, e + 1))
        strand = p[0][3]
        models[tx] = (p[0][0], np.array(pos[::-1] if strand == '-' else pos), strand)
    return models


def anchored(path, anchor, other_prefix):
    """Map an alignment onto ungapped anchor positions: (anchor seq, other base per position)."""
    s = mk.read_fasta(path)
    a = s[anchor]
    oth = [v for k, v in s.items() if k.lower().startswith(other_prefix)]
    assert len(oth) == 1, (path, list(s))
    keep = [i for i, c in enumerate(a) if c != '-']
    return ''.join(a[i] for i in keep), ''.join(oth[0][i] for i in keep)


def anchored_ix(path, anchor, other_prefix):
    """As anchored(), plus the ungapped other-sequence index per anchor position (-1 = gap) and
    the ungapped other sequence (EXCLUDE_SHARED only)."""
    s = mk.read_fasta(path)
    a = s[anchor]
    oth = [v for k, v in s.items() if k.lower().startswith(other_prefix)]
    assert len(oth) == 1, (path, list(s))
    o = oth[0]
    keep = [i for i, c in enumerate(a) if c != '-']
    idx, j = [], -1
    for ca, co in zip(a, o):
        if co != '-':
            j += 1
        if ca != '-':
            idx.append(j if co != '-' else -1)
    return ''.join(a[i] for i in keep), ''.join(o[i] for i in keep), np.array(idx), o.replace('-', '').upper()


class OtherPools:
    """Allele-filter flags of the other species' pools, one byte per sister-genome position:
    bits 0-3 = base A, C, G, T passes the allele filter; bit 4 = callable (depth 20..P99 of that
    pool, outside its indel mask). Loaded per chromosome on first use."""

    def __init__(self, thr):
        self.thr, self.cache = thr, {}

    def flags(self, pool, chrom):
        if (pool, chrom) not in self.cache:
            cnt, indel = load_chrom(pool, chrom)
            tot = cnt[:, 0::2].astype(np.int32) + cnt[:, 1::2]
            depth = tot.sum(1)
            dmax = int(self.thr.loc[(pool, 'X' if chrom in cs.XCHROM else 'autosome'), 'p99'])
            passed = (tot >= MIN_ALT) & (cnt[:, 0::2] > 0) & (cnt[:, 1::2] > 0) & \
                     (tot > MIN_FRAC * np.maximum(depth, 1)[:, None])
            f = (passed.astype(np.uint8) << np.arange(4, dtype=np.uint8)).sum(1).astype(np.uint8)
            ok = (depth >= MIN_DEPTH) & (depth <= dmax) & ~indel
            self.cache[(pool, chrom)] = f | (ok.astype(np.uint8) << 4)
        return self.cache[(pool, chrom)]


def load_chrom(sample, chrom):
    """Per-site counts (L, 8) and indel mask for a whole chromosome; index = position - 1."""
    files = sorted(glob.glob(f'{QC}/counts/{sample}/{chrom}_*.npz'),
                   key=lambda f: int(f.split('_')[-2]))
    cnt, ind = [], []
    for f in files:
        d = np.load(f)
        start, end = int(d['start']), int(d['end'])
        m = np.zeros(end - start + 2, dtype=np.int32)
        for s, e in d['indel_mask']:
            s, e = max(s, start), min(e, end)
            if s <= e:
                m[s - start] += 1
                m[e - start + 1] -= 1
        cnt.append(d['counts'])
        ind.append(np.cumsum(m)[:-1] > 0)
    return np.concatenate(cnt), np.concatenate(ind)


def chrom_seq(fna, chrom):
    out, keep = [], False
    with open(fna) as f:
        for line in f:
            if line.startswith('>'):
                if keep:
                    break
                keep = line[1:].split()[0] == chrom
            elif keep:
                out.append(line.strip().upper())
    return ''.join(out)


def main():
    focal = sys.argv[1]
    sp = 'Dsub' if focal.startswith('Dsub') else 'Dsuz'
    if sp == 'Dsub':
        sis_dir, out_dir = f'{MK}/prank_out/Dsub_vs_Dsuz', f'{MK}/prank_out/Dsub_vs_Dbia'
        sis_prefix, anchor = 'dsuz', 'Dsub_0'
        dmel_dir = f'{MK}/prank_out/Dsub_vs_Dmel'
    else:
        sis_dir, out_dir = f'{MK}/prank_out/Dsuz/{focal}/Dsuz_vs_Dsub', f'{MK}/prank_out/Dsuz/{focal}/Dsuz_vs_Dbia'
        sis_prefix, anchor, dmel_dir = 'dsub', f'{focal}_0', None
    regions = f'{A}/population_polymorphism/reanalysis_20260821/config/regions'
    major = {l.split()[0] for cls in ('autosome', 'X') for l in open(f'{regions}/{sp}.{cls}.bed')}
    tri, n_pairs = orthologs()
    tri = tri[sp]
    models = cds_models(GFF[sp], major)
    thr = pd.read_csv(f'{PROJ}/results/depth_thresholds.tsv', sep='\t')
    thr = thr[thr['sample'] == focal].set_index('cls')
    qc = defaultdict(int)
    qc['sister_pairs'] = n_pairs
    qc['closed_triangles'] = len(tri)
    if EXCL:
        sis_sp = 'Dsuz' if sp == 'Dsub' else 'Dsub'
        major_s = {l.split()[0] for cls in ('autosome', 'X') for l in open(f'{regions}/{sis_sp}.{cls}.bed')}
        smodels = cds_models(GFF[sis_sp], major_s)
        sgen = {}
        others = OtherPools(pd.read_csv(f'{PROJ}/results/depth_thresholds.tsv', sep='\t').set_index(['sample', 'cls']))
        opools = OTHER_POOLS[sp]
    genes, polys, contam = [], [], []
    by_chrom = defaultdict(list)
    for tx in tri:
        if tx in models:
            by_chrom[models[tx][0]].append(tx)
        else:
            qc['skip_not_major_chromosome'] += 1
    test_chrom = os.environ.get('TEST_CHROM')          # test mode: one chromosome only
    for chrom, txs in sorted(by_chrom.items()):
        if test_chrom and chrom != test_chrom:
            continue
        cls = 'X' if chrom in cs.XCHROM else 'autosome'
        dmax = int(thr.loc[cls, 'p99'])
        cnt, indel = load_chrom(focal, chrom)
        gseq = chrom_seq(FNA[sp], chrom)
        assert len(gseq) == len(cnt), (chrom, len(gseq), len(cnt))
        for tx in txs:
            f_s, f_o = f'{sis_dir}/rna-{tx}.best.fas', f'{out_dir}/rna-{tx}.best.fas'
            if not (os.path.exists(f_s) and os.path.exists(f_o)):
                qc['skip_alignment_missing'] += 1
                continue
            if EXCL:
                ref1, sis, s_ix, s_u = anchored_ix(f_s, anchor, sis_prefix)
            else:
                ref1, sis = anchored(f_s, anchor, sis_prefix)
            ref2, outg = anchored(f_o, anchor, 'dbia')
            dmel = None
            if dmel_dir and os.path.exists(f'{dmel_dir}/rna-{tx}.best.fas'):
                ref3, dmel = anchored(f'{dmel_dir}/rna-{tx}.best.fas', anchor, 'dmel')
                if ref3 != ref1:
                    dmel = None
            _, pos, strand = models[tx]
            if len(pos) == len(ref1) + 3:                 # alignments omit the stop codon
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
            s_ok = False
            if EXCL:                                      # sister genome positions of the alignment
                stx = tri[tx][0]
                if stx in smodels:
                    schrom, spos, sstrand = smodels[stx]
                    if schrom not in sgen:
                        sgen[schrom] = chrom_seq(FNA[sis_sp], schrom)
                    gs_ = ''.join(sgen[schrom][p - 1] for p in spos)
                    if sstrand == '-':
                        gs_ = ''.join(COMP.get(b, 'N') for b in gs_)
                    su = s_u[:-3] if len(s_u) == len(gs_) + 3 else s_u
                    s_ok = su == gs_[:len(su)] and len(su) in (len(gs_), len(gs_) - 3)
                qc['excl_genes_sister_mapped' if s_ok else 'excl_genes_sister_not_mapped'] += 1
                if s_ok:
                    fl = [others.flags(p_, schrom) for p_ in opools]
            c = cnt[pos - 1].astype(np.int64)             # (L, 8): A+ A- C+ C- G+ G- T+ T-
            if strand == '-':                             # complement bases, swap strands
                c = c[:, [7, 6, 5, 4, 3, 2, 1, 0]]
            tot = c[:, 0::2] + c[:, 1::2]                 # (L, 4) in ACGT order
            depth = tot.sum(1)
            ok_site = (depth >= MIN_DEPTH) & (depth <= dmax) & ~indel[pos - 1]
            passed = (tot >= MIN_ALT) & (c[:, 0::2] > 0) & (c[:, 1::2] > 0) & \
                     (tot > MIN_FRAC * np.maximum(depth, 1)[:, None])
            g = dict(gene=tx, chrom=chrom, cls=cls, codons=len(pos) // 3, callable=0,
                     Pn=0.0, Ps=0.0, Dn_focal=0.0, Ds_focal=0.0, Dn_sister=0.0, Ds_sister=0.0,
                     unpolarized_codons=0)
            if EXCL:
                g.update(Dn_focal_all=0.0, Ds_focal_all=0.0, fixed_positions=0, fixed_testable=0, fixed_shared=0)
                if len(opools) > 1:
                    for p_ in opools:
                        g[f'Dn_focal_only_{p_}'] = 0.0
                        g[f'Ds_focal_only_{p_}'] = 0.0
            for k in range(0, len(pos), 3):
                idx = range(k, k + 3)
                if not ok_site[k:k + 3].all():
                    continue
                if any(sis[i] not in BASES or outg[i] not in BASES for i in idx):
                    continue
                major = ''.join(BASES[int(np.argmax(tot[i]))] for i in idx)
                if mk.CODE[major] == '*' or mk.CODE[sis[k:k + 3]] == '*':
                    continue
                alleles = [[BASES[j] for j in np.argsort(-tot[i]) if passed[i, j]][:2] for i in idx]
                if any(len(a) == 0 for a in alleles):     # major allele fails the strand filter
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
                    g['unpolarized_codons'] += 1
                    continue
                anc = ''.join(anc)
                if mk.CODE[anc] == '*':
                    continue
                g['callable'] += 1
                poly_pos = [q for q in range(3) if len(alleles[q]) == 2]
                for q in poly_pos:
                    a1, a2 = alleles[q]
                    var = major[:q] + (a2 if a1 == major[q] else a1) + major[q + 1:]
                    if mk.CODE[var] == '*':
                        continue
                    s_, n_ = mk.path_counts(major, var)
                    i = k + q
                    d = depth[i]
                    fa = {b: tot[i, BASES.index(b)] / d for b in (a1, a2)}
                    y, o = sis[i], outg[i]
                    daf = 1 - fa[y] if (y == o and y in fa) else np.nan
                    polys.append((tx, 'n' if n_ > 0 else 's', round(daf, 4) if daf == daf else '',
                                  round(min(fa.values()), 4)))
                    g['Pn'] += n_
                    g['Ps'] += s_
                fix_f = ''.join(major[q] if (q not in poly_pos and major[q] != anc[q]) else anc[q]
                                for q in range(3))
                if EXCL:                                  # drop positions polymorphic in the other species
                    shared = {p_: set() for p_ in opools}
                    for q in range(3):
                        if fix_f[q] == anc[q]:
                            continue
                        g['fixed_positions'] += 1
                        j = s_ix[k + q] if s_ok else -1
                        if j < 0 or j >= len(spos):
                            continue
                        b = BASES.index(fix_f[q])
                        b = 3 - b if sstrand == '-' else b   # genome strand of the sister base
                        tested = False
                        for p_, f in zip(opools, fl):
                            v = int(f[spos[j] - 1])
                            if v & 16:
                                tested = True
                                if (v >> b) & 1:
                                    shared[p_].add(q)
                        g['fixed_testable'] += int(tested)
                    any_sh = set().union(*shared.values())
                    g['fixed_shared'] += len(any_sh)
                    s_, n_ = mk.path_counts(anc, fix_f)
                    g['Ds_focal_all'] += s_
                    g['Dn_focal_all'] += n_
                    if len(opools) > 1:
                        for p_ in opools:
                            s_, n_ = mk.path_counts(anc, ''.join(anc[q] if q in shared[p_] else fix_f[q]
                                                                 for q in range(3)))
                            g[f'Ds_focal_only_{p_}'] += s_
                            g[f'Dn_focal_only_{p_}'] += n_
                    fix_f = ''.join(anc[q] if q in any_sh else fix_f[q] for q in range(3))
                s_, n_ = mk.path_counts(anc, fix_f)
                g['Ds_focal'] += s_
                g['Dn_focal'] += n_
                fix_s = ''.join(sis[k + q] if (q not in poly_pos and sis[k + q] != anc[q]) else anc[q]
                                for q in range(3))
                s_, n_ = mk.path_counts(anc, fix_s)
                g['Ds_sister'] += s_
                g['Dn_sister'] += n_
            genes.append(g)
            if sp == 'Dsub' and not EXCL:                 # contamination evidence
                # species-specific derived bases: D. suzukii base differs from D. subpulchrella
                # while D. biarmipes matches D. subpulchrella; D. biarmipes-specific: D. suzukii
                # matches D. subpulchrella; D. melanogaster-specific: both others match.
                r = np.frombuffer(ref1.upper().encode(), dtype='S1')
                y = np.frombuffer(sis.upper().encode(), dtype='S1')
                o = np.frombuffer(outg.upper().encode(), dtype='S1')
                m = np.frombuffer(dmel.upper().encode(), dtype='S1') if dmel else None
                isb = lambda v: np.isin(v, [b'A', b'C', b'G', b'T'])
                base_ok = ok_site & isb(r) & isb(y) & isb(o)
                tests = [('dsuz', y, base_ok & (y != r) & (o == r)),
                         ('dbia', o, base_ok & (o != r) & (y == r))]
                if m is not None:
                    tests.append(('dmel', m, base_ok & isb(m) & (m != r) & (y == r) & (o == r)))
                code = {b'A': 0, b'C': 1, b'G': 2, b'T': 3}
                for label, other, sel in tests:
                    for i in np.flatnonzero(sel):
                        ri, oi = code[r[i]], code[other[i]]
                        contam.append((label, chrom, int(pos[i]), int(depth[i]), int(tot[i, ri]),
                                       int(tot[i, oi]), int(depth[i] - tot[i, ri] - tot[i, oi])))
        print(f'{focal} {chrom}: genes so far {len(genes)}', flush=True)
    out = OUT + ('_test' if test_chrom else '')
    os.makedirs(out, exist_ok=True)
    pd.DataFrame(genes).to_csv(f'{out}/{focal}.genes.tsv', sep='\t', index=False)
    pd.DataFrame(polys, columns=['gene', 'type', 'daf', 'maf']).to_csv(
        f'{out}/{focal}.polymorphisms.tsv.gz', sep='\t', index=False)
    qc['genes_written'] = len(genes)
    if EXCL and genes:
        gd = pd.DataFrame(genes)
        for k_ in ('fixed_positions', 'fixed_testable', 'fixed_shared'):
            qc[f'excl_{k_}'] = int(gd[k_].sum())
    pd.Series(qc).to_csv(f'{out}/{focal}.qc.tsv', sep='\t', header=False)
    if contam:
        pd.DataFrame(contam, columns=['species', 'chrom', 'pos', 'depth', 'n_focal_allele',
                                      'n_species_allele', 'n_other_two_bases']).to_csv(
            f'{out}/{focal}.contamination_sites.tsv.gz', sep='\t', index=False)
    print(dict(qc))


if __name__ == '__main__':
    main()
