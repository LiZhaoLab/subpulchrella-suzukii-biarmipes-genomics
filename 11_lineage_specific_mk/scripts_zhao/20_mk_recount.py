#!/usr/bin/env python3
"""Independent McDonald-Kreitman recount from Airi's PRANK codon alignments.

For each alignment (ingroup sequences named <ingroup>_0.._k, one outgroup sequence),
count per codon column (columns with a gap or non-ACGT base in any sequence are skipped):
  polymorphic sites: positions that vary among ingroup sequences;
  fixed differences: positions where all ingroup sequences agree and the outgroup differs.
Codons with more than one change are resolved by averaging over all mutational
pathways that avoid stop codons (Nei-Gojobori style). Polymorphic codons with three
distinct ingroup codons are connected through the minimum spanning path.
P value: two-sided Fisher's exact test on [[Pn, Ps], [Dn, Ds]]; alpha = 1 - Ds*Pn/(Dn*Ps).
Environment options: KEEP_REF=1 adds the reference-strain sequence to the ingroup sample;
THIRD_ALLELE=0 counts a polymorphic position only as polymorphic even when the outgroup
carries a third base (default 1, as cmktest does).

Usage: 20_mk_recount.py ALIGN_DIR INGROUP_PREFIX OUT_TSV
"""
import glob
import itertools
import os
import sys

from scipy.stats import fisher_exact

BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AA[16 * i + 4 * j + k]
        for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
_cache = {}


def path_counts(c1, c2):
    """(syn, nonsyn) changes between two codons, averaged over stop-free pathways."""
    key = (c1, c2)
    if key in _cache:
        return _cache[key]
    diff = [i for i in range(3) if c1[i] != c2[i]]
    if not diff:
        r = (0.0, 0.0)
    elif len(diff) == 1:
        r = (1.0, 0.0) if CODE[c1] == CODE[c2] else (0.0, 1.0)
    else:
        tot_s = tot_n = 0.0
        npath = 0
        for order in itertools.permutations(diff):
            cur = list(c1)
            s = n = 0
            ok = True
            for step, pos in enumerate(order):
                prev = ''.join(cur)
                cur[pos] = c2[pos]
                nxt = ''.join(cur)
                if step < len(order) - 1 and CODE[nxt] == '*':
                    ok = False
                    break
                if CODE[prev] == CODE[nxt]:
                    s += 1
                else:
                    n += 1
            if ok:
                tot_s += s
                tot_n += n
                npath += 1
        if npath == 0:          # every pathway passes a stop codon: count all as nonsynonymous
            r = (0.0, float(len(diff)))
        else:
            r = (tot_s / npath, tot_n / npath)
    _cache[key] = r
    return r


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def mk_counts(seqs, ingroup):
    # cmktest reports Npopulation = 2: the pool haplotypes <ingroup>_1 and _2 form the
    # population sample; <ingroup>_0 (the reference strain) is excluded unless KEEP_REF=1.
    keep = ('_0', '_1', '_2') if os.environ.get('KEEP_REF') == '1' else ('_1', '_2')
    ing = [v for k, v in seqs.items() if k.startswith(ingroup + '_') and k.endswith(keep)]
    out = [v for k, v in seqs.items() if not k.startswith(ingroup + '_')]
    assert sum(k.startswith(ingroup + '_') for k in seqs) == 3, list(seqs)
    assert len(out) == 1 and len(ing) >= 2, (list(seqs), ingroup)
    out = out[0]
    L = len(out)
    assert all(len(s) == L for s in ing) and L % 3 == 0
    Pn = Ps = Dn = Ds = 0.0
    ncod = 0
    for i in range(0, L, 3):
        cods = [s[i:i + 3] for s in ing]
        oc = out[i:i + 3]
        if any(set(c) - set('ACGT') for c in cods + [oc]):
            continue
        if any(CODE[c] == '*' for c in cods + [oc]) and i + 3 < L:
            continue
        ncod += 1
        distinct = sorted(set(cods))
        poly_pos = {p for p in range(3) if len({c[p] for c in cods}) > 1}
        if len(distinct) == 2:
            s, n = path_counts(distinct[0], distinct[1])
            Ps += s
            Pn += n
        elif len(distinct) == 3:
            pairs = sorted(itertools.combinations(distinct, 2),
                           key=lambda ab: sum(x != y for x, y in zip(*ab)))
            for a, b in pairs[:2]:
                s, n = path_counts(a, b)
                Ps += s
                Pn += n
        ref = cods[0]
        # THIRD_ALLELE=1 (cmktest convention): a polymorphic position where the outgroup
        # base is absent from every ingroup sequence also counts as a fixed difference.
        third = os.environ.get('THIRD_ALLELE', '1') == '1'
        fixed = [p for p in range(3) if oc[p] != ref[p] and
                 (p not in poly_pos or (third and all(c[p] != oc[p] for c in cods)))]
        if fixed:
            tgt = ''.join(oc[p] if p in fixed else ref[p] for p in range(3))
            s, n = path_counts(ref, tgt)
            Ds += s
            Dn += n
    return ncod, Pn, Ps, Dn, Ds


def main():
    adir, ingroup, out_tsv = sys.argv[1:4]
    files = sorted(glob.glob(os.path.join(adir, '*.best.fas')))
    assert files, adir
    with open(out_tsv, 'w') as fh:
        fh.write('id\tNcodon\tPn\tPs\tDn\tDs\talpha\tpvalue\n')
        for f in files:
            gid = os.path.basename(f).replace('.best.fas', '').replace('rna-', '')
            ncod, Pn, Ps, Dn, Ds = mk_counts(read_fasta(f), ingroup)
            r = [round(x) for x in (Pn, Ps, Dn, Ds)]
            p = fisher_exact([[r[0], r[1]], [r[2], r[3]]])[1]
            alpha = 1 - (Ds * Pn) / (Dn * Ps) if Dn > 0 and Ps > 0 else float('nan')
            fh.write(f'{gid}\t{ncod}\t{Pn:.3f}\t{Ps:.3f}\t{Dn:.3f}\t{Ds:.3f}\t{alpha:.6f}\t{p:.6g}\n')
    print(out_tsv, len(files))


if __name__ == '__main__':
    main()
