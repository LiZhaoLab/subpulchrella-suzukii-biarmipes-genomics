"""Estimate of gene copies hidden by the locus clustering of scripts/137.
Clustering merges an alignment into a locus when the spans overlap by >= 50% of the shorter span, with no exon test,
so a copy inside an intron of a larger same-strand gene is absorbed. Here, per genome: alignments that pass the
scripts/137 filters (cov >= 0.5, query set) with no locus of the same paralog group overlapping them on the same
strand, lying >= 50% inside the span of a locus of another paralog group on the same strand, and with no exon overlap
with that locus's seed alignment. Overlapping such alignments are merged into regions. No random numbers."""
import importlib.util
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('do', f'{PROJ}/scripts/144_duplication_origin.py')
do = importlib.util.module_from_spec(spec)
spec.loader.exec_module(do)
ev = do.ev
OUTF = sys.argv[1]

LO = pd.read_csv(f'{ev.OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
E = pd.read_csv(f'{ev.OUT}/C2_per_gene.tsv.gz', sep='\t', index_col=0, low_memory=False)
LO['pgroup'] = LO.dmel.map(E.pgroup)
dmel_of, full = ev.query_tables()
TCONS = {'M': 3, 'D': 3, 'I': 0, 'F': 1, 'G': 1, 'N': 1, 'U': 1, 'V': 1}


def blocks(cg, ts, te, strand):
    """Exon blocks (target coordinates) of a miniprot alignment; N/U/V = introns."""
    ops = [(int(k), o) for k, o in re.findall(r'(\d+)([MIDFGNUV])', cg)]
    tot = sum(k * TCONS[o] for k, o in ops)
    out, pos, cur = [], 0, 0
    for k, o in ops:
        L = k * TCONS[o]
        if o in 'NUV':
            if cur < pos:
                out.append((cur, pos))
            pos += L
            cur = pos
        else:
            pos += L
    if cur < pos:
        out.append((cur, pos))
    if strand == '+':
        return [(ts + a, ts + b) for a, b in out], tot - (te - ts)
    return [(te - b, te - a) for a, b in out], tot - (te - ts)


def ov_blocks(x, y):
    return any(min(b, d) - max(a, c) > 0 for a, b in x for c, d in y)


rows, diffs = [], []
for g in ['Dbia', 'Dsub', 'DsuzI', 'DsuzR']:
    A = ev.read_alignments(g, dmel_of, full)
    A['pgroup'] = ('gene-Dmel_' + A.dmel.astype(str).str.replace('gene-Dmel_', '')).map(E.pgroup)
    L = LO[LO.genome == g].copy()
    # seed exon blocks
    key = A.set_index(['q', 'tname', 'ts', 'te']).sort_index()
    seed_blocks = {}
    for r in L.itertuples():
        try:
            s = key.loc[(r.seed_q, r.tname, r.ts, r.te)]
        except KeyError:
            continue
        s = s.iloc[0] if isinstance(s, pd.DataFrame) else s
        b, d = blocks(s.cg, r.ts, r.te, r.strand)
        diffs.append(d)
        seed_blocks[r.Index] = b
    byseq = {k: x for k, x in L.groupby(['tname', 'strand'])}
    hid = []
    for a in A.itertuples():
        x = byseq.get((a.tname, a.strand))
        if x is None:
            continue
        o = np.minimum(x.te, a.te) - np.maximum(x.ts, a.ts)
        near = x[o > 0]
        if not len(near) or (near.pgroup == a.pgroup).any():
            continue  # represented by a locus of its own group (or no locus at all: not this problem)
        o2 = np.minimum(near.te, a.te) - np.maximum(near.ts, a.ts)
        host = near[o2 >= 0.5 * np.minimum(near.te - near.ts, a.te - a.ts)]
        if not len(host):
            continue
        ab, _ = blocks(a.cg, a.ts, a.te, a.strand)
        if any(h in seed_blocks and ov_blocks(ab, seed_blocks[h]) for h in host.index):
            continue
        h = host.iloc[0]
        hid.append(dict(genome=g, tname=a.tname, strand=a.strand, ts=a.ts, te=a.te, q=a.q, dmel=a.dmel, pgroup=a.pgroup,
                        AS=a.AS, cov=a.cov, fs=a.fs, st=a.st, introns=sum(o in 'NUV' for o in re.findall(r'\d+([MIDFGNUV])', a.cg)),
                        host_dmel=h.dmel, host_pgroup=h.pgroup, host_span=h.te - h.ts))
    H = pd.DataFrame(hid)
    if not len(H):
        continue
    # merge overlapping hidden alignments into regions; keep best-scoring alignment per region
    H = H.sort_values(['tname', 'strand', 'ts'])
    reg, cur_end, cur_key, rid = [], -1, None, -1
    for r in H.itertuples():
        if (r.tname, r.strand) != cur_key or r.ts >= cur_end:
            rid += 1
            cur_key, cur_end = (r.tname, r.strand), r.te
        cur_end = max(cur_end, r.te)
        reg.append(rid)
    H['region'] = reg
    H['intact'] = (H['cov'] >= 0.8) & (H.fs == 0) & (H.st == 0)
    best = H.sort_values('AS', ascending=False).drop_duplicates('region')
    best = best.assign(any_intact=best.region.map(H.groupby('region').intact.any()))
    rows.append(best)
R = pd.concat(rows)
R.to_csv(OUTF, sep='\t', index=False)
d = pd.Series(diffs)
print(f'CIGAR length check (target consumed - span) over seeds: exact {np.mean(d == 0):.3f}, |diff|<=3 {np.mean(d.abs() <= 3):.3f}')
print('hidden regions per genome:', R.genome.value_counts().to_dict())
print('  with an intact alignment:', R[R.any_intact].genome.value_counts().to_dict())
print('  intronless (best alignment):', R[R.introns == 0].genome.value_counts().to_dict())
print('  distinct paralog groups:', R.groupby('genome').pgroup.nunique().to_dict())
