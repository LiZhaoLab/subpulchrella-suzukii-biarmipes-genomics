#!/usr/bin/env python3
"""Inversion analysis: do the robust pool inversion events coincide with breakpoints of the
D. suzukii - D. subpulchrella synteny, and does reversing the segment restore D. suzukii order?

For each robust event (results/INV_robust_events.tsv) whose two ends both lie within TOL of a
boundary between collinear synteny blocks (results/synteny_segments.tsv, 'main',
Dsuz_vs_Dsub; D. subpulchrella = query), the four blocks around the breakpoints are found: L (ends at
the left breakpoint), A (first block inside), B (last block inside), R (starts at the right breakpoint).
Reversing the segment creates the adjacencies L-B' and A'-R. An adjacency is collinear with D. suzukii
when both blocks lie on the same D. suzukii chromosome, have the same orientation after the reversal,
and their facing ends are within GAP in D. suzukii coordinates. The same test is applied to the current
(reference) adjacencies L-A and B-R. Ancestral state: the same adjacency test in the D. subpulchrella -
D. biarmipes synteny (Dsub_vs_Dbia; D. subpulchrella = ref in that comparison) tells whether the
D. biarmipes order agrees with the reference or with the inverted arrangement.
Output: results/INV_synteny_check.tsv and .log
"""
import numpy as np
import pandas as pd

RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
TOL, GAP = 20_000, 300_000


def blocks(comp, dsub_is):
    s = pd.read_csv(f'{RES}/synteny_segments.tsv', sep='\t')
    s = s[(s['params'] == 'main') & (s['comp'] == comp)].copy()
    if dsub_is == 'query':
        s = s.rename(columns={'query': 'dchr', 'qs': 'ds', 'qe': 'de', 'ref': 'ochr', 'rs': 'os', 're': 'oe'})
    else:
        s = s.rename(columns={'ref': 'dchr', 'rs': 'ds', 're': 'de', 'query': 'ochr', 'qs': 'os', 'qe': 'oe'})
    return s[['dchr', 'ds', 'de', 'ochr', 'os', 'oe', 'ori']]


def end_coord(b, side, ori):
    """Other-species coordinate at the left/right end (in D. subpulchrella order) of block b with orientation ori."""
    if side == 'right':
        return b.oe if ori == 1 else b.os
    return b.os if ori == 1 else b.oe


def collinear(b1, o1, b2, o2):
    if b1.ochr != b2.ochr or o1 != o2:
        return False
    return abs(end_coord(b1, 'right', o1) - end_coord(b2, 'left', o2)) <= GAP


def check(ev, bl):
    d = bl[bl.dchr == ev.chrom]
    L = d[(d.de - ev.start).abs() <= TOL]
    A = d[(d.ds - ev.start).abs() <= TOL]
    B = d[(d.de - ev.end).abs() <= TOL]
    R = d[(d.ds - ev.end).abs() <= TOL]
    if min(len(L), len(A), len(B), len(R)) == 0:
        return None
    L, A, B, R = L.iloc[0], A.iloc[0], B.iloc[0], R.iloc[0]
    ref_ok = collinear(L, L.ori, A, A.ori) or collinear(B, B.ori, R, R.ori)
    inv_l = collinear(L, L.ori, B, -B.ori)
    inv_r = collinear(A, -A.ori, R, R.ori)
    return ref_ok, inv_l, inv_r


def main():
    ev = pd.read_csv(f'{RES}/INV_robust_events.tsv', sep='\t')
    suz = blocks('Dsuz_vs_Dsub', 'query')
    bia = blocks('Dsub_vs_Dbia', 'ref')
    rows = []
    for e in ev.itertuples():
        cs, cb = check(e, suz), check(e, bia)
        if cs is None and cb is None:
            continue
        rows.append(dict(arm=e.arm, chrom=e.chrom, start=e.start, end=e.end, L_mb=e.L_mb, pool_pairs=e.pool_pairs,
                         f_pairs=e.f_pairs, pool_split=e.pool_split, f_split=e.f_split,
                         dsuz_blocks_at_both_ends=cs is not None,
                         dsuz_reference_adjacency_collinear=None if cs is None else cs[0],
                         dsuz_inverted_left_collinear=None if cs is None else cs[1],
                         dsuz_inverted_right_collinear=None if cs is None else cs[2],
                         dbia_blocks_at_both_ends=cb is not None,
                         dbia_reference_adjacency_collinear=None if cb is None else cb[0],
                         dbia_inverted_left_collinear=None if cb is None else cb[1],
                         dbia_inverted_right_collinear=None if cb is None else cb[2]))
    r = pd.DataFrame(rows)
    r.to_csv(f'{RES}/INV_synteny_check.tsv', sep='\t', index=False)
    log = [f'robust events {len(ev)}; with synteny block boundaries (within {TOL // 1000} kb) at both ends in D. suzukii '
           f'or D. biarmipes comparison: {len(r)}']
    for x in r.itertuples():
        log.append(f'{x.arm} {x.chrom}:{x.start}-{x.end} ({x.L_mb} Mb) pool pairs {x.pool_pairs} f {x.f_pairs} | split {x.pool_split} '
                   f'f {x.f_split} | Dsuz: blocks {x.dsuz_blocks_at_both_ends}, reference collinear '
                   f'{x.dsuz_reference_adjacency_collinear}, inverted collinear left {x.dsuz_inverted_left_collinear} right '
                   f'{x.dsuz_inverted_right_collinear} | Dbia: blocks {x.dbia_blocks_at_both_ends}, reference collinear '
                   f'{x.dbia_reference_adjacency_collinear}, inverted collinear left {x.dbia_inverted_left_collinear} right '
                   f'{x.dbia_inverted_right_collinear}')
    open(f'{RES}/INV_synteny_check.log', 'w').write('\n'.join(log) + '\n')
    print('\n'.join(log))


if __name__ == '__main__':
    main()
