#!/usr/bin/env python3
"""Origin of gained gene copies: tandem duplication, dispersed (DNA-based) duplication or retroposition.

Question: for each solid copy-number gain (scripts/139), how did the additional copy arise?
Units: rows of the solid table with direction 'gain' (or 'more copies in D. subpulchrella', lineage unresolved).
Genomes: the gaining genome(s) X and a reference genome R that has the ancestral copy number:
  D. subpulchrella gain  X = Dsub;                 R = DsuzI (sister species; same count as Dbia)
  D. suzukii gain        X = DsuzI and DsuzR;      R = Dsub
  shared-branch gain     X = Dsub, DsuzI, DsuzR;   R = Dbia
  unresolved (Dsub > Dsuz, Dbia differs) X = Dsub; R = DsuzI (copies additional relative to D. suzukii)
Loci: miniprot loci of the paralog group (scripts/137), in X and R.
1. New copies. Synteny score of a locus = number of D. melanogaster genes shared between its flanks (the 5 nearest
   loci on each side on the same sequence, loci of the same group excluded) and the flanks of a locus in R. X loci
   are matched one-to-one to R loci in order of decreasing score (score >= 2). The unmatched X loci are the
   candidate new copies; when there are more than g = (loci in X) - (loci in R), the g with the lowest best score
   are taken and the identification is marked uncertain. (Tandem copies share flanks with their parent; one of
   them is matched and the other is left as the new copy, which is correct for the classification below.)
2. Clusters. New copies on the same sequence with gaps <= 50 kb form a cluster. Parent: the ancestral (matched)
   X locus of the same D. melanogaster gene, nearest on the same sequence (otherwise any ancestral locus).
3. Introns. For each new copy and the parent, the alignments of one common query (the D. melanogaster protein of
   the new copy's gene if both loci have it; otherwise the common query with the highest lower score) are walked
   (miniprot CIGAR; N, U, V = introns), and introns are counted within the query range covered by both
   alignments (comparisons covering < 50% of the query are not used).
4. Mechanism of each cluster:
   - tandem: <= 50 kb (edge to edge) from an ancestral copy on the same sequence; 'tandem or retroposition
     (adjacent copies differ in introns)' when the new copies are intronless and the parent is not (which of two
     adjacent copies is the new one cannot be decided from synteny);
   - retroposition: not adjacent to an ancestral copy, all new copies intronless and the parent with >= 1 intron;
   - dispersed (DNA-based): not adjacent, introns retained, on an assembled chromosome (RefSeq NC_ sequence) as is
     the parent; 'segmental, k genes' when the cluster contains copies of k different D. melanogaster genes
     (duplicated together); 'partial intron loss' when a copy has fewer introns than the parent but >= 1;
   - '+ tandem': the cluster has >= 2 copies of the same gene away from the ancestral copy (a dispersed or
     retroposed founding copy, then tandem duplication);
   - not determined (unplaced scaffold): not adjacent, introns retained, and the copy or its parent lies on an
     unplaced scaffold (NW_), so its position relative to the parent is unknown;
   - no parent copy in this genome: the group has no ancestral locus in X;
   - '[parent intronless: retroposition not testable]' is added when the parent has no intron.
   The summary column gives the calls per genome (with the number of copies) when the genomes differ.
No random numbers. Run as a module from scripts/139 (annotate()), or as a script on the solid table for a check.
"""
import importlib.util
import re
from collections import Counter

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('ev', f'{PROJ}/scripts/142_pseudogene_event_support.py')
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)
ev.OWN['DsuzR'] = 'DsuzI'  # as scripts/137: DsuzI proteins are not used on either D. suzukii assembly

K_FLANK, MIN_SYN, TANDEM_BP, MIN_SHARED_Q = 5, 2, 50000, 0.5
SHARED = 'shared (D. suzukii + D. subpulchrella)'
SETUP = {'D. subpulchrella': (['Dsub'], 'DsuzI'), 'D. suzukii': (['DsuzI', 'DsuzR'], 'Dsub'),
         SHARED: (['Dsub', 'DsuzI', 'DsuzR'], 'Dbia'), 'unresolved': (['Dsub'], 'DsuzI')}


def introns_q(cg, qs):
    """Query positions (0-based aa) of the introns of an alignment."""
    q, out = qs, []
    for k, o in re.findall(r'(\d+)([MIDFGNUV])', cg):
        k = int(k)
        if o in 'MI':
            q += k
        elif o == 'G':
            q += 1
        elif o in 'NUV':
            out.append(q)
            q += 0 if o == 'N' else 1
    return out


class Genomes:
    def __init__(self, LO):
        self.LO = LO.sort_values(['genome', 'tname', 'ts']).reset_index(drop=True)
        dmel_of, full = ev.query_tables()
        self.A = {g: ev.read_alignments(g, dmel_of, full) for g in ['Dbia', 'Dsub', 'DsuzI', 'DsuzR']}
        self.pos = {}  # genome -> {tname: loci frame sorted by ts}
        for (g, t), x in self.LO.groupby(['genome', 'tname']):
            self.pos.setdefault(g, {})[t] = x

    def flanks(self, L):
        x = self.pos[L.genome][L.tname]
        x = x[x.pgroup != L.pgroup]
        up = x[x.ts < L.ts].tail(K_FLANK)
        dn = x[x.ts > L.ts].head(K_FLANK)
        return set(pd.concat([up, dn]).dmel.dropna())

    def between(self, L, M):
        """Loci of other genes between two loci on the same sequence."""
        x = self.pos[L.genome][L.tname]
        a, b = min(L.te, M.te), max(L.ts, M.ts)
        return int(((x.ts >= a) & (x.te <= b) & (x.pgroup != L.pgroup)).sum())

    def alns(self, L):
        A = self.A[L.genome]
        x = A[(A.tname == L.tname) & (A.strand == L.strand)]
        ov = np.minimum(x.te, L.te) - np.maximum(x.ts, L.ts)
        return x[(ov > 0) & (ov >= 0.5 * np.minimum(x.te - x.ts, L.te - L.ts))]

    def intron_compare(self, N, P):
        an, ap = self.alns(N), self.alns(P)
        common = set(an.q) & set(ap.q)
        if not common:
            return None
        pref = 'Dmel|' + str(N.dmel)
        if pref in common:
            q = pref
        else:
            q = max(common, key=lambda q_: min(an.AS[an.q == q_].max(), ap.AS[ap.q == q_].max()))
        a1 = an[an.q == q].sort_values('AS').iloc[-1]
        a2 = ap[ap.q == q].sort_values('AS').iloc[-1]
        lo, hi = max(a1.qs, a2.qs), min(a1.qe, a2.qe)
        if hi <= lo:
            return None
        n1 = sum(lo <= p < hi for p in introns_q(a1.cg, a1.qs))
        n2 = sum(lo <= p < hi for p in introns_q(a2.cg, a2.qs))
        return dict(query=q, new=n1, parent=n2, shared_frac=(hi - lo) / a1.qlen)


def gap(a, b):
    return int(max(0, max(a.ts, b.ts) - min(a.te, b.te)))


def classify_genome(G, pg, x_gen, r_gen):
    """Mechanism calls for the new copies of paralog group pg in genome x_gen (reference r_gen).
    Returns (list of (call, detail) per new-copy cluster, note)."""
    LX = G.LO[(G.LO.genome == x_gen) & (G.LO.pgroup == pg)].reset_index(drop=True)
    LR = G.LO[(G.LO.genome == r_gen) & (G.LO.pgroup == pg)].reset_index(drop=True)
    g = len(LX) - len(LR)
    if g <= 0:
        return [], f'{x_gen}: {len(LX)} loci vs {len(LR)} in {r_gen}; no additional locus to classify'
    fx = [G.flanks(L) for _, L in LX.iterrows()]
    fr = [G.flanks(L) for _, L in LR.iterrows()]
    S = np.array([[len(a & b) for b in fr] for a in fx]) if len(fr) else np.zeros((len(LX), 0), int)
    matched, used = {}, set()
    for i, j in sorted(((i, j) for i in range(len(LX)) for j in range(len(LR))), key=lambda ij: -S[ij]):
        if S[i, j] < MIN_SYN:
            break
        if i in matched or j in used:
            continue
        matched[i] = j
        used.add(j)
    best = S.max(axis=1) if S.shape[1] else np.zeros(len(LX), int)
    un = [i for i in range(len(LX)) if i not in matched]
    uncertain = len(un) > g
    new = sorted(un, key=lambda i: best[i])[:g]
    anc = LX.drop(index=new)
    note = f'{x_gen} vs {r_gen}: {len(LX)} vs {len(LR)} loci, {g} new' + (
        '; new-copy identification uncertain (more loci without synteny than new copies)' if uncertain else '')
    if not len(anc):
        return [('no parent copy in this genome', f'no locus of the group in {r_gen}', g)], note
    # clusters of new copies: same sequence, gaps <= TANDEM_BP (single linkage)
    NX = LX.loc[new].sort_values(['tname', 'ts'])
    clusters, cur = [], []
    for _, N in NX.iterrows():
        if cur and (N.tname != cur[-1].tname or gap(N, cur[-1]) > TANDEM_BP):
            clusters.append(cur)
            cur = []
        cur.append(N)
    clusters.append(cur)
    calls = []
    for cl in clusters:
        same = anc[anc.tname == cl[0].tname]
        d_anc = min((gap(N, M) for N in cl for _, M in same.iterrows()), default=None)
        near_anc = None
        if d_anc is not None:
            near_anc = min(((gap(N, M), M) for N in cl for _, M in same.iterrows()), key=lambda x: x[0])[1]
        # parent: ancestral locus of the same D. melanogaster gene, nearest on the same sequence
        N0 = cl[0]
        pool = anc[anc.dmel == N0.dmel] if (anc.dmel == N0.dmel).any() else anc
        P = sorted(pool.itertuples(), key=lambda r: (r.tname != N0.tname, abs(r.ts - N0.ts)))[0]
        P = anc.loc[P.Index]
        ics = [G.intron_compare(N, P) for N in cl]
        ok_ic = [ic for ic in ics if ic is not None and ic['shared_frac'] >= MIN_SHARED_Q]
        intronless = bool(ok_ic) and all(ic['new'] == 0 for ic in ok_ic) and any(ic['parent'] >= 1 for ic in ok_ic)
        partial = bool(ok_ic) and not intronless and any(0 < ic['new'] < ic['parent'] for ic in ok_ic)
        parent_intronless = bool(ok_ic) and all(ic['parent'] == 0 for ic in ok_ic)
        genes = sorted({str(N.dmel).replace('gene-Dmel_', '') for N in cl})
        placed = cl[0].tname.startswith('NC_')
        if d_anc is not None and d_anc <= TANDEM_BP:
            m = 'tandem or retroposition (adjacent copies differ in introns)' if intronless else 'tandem'
        elif intronless:
            m = 'retroposition'
        elif not placed or not P.tname.startswith('NC_'):
            m = 'not determined (unplaced scaffold)'
        elif len(genes) > 1:
            m = f'dispersed (segmental, {len(genes)} genes)'
        else:
            m = 'dispersed'
        if m == 'dispersed' and partial:
            m = 'dispersed (partial intron loss)'
        if len(cl) > 1 and (d_anc is None or d_anc > TANDEM_BP) and len(genes) == 1:
            m += ' + tandem'
        if parent_intronless and m in ('tandem', 'dispersed', 'dispersed + tandem'):
            m += ' [parent intronless: retroposition not testable]'
        pos = (f'{d_anc / 1000:.1f} kb from the nearest ancestral copy ({near_anc.tname}:{near_anc.ts + 1}-'
               f'{near_anc.te}{near_anc.strand}; {"same" if near_anc.strand == cl[0].strand else "opposite"} '
               f'orientation, {G.between(cl[0], near_anc)} other loci between)'
               if d_anc is not None else 'no ancestral copy on this sequence')
        locs = ','.join(f'{N.tname}:{N.ts + 1}-{N.te}{N.strand}' for N in cl)
        it = '; '.join(f'introns new/parent {ic["new"]}/{ic["parent"]} over {ic["shared_frac"]:.0%} of {ic["query"]}'
                       for ic in ics if ic is not None) or 'introns not compared (no common query)'
        calls.append((m, f'{locs} [{m}; {pos}; parent {P.tname}:{P.ts + 1}-{P.te}; {it}; synteny score with '
                         f'{r_gen} {",".join(str(best[i]) for i in LX.index if i in [N.name for N in cl])}]', len(cl)))
    return calls, note


def annotate(T, LO, mask=None):
    """Adds gain_origin and gain_origin_detail to table T (rows of scripts/139); rows outside mask are 'not
    assessed (not solid)'."""
    G = Genomes(LO)
    org, det = [], []
    for k, r in enumerate(T.itertuples()):
        gain = r.change_type == 'copy number' and (r.direction == 'gain' or r.branch == 'unresolved')
        if mask is not None and not mask[k]:
            org.append('not assessed (not solid)')
            det.append('')
            continue
        if not gain:
            org.append('not applicable' if r.change_type != 'copy number' else 'not applicable (loss)')
            det.append('')
            continue
        xs, rg = SETUP[r.branch]
        per, txt = {}, []
        for xg in xs:
            calls, note = classify_genome(G, r.paralog_group, xg, rg)
            c = Counter()
            for m, _, n_ in calls:
                c[m.split(' [parent')[0]] += n_
            per[xg] = c
            txt.append(note + (': ' + ' | '.join(d for _, d, _ in calls) if calls else ''))
        fmt = {xg: ', '.join(m + (f' ({n} copies)' if n > 1 else '') for m, n in sorted(c.items()))
               for xg, c in per.items() if c}
        if len(set(fmt.values())) == 1:
            org.append(next(iter(fmt.values())))
        elif fmt:
            org.append('; '.join(f'{xg}: {v}' for xg, v in fmt.items()))
        else:
            org.append('not determined')
        det.append(' || '.join(txt))
    T = T.copy()
    i = list(T.columns).index('size_copies') + 1
    T.insert(i, 'gain_origin', org)
    T.insert(i + 1, 'gain_origin_detail', det)
    return T


if __name__ == '__main__':
    import sys
    LO = pd.read_csv(f'{ev.OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
    E = pd.read_csv(f'{ev.OUT}/C2_per_gene.tsv.gz', sep='\t', index_col=0, low_memory=False)
    LO['pgroup'] = LO.dmel.map(E.pgroup)
    T = pd.read_csv(sys.argv[1], sep='\t')
    T = T.drop(columns=[c for c in ['gain_origin', 'gain_origin_detail'] if c in T])
    O = annotate(T, LO)
    O.to_csv(sys.argv[2], sep='\t', index=False)
    print(O.gain_origin.value_counts().to_string())
