#!/usr/bin/env python3
"""Table of solid gene copy-number changes between D. subpulchrella and D. suzukii, by lineage.

Inputs: scripts/137 outputs (paralog-group calls results/C2_paralog_groups.tsv; per-gene and locus tables
in Data_scratch/.../gene_family_cn) and the read-check functions of scripts/138.

Estimators of copy number. Loci: miniprot loci in Dbia, Dsub, DsuzI (IsoJpt1.0), DsuzR (RU_1.0). Depth: sum of
locus read-depth weights (1 = one copy per haploid genome) in the Dsub pool, the Dsub reference strain, each of the
four D. suzukii pools (on IsoJpt1.0) and their median. Dmel: annotated genes.

A. Copy-number differences between the species (units: paralog groups, >= 40% identity; scripts/137 call
   'supported'). Solid if every estimator gives the same direction: locus count Dsub vs IsoJpt1.0 and vs RU_1.0
   (non-zero), depth Dsub pool vs Dsuz median and Dsub reference strain vs Dsuz median (>= 0.5 copies), and
   Dsub pool vs each Dsuz pool (>= 0.5 copies, >= 3 of 4 pools). Size = smallest absolute difference over these
   estimators (depth values rounded). Lineage: the species whose count differs from Dbia. Because Dsub has one
   assembly and absent loci cannot be tested with reads, a difference with fewer copies in Dsub is solid only when
   Dbia has the Dsub count (D. suzukii gain).
B. Changes on the branch shared by the two species (units: paralog groups without a Dsub-vs-Dsuz call). Solid if
   Dsub, IsoJpt1.0 and RU_1.0 have the same locus count n, every depth estimate is within max(0.5, 2 SD sqrt(n)) of
   n (SD = per-copy robust SD of that sample), and the Dbia locus count equals the D. melanogaster gene count and
   differs from n.
C. Loss of intact copies (pseudogenisation; units: single D. melanogaster genes with the same locus count n in all
   four genomes, and depth within tolerance of n in Dsub and Dsuz). Between species: one species has fewer intact
   loci than the other and than Dbia (both D. suzukii assemblies agree for D. suzukii). Shared: Dsub and both Dsuz
   assemblies have fewer intact loci than Dbia, which has n intact loci. Each disrupted locus of the losing
   species is tested at the level of its disrupting mutations (v3; scripts/142): every frameshift and
   premature stop codon of the seed alignment is located in the genome; it is excluded as a possible artefact if it
   lies within 10 nt of an intron boundary, is compensated by a nearby frameshift, lies in the first or last 10%
   of the protein, is avoided by another query's alignment with >= 0.9 of the seed score, or is not confirmed by
   DNA comparison with the intact copies of D. biarmipes and the sister species (the intact copy must lack it);
   it is population-supported if its allele frequency is >= 0.8 (frameshift) or >= 0.9 (stop) in the Dsub pool
   (with the reference-strain reads agreeing with the assembly) or in >= 3 of the 4 D. suzukii pools. A locus is
   validated if >= 1 mutation passes both. Annotation control: the NCBI model at the locus is full length and of
   normal quality -> not validated. Solid if the validated disrupted loci are >= the difference in intact copies
   (shared: in both species).
   (v2, superseded: locus-level read check without locating the mutation; results/*_SUPERSEDED_v2.*)
Origin of gained copies (tandem, dispersed, retroposition): scripts/144, columns gain_origin and gain_origin_detail.
Outputs: results/C2_solid_copy_number_changes.tsv (solid only),
         results/C2_copy_number_changes_all_tiers.tsv (all candidates with tier and reason),
         results/C2_solid_readcheck_loci.tsv (per disrupted locus), results/C2_pseudogene_events.tsv
         (per mutation), results/C2_frameshift_estimator_calibration.tsv, results/C2_solid_changes.log
Run with envs/popgen/bin/python. No random numbers are used.
"""
import importlib.util
from collections import Counter

import numpy as np
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
OUT = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn'
POOLS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
MAX_LIST = 30
MIN_D = 0.5
SHARED = 'shared (D. suzukii + D. subpulchrella)'
FOCUS_CATS = ['Or', 'Gr_sugar', 'Gr_other', 'Ir', 'Obp', 'ppk', 'mechanosensory_channel', 'Cyp', 'Gst', 'Ugt',
              'carboxylesterase', 'cuticle_structural', 'chitin_metabolism', 'pigmentation_sclerotization',
              'genitalia_development', 'eggshell', 'immune_GO0006955']

spec = importlib.util.spec_from_file_location('rc', f'{PROJ}/scripts/138_disruption_read_check.py')
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)
log = open(f'{RES}/C2_solid_changes.log', 'w')


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + '\n')


P = pd.read_csv(f'{RES}/C2_paralog_groups.tsv', sep='\t', index_col=0)
E = pd.read_csv(f'{OUT}/C2_per_gene.tsv.gz', sep='\t', index_col=0, low_memory=False)
LO = rc.load_loci()
dmg = pd.read_csv(f'{OUT}/Dmel_genes.tsv.gz', sep='\t', low_memory=False).set_index('gene')
LO['pgroup'] = LO.dmel.map(E.pgroup)


# TE association: fraction of each locus' exon bp inside RepeatMasker annotations of the project TE library
# (te/lib/rm_{Dsub,Dsuz}/*.fa.out; EDTA library, unclassified, run on RU_Dsub_v1.1 and IsoJpt1.0; analysis E).
# A unit is TE-associated if >= 50% of its Dsub and IsoJpt1.0 loci have >= 50% of exon bp in TE annotation:
# its "copies" are transposable-element sequences recognised by a domesticated or TE-derived protein (e.g. Iris).
def te_fraction():
    T = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/te/lib'
    res = {}
    for g, f in [('Dsub', 'rm_Dsub/Dsub.fa.out'), ('DsuzI', 'rm_Dsuz/Dsuz.fa.out')]:
        rm = pd.read_csv(f'{T}/{f}', sep=r'\s+', skiprows=3, header=None, usecols=[4, 5, 6],
                         names=['seqid', 's', 'e'], engine='python')
        iv = {}
        for sq, x in rm.groupby('seqid'):
            st, en = [], []
            for s_, e_ in zip(x.sort_values('s').s.values - 1, x.sort_values('s').e.values):
                if st and s_ <= en[-1]:
                    en[-1] = max(en[-1], e_)
                else:
                    st.append(s_)
                    en.append(e_)
            iv[sq] = (np.array(st), np.array(en))
        bed = pd.read_csv(f'{OUT}/mp/{g}_loci_exons.bed', sep='\t', header=None, names=['seqid', 'start', 'end', 'name'])
        ov = []
        for sq, a, b in zip(bed.seqid, bed.start, bed.end):
            if sq not in iv:
                ov.append(0)
                continue
            st, en = iv[sq]
            lo, hi = np.searchsorted(en, a, side='right'), np.searchsorted(st, b)
            ov.append(sum(max(0, min(b, en[k]) - max(a, st[k])) for k in range(lo, hi)))
        bed['ov'], bed['len'] = ov, bed.end - bed.start
        f_ = bed.groupby('name').agg(ov=('ov', 'sum'), len=('len', 'sum'))
        res[g] = f_.ov / f_.len
    out = pd.Series(np.nan, index=LO.index)
    for g in res:
        m = LO.genome == g
        out[m] = LO.loc[m, 'name'].map(res[g]).values
    return out


LO['te_frac'] = te_fraction()
one_ = E.index[(E[['L_Dbia', 'L_Dsub', 'L_DsuzI', 'L_DsuzR']] == 1).all(axis=1)]
xx = LO[LO.genome.isin(['Dsub', 'DsuzI'])]
say(f'loci with >= 50% of exon bp in TE annotation: all {(xx.te_frac >= 0.5).mean():.4f}; loci of genes with one locus '
    f'in every genome {(xx[xx.dmel.isin(one_)].te_frac >= 0.5).mean():.4f}')


def te_assoc(locsel):
    """(fraction of Dsub + IsoJpt1.0 loci with >= 50% of exon bp in TE annotation, largest number of different
    sequences carrying loci in one of these assemblies). TE annotation alone is not decisive: the de novo library
    also covers two tandem host-gene clusters (Or67a, cecropins); TE copies are in addition dispersed."""
    y = locsel[locsel.genome.isin(['Dsub', 'DsuzI'])]
    if not len(y):
        return 0.0, 0
    return float((y.te_frac >= 0.5).mean()), int(y.groupby('genome').tname.nunique().max())


for c in ['W_Dsub_ref'] + [f'W_{p}' for p in POOLS]:
    P[c] = E.groupby('pgroup')[c].sum(min_count=1).reindex(P.index)
say('C2 solid copy-number changes')
say(f'paralog groups (>= 40% identity, >= 60% coverage) with >= 1 locus: {len(P)}')

# per-copy robust SD per sample, from loci of genes with one intact locus in all four genomes (placed chromosomes)
one = E.index[(E[['L_Dbia', 'L_Dsub', 'L_DsuzI', 'L_DsuzR']] == 1).all(axis=1) &
              (E[['Lint_Dbia', 'Lint_Dsub', 'Lint_DsuzI', 'Lint_DsuzR']] == 1).all(axis=1)]
ref_loci = LO[LO.dmel.isin(one) & LO.tname.str.startswith('NC_')]
SD = {}
for col, g in [('w', 'Dsub'), ('w_ref', 'Dsub'), ('w', 'DsuzI')] + [(f'w_{p}', 'DsuzI') for p in POOLS]:
    v = ref_loci.loc[ref_loci.genome == g, col].dropna()
    SD[(col, g)] = (v.quantile(0.9) - v.quantile(0.1)) / 2.563
say('per-copy robust SD (10-90% range / 2.563) over loci of genes with one intact locus in all genomes: ' +
    ', '.join(f'{g} {c}: {s:.3f}' for (c, g), s in SD.items()))
SD_SUB, SD_REF, SD_SUZ = SD[('w', 'Dsub')], SD[('w_ref', 'Dsub')], SD[('w', 'DsuzI')]


def tol(n, sd):
    return max(MIN_D, 2 * sd * np.sqrt(max(n, 1)))


# ------------------------------------------------------------------ A. copy-number differences
rows = []
for i in P.index[P.call_tot == 'supported']:
    r = P.loc[i]
    dW = r.W_Dsub - r.W_DsuzI
    s = np.sign(dW)
    dL_I, dL_R, dW_ref = r.L_Dsub - r.L_DsuzI, r.L_Dsub - r.L_DsuzR, r.W_Dsub_ref - r.W_DsuzI
    dW_p = [r.W_Dsub - r[f'W_{p}'] for p in POOLS]
    n_pool = sum(np.sign(x) == s and abs(x) >= MIN_D for x in dW_p if pd.notna(x))
    reasons = []
    if np.sign(dL_I) != s or np.sign(dL_R) != s:
        reasons.append(f'locus counts do not agree in direction (Dsub {r.L_Dsub}, IsoJpt1.0 {r.L_DsuzI}, RU_1.0 {r.L_DsuzR})')
    if not (np.sign(dW_ref) == s and abs(dW_ref) >= MIN_D):
        reasons.append(f'reference strain ({r.W_Dsub_ref:.2f}) does not confirm the pool ({r.W_Dsub:.2f})')
    if n_pool < 3:
        reasons.append(f'difference in only {n_pool} of 4 D. suzukii pools')
    size = int(min(abs(dL_I), abs(dL_R), round(abs(dW)), round(abs(dW_ref)) if pd.notna(dW_ref) else 99))
    n_sub = r.L_Dsub if abs(r.W_Dsub - r.L_Dsub) < tol(r.L_Dsub, SD_SUB) else round(r.W_Dsub)
    n_suz = r.L_DsuzR if r.L_DsuzI == r.L_DsuzR else round(r.W_DsuzI)
    b = r.L_Dbia
    if b == n_suz:
        lin, dirn = 'D. subpulchrella', ('gain' if n_sub > b else 'loss')
    elif b == n_sub:
        lin, dirn = 'D. suzukii', ('gain' if n_suz > b else 'loss')
    else:
        lin, dirn = 'unresolved', ('more copies in D. subpulchrella' if s > 0 else 'more copies in D. suzukii')
    if s < 0 and lin != 'D. suzukii':
        reasons.append('fewer copies in D. subpulchrella and Dbia differs from D. subpulchrella '
                       '(single Dsub assembly; absent loci cannot be tested)')
    if size < 1:
        reasons.append('smallest estimated difference rounds to 0')
    rows.append(dict(unit='paralog group', pgroup=i, dmel=None, change_type='copy number', branch=lin,
                     direction=dirn, size=max(size, 0), solid=not reasons, reason='; '.join(reasons)))
A = pd.DataFrame(rows)
say(f'\nA. supported Dsub-vs-Dsuz copy-number differences (scripts/137): {len(A)}; solid {int(A.solid.sum())}')
say('   solid by branch/direction: ' + str(dict(Counter(zip(A.branch[A.solid], A.direction[A.solid])))))
rc_ = Counter(x.split(' (')[0] for v in A.reason[~A.solid] for x in v.split('; '))
say('   reasons (a difference can fail several): ' + str(dict(rc_)))

# ------------------------------------------------------------------ B. shared-branch copy-number changes
rowsB = []
nocall = P.call_tot.isna() | P.call_tot.isin(['', 'small'])
for i in P.index[nocall & (P.L_Dsub == P.L_DsuzI) & (P.L_DsuzI == P.L_DsuzR) & (P.L_Dbia == P.ann_Dmel) &
                 (P.L_Dsub != P.L_Dbia) & ~P.histone.astype(bool)]:
    r = P.loc[i]
    n = r.L_Dsub
    ok = (abs(r.W_Dsub - n) <= tol(n, SD_SUB) and pd.notna(r.W_Dsub_ref) and abs(r.W_Dsub_ref - n) <= tol(n, SD_REF)
          and abs(r.W_DsuzI - n) <= tol(n, SD_SUZ))
    rowsB.append(dict(unit='paralog group', pgroup=i, dmel=None, change_type='copy number', branch=SHARED,
                      direction='gain' if n > r.L_Dbia else 'loss', size=int(abs(n - r.L_Dbia)), solid=bool(ok),
                      reason='' if ok else 'read depth not consistent with the locus count in Dsub or Dsuz'))
B = pd.DataFrame(rowsB)
say(f'\nB. shared-branch candidates (equal loci in Dsub and both Dsuz assemblies; Dbia = Dmel, different): {len(B)}; '
    f'solid {int(B.solid.sum())} ({dict(Counter(B.direction[B.solid]))})')

# ------------------------------------------------------------------ C. intact copies (per gene)
Lc = ['L_Dbia', 'L_Dsub', 'L_DsuzI', 'L_DsuzR']
eq = E[(E[Lc].nunique(axis=1) == 1) & (E.L_Dsub >= 1)].copy()
eq['n'] = eq.L_Dsub
eq = eq[[abs(w - n) <= tol(n, SD_SUB) and abs(z - n) <= tol(n, SD_SUZ) for w, z, n in zip(eq.W_Dsub, eq.W_DsuzI, eq.n)]]
ib, is_, iz, ir = eq.Lint_Dbia, eq.Lint_Dsub, eq.Lint_DsuzI, eq.Lint_DsuzR
cand = {}
for g in eq.index[(is_ < iz) & (iz == ir) & (iz == ib)]:
    cand[g] = ('D. subpulchrella', ['Dsub'], int(iz[g] - is_[g]))
for g in eq.index[(iz == ir) & (iz < is_) & (is_ == ib)]:
    cand[g] = ('D. suzukii', ['DsuzI'], int(is_[g] - iz[g]))
for g in eq.index[(is_ == iz) & (iz == ir) & (ib == eq.n) & (is_ < ib)]:
    cand[g] = (SHARED, ['Dsub', 'DsuzI'], int(ib[g] - is_[g]))
sel = LO[~LO.intact & LO.genome.isin(['Dsub', 'DsuzI'])]
sel = sel[[g in cand and gn in cand[g][1] for g, gn in zip(sel.dmel, sel.genome)]].copy()
say(f'\nC. intact-copy candidates (genes with equal locus counts in all four genomes): Dsub '
    f'{sum(v[0] == "D. subpulchrella" for v in cand.values())}, Dsuz {sum(v[0] == "D. suzukii" for v in cand.values())}, '
    f'shared {sum(v[0] == SHARED for v in cand.values())}; disrupted loci to read-check {len(sel)}')
# each disrupting mutation is located and tested (scripts/142): artefact filters (intron
# boundary, compensating frameshift, first or last 10% of the protein, a near-equal alternative alignment without
# the event, DNA comparison with the intact copies of D. biarmipes and the sister species) and allele frequency in
# the population pools. A locus is validated if >= 1 of its mutations passes all filters and is at high frequency in
# the pool(s). This replaces the locus-level read check of v2 (results/*_SUPERSEDED_v2.*), which did not locate the
# mutation and counted a disruption as confirmed whenever the reads did not contradict the assembly.
spec_ev = importlib.util.spec_from_file_location('ev', f'{PROJ}/scripts/142_pseudogene_event_support.py')
ev = importlib.util.module_from_spec(spec_ev)
spec_ev.loader.exec_module(ev)
CAL = ev.calibrate()
CAL.to_csv(f'{RES}/C2_frameshift_estimator_calibration.tsv', sep='\t', index=False, float_format='%.3f')
for smp in ['Dsub_L5_7', 'refstrain'] + POOLS:
    x_ = CAL[f'freq_{smp}'].dropna()
    say(f'   calibration {smp:10s}: {len(x_)} control sites; frameshift-allele estimate median {x_.median():.3f}, '
        f'< 0.9 at {(x_ < 0.9).mean():.1%}, < 0.8 at {(x_ < 0.8).mean():.1%}')
EV = ev.run(list(sel.name))
EV.to_csv(f'{RES}/C2_pseudogene_events.tsv', sep='\t', index=False, float_format='%.3f')
say(f'   events (frameshifts and stop codons in the seed alignments): {len(EV)} in {EV.locus.nunique()} of {len(sel)} '
    f'loci; passing the artefact filters {int((EV.artefact_filter == "").sum())}; '
    f'population-supported {int((EV.population == "supported").sum())}; both {int(EV.supported.sum())}')
say('   artefact filters (an event can fail several): ' + str(dict(Counter(
    x for v in EV.artefact_filter.fillna('') for x in v.split('; ') if x))))
say('   population: ' + str(dict(Counter(EV.population))))


def freq_txt(r):
    if r.genome == 'Dsub':
        return f'pool {r.freq_Dsub_L5_7:.2f}, reference strain {r.freq_refstrain:.2f}'
    return 'pools ' + ','.join(f'{r[f"freq_{p}"]:.2f}' for p in POOLS)


def event_txt(r):
    head = f'{r.type} at codon {r.qpos}/{r.qlen} ({r.seqid}:{r.g0})'
    if r.supported:
        return f'{head}: supported; {freq_txt(r)}; DNA: {r.dna_test}'
    why = r.artefact_filter if isinstance(r.artefact_filter, str) and r.artefact_filter else r.population
    return f'{head}: not supported ({why}; {freq_txt(r)})'


S = sel[['name', 'genome', 'dmel', 'pgroup', 'symbol', 'seed_fs', 'seed_st', 'seed_cov']].rename(columns={'name': 'locus'})
S['n_events'] = S.locus.map(EV.groupby('locus').size()).fillna(0).astype(int)
S['n_supported_events'] = S.locus.map(EV[EV.supported].groupby('locus').size()).fillna(0).astype(int)
S['shared_supported_event'] = S.locus.map(EV[EV.supported].groupby('locus').shared_with_sister.any()).fillna(False)
S['validated'] = S.n_supported_events >= 1
evt = EV.groupby('locus').apply(lambda x: ' | '.join(event_txt(r) for _, r in x.iterrows()), include_groups=False)
S['validation'] = [evt.get(l_, 'no frameshift or stop codon in the alignment (truncated; not testable)')
                   for l_ in S.locus]
S.loc[S.validated, 'validation'] = [' | '.join(event_txt(r) for _, r in EV[EV.supported & (EV.locus == l_)].iterrows())
                                    for l_ in S.locus[S.validated]]
S.to_csv(f'{RES}/C2_solid_readcheck_loci.tsv', sep='\t', index=False, float_format='%.3f')
say(f'   loci validated (>= 1 supported mutation): {int(S.validated.sum())} of {len(S)} '
    f'({dict(Counter(S.genome[S.validated]))})')
# Annotation control. The read check tests the assembly sequence, not the miniprot exon structure: a misplaced
# splice junction in the alignment creates a frameshift or a stop codon that reads confirm. If the NCBI annotation
# of the same genome (built with RNA evidence) reports a protein of normal quality (no "LOW QUALITY PROTEIN" tag)
# covering >= 80% of the D. melanogaster protein at the disrupted locus, the disruption is not counted as solid.
asg = pd.read_csv(f'{OUT}/C2_gene_assignments.tsv.gz', sep='\t', low_memory=False)
hits = []
for sp in ['Dsub', 'DsuzI']:
    h = pd.read_csv(f'{OUT}/{sp}_vs_Dmel.tsv', sep='\t',
                    names=['q', 's', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore'])
    h = h[h.evalue <= 1e-5].sort_values(['q', 'bitscore'], ascending=[True, False]).drop_duplicates('q')
    hits.append(pd.DataFrame({'species': sp, 'gene': h.q.str.split('|').str[1], 'scov': h.length / h.slen}))
asg = asg.merge(pd.concat(hits), on=['species', 'gene'], how='left')
asg['normal'] = ~asg.name.fillna('').str.startswith('LOW QUALITY PROTEIN') & (asg.scov >= 0.8)
normal = set(zip(asg.species[asg.normal], asg.gene[asg.normal], asg.dmel[asg.normal]))
lo_ann = LO.set_index('name').annotated_cds


def annotated_intact(row):
    genes_ = [g for g in str(lo_ann.get(row.locus, '') or '').split(',') if g and g != 'nan']
    return any((row.genome, g, row.dmel) in normal for g in genes_)


S['annotation_reports_intact'] = S.apply(annotated_intact, axis=1)
S.loc[S.validated & S.annotation_reports_intact, 'validation'] += \
    '; NCBI model at this locus is full length and of normal quality (possible alignment artefact)'
S['validated'] = S.validated & ~S.annotation_reports_intact
S.to_csv(f'{RES}/C2_solid_readcheck_loci.tsv', sep='\t', index=False, float_format='%.3f')
say(f'   annotation control: {int(S.annotation_reports_intact.sum())} of {len(S)} disrupted loci carry a full-length NCBI '
    f'model of normal quality; validated after the control: {int(S.validated.sum())}')
nval = S[S.validated].groupby(['genome', 'dmel']).size()
note = S.groupby(['genome', 'dmel']).validation.agg(lambda x: ' | '.join(x))
rowsC = []
for g, (lin, gens, d) in cand.items():
    nv = [int(nval.get((gn, g), 0)) for gn in gens]
    ok = min(nv) >= d
    txt = '; '.join(f'{gn}: {note.get((gn, g), "no disrupted locus checked")}' for gn in gens)
    if lin == SHARED and ok:
        sh_ = bool(S.shared_supported_event[S.validated & (S.dmel == g)].any())
        txt += ('; the same mutation is present in both species' if sh_ else
                '; the supported mutations differ between the species')
    rowsC.append(dict(unit='gene', pgroup=E.pgroup[g], dmel=g, change_type='intact copies (pseudogenisation)',
                      branch=lin, direction='loss of intact copy', size=d, solid=ok,
                      reason='' if ok else f'validated disrupted loci {nv} < {d}', readcheck=txt,
                      validated_disrupted_loci=min(nv)))
C = pd.DataFrame(rowsC)
say(f'   solid: {dict(Counter(C.branch[C.solid]))} of {dict(Counter(C.branch))}')


# ------------------------------------------------------------------ assemble table
def loci_list(sel_):
    x = sel_.sort_values(['tname', 'ts'])
    items = [f'{t}:{a + 1}-{b}{s}{"" if it else "[disrupted]"}' for t, a, b, s, it in
             zip(x.tname, x.ts, x.te, x.strand, x.intact)]
    return ';'.join(items[:MAX_LIST]) + (f';...(+{len(items) - MAX_LIST})' if len(items) > MAX_LIST else '')


def gene_ids(sel_):
    ids = sorted({y.replace('gene-', '') for v in sel_.annotated_cds.dropna() for y in str(v).split(',') if y})
    return ','.join(ids[:MAX_LIST]) + (f',...(+{len(ids) - MAX_LIST})' if len(ids) > MAX_LIST else '')


T = pd.concat([A, B, C], ignore_index=True)
out = []
for r in T.itertuples():
    if r.unit == 'gene':
        mem = [r.dmel]
        q = P.loc[r.pgroup] if r.pgroup in P.index else None
        src = E.loc[r.dmel]
        locsel = LO[LO.dmel == r.dmel]
    else:
        mem = list(E.index[E.pgroup == r.pgroup])
        src = P.loc[r.pgroup]
        locsel = LO[LO.pgroup == r.pgroup]
    sy = [str(E.symbol[g]) for g in mem]
    changed = [f'{E.symbol[g]}:{E.L_Dbia[g]}/{E.L_Dsub[g]}/{E.L_DsuzI[g]}/{E.L_DsuzR[g]}' for g in mem
               if len({E.L_Dbia[g], E.L_Dsub[g], E.L_DsuzI[g], E.L_DsuzR[g]}) > 1
               or len({E.Lint_Dbia[g], E.Lint_Dsub[g], E.Lint_DsuzI[g], E.Lint_DsuzR[g]}) > 1]
    cats = sorted({c for g in mem for c in str(E.categories[g] if pd.notna(E.categories[g]) else '').split(',') if c})
    tef, nseq = te_assoc(locsel)
    solid, reason = bool(r.solid), r.reason
    if tef >= 0.5 and nseq >= 5:
        solid = False
        reason = '; '.join(x for x in [reason, f'TE-associated ({tef:.0%} of Dsub/Dsuz loci inside TE annotation, '
                                                f'loci on {nseq} sequences)'] if x)
    rec = dict(change_id='', change_type=r.change_type, branch=r.branch, direction=r.direction, size_copies=r.size,
               solid=solid, reason_not_solid=reason, unit=r.unit, te_loci_fraction=round(tef, 2),
               max_sequences_with_loci=nseq,
               dmel_genes=','.join(sy[:MAX_LIST]) + (f',...(+{len(sy) - MAX_LIST})' if len(sy) > MAX_LIST else ''),
               n_dmel_genes=len(mem),
               members_with_differences_Dbia_Dsub_DsuzI_DsuzR=','.join(changed[:MAX_LIST]),
               dmel_CG=','.join(g.replace('gene-Dmel_', '') for g in mem[:MAX_LIST]),
               dmel_FBgn=','.join(str(dmg.fbgn.get(g, '')) for g in mem[:MAX_LIST] if pd.notna(dmg.fbgn.get(g))),
               categories=','.join(cats),
               dmel_annotated=int(src.ann_Dmel),
               loci_Dbia=int(src.L_Dbia), loci_Dsub=int(src.L_Dsub), loci_DsuzI=int(src.L_DsuzI),
               loci_DsuzR=int(src.L_DsuzR), intact_Dbia=int(src.Lint_Dbia), intact_Dsub=int(src.Lint_Dsub),
               intact_DsuzI=int(src.Lint_DsuzI), intact_DsuzR=int(src.Lint_DsuzR),
               cn_Dsub_pool=round(src.W_Dsub, 2), cn_Dsub_refstrain=round(src.W_Dsub_ref, 2),
               cn_Dsuz_median4pools=round(src.W_DsuzI, 2), **{f'cn_{p}': round(src[f'W_{p}'], 2) for p in POOLS},
               intact_cn_Dsub=round(src.Wint_Dsub, 2), intact_cn_Dsuz=round(src.Wint_DsuzI, 2),
               pseudogene_evidence=getattr(r, 'readcheck', '') if isinstance(getattr(r, 'readcheck', ''), str) else '',
               loci_coordinates_Dbia=loci_list(locsel[locsel.genome == 'Dbia']),
               loci_coordinates_Dsub=loci_list(locsel[locsel.genome == 'Dsub']),
               loci_coordinates_DsuzI=loci_list(locsel[locsel.genome == 'DsuzI']),
               loci_coordinates_DsuzR=loci_list(locsel[locsel.genome == 'DsuzR']),
               annotated_genes_Dbia=gene_ids(locsel[locsel.genome == 'Dbia']),
               annotated_genes_Dsub=gene_ids(locsel[locsel.genome == 'Dsub']),
               annotated_genes_DsuzI=gene_ids(locsel[locsel.genome == 'DsuzI']),
               annotated_genes_DsuzR=gene_ids(locsel[locsel.genome == 'DsuzR']),
               paralog_group=int(r.pgroup))
    out.append(rec)
O = pd.DataFrame(out)
order = {'D. subpulchrella': 0, 'D. suzukii': 1, SHARED: 2, 'unresolved': 3}
O['_o'] = O.branch.map(order)
O['_f'] = O.categories.apply(lambda c: 0 if any(x in FOCUS_CATS for x in c.split(',')) else 1)
O = O.sort_values(['solid', 'change_type', '_o', 'direction', '_f', 'dmel_genes'],
                  ascending=[False, True, True, True, True, True]).drop(columns=['_o', '_f'])
O['change_id'] = [f'C2-{k + 1:04d}' for k in range(len(O))]
# origin of gained copies (tandem, dispersed, retroposition) for solid gains; scripts/144
spec_og = importlib.util.spec_from_file_location('og', f'{PROJ}/scripts/144_duplication_origin.py')
og = importlib.util.module_from_spec(spec_og)
spec_og.loader.exec_module(og)
O = og.annotate(O, LO, mask=O.solid.to_numpy())
O.to_csv(f'{RES}/C2_copy_number_changes_all_tiers.tsv', sep='\t', index=False)
SOL = O[O.solid].drop(columns=['solid', 'reason_not_solid'])
SOL.to_csv(f'{RES}/C2_solid_copy_number_changes.tsv', sep='\t', index=False)
say(f'\nSolid changes: {len(SOL)} (candidates of all tiers {len(O)})')
say('Origin of solid gains (scripts/144): ' + str(dict(Counter(SOL.gain_origin))))
say(SOL.groupby(['change_type', 'branch', 'direction']).agg(n=('size_copies', 'size'), copies=('size_copies', 'sum')).to_string())
foc = SOL[SOL.categories.apply(lambda c: any(x in FOCUS_CATS for x in c.split(',')))]
say(f'\nSolid changes in the pre-specified categories: {len(foc)}')
for r in foc.itertuples():
    say(f'  {r.change_id} {r.change_type[:13]:13s} {r.branch[:16]:16s} {r.direction[:20]:20s} {r.size_copies} | '
        f'{r.members_with_differences_Dbia_Dsub_DsuzI_DsuzR[:70]:70s} | intact {r.intact_Dbia}/{r.intact_Dsub}/'
        f'{r.intact_DsuzI}/{r.intact_DsuzR} | cn {r.cn_Dsub_pool}/{r.cn_Dsub_refstrain}/{r.cn_Dsuz_median4pools} | {r.categories}')
log.close()

# ------------------------------------------------------------------ category test on solid copy-number changes
# Units: paralog groups with >= 1 locus. A group belongs to a category if any member does. Tests: (1) Fisher exact
# against all other groups; (2) Cochran-Mantel-Haenszel test stratified by the number of D. melanogaster genes in the
# group (1, 2-3, >= 4), which is fixed independently of the two ingroup species. (A background defined by >= 2 loci
# in any genome would be partly defined by the outcome, because a gained copy makes a group multi-copy.)
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests
log = open(f'{RES}/C2_solid_changes.log', 'a')
gcat = E.groupby('pgroup').categories.agg(lambda x: {c for v in x.fillna('') for c in v.split(',') if c})
gcat = gcat.reindex(P.index)
stratum = np.where(P.n_dmel_genes == 1, 0, np.where(P.n_dmel_genes <= 3, 1, 2))
solid_cn = set(SOL.paralog_group[SOL.change_type == 'copy number'])
has = P.index.isin(list(solid_cn))
rows = []
for c in FOCUS_CATS:
    mk = np.array([isinstance(s_, set) and c in s_ for s_ in gcat])
    a_, b_ = int((mk & has).sum()), int((mk & ~has).sum())
    c_, d_ = int((~mk & has).sum()), int((~mk & ~has).sum())
    orr, pf = fisher_exact([[a_, b_], [c_, d_]])
    tabs = []
    for k in range(3):
        m = stratum == k
        t = np.array([[(mk & has & m).sum(), (mk & ~has & m).sum()], [(~mk & has & m).sum(), (~mk & ~has & m).sum()]])
        if t[0].sum() > 0:  # strata in which the category is present; no continuity correction of the cells
            tabs.append(t.astype(float))
    if a_ > 0:
        st_ = StratifiedTable(tabs)
        or_mh, p_mh = st_.oddsratio_pooled, st_.test_null_odds(correction=True).pvalue
    elif tabs:  # no changed group in the category: pooled OR 0; P from the same test
        st_ = StratifiedTable(tabs)
        or_mh, p_mh = 0.0, st_.test_null_odds(correction=True).pvalue
    else:
        or_mh, p_mh = np.nan, 1.0
    rows.append(dict(category=c, changed=a_, groups=a_ + b_, bg_changed=c_, bg_groups=c_ + d_, OR_all=orr, P_all=pf,
                     OR_MH_by_Dmel_family_size=or_mh, P_MH=p_mh))
CT = pd.DataFrame(rows)
CT['BH_all'] = multipletests(CT.P_all, method='fdr_bh')[1]
CT['BH_MH'] = multipletests(CT.P_MH.fillna(1), method='fdr_bh')[1]
CT.to_csv(f'{RES}/C2_solid_categories.tsv', sep='\t', index=False, float_format='%.4g')
say('\nCategory test, solid copy-number changes per paralog group: Fisher vs all other groups; Mantel-Haenszel '
    'stratified by D. melanogaster family size (1, 2-3, >= 4 genes); BH over categories')
for r in CT.itertuples():
    say(f'  {r.category:28s} {r.changed}/{r.groups} vs {r.bg_changed}/{r.bg_groups}; all: OR {r.OR_all:.2f} '
        f'BH {r.BH_all:.3g}; stratified: OR {r.OR_MH_by_Dmel_family_size:.2f} P {r.P_MH:.3g} BH {r.BH_MH:.3g}')
log.close()
