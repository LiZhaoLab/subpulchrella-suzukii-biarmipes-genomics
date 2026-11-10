#!/usr/bin/env python3
"""C2 second pass: gene copy number by lineage from miniprot loci weighted by read depth.

Question: which differences in gene copy number between D. subpulchrella (Dsub) and D. suzukii (Dsuz)
are real, on which lineage did they arise, and do any concern chemoreception, detoxification,
mechanosensation, cuticle or eggshell genes (pest status, enlarged ovipositor, substrate preference)?

Inputs: scripts/135 tables (NCBI annotations, DIAMOND best D. melanogaster hit per protein) and the outputs
of scripts/136_copy_number_validation.sbatch (bedcov per gene CDS; miniprot 0.13-r248 GFF/PAF).

Steps
 1. Read depth. Mean depth over the CDS of each annotated gene / median of strict single-copy genes on the
    same chromosome class (A or X) -> normalised depth nd (1 = one copy per haploid genome). Dsub: pool
    Dsub_L5_7 (and the reference strain); Dsuz (IsoJpt1.0): median over JP-Kan, JP-Sap, JP-Sap2, JP-Tok.
 2. miniprot loci. Alignments of proteins of the other species (never the genome's own annotation; D. mel.
    queries plus full-length proteins only, i.e. covering >= 50% of their D. melanogaster best hit) with
    query coverage >= 0.5. Greedy clustering: alignments in decreasing score; an alignment joins an
    existing locus if it overlaps that locus' seed alignment by >= 50% of the shorter span, otherwise it
    seeds a new locus. The locus is assigned to the D. melanogaster gene of its seed (best) alignment
    (D. melanogaster query: itself; other query: its DIAMOND best hit). Intact = some alignment of a query
    assigned to the same D. melanogaster gene covers >= 80% of the query without frameshift or in-frame
    stop codon; otherwise disrupted.
 3. Copy number per D. melanogaster gene and genome = sum over its loci of a depth weight w: nd of the
    annotated gene whose CDS the locus overlaps (same D. melanogaster assignment preferred), w = 1 if the
    locus overlaps no annotated CDS or in Dbia / DsuzR (no reads). An assembly duplicate contributes 0.5 +
    0.5, a duplicate collapsed into one assembled copy about 2. Totals (all loci) and intact copies.
 4. Differences Dsub vs Dsuz: rounded copy numbers differ, |difference| >= 0.5, the two D. suzukii
    assemblies agree (rounded IsoJpt1.0 value = RU_1.0 locus count), no fractional part in 0.3-0.7
    (otherwise 'intermediate': possible polymorphic CNV or mapping noise).
 5. Lineage: Wagner parsimony on the quartet ((Dsub, Dsuz), Dbia, Dmel) with integer copy numbers (Dbia:
    loci; Dmel: annotated protein-coding genes); branch changes reported only where all most-parsimonious
    reconstructions agree.
Outputs: results/C2_{validation.log, calls.tsv, categories.tsv, positive_controls.tsv,
         annotation_candidates_outcome.tsv}; Data_scratch/.../gene_family_cn/C2_{loci,per_gene}.tsv.gz
Run with envs/popgen/bin/python. No random numbers are used.
"""
import gzip
import importlib.util
import itertools
import os
import re
from collections import Counter, defaultdict
from urllib.parse import unquote

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
OUT = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn'
G = '/ru-auth/local/home/asato/scratch/genome'
SPECIES = ['Dmel', 'Dbia', 'Dsub', 'DsuzI', 'DsuzR']
GENOMES = ['Dbia', 'Dsub', 'DsuzI', 'DsuzR']
DSUZ_POOLS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
OWN = {'Dbia': 'Dbia', 'Dsub': 'Dsub', 'DsuzI': 'DsuzI', 'DsuzR': 'DsuzI'}  # query set not used on this genome
GFF = {'Dbia': f'{G}/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.gff',
       'DsuzR': f'{G}/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.gff'}
SCOV_FULL = 0.5
MIN_COV = 0.5
INTACT_COV = 0.8
MIN_DIFF = 0.5
MAX_LEN_RATIO = 1.5
HISTONE_ARRAY = re.compile(r'^His(1|2A|2B|3|4):CG\d+$')

spec = importlib.util.spec_from_file_location('c2', f'{PROJ}/scripts/135_gene_family_copy_number.py')
c2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c2)
log = open(f'{RES}/C2_validation.log', 'w')


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    log.write(s + '\n')


say('C2 second pass: copy number from miniprot loci weighted by read depth')
# ------------------------------------------------------------------ scripts/135 tables
genes = {sp: pd.read_csv(f'{OUT}/{sp}_genes.tsv.gz', sep='\t', low_memory=False) for sp in SPECIES}
asg = pd.read_csv(f'{OUT}/C2_gene_assignments.tsv.gz', sep='\t', low_memory=False)
cols = ['q', 's', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']
scov, lratio = {}, {}
for sp in SPECIES:
    h = pd.read_csv(f'{OUT}/{sp}_vs_Dmel.tsv', sep='\t', names=cols)
    h = h[h.evalue <= 1e-5].sort_values(['q', 'bitscore'], ascending=[True, False]).drop_duplicates('q')
    scov.update(dict(zip(h.q, h.length / h.slen)))
    lratio.update(dict(zip(h.q, h.qlen / h.slen)))
asg['key'] = asg.species + '|' + asg.gene
asg['scov'] = asg.key.map(scov)
asg['len_ratio'] = asg.key.map(lratio)
assert asg.scov.notna().all()
# full-length query: covers >= 50% of its D. melanogaster best hit and is not > 1.5x longer (fused gene models
# align across tandem paralogs and merge them into one locus; found at the alpha-Est cluster)
asg['full'] = (asg.scov >= SCOV_FULL) & (asg.len_ratio <= MAX_LEN_RATIO)
say(f'annotated proteins > {MAX_LEN_RATIO}x the length of their D. melanogaster best hit (possible fused models, '
    f'not used as miniprot queries): {((asg.len_ratio > MAX_LEN_RATIO) & (asg.species != "Dmel")).groupby(asg.species).sum().to_dict()}')
lowq = {}
for sp in SPECIES:
    t = genes[sp]
    lowq.update(dict(zip(sp + '|' + t.gene, t.name.fillna('').str.startswith('LOW QUALITY PROTEIN'))))
asg['low_quality'] = asg.key.map(lowq).fillna(False).astype(bool)
dmel_of = dict(zip(asg.key, asg.dmel))
full_key = set(asg.key[asg.full])
dm = genes['Dmel'][genes['Dmel'].coding]
sym = dict(zip(dm.gene, dm.symbol))
fbgn = dict(zip(dm.gene, dm.fbgn))
aCN = asg.pivot_table(index='dmel', columns='species', values='gene', aggfunc='count', fill_value=0)
aCN = aCN.reindex(index=dm.gene, columns=SPECIES, fill_value=0)
first = pd.read_csv(f'{RES}/C2_dmel_gene_copy_number.tsv.gz', sep='\t', index_col=0)
assert (first[SPECIES] == aCN[SPECIES]).all().all(), 'annotation counts do not reproduce scripts/135'
say(f'annotation counts reproduce the first-pass table ({len(aCN)} D. melanogaster genes)')

# ------------------------------------------------------------------ 1. read depth
cls = {}
for sp in ['Dsub', 'DsuzI']:
    cc = pd.read_csv(f'{OUT}/{sp}_chrom_class.tsv', sep='\t', index_col=0).chromosome
    cls[sp] = cc.map(lambda c: 'X' if c == 'X' else ('Y' if c == 'Y' else 'A'))
single = aCN.index[(aCN[SPECIES] == 1).all(axis=1)]
SAMPLES = {'Dsub_L5_7': 'Dsub', 'refstrain': 'Dsub', **{p: 'DsuzI' for p in DSUZ_POOLS}}
nd = {}
say('\n1. Read depth over CDS / median of strict single-copy genes (annotation count 1 in all five annotations, '
    'placed, no near-identical paralog, >= 80% of the D. melanogaster protein)')
for smp, sp in SAMPLES.items():
    b = pd.read_csv(f'{OUT}/depth/{smp}.Q0.bedcov', sep='\t', header=None,
                    names=['seqid', 'start', 'end', 'gene', 'sum', 'n'])
    b['len'] = b.end - b.start
    g = b.groupby('gene').agg(seqid=('seqid', 'first'), sum=('sum', 'sum'), len=('len', 'sum'))
    g['mean'] = g['sum'] / g['len']
    g['cls'] = g.seqid.map(cls[sp])
    a = asg[asg.species == sp].set_index('gene')
    ref = a.index[a.dmel.isin(single) & a.placed & ~a.near_identical_paralog & (a.scov >= 0.8)]
    norm = {}
    for c in ('A', 'X'):
        sel = g.loc[g.index.isin(ref) & (g.cls == c), 'mean']
        norm[c] = sel.median()
        r = sel / norm[c]
        say(f'  {smp:10s} {c}: n {len(sel):5d}; median depth {norm[c]:6.1f}; single-copy genes 2.5-97.5% '
            f'{r.quantile(0.025):.2f}-{r.quantile(0.975):.2f}, 10-90% {r.quantile(0.1):.2f}-{r.quantile(0.9):.2f}')
    nd[smp] = (g['mean'] / g.cls.map(norm)).to_dict()  # Y-linked genes -> NaN
nd_suz = pd.DataFrame({p: pd.Series(nd[p]) for p in DSUZ_POOLS})
nd['Dsuz_median'] = nd_suz.median(axis=1).to_dict()

# ------------------------------------------------------------------ CDS intervals of annotated genes
def attrs(s):
    return {k: unquote(v) for k, v in (x.split('=', 1) for x in s.rstrip(';').split(';') if '=' in x)}


def cds_bed_other(sp):
    """CDS of the chosen protein (Dbia) or chosen transcript (DsuzR) per gene, cached as BED."""
    path = f'{OUT}/{sp}_cds.bed'
    if os.path.exists(path):
        return
    t = genes[sp][genes[sp].coding]
    want = dict(zip(t.protein, t.gene))
    iv = defaultdict(list)
    with open(GFF[sp]) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'CDS':
                continue
            a = attrs(p[8])
            unit = a.get('protein_id') if sp == 'Dbia' else a.get('Parent')
            if unit in want:
                iv[want[unit]].append((p[0], int(p[3]) - 1, int(p[4])))
    assert len(iv) == len(want), (sp, len(iv), len(want))
    rows = sorted((s, a, b, g) for g, x in iv.items() for s, a, b in x)
    pd.DataFrame(rows).to_csv(path, sep='\t', header=False, index=False)


for sp in GFF:
    cds_bed_other(sp)
CDS = {}
for gnm in GENOMES:
    b = pd.read_csv(f'{OUT}/{gnm}_cds.bed', sep='\t', header=None, names=['seqid', 'start', 'end', 'gene'])
    CDS[gnm] = {s: x.sort_values('start').reset_index(drop=True) for s, x in b.groupby('seqid')}


def overlapping_genes(gnm, seqid, ts, te):
    x = CDS[gnm].get(seqid)
    if x is None:
        return []
    st = x.start.to_numpy()
    lo = np.searchsorted(st, ts - 50000)
    hi = np.searchsorted(st, te)
    y = x.iloc[lo:hi]
    return sorted(set(y.gene[y.end.to_numpy() > ts]))


# ------------------------------------------------------------------ 2. miniprot loci
def read_paf(path):
    """PAF records of a miniprot --gff file, plus the CDS exons (0-based half-open) of each alignment."""
    rows, exons = [], []
    with open(path) as f:
        for line in f:
            if line.startswith('##PAF'):
                p = line.rstrip('\n').split('\t')[1:]
                tags = dict((t[:2], t[5:]) for t in p[12:] if t[2:5] == ':i:')
                rows.append((p[0], int(p[1]), int(p[2]), int(p[3]), p[4], p[5], int(p[7]), int(p[8]),
                             int(p[9]), int(p[10]), int(tags.get('AS', 0)), int(tags.get('fs', 0)),
                             int(tags.get('st', 0))))
                exons.append([])
            elif not line.startswith('#'):
                p = line.split('\t', 5)
                if len(p) > 4 and p[2] == 'CDS':
                    exons[-1].append((int(p[3]) - 1, int(p[4])))
    d = pd.DataFrame(rows, columns=['q', 'qlen', 'qs', 'qe', 'strand', 'tname', 'ts', 'te', 'nmatch', 'blen',
                                    'AS', 'fs', 'st'])
    d['aid'] = np.arange(len(d))
    d['qsp'] = d.q.str.split('|').str[0]
    d['cov'] = (d.qe - d.qs) / d.qlen
    d['ident'] = d.nmatch / d.blen
    return d, exons


BIN = 100000
say(f'\n2. miniprot loci (other-species queries: D. melanogaster + full-length proteins; query coverage >= {MIN_COV}; '
    f'greedy clustering with >= 50% overlap of the seed alignment)')
loci = []
for gnm in GENOMES:
    d, exons = read_paf(f'{OUT}/mp/{gnm}.gff')
    n_all = len(d)
    keep = (d.qsp != OWN[gnm]) & (d['cov'] >= MIN_COV) & ((d.qsp == 'Dmel') | d.q.isin(full_key))
    d = d[keep].copy()
    d['dmel'] = np.where(d.qsp == 'Dmel', d.q.str.split('|').str[1], d.q.map(dmel_of))
    d['clean'] = (d['cov'] >= INTACT_COV) & (d.fs == 0) & (d.st == 0)
    # chimeric alignments: the span contains >= 2 non-overlapping other alignments that each cover >= 80% of
    # their query (the alignment chains exons of tandem paralogs); removed before clustering
    chim = np.zeros(len(d), dtype=bool)
    d = d.reset_index(drop=True)
    for (t, s), x in d.groupby(['tname', 'strand']):
        f = x[x['cov'] >= INTACT_COV].sort_values('ts')
        fts, fte, fidx = f.ts.to_numpy(), f.te.to_numpy(), f.index.to_numpy()
        for i, a, b in zip(x.index, x.ts, x.te):
            lo, hi = np.searchsorted(fts, a), np.searchsorted(fts, b)
            if hi - lo < 2:
                continue
            sel = (fte[lo:hi] <= b) & (fidx[lo:hi] != i)
            if sel.sum() < 2:
                continue
            order = np.argsort(fte[lo:hi][sel])
            st_, en_ = fts[lo:hi][sel][order], fte[lo:hi][sel][order]
            n, last = 0, -1
            for u, v in zip(st_, en_):
                if u >= last:
                    n += 1
                    last = v
                    if n >= 2:
                        chim[i] = True
                        break
    n_chim = int(chim.sum())
    d = d[~chim]
    d = d.sort_values('AS', ascending=False).reset_index(drop=True)
    seeds = []  # (tname, strand, ts, te)
    bins = defaultdict(list)
    lid = np.empty(len(d), dtype=int)
    for i, (t, s, a, b) in enumerate(zip(d.tname, d.strand, d.ts, d.te)):
        found = -1
        for k in range(a // BIN, b // BIN + 1):
            for j in bins[(t, s, k)]:
                _, _, sa, sb = seeds[j]
                ov = min(b, sb) - max(a, sa)
                if ov > 0 and ov >= 0.5 * min(b - a, sb - sa):
                    found = j
                    break
            if found >= 0:
                break
        if found < 0:
            found = len(seeds)
            seeds.append((t, s, a, b))
            for k in range(a // BIN, b // BIN + 1):
                bins[(t, s, k)].append(found)
        lid[i] = found
    d['locus'] = lid
    seed = d.drop_duplicates('locus').set_index('locus')  # first row = highest score
    d['seed_dmel'] = d.locus.map(seed.dmel)
    intact = d[d.dmel == d.seed_dmel].groupby('locus').clean.any()
    L = pd.DataFrame({'genome': gnm, 'tname': seed.tname, 'strand': seed.strand, 'ts': seed.ts, 'te': seed.te,
                      'dmel': seed.dmel, 'seed_q': seed.q, 'seed_ident': seed.ident, 'seed_cov': seed['cov'],
                      'seed_fs': seed.fs, 'seed_st': seed.st, 'n_aln': d.groupby('locus').size()})
    L['intact'] = intact.reindex(L.index).fillna(False).astype(bool)
    ov = [overlapping_genes(gnm, t, a, b) for t, a, b in zip(L.tname, L.ts, L.te)]
    own = 'DsuzR' if gnm == 'DsuzR' else OWN[gnm]
    L['annotated_cds'] = [','.join(x) for x in ov]
    L['annotated_same_dmel'] = [','.join(g for g in x if dmel_of.get(f'{own}|{g}') == dm_) for x, dm_ in zip(ov, L.dmel)]
    L['unannotated'] = [len(x) == 0 for x in ov]
    L['name'] = [f'{gnm}|{l}' for l in L.index]
    if gnm in ('Dsub', 'DsuzI'):  # exons of the seed alignment, for read depth of every locus
        rows_ = sorted((t, a, b, nm) for nm, t, aid in zip(L.name, L.tname, seed.aid) for a, b in exons[aid])
        txt = ''.join(f'{t}\t{a}\t{b}\t{nm}\n' for t, a, b, nm in rows_)  # sorted by position: faster bedcov
        bp = f'{OUT}/mp/{gnm}_loci_exons.bed'
        if not (os.path.exists(bp) and open(bp).read() == txt):  # rewrite only if changed (a running bedcov
            with open(bp + '.tmp', 'w') as fo:                   # may be reading it)
                fo.write(txt)
            os.replace(bp + '.tmp', bp)
    loci.append(L.reset_index().rename(columns={'index': 'locus'}))
    say(f'  {gnm:5s}: alignments {n_all}; used {keep.sum()}; chimeric removed {n_chim}; loci {len(L)}; intact {L.intact.sum()}; disrupted '
        f'{(~L.intact).sum()}; loci overlapping no annotated CDS {L.unannotated.sum()} '
        f'(intact {(L.unannotated & L.intact).sum()})')
LO = pd.concat(loci, ignore_index=True)
need = [f'{OUT}/depth/loci_{s}.Q0.bedcov' for s in SAMPLES]
stale = False
for s, sp in SAMPLES.items():  # depth files must match the current locus BED exactly
    p = f'{OUT}/depth/loci_{s}.Q0.bedcov'
    if os.path.exists(p):
        a_ = pd.read_csv(p, sep='\t', header=None, usecols=[0, 1, 2, 3])
        b_ = pd.read_csv(f'{OUT}/mp/{sp}_loci_exons.bed', sep='\t', header=None)
        if a_.shape != b_.shape or not (a_.values == b_.values).all():
            say(f'stale locus depth file for {s} (does not match the current BED)')
            stale = True
LO.to_csv(f'{OUT}/C2_loci_predepth.tsv.gz', sep='\t', index=False)
if stale or not all(os.path.exists(p) for p in need):
    say('locus exon BED files written (mp/{Dsub,DsuzI}_loci_exons.bed); run bedcov on them '
        '(scripts/137_loci_depth.sbatch), then rerun this script')
    raise SystemExit(0)

# depth of every locus from its own exons; normaliser = loci of strict single-copy genes (one locus per genome)
say('\n   locus read depth (own exons) / median over loci of strict single-copy genes (one locus in the genome, '
    'intact, placed chromosome)')
LO['w'] = 1.0
LO['w_ref'] = np.nan
for p in DSUZ_POOLS:
    LO[f'w_{p}'] = np.nan
SD = {}
for smp, sp in SAMPLES.items():
    b = pd.read_csv(f'{OUT}/depth/loci_{smp}.Q0.bedcov', sep='\t', header=None,
                    names=['seqid', 'start', 'end', 'name', 'sum', 'n'])
    b['len'] = b.end - b.start
    g = b.groupby('name').agg(sum=('sum', 'sum'), len=('len', 'sum'))
    m = LO.genome == sp
    x = LO[m].set_index('name')
    mean = (g['sum'] / g['len']).reindex(x.index)
    c = x.tname.map(cls[sp])
    n_loc = x.groupby('dmel').dmel.transform('size')
    ref = x.dmel.isin(single) & (n_loc == 1) & x.intact & x.tname.str.startswith('NC_')
    v = pd.Series(np.nan, index=x.index)
    for k in ('A', 'X'):
        med = mean[ref & (c == k)].median()
        v[c == k] = mean[c == k] / med
        r = v[ref & (c == k)]
        SD[(smp, k)] = (r.quantile(0.9) - r.quantile(0.1)) / 2.563
        say(f'   {smp:10s} {k}: reference loci {int((ref & (c == k)).sum()):5d}; median depth {med:6.1f}; '
            f'2.5-97.5% {r.quantile(0.025):.2f}-{r.quantile(0.975):.2f}; robust SD per copy {SD[(smp, k)]:.3f}')
    col = {'Dsub_L5_7': 'w', 'refstrain': 'w_ref'}.get(smp, f'w_{smp}')
    LO.loc[m, col] = v.values
dz = LO.genome == 'DsuzI'
LO.loc[dz, 'w'] = LO.loc[dz, [f'w_{p}' for p in DSUZ_POOLS]].median(axis=1)
LO['has_depth'] = LO.genome.isin(['Dsub', 'DsuzI']) & LO.w.notna()
LO['w'] = LO.w.fillna(1.0)  # Y-linked IsoJpt1.0 loci and genomes without reads
# per-copy noise of the Dsuz median: from the reference loci
xz = LO[dz].copy()
n_loc = xz.groupby('dmel').dmel.transform('size')
r = xz.w[xz.dmel.isin(single) & (n_loc == 1) & xz.intact & xz.tname.str.startswith('NC_')]
SD_SUB = max(SD[('Dsub_L5_7', 'A')], SD[('Dsub_L5_7', 'X')])
SD_SUZ = (r.quantile(0.9) - r.quantile(0.1)) / 2.563
say(f'   per-copy SD used in the difference test: Dsub {SD_SUB:.3f}; Dsuz (median of four pools) {SD_SUZ:.3f}')
LO.to_csv(f'{OUT}/C2_loci.tsv.gz', sep='\t', index=False)

# ------------------------------------------------------------------ 3. copy number per D. melanogaster gene
E = pd.DataFrame(index=dm.gene)
E['symbol'] = E.index.map(sym)
for sp in SPECIES:
    E[f'ann_{sp}'] = aCN[sp]
for gnm in GENOMES:
    x = LO[LO.genome == gnm]
    gb = x.groupby('dmel')
    E[f'L_{gnm}'] = gb.size().reindex(E.index).fillna(0).astype(int)
    E[f'Lint_{gnm}'] = x[x.intact].groupby('dmel').size().reindex(E.index).fillna(0).astype(int)
    E[f'Lun_{gnm}'] = x[x.unannotated].groupby('dmel').size().reindex(E.index).fillna(0).astype(int)
    E[f'W_{gnm}'] = gb.w.sum().reindex(E.index).fillna(0)
    E[f'Wint_{gnm}'] = x[x.intact].groupby('dmel').w.sum().reindex(E.index).fillna(0)
E['W_Dsub_ref'] = LO[LO.genome == 'Dsub'].groupby('dmel').w_ref.sum(min_count=1).reindex(E.index)
for p in DSUZ_POOLS:
    E[f'W_{p}'] = LO[LO.genome == 'DsuzI'].groupby('dmel')[f'w_{p}'].sum(min_count=1).reindex(E.index)
E['nodepth_Dsub'] = LO[(LO.genome == 'Dsub') & ~LO.has_depth].groupby('dmel').size().reindex(E.index).fillna(0)
E['nodepth_DsuzI'] = LO[(LO.genome == 'DsuzI') & ~LO.has_depth].groupby('dmel').size().reindex(E.index).fillna(0)
for sp, tag in (('Dsub', 'Dsub'), ('DsuzI', 'DsuzI')):
    a = asg[asg.species == sp]
    E[f'unplaced_{tag}'] = a[~a.placed].groupby('dmel').size().reindex(E.index).fillna(0).astype(int)
    E[f'lowq_{tag}'] = a[a.low_quality].groupby('dmel').size().reindex(E.index).fillna(0).astype(int)


def call(r, kind):
    """Dsub vs Dsuz for all copies ('tot') or intact copies ('int').

    W = depth-weighted copy number, L = locus count. Supported: the two D. suzukii assemblies agree (equal
    locus counts, or the depth-weighted IsoJpt1.0 value rounds to the RU_1.0 count); |W_Dsub - W_Dsuz| >=
    max(0.5, 2 SE) with SE from the per-copy SD of single-copy loci; locus and depth differences in the same
    direction; rounded values differ ('intermediate' otherwise). Equal locus counts with a depth difference =
    'depth_only' (copies collapsed in one assembly, or a population CNV)."""
    p = '' if kind == 'tot' else 'int'
    Ws, Wz = r[f'W{p}_Dsub'], r[f'W{p}_DsuzI']
    Ls, Lz, LzR = r[f'L{p}_Dsub'], r[f'L{p}_DsuzI'], r[f'L{p}_DsuzR']
    rs, rz = int(round(Ws)), int(round(Wz))
    R, D = Ls - Lz, Ws - Wz
    if R == 0 and rs == rz:
        return ''
    if not (Lz == LzR or rz == LzR):
        return 'dsuz_assemblies_disagree'
    se = np.sqrt(max(Ls, 1) * SD_SUB ** 2 + max(Lz, 1) * SD_SUZ ** 2)
    if abs(D) < max(MIN_DIFF, 2 * se):
        if R > 0:
            return 'assembly_duplicate_Dsub'
        if R < 0:
            return 'assembly_duplicate_DsuzI'
        return 'small'
    if R != 0 and np.sign(R) != np.sign(D):
        return 'depth_opposite'
    if rs == rz:
        return 'intermediate'
    return 'depth_only' if R == 0 else 'supported'


E['call_tot'] = E.apply(lambda r: call(r, 'tot'), axis=1)
E['call_int'] = E.apply(lambda r: call(r, 'int'), axis=1)
hist = E.symbol.astype(str).str.match(HISTONE_ARRAY)
E.loc[hist & (E.call_tot != ''), 'call_tot'] = 'histone_array_excluded'
E.loc[hist & (E.call_int != ''), 'call_int'] = 'histone_array_excluded'


def parsimony(sub, suz, bia, mel):
    """Wagner parsimony on ((sub, suz)X, bia, mel)Y; returns branch changes (None = ambiguous)."""
    vals = [sub, suz, bia, mel]
    rng = range(min(vals), max(vals) + 1)
    best, sols = None, []
    for X, Y in itertools.product(rng, rng):
        c = abs(sub - X) + abs(suz - X) + abs(X - Y) + abs(bia - Y) + abs(mel - Y)
        if best is None or c < best:
            best, sols = c, [(X, Y)]
        elif c == best:
            sols.append((X, Y))
    out = {}
    for name, f in (('Dsub', lambda X, Y: sub - X), ('Dsuz', lambda X, Y: suz - X),
                    ('shared', lambda X, Y: X - Y), ('Dbia', lambda X, Y: bia - Y)):
        v = {f(X, Y) for X, Y in sols}
        out[name] = float(v.pop()) if len(v) == 1 else np.nan
    out['anc_ss'] = sols[0][0] if len({X for X, _ in sols}) == 1 else None
    return out


for kind, wsub, wsuz, lb in (('tot', 'W_Dsub', 'W_DsuzI', 'L_Dbia'), ('int', 'Wint_Dsub', 'Wint_DsuzI', 'Lint_Dbia')):
    res = [parsimony(int(round(a)), int(round(b)), int(c), int(m))
           for a, b, c, m in zip(E[wsub], E[wsuz], E[lb], E.ann_Dmel)]
    for br in ('Dsub', 'Dsuz', 'shared', 'Dbia'):
        E[f'{kind}_{br}'] = [x[br] for x in res]
E.to_csv(f'{OUT}/C2_per_gene.tsv.gz', sep='\t')

# ------------------------------------------------------------------ genome-wide summary
present = E[E[[f'L_{g}' for g in GENOMES]].sum(axis=1) > 0].copy()
say(f'\n3. D. melanogaster genes with >= 1 locus in Dbia/Dsub/Dsuz: {len(present)}')
for kind, lab in (('tot', 'all copies'), ('int', 'intact copies')):
    c = Counter(present[f'call_{kind}'])
    c.pop('', None)
    say(f'  Dsub vs Dsuz, {lab}: {dict(c)}')
sup = present[(present.call_tot == 'supported') | (present.call_int == 'supported')].copy()
sup['type'] = np.where(sup.call_tot == 'supported',
                       np.where(sup.W_Dsub > sup.W_DsuzI, 'more_copies_Dsub', 'more_copies_Dsuz'),
                       np.where(sup.Wint_Dsub > sup.Wint_DsuzI, 'disrupted_in_Dsuz', 'disrupted_in_Dsub'))
say(f'  supported differences (all or intact copies): {len(sup)}; {dict(Counter(sup.type))}')
say(f'  of these, the extra copies or the intact copies lack read depth (locus without depth) in '
    f'{int(((sup.nodepth_Dsub > 0) | (sup.nodepth_DsuzI > 0)).sum())}')


def branch_events(df, kind):
    ev = Counter()
    for br in ('Dsub', 'Dsuz', 'shared'):
        v = df[f'{kind}_{br}']
        ev[f'{br}_gain'] = int((v > 0).sum())
        ev[f'{br}_loss'] = int((v < 0).sum())
        ev[f'{br}_ambiguous'] = int(v.isna().sum())
    return ev


st = sup[sup.call_tot == 'supported']
si = sup[(sup.call_tot != 'supported') & (sup.call_int == 'supported')]
say(f'  lineage (parsimony), copy-number differences: {dict(branch_events(st, "tot"))}')
say(f'  lineage (parsimony), intact-copy differences without copy-number difference (pseudogenes): '
    f'{dict(branch_events(si, "int"))}')
# shared branch: Dsub and Dsuz agree (rounded depth-weighted values equal, RU_1.0 agrees), differ from outgroups
def shared_ok(df, p):
    rs, rz = df[f'W{p}_Dsub'].round(), df[f'W{p}_DsuzI'].round()
    return (rs == rz) & ((df[f'L{p}_DsuzI'] == df[f'L{p}_DsuzR']) | (rz == df[f'L{p}_DsuzR']))


hist_p = present.symbol.astype(str).str.match(HISTONE_ARRAY)
shd = present[shared_ok(present, '') & present.tot_shared.notna() & (present.tot_shared != 0) & ~hist_p]
shi = present[shared_ok(present, 'int') & present.int_shared.notna() & (present.int_shared != 0) & ~hist_p
              & ~present.index.isin(shd.index)]
say(f'  changes on the branch shared by Dsub and Dsuz (Dsub = Dsuz, both D. suzukii assemblies agree): '
    f'copies gained {(shd.tot_shared > 0).sum()}, lost {(shd.tot_shared < 0).sum()}; intact copies only (copy '
    f'present, disrupted in both species): gained {(shi.int_shared > 0).sum()}, lost {(shi.int_shared < 0).sum()}')
for kind in ('tot', 'int'):
    c = Counter(present[f'call_{kind}'])
    say(f'  depth_only ({kind}; equal locus counts, depth differs: collapsed copies or population CNV): '
        f'{c.get("depth_only", 0)}')

# ------------------------------------------------------------------ annotation candidates of the first pass
lc = pd.read_csv(f'{OUT}/C2_lineage_calls.tsv.gz', sep='\t', index_col=0)
ac = lc[lc.sub_ne_suz].index
A = E.loc[ac].copy()
A['ann_diff'] = A.ann_Dsub - A.ann_DsuzI


def outcome(r):
    if HISTONE_ARRAY.match(str(r.symbol)):
        return 'histone_array'
    if r.call_tot in ('supported', 'depth_only') and np.sign(r.W_Dsub - r.W_DsuzI) == np.sign(r.ann_diff):
        return 'confirmed'
    if r.call_int == 'supported' and np.sign(r.Wint_Dsub - r.Wint_DsuzI) == np.sign(r.ann_diff):
        return 'confirmed_intact_only'
    if r.call_tot.startswith('assembly_duplicate'):
        return r.call_tot
    if r.L_Dsub == r.L_DsuzI and r.call_tot in ('', 'small'):
        return 'annotation_only (same number of loci)'
    return r.call_tot or 'other'


A['outcome'] = A.apply(outcome, axis=1)
A[['symbol', 'ann_Dbia', 'ann_Dsub', 'ann_DsuzI', 'ann_DsuzR', 'L_Dsub', 'L_DsuzI', 'L_DsuzR', 'W_Dsub',
   'W_DsuzI', 'outcome']].rename_axis('dmel_gene').to_csv(
    f'{RES}/C2_annotation_candidates_outcome.tsv', sep='\t', float_format='%.2f')
say(f'\n4. The {len(A)} annotation-based differences of the first pass (scripts/135): {dict(Counter(A.outcome))}')
dsub_gain = A[A.ann_diff > 0]
say(f'   of these, {len(dsub_gain)} with more annotated copies in Dsub: {dict(Counter(dsub_gain.outcome))}')

# ------------------------------------------------------------------ categories
go_cache = {}


def go_set(root):
    if root not in go_cache:
        terms = c2.go_descendants(root)
        s = set()
        with gzip.open(f'{c2.GO}/fb.gaf.gz', 'rt') as f:
            for line in f:
                if line.startswith('!'):
                    continue
                t = line.split('\t')
                if 'NOT' not in t[3] and t[4] in terms:
                    s.add(t[1])
        go_cache[root] = s
    return go_cache[root]


SUGAR_GR = {'Gr5a', 'Gr43a', 'Gr61a', 'Gr64a', 'Gr64b', 'Gr64c', 'Gr64d', 'Gr64e', 'Gr64f'}
CATS = {
    'Or': lambda s, f: bool(re.match(r'^Or\d', s)),
    'Gr_sugar': lambda s, f: s in SUGAR_GR,
    'Gr_other': lambda s, f: bool(re.match(r'^Gr\d', s)) and s not in SUGAR_GR,
    'Ir': lambda s, f: bool(re.match(r'^Ir\d', s)),
    'Obp': lambda s, f: s.startswith('Obp'),
    'ppk': lambda s, f: s.startswith('ppk'),
    'mechanosensory_channel': lambda s, f: s in c2.MECHANO or s == 'nompA',
    'Cyp': lambda s, f: s.startswith('Cyp'),
    'Gst': lambda s, f: s.startswith('Gst'),
    'Ugt': lambda s, f: s.startswith('Ugt'),
    'carboxylesterase': lambda s, f: bool(re.match(r'^(alpha-Est|Est-|Jhe)', s)),
    'cuticle_structural': lambda s, f: f in go_set('GO:0042302') or bool(re.match(r'^(Cpr|Twdl|Lcp|Edg)', s)),
    'chitin_metabolism': lambda s, f: f in go_set('GO:0006030'),
    'pigmentation_sclerotization': lambda s, f: f in go_set('GO:0048067') or
        s in {'y', 'e', 't', 'ple', 'Ddc', 'Lac', 'Pu', 'Dat', 'black'} or s.startswith('yellow'),
    'genitalia_development': lambda s, f: f in go_set('GO:0030540') or f in go_set('GO:0035215'),
    'eggshell': lambda s, f: f in go_set('GO:0007304'),
    'immune_GO0006955': lambda s, f: f in go_set('GO:0006955'),
}
for c, fn in CATS.items():
    E[f'cat_{c}'] = [fn(str(s), fb) for s, fb in zip(E.symbol, E.index.map(fbgn))]
E['categories'] = E[[f'cat_{c}' for c in CATS]].apply(lambda x: ','.join(c for c, v in zip(CATS, x) if v), axis=1)
present = E.loc[present.index]
issup = (present.call_tot == 'supported') | (present.call_int == 'supported')
isshared = present.index.isin(shd.index) | present.index.isin(shi.index)
present['shared_change'] = np.where(present.index.isin(shd.index), present.tot_shared,
                                    np.where(present.index.isin(shi.index), present.int_shared, 0))
say('\n5. Pre-specified categories: genes present; supported Dsub-vs-Dsuz differences (all or intact copies); '
    'shared-branch changes; Fisher test of the supported-difference rate vs all other genes; BH FDR over categories')
rows = []
for c in CATS:
    mk = present[f'cat_{c}'].to_numpy()
    s = issup.to_numpy()
    a_, b_, c_, d_ = int((mk & s).sum()), int((mk & ~s).sum()), int((~mk & s).sum()), int((~mk & ~s).sum())
    orr, p = fisher_exact([[a_, b_], [c_, d_]])
    x = present[mk & s]
    ev = Counter()
    for _, r in x.iterrows():
        kind = 'tot' if r.call_tot == 'supported' else 'int'
        for br in ('Dsub', 'Dsuz'):
            v = r[f'{kind}_{br}']
            if v is not None and pd.notna(v) and v != 0:
                ev[f'{br}_{"gain" if v > 0 else "loss"}{"" if kind == "tot" else "_intact"}'] += 1
        if all(pd.isna(r[f'{kind}_{br}']) for br in ('Dsub', 'Dsuz')):
            ev['ambiguous'] += 1
    xs = present[mk & isshared]
    rows.append(dict(category=c, n_genes=int(mk.sum()), supported=a_, OR=orr, P=p,
                     events=';'.join(f'{k} {v}' for k, v in sorted(ev.items())),
                     shared_gain=int((xs.shared_change > 0).sum()), shared_loss=int((xs.shared_change < 0).sum()),
                     genes=','.join(sorted(x.symbol.astype(str))),
                     shared_genes=','.join(f'{s}{"+" if v > 0 else "-"}' for s, v in
                                           sorted(zip(xs.symbol.astype(str), xs.shared_change)))))
CR = pd.DataFrame(rows)
CR['BH_q'] = multipletests(CR.P, method='fdr_bh')[1]
CR.to_csv(f'{RES}/C2_categories.tsv', sep='\t', index=False)
for r in CR.itertuples():
    say(f'  {r.category:28s} n {r.n_genes:4d} supported {r.supported:3d} OR {r.OR:5.2f} P {r.P:.3g} BH {r.BH_q:.3g}; '
        f'{r.events}; shared +{r.shared_gain}/-{r.shared_loss}')
    if r.genes:
        say(f'      Dsub vs Dsuz: {r.genes}')
    if r.shared_genes:
        say(f'      shared branch: {r.shared_genes}')

# ------------------------------------------------------------------ gene clusters (tandem paralogs)
def stem(s):
    s = str(s)
    for rx in (r'^((?:Or|Gr|Ir)\d+)[a-z]?$', r'^(Obp\d+)[a-z]$', r'^(Gst[A-Z])\d+$', r'^(Cyp\d+[a-z])\d+(?:-p)?$',
               r'^(Ugt\d+[A-Z])\d+$', r'^(alpha-Est)\d+$', r'^(Cpr\d+[A-Z])[a-z]?\d*$', r'^(lectin-\d+[A-Z])[a-z]?$',
               r'^(Jon\d+[A-Z])[a-z]+$', r'^(Twdl)[A-Za-z]+$', r'^(Lcp65A)[a-z]+\d*$', r'^(Idgf)\d$',
               r'^(Cht)\d+$', r'^(ppk)\d+$'):
        m = re.match(rx, s)
        if m:
            return m.group(1)
    return s


E['cluster'] = E.symbol.map(stem)
num = [c for c in E.columns if re.match(r'^(L|Lint|W|Wint)_(Dbia|Dsub|DsuzI|DsuzR)$', c)] + ['ann_Dmel']
K = E.groupby('cluster')[num].sum()
K['n_dmel_genes'] = E.groupby('cluster').size()
K['categories'] = E.groupby('cluster').categories.agg(lambda x: ','.join(sorted({c for v in x for c in v.split(',') if c})))
K = K[K.n_dmel_genes > 1]
K['call_tot'] = K.apply(lambda r: call(r, 'tot'), axis=1)
K['call_int'] = K.apply(lambda r: call(r, 'int'), axis=1)
for kind, wsub, wsuz, lb in (('tot', 'W_Dsub', 'W_DsuzI', 'L_Dbia'), ('int', 'Wint_Dsub', 'Wint_DsuzI', 'Lint_Dbia')):
    res = [parsimony(int(round(a)), int(round(b)), int(c), int(m))
           for a, b, c, m in zip(K[wsub], K[wsuz], K[lb], K.ann_Dmel)]
    for br in ('Dsub', 'Dsuz', 'shared', 'Dbia'):
        K[f'{kind}_{br}'] = [x[br] for x in res]
K.to_csv(f'{RES}/C2_clusters.tsv', sep='\t', float_format='%.2f')
say(f'\n6. Gene clusters (D. melanogaster genes with a common symbol stem, e.g. Gr36a-c, GstD1-10; {len(K)} clusters '
    f'of >= 2 genes). Clusters in the pre-specified categories with a Dsub-vs-Dsuz call or a shared-branch change:')
say('   cluster (n genes): Dmel / loci Dbia, Dsub, DsuzI, DsuzR / intact / depth-weighted Dsub, DsuzI (intact) | '
    'calls tot/int | parsimony tot Dsub, Dsuz, shared | int Dsub, Dsuz, shared')
foc = K[(K.categories != '') & (K.categories != 'immune_GO0006955')]
for cl, r in foc.iterrows():
    if (r.call_tot in ('', 'small') and r.call_int in ('', 'small') and (r.tot_shared in (0, 0.0))
            and (r.int_shared in (0, 0.0))):
        continue
    say(f'   {cl:10s} ({r.n_dmel_genes:2d}) {int(r.ann_Dmel):2d} / {int(r.L_Dbia)}, {int(r.L_Dsub)}, {int(r.L_DsuzI)}, '
        f'{int(r.L_DsuzR)} / {int(r.Lint_Dbia)}, {int(r.Lint_Dsub)}, {int(r.Lint_DsuzI)}, {int(r.Lint_DsuzR)} / '
        f'{r.W_Dsub:.1f}, {r.W_DsuzI:.1f} ({r.Wint_Dsub:.1f}, {r.Wint_DsuzI:.1f}) | {r.call_tot or "-"}/{r.call_int or "-"} | '
        f'{r.tot_Dsub}, {r.tot_Dsuz}, {r.tot_shared} | {r.int_Dsub}, {r.int_Dsuz}, {r.int_shared}')

# ------------------------------------------------------------------ paralog groups (genome-wide counts)
# A copy assigned to paralog X in one species and to paralog Y in the other creates two opposite per-gene
# differences. Genome-wide counts are therefore made on paralog groups: connected components of D. melanogaster
# genes linked by a DIAMOND hit with >= 40% amino-acid identity covering >= 60% of both proteins.
# The first definition (E <= 1e-30, >= 60% coverage, plus common symbol stems) chained
# diverged paralogs into large groups (48 odorant receptors in one group; Gr36 with Gr59) whose summed read depth is
# dominated by noise. Copies are reassigned only between close paralogs, so the linkage is now identity-based.
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
PG_MIN_IDENT, PG_MIN_COV = 40, 0.6
h = pd.read_csv(f'{OUT}/Dmel_vs_Dmel.tsv', sep='\t', names=cols)
h['q'] = h.q.str.split('|').str[1]
h['s'] = h.s.str.split('|').str[1]
h = h[(h.q != h.s) & (h.pident >= PG_MIN_IDENT) & (h.length >= PG_MIN_COV * h.qlen) & (h.length >= PG_MIN_COV * h.slen)]
ids = pd.Index(E.index)
ei, ej = np.array(ids.get_indexer(h.q)), np.array(ids.get_indexer(h.s))
ok = (ei >= 0) & (ej >= 0)
ncomp, lab = connected_components(coo_matrix((np.ones(ok.sum()), (ei[ok], ej[ok])), shape=(len(ids), len(ids))),
                                  directed=False)
E['pgroup'] = lab
P = E.groupby('pgroup')[num].sum()
P['n_dmel_genes'] = E.groupby('pgroup').size()
P['members'] = E.groupby('pgroup').symbol.agg(lambda x: ','.join(sorted(map(str, x))[:6]) + (',...' if len(x) > 6 else ''))
P['histone'] = E.groupby('pgroup').symbol.agg(lambda x: any(HISTONE_ARRAY.match(str(s)) for s in x))
P = P[(P[[f'L_{g}' for g in GENOMES]].sum(axis=1) > 0) & ~P.histone]
P['call_tot'] = P.apply(lambda r: call(r, 'tot'), axis=1)
P['call_int'] = P.apply(lambda r: call(r, 'int'), axis=1)
for kind, wsub, wsuz, lb in (('tot', 'W_Dsub', 'W_DsuzI', 'L_Dbia'), ('int', 'Wint_Dsub', 'Wint_DsuzI', 'Lint_Dbia')):
    res = [parsimony(int(round(a)), int(round(b)), int(c), int(m))
           for a, b, c, m in zip(P[wsub], P[wsuz], P[lb], P.ann_Dmel)]
    for br in ('Dsub', 'Dsuz', 'shared', 'Dbia'):
        P[f'{kind}_{br}'] = [x[br] for x in res]
P.to_csv(f'{RES}/C2_paralog_groups.tsv', sep='\t', float_format='%.2f')
say(f'\n7. Paralog groups (>= {PG_MIN_IDENT}% identity over >= {PG_MIN_COV:.0%} of both proteins): {len(P)} units with >= 1 locus (histone cluster excluded); multi-gene groups '
    f'{(P.n_dmel_genes > 1).sum()} (largest {P.n_dmel_genes.max()} genes)')
for kind, lab_ in (('tot', 'all copies'), ('int', 'intact copies')):
    c = Counter(P[f'call_{kind}'])
    c.pop('', None)
    say(f'   Dsub vs Dsuz, {lab_}: {dict(c)}')
# technical baseline: raw locus counts, two assemblies of one species vs the two species
say(f'   raw locus counts differ: DsuzI vs DsuzR (same species) {(P.L_DsuzI != P.L_DsuzR).sum()}; Dsub vs DsuzI '
    f'{(P.L_Dsub != P.L_DsuzI).sum()}; Dsub vs DsuzR {(P.L_Dsub != P.L_DsuzR).sum()}')
ps = P[P.call_tot == 'supported']
pi = P[(P.call_tot != 'supported') & (P.call_int == 'supported')]
say(f'   supported copy-number differences {len(ps)}: more copies in Dsub {(ps.W_Dsub > ps.W_DsuzI).sum()}, '
    f'in Dsuz {(ps.W_Dsub < ps.W_DsuzI).sum()}; lineage {dict(branch_events(ps, "tot"))}')
say(f'   supported intact-copy differences without copy-number difference {len(pi)}: fewer intact in Dsub '
    f'{(pi.Wint_Dsub < pi.Wint_DsuzI).sum()}, in Dsuz {(pi.Wint_Dsub > pi.Wint_DsuzI).sum()}; lineage '
    f'{dict(branch_events(pi, "int"))}')
pshd = P[shared_ok(P, '') & P.tot_shared.notna() & (P.tot_shared != 0)]
say(f'   shared Dsub + Dsuz branch: copy gains {(pshd.tot_shared > 0).sum()}, losses {(pshd.tot_shared < 0).sum()}')


def lineage_bia(df, p=''):
    """Polarisation with D. biarmipes only: the species that differs from Dbia carries the change."""
    rs, rz, rb = df[f'W{p}_Dsub'].round(), df[f'W{p}_DsuzI'].round(), df[f'L{p}_Dbia']
    lin = np.where(rb == rz, 'Dsub', np.where(rb == rs, 'Dsuz', 'unresolved'))
    dirn = np.where(lin == 'Dsub', np.where(rs > rb, 'gain', 'loss'),
                    np.where(lin == 'Dsuz', np.where(rz > rb, 'gain', 'loss'), ''))
    return pd.Series([f'{a} {b}'.strip() for a, b in zip(lin, dirn)], index=df.index)


P['lineage_bia_tot'] = lineage_bia(P)
P['lineage_bia_int'] = lineage_bia(P, 'int')
say(f'   polarised with Dbia only, copy-number differences: {dict(Counter(P.lineage_bia_tot[P.call_tot == "supported"]))}')
say(f'   polarised with Dbia only, intact-only differences: '
    f'{dict(Counter(P.lineage_bia_int[(P.call_tot != "supported") & (P.call_int == "supported")]))}')
# per-assembly rate of missing loci, from the two D. suzukii assemblies (IsoJpt1.0 depth-checked: its locus count
# equals its rounded depth-weighted count, so the IsoJpt1.0 loci are not assembly duplicates)
okI = P.W_DsuzI.round() == P.L_DsuzI
say(f'   loci present in one D. suzukii assembly only (groups): RU_1.0 has fewer {(okI & (P.L_DsuzR < P.L_DsuzI)).sum()}, '
    f'IsoJpt1.0 has fewer {(P.L_DsuzI < P.L_DsuzR).sum()} (of these, IsoJpt1.0 has no locus at all '
    f'{((P.L_DsuzI == 0) & (P.L_DsuzR > 0)).sum()}; RU_1.0 none {((P.L_DsuzR == 0) & (P.L_DsuzI > 0) & okI).sum()})')
P.to_csv(f'{RES}/C2_paralog_groups.tsv', sep='\t', float_format='%.2f')
# category test against a family-matched background: genes in multi-gene paralog groups outside the category
# (tandem families carry more copy-number change, real and technical, than single-copy genes)
E['multi_group'] = E.pgroup.map(E.groupby('pgroup').size()) > 1
pr = E.loc[present.index]
s_ = ((pr.call_tot == 'supported') | (pr.call_int == 'supported')).to_numpy()
say('\n   categories vs genes of other multi-gene paralog groups (family-matched background):')
say(f'   supported-difference rate: all genes {s_.mean():.4f}; multi-gene groups {s_[pr.multi_group.to_numpy()].mean():.4f}; '
    f'single-gene groups {s_[~pr.multi_group.to_numpy()].mean():.4f}')
rows = []
for c in CATS:
    mk = pr[f'cat_{c}'].to_numpy(bool)
    bg = pr.multi_group.to_numpy() & ~mk
    a_, b_, c_, d_ = int((mk & s_).sum()), int((mk & ~s_).sum()), int((bg & s_).sum()), int((bg & ~s_).sum())
    orr, p = fisher_exact([[a_, b_], [c_, d_]])
    rows.append(dict(category=c, supported=a_, n=a_ + b_, bg_supported=c_, bg_n=c_ + d_, OR_matched=orr, P_matched=p))
CM = pd.DataFrame(rows)
CM['BH_matched'] = multipletests(CM.P_matched, method='fdr_bh')[1]
CM.to_csv(f'{RES}/C2_categories_matched_background.tsv', sep='\t', index=False)
for r in CM.itertuples():
    say(f'   {r.category:28s} {r.supported}/{r.n} vs {r.bg_supported}/{r.bg_n}; OR {r.OR_matched:.2f} P {r.P_matched:.3g} '
        f'BH {r.BH_matched:.3g}')
E['pgroup_members'] = E.pgroup.map(P.members)
E['pgroup_call_tot'] = E.pgroup.map(P.call_tot)

# ------------------------------------------------------------------ outputs
keep = ['symbol', 'categories', 'call_tot', 'call_int', 'tot_Dsub', 'tot_Dsuz', 'tot_shared', 'tot_Dbia',
        'int_Dsub', 'int_Dsuz', 'int_shared', 'ann_Dmel', 'ann_Dbia', 'ann_Dsub', 'ann_DsuzI', 'ann_DsuzR',
        'L_Dbia', 'L_Dsub', 'L_DsuzI', 'L_DsuzR', 'Lint_Dbia', 'Lint_Dsub', 'Lint_DsuzI', 'Lint_DsuzR',
        'Lun_Dsub', 'Lun_DsuzI', 'W_Dsub', 'W_Dsub_ref', 'W_DsuzI', 'Wint_Dsub', 'Wint_DsuzI'] + \
       [f'W_{p}' for p in DSUZ_POOLS] + ['nodepth_Dsub', 'nodepth_DsuzI', 'unplaced_Dsub', 'unplaced_DsuzI',
                                         'lowq_Dsub', 'lowq_DsuzI']
keep.insert(4, 'shared_change')
present['cluster'] = E.cluster
for c in ('pgroup', 'pgroup_members', 'pgroup_call_tot'):
    present[c] = E[c]
keep[1:1] = ['cluster', 'pgroup', 'pgroup_members', 'pgroup_call_tot']
out = present[(present.call_tot != '') | (present.call_int != '') | isshared][keep]
out.rename_axis('dmel_gene').to_csv(f'{RES}/C2_calls.tsv', sep='\t', float_format='%.2f')
say(f'\nwrote {len(out)} genes to results/C2_calls.tsv')

PC = ['Or23a', 'Or67a', 'Or74a', 'Or85a', 'Or98b', 'Gr36a', 'Gr36b', 'Gr36c', 'Gr59a', 'Gr59b', 'Gr59c',
      'Gr59d', 'Gr93a', 'Gr98b', 'Gr98d']
pc = present[present.symbol.isin(PC)].sort_values('symbol')
pc[keep].rename_axis('dmel_gene').to_csv(f'{RES}/C2_positive_controls.tsv', sep='\t', float_format='%.2f')
say('\nPositive controls (Hickner et al. 2016: D. suzukii Or23a and Or67a expanded; Or74a, Or85a, Or98b '
    'pseudogenes; Gr36a-c, Gr59a/b, Gr59c/d, Gr93a, Gr98b, Gr98d expanded)')
say('  symbol | loci Dbia/Dsub/DsuzI/DsuzR | intact | depth-weighted Dsub, DsuzI (intact) | calls tot/int | '
    'parsimony Dsub,Dsuz,shared (tot)')
for r in pc.itertuples():
    say(f'  {r.symbol:6s} | {r.L_Dbia}/{r.L_Dsub}/{r.L_DsuzI}/{r.L_DsuzR} | {r.Lint_Dbia}/{r.Lint_Dsub}/'
        f'{r.Lint_DsuzI}/{r.Lint_DsuzR} | {r.W_Dsub:.2f}, {r.W_DsuzI:.2f} ({r.Wint_Dsub:.2f}, {r.Wint_DsuzI:.2f}) | '
        f'{r.call_tot or "-"}/{r.call_int or "-"} | {r.tot_Dsub},{r.tot_Dsuz},{r.tot_shared} | int '
        f'{r.int_Dsub},{r.int_Dsuz},{r.int_shared}')
E.to_csv(f'{OUT}/C2_per_gene.tsv.gz', sep='\t')
log.close()
