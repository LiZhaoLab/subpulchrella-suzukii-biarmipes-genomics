#!/usr/bin/env python3
"""Gene-family copy number by lineage (first pass of plan item C2).

Question: do the annotated copy numbers of the C1 gene families (chemoreceptors, detoxification,
cuticle, immunity, mechanosensory channels) differ between D. subpulchrella and D. suzukii by more
than the technical difference between two assemblies of D. suzukii?

Annotations (all NCBI RefSeq GFF3 + protein FASTA; Airi's genome folder):
  Dmel  GCF_000001215.4 Release 6 plus ISO1 MT, FlyBase Release 6.54 (reference for families)
  Dbia  GCF_025231255.1 RU_DBia_V1.1, NCBI Annotation Release 103 (outgroup)
  Dsub  GCF_014743375.2 RU_Dsub_v1.1, NCBI Annotation Release 100
  DsuzI GCF_043229965.1 CBGP_Dsuzu_IsoJpt1.0, RefSeq RS_2025_01
  DsuzR GCF_037355615.1 Dsuz_RU_1.0, RefSeq RS_2024_06 (no protein FASTA: translated with gffread 0.9.9)
One protein per protein-coding gene (longest isoform; ties: first in file).

Mode 'prep' (login node, parsing only): writes per-species longest-protein FASTA, gene tables
(location, placed NC_ vs unplaced NW_/NT_ sequence, NCBI product name), and name-based family
calls. Mode 'analyze' (after 135_gene_family_diamond.sbatch): best D. melanogaster hit per gene
(DIAMOND blastp, highest bitscore, E <= 1e-5) defines the family; Dmel family membership uses the
C1 rules of scripts/106_gene_families.py (symbol regex; mechanosensory list; GO:0006955 and
descendants). Outputs per-family and per-Dmel-gene copy numbers, lineage assignment by parsimony with
Dbia, and flags for assembly duplication (near-identical paralog: >= 98% identity over >= 90% of
both proteins; gene on an unplaced sequence).
Outputs: results/C2_family_copy_number.tsv, results/C2_dmel_gene_copy_number.tsv.gz,
         results/C2_family_copy_number.log
"""
import gzip
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from urllib.parse import unquote

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
GO = f'{PROJ}/inputs/GO_20261022'
OUT = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn'
G = '/ru-auth/local/home/asato/scratch/genome'
GFFREAD = '/ru-auth/local/home/lzhao/bin/gffread/gffread'
SP = {
    'Dmel': (f'{G}/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff',
             f'{G}/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa'),
    'Dbia': (f'{G}/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.gff',
             f'{G}/Dbia/GCF_025231255.1_RU_DBia_V1.1_protein.faa'),
    'Dsub': (f'{G}/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.gff',
             f'{G}/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_protein.faa'),
    'DsuzI': (f'{G}/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.gff',
              f'{G}/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_protein.faa'),
    'DsuzR': (f'{G}/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.gff', None),
}
DSUZR_FNA = f'{G}/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.fna'
# C1 family rules (copied from scripts/106_gene_families.py)
MECHANO = ['nompC', 'iav', 'nan', 'pain', 'pyx', 'TrpA1', 'trp', 'trpl', 'Trpgamma', 'wtrw', 'Trpm',
           'Trpml', 'Pkd2', 'brv1', 'brv2', 'brv3', 'Piezo', 'Tmc']
PATTERNS = {'Or': r'^Or\d', 'Gr': r'^Gr\d', 'Ir': r'^Ir\d', 'Obp': r'^Obp', 'ppk': r'^ppk',
            'Cpr': r'^Cpr', 'Cyp': r'^Cyp', 'Gst': r'^Gst', 'Ugt': r'^Ugt'}
# name-based calls (secondary; naming differs between annotation releases)
NAME_RX = {'Or': r'(?i)\bodorant receptor\b(?! coreceptor)', 'Gr': r'(?i)gustatory receptor',
           'Ir': r'(?i)^(?:putative )?ionotropic receptor', 'Obp': r'(?i)odorant[- ]binding protein',
           'ppk': r'(?i)pickpocket', 'Cyp': r'(?i)cytochrome P450 (?!reductase)',
           'Gst': r'(?i)glutathione S-transferase', 'Ugt': r'(?i)UDP-gl[uy]c\w*transferase'}
FAMS = list(PATTERNS) + ['mechanosensory', 'immune_GO0006955']


def attrs(s):
    return {k: unquote(v) for k, v in (x.split('=', 1) for x in s.rstrip(';').split(';') if '=' in x)}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def parse_gff(path):
    genes, tx2gene, prot2gene, tx_prod, prot_prod = {}, {}, {}, {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            t = p[2]
            if t not in ('gene', 'pseudogene', 'mRNA', 'CDS', 'transcript'):
                continue
            a = attrs(p[8])
            if t in ('gene', 'pseudogene'):
                fb = re.search(r'FLYBASE:(FBgn\d+)', a.get('Dbxref', ''))
                genes[a['ID']] = dict(seqid=p[0], start=int(p[3]), end=int(p[4]), strand=p[6],
                                      biotype=a.get('gene_biotype', ''), symbol=a.get('gene', a.get('Name', '')),
                                      description=a.get('description', ''), fbgn=fb.group(1) if fb else '')
            elif t in ('mRNA', 'transcript'):
                tx2gene[a['ID']] = a.get('Parent', '')
                tx_prod[a['ID']] = a.get('product', '')
            elif t == 'CDS' and 'protein_id' in a:
                par = a.get('Parent', '')
                prot2gene[a['protein_id']] = tx2gene.get(par, par)
                prot_prod[a['protein_id']] = a.get('product', '')
    return genes, tx2gene, prot2gene, tx_prod, prot_prod


def prep():
    os.makedirs(OUT, exist_ok=True)
    for sp, (gff, faa) in SP.items():
        genes, tx2gene, prot2gene, tx_prod, prot_prod = parse_gff(gff)
        if faa is None:  # translate RU_1.0 transcripts
            tfaa = f'{OUT}/{sp}_gffread_proteins.faa'
            if not os.path.exists(tfaa):
                subprocess.run([GFFREAD, gff, '-g', DSUZR_FNA, '-y', tfaa], check=True)
            seqs = read_fasta(tfaa)
            unit2gene, unit_prod = tx2gene, tx_prod
        else:
            seqs = read_fasta(faa)
            unit2gene, unit_prod = prot2gene, prot_prod
        best = {}
        n_unmapped = 0
        for pid, s in seqs.items():
            g = unit2gene.get(pid)
            if g is None or g not in genes:
                n_unmapped += 1
                continue
            s = s.rstrip('*.')
            if '*' in s or '.' in s:  # internal stop (gffread on frameshifted models)
                s = s.replace('*', 'X').replace('.', 'X')
            if g not in best or len(s) > len(best[g][1]):
                best[g] = (pid, s)
        rows = []
        with open(f'{OUT}/{sp}_longest.faa', 'w') as fo:
            for g, (pid, s) in best.items():
                fo.write(f'>{sp}|{g}\n{s}\n')
        for gid, d in genes.items():
            pc = gid in best
            name = unit_prod.get(best[gid][0], '') if pc else ''
            name = re.sub(r' isoform \S+$', '', name) or d['description']
            rows.append(dict(species=sp, gene=gid, protein=best[gid][0] if pc else '',
                             length=len(best[gid][1]) if pc else 0, coding=pc, name=name, **d))
        df = pd.DataFrame(rows)
        df['placed'] = ~df.seqid.str.startswith('NW_')  # NC_ and NT_ (Dmel arms) are chromosomes
        for fam, rx in NAME_RX.items():
            df[f'name_{fam}'] = df.name.fillna('').str.contains(rx, regex=True)
        df.to_csv(f'{OUT}/{sp}_genes.tsv.gz', sep='\t', index=False)
        pcg = df[df.biotype == 'protein_coding']
        print(f'{sp}: genes {len(genes)}; protein-coding {len(pcg)}, with protein {pcg.coding.sum()}; '
              f'pseudogenes {(df.biotype == "pseudogene").sum()}; protein records not mapped to a gene {n_unmapped}; '
              f'protein-coding on unplaced sequences {(~pcg.placed).sum()}', flush=True)


def go_descendants(root):
    children = defaultdict(set)
    cur = None
    with open(f'{GO}/go-basic.obo') as f:
        for line in f:
            line = line.strip()
            if line == '[Term]':
                cur = None
            elif line.startswith('id: GO:'):
                cur = line[4:]
            elif cur and line.startswith('is_a: '):
                children[line[6:16]].add(cur)
            elif cur and line.startswith('relationship: part_of '):
                children[line[22:32]].add(cur)
    out, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def dmel_families(dm):
    imm_terms = go_descendants('GO:0006955')
    imm = set()
    with gzip.open(f'{GO}/fb.gaf.gz', 'rt') as f:
        for line in f:
            if line.startswith('!'):
                continue
            t = line.split('\t')
            if 'NOT' in t[3]:
                continue
            if t[4] in imm_terms:
                imm.add(t[1])
    fam = {}
    for _, r in dm.iterrows():
        sym = str(r.symbol)
        hit = [n for n, pat in PATTERNS.items() if re.search(pat, sym)]
        if sym in MECHANO:
            hit.append('mechanosensory')
        if r.fbgn in imm:
            hit.append('immune_GO0006955')
        fam[r.gene] = hit
    return fam


def analyze():
    log = open(f'{RES}/C2_family_copy_number.log', 'w')

    def say(*a):
        s = ' '.join(str(x) for x in a)
        print(s, flush=True)
        log.write(s + '\n')

    say('C2 first pass: gene-family copy number from NCBI annotations; DIAMOND best hit to D. melanogaster')
    tabs = {sp: pd.read_csv(f'{OUT}/{sp}_genes.tsv.gz', sep='\t', low_memory=False) for sp in SP}
    dm = tabs['Dmel'][tabs['Dmel'].coding]
    fam_of = dmel_families(dm)
    sym = dict(zip(dm.gene, dm.symbol))
    cols = ['q', 's', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']
    assign = {}
    for sp in SP:
        h = pd.read_csv(f'{OUT}/{sp}_vs_Dmel.tsv', sep='\t', names=cols)
        h = h[h.evalue <= 1e-5].sort_values(['q', 'bitscore'], ascending=[True, False]).drop_duplicates('q')
        h['gene'] = h.q.str.split('|').str[1]
        h['dmel'] = h.s.str.split('|').str[1]
        h['qcov'] = h.length / h.qlen
        assign[sp] = h.set_index('gene')[['dmel', 'pident', 'qcov', 'bitscore']]
        t = tabs[sp]
        n = t.coding.sum()
        say(f'{sp}: protein-coding genes with protein {n}; with a D. melanogaster hit (E <= 1e-5) {len(h)} '
            f'({len(h) / n:.1%})')
    # sanity: Dmel genes should hit themselves
    a = assign['Dmel']
    self_hit = (a.index == a.dmel).mean()
    say(f'sanity: Dmel genes whose best hit is the gene itself: {self_hit:.4f}')
    assert self_hit > 0.97, self_hit

    # near-identical paralogs within each genome (assembly duplication flag)
    dup = {}
    for sp in SP:
        h = pd.read_csv(f'{OUT}/{sp}_self.tsv', sep='\t', names=cols)
        h = h[(h.q != h.s) & (h.pident >= 98) & (h.length >= 0.9 * h.qlen) & (h.length >= 0.9 * h.slen)]
        dup[sp] = set(h.q.str.split('|').str[1])
        say(f'{sp}: protein-coding genes with a near-identical paralog (>= 98% identity, >= 90% of both): '
            f'{len(dup[sp])}')

    # per-gene table
    rows = []
    for sp in SP:
        t = tabs[sp][tabs[sp].coding].set_index('gene')
        a = assign[sp]
        for g, r in a.iterrows():
            fams = fam_of.get(r.dmel, [])
            rows.append(dict(species=sp, gene=g, dmel=r.dmel, dmel_symbol=sym.get(r.dmel, ''),
                             pident=r.pident, qcov=r.qcov, families=','.join(fams),
                             placed=bool(t.at[g, 'placed']), near_identical_paralog=g in dup[sp],
                             name=t.at[g, 'name']))
    pg = pd.DataFrame(rows)
    pg.to_csv(f'{OUT}/C2_gene_assignments.tsv.gz', sep='\t', index=False)

    # copy number per Dmel gene
    cn = pg.pivot_table(index='dmel', columns='species', values='gene', aggfunc='count', fill_value=0)
    cn = cn.reindex(dm.gene, fill_value=0)[list(SP)]
    cn['symbol'] = cn.index.map(sym)
    cn['families'] = cn.index.map(lambda g: ','.join(fam_of.get(g, [])))
    cn.to_csv(f'{RES}/C2_dmel_gene_copy_number.tsv.gz', sep='\t')

    # family totals
    out = []
    say('\nFamily copy number (genes whose best D. melanogaster hit is a family member; '
        'in brackets: on unplaced sequences / with a near-identical paralog)')
    say(f'{"family":18s} ' + ' '.join(f'{s:>14s}' for s in SP))
    for fam in FAMS:
        r = dict(family=fam)
        cells = []
        for sp in SP:
            sub = pg[(pg.species == sp) & pg.families.str.split(',').apply(lambda x: fam in x)]
            r[sp] = len(sub)
            r[f'{sp}_unplaced'] = int((~sub.placed).sum())
            r[f'{sp}_nearident'] = int(sub.near_identical_paralog.sum())
            cells.append(f'{len(sub):>4d} [{r[f"{sp}_unplaced"]:>2d}/{r[f"{sp}_nearident"]:>2d}]')
        # name-based counts (protein-coding genes with a matching NCBI product name)
        if fam in NAME_RX:
            for sp in SP:
                t = tabs[sp]
                r[f'{sp}_by_name'] = int((t.coding & t[f'name_{fam}']).sum())
        out.append(r)
        say(f'{fam:18s} ' + ' '.join(f'{c:>14s}' for c in cells))
    fam_df = pd.DataFrame(out)
    say('\nName-based counts (NCBI product names; not comparable between annotation releases):')
    for _, r in fam_df.iterrows():
        if f'Dsub_by_name' in r and not pd.isna(r.get('Dsub_by_name')):
            say(f'{r.family:18s} ' + ' '.join(f'{sp} {int(r[f"{sp}_by_name"])}' for sp in SP))

    # lineage assignment per Dmel gene, genome-wide and per family
    c = cn[cn[['Dbia', 'Dsub', 'DsuzI', 'DsuzR']].sum(axis=1) > 0].copy()
    c['suz_agree'] = c.DsuzI == c.DsuzR
    c['sub_ne_suz'] = c.suz_agree & (c.Dsub != c.DsuzI)
    c['tech'] = ~c.suz_agree
    c['lineage'] = np.where(~c.sub_ne_suz, '',
                            np.where(c.Dbia == c.DsuzI, 'Dsub', np.where(c.Dbia == c.Dsub, 'Dsuz', 'unresolved')))
    c['direction'] = np.where(~c.sub_ne_suz, '',
                              np.where(c.lineage == 'Dsub', np.where(c.Dsub > c.Dbia, 'gain', 'loss'),
                                       np.where(c.lineage == 'Dsuz', np.where(c.DsuzI > c.Dbia, 'gain', 'loss'),
                                                '')))
    say(f'\nD. melanogaster genes with >= 1 hit in Dbia/Dsub/Dsuz: {len(c)}')
    say(f'  two D. suzukii assemblies disagree (technical): {c.tech.sum()} ({c.tech.mean():.2%})')
    say(f'  D. suzukii assemblies agree and Dsub differs: {c.sub_ne_suz.sum()} ({c.sub_ne_suz.mean():.2%}); '
        f'lineage: {dict(Counter(c.lineage[c.sub_ne_suz]))}')
    say(f'  direction: {dict(Counter(zip(c.lineage[c.sub_ne_suz], c.direction[c.sub_ne_suz])))}')
    say('\nPer family: Dmel genes; technical disagreement; Dsub != Dsuz (both assemblies); '
        'Fisher test of Dsub != Dsuz rate vs all other Dmel genes')
    pv = []
    for fam in FAMS:
        m = c.families.str.split(',').apply(lambda x: fam in x)
        a_, b_ = int(c.sub_ne_suz[m].sum()), int((~c.sub_ne_suz[m]).sum())
        c_, d_ = int(c.sub_ne_suz[~m].sum()), int((~c.sub_ne_suz[~m]).sum())
        orr, p = fisher_exact([[a_, b_], [c_, d_]])
        pv.append(p)
        lin = Counter(zip(c.lineage[m & c.sub_ne_suz], c.direction[m & c.sub_ne_suz]))
        say(f'{fam:18s} n {m.sum():4d}; tech {c.tech[m].sum():3d} ({c.tech[m].mean():.1%}); '
            f'Dsub!=Dsuz {a_:3d} ({a_ / max(m.sum(), 1):.1%}); OR {orr:.2f} P {p:.3g}; {dict(lin)}')
    from statsmodels.stats.multitest import multipletests
    q = multipletests(pv, method='fdr_bh')[1]
    say('BH FDR over families: ' + ', '.join(f'{f} {x:.3g}' for f, x in zip(FAMS, q)))
    fam_df['BH_q_sub_ne_suz'] = q
    fam_df.to_csv(f'{RES}/C2_family_copy_number.tsv', sep='\t', index=False)
    c.to_csv(f'{OUT}/C2_lineage_calls.tsv.gz', sep='\t')

    # list family genes with lineage calls
    say('\nFamily genes with Dsub != Dsuz (both D. suzukii assemblies agree): Dmel symbol, copies '
        'Dmel/Dbia/Dsub/DsuzI/DsuzR, lineage, direction')
    fm = c[c.sub_ne_suz & (c.families != '')]
    for g, r in fm.sort_values(['families', 'symbol']).iterrows():
        say(f'  {r.symbol:12s} {r.families:28s} {r.Dmel}/{r.Dbia}/{r.Dsub}/{r.DsuzI}/{r.DsuzR} {r.lineage} {r.direction}')
    log.close()


if __name__ == '__main__':
    {'prep': prep, 'analyze': analyze}[sys.argv[1]]()
