#!/usr/bin/env python3
"""Per-site divergence and polymorphism rates behind the polarized MK counts (Table 7).

Table 7 reports Dn/Ds and Pn/Ps as ratios of COUNTS (pathway-averaged changes, script 60).
There are about three nonsynonymous sites per synonymous site, so a count ratio is not dN/dS.
This script:
  1. counts synonymous (Ls) and nonsynonymous (Ln) sites per gene on the focal reference CDS
     (Nei-Gojobori 1986; changes to stop codons excluded, Ls + Ln = 3 per codon), scaled to
     the callable codons of each gene (Ls_gene * callable / codons; approximation, because
     the callable-codon list is not saved);
  2. reports Ka = Dn/Ln, Ks = Ds/Ls (raw and Jukes-Cantor corrected), pN = Pn/Ln, pS = Ps/Ls
     (DAF >= 15%) for (a) the tested genes of Table 7 and (b) all genes with callable codons;
  3. tests three explanations for a high Ka/Ks:
     - gene filter: Table 7 keeps only genes with no empty row or column, which removes
       genes without nonsynonymous changes;
     - polarization: codons where the outgroup differs from both ingroup states are skipped;
       compares both-lineage polarized sums with my unpolarized recount (script 20);
     - a few genes with very many nonsynonymous differences (annotation or alignment errors);
  4. compares count and per-site ratios over increasing divergence (Dsub vs Dsuz, Dbia, Dmel;
     my unpolarized recount from Airi's alignments) to test whether Ka/Ks is higher on short
     branches.
Outputs: results/mk_per_site_rates.tsv, results/mk_per_site_rates.log
(with PMK_VARIANT=excl: the same names in results/polarized_mk_excl/summary/).
"""
import importlib.util
import math
import os

import numpy as np
import pandas as pd

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
VARIANT = os.environ.get('PMK_VARIANT', '')  # 'excl': run (scripts/114), shared sites excluded
R = f'{P}/results/polarized_mk' + (f'_{VARIANT}' if VARIANT else '')
OUTD = f'{R}/summary' if VARIANT else f'{P}/results'
POPS = ['JP-Kan', 'JP-Sap', 'JP-Sap2', 'JP-Tok']
spec = importlib.util.spec_from_file_location('pm', f'{P}/scripts/60_polarized_mk.py')
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)
CODE = pm.mk.CODE
log = []


def say(s=''):
    print(s, flush=True)
    log.append(str(s))


def codon_sites():
    """Nei-Gojobori (S, N) per sense codon, changes to stop codons excluded."""
    out = {}
    for c, aa in CODE.items():
        if aa == '*':
            continue
        s = n = 0.0
        for i in range(3):
            syn = non = 0
            for b in 'ACGT':
                if b == c[i]:
                    continue
                m = c[:i] + b + c[i + 1:]
                if CODE[m] == '*':
                    continue
                if CODE[m] == aa:
                    syn += 1
                else:
                    non += 1
            s += syn / (syn + non)
            n += non / (syn + non)
        out[c] = (s, n)
    return out


SITES = codon_sites()
assert abs(SITES['TTT'][0] - 1 / 3) < 1e-9 and abs(sum(SITES['CTG']) - 3) < 1e-9


def gene_sites(sp, genes):
    """Return DataFrame gene -> codons, Ls, Ln on the reference CDS (stop codons skipped)."""
    regions = f'{pm.A}/population_polymorphism/reanalysis_20260821/config/regions'
    major = {l.split()[0] for cls in ('autosome', 'X') for l in open(f'{regions}/{sp}.{cls}.bed')}
    models = pm.cds_models(pm.GFF[sp], major)
    need = {models[g][0] for g in genes if g in models}
    seqs, cur, keep = {}, None, False
    with open(pm.FNA[sp]) as f:
        for line in f:
            if line.startswith('>'):
                cur = line[1:].split()[0]
                keep = cur in need
                if keep:
                    seqs[cur] = []
            elif keep:
                seqs[cur].append(line.strip().upper())
    seqs = {k: ''.join(v) for k, v in seqs.items()}
    rows = []
    for g in genes:
        chrom, pos, strand = models[g]
        s = ''.join(seqs[chrom][p - 1] for p in pos)
        if strand == '-':
            s = ''.join(pm.COMP.get(b, 'N') for b in s)
        ls = ln = 0.0
        nc = 0
        for k in range(0, len(s) - len(s) % 3, 3):
            c = s[k:k + 3]
            if c in SITES:
                ls += SITES[c][0]
                ln += SITES[c][1]
                nc += 1
        rows.append(dict(gene=g, ref_codons=nc, Ls_ref=ls, Ln_ref=ln))
    return pd.DataFrame(rows).set_index('gene')


def jc(p):
    return -0.75 * math.log(1 - 4 * p / 3) if p < 0.75 else float('nan')


def poly_counts(focal, genes, thr=0.15):
    p = pd.read_csv(f'{R}/{focal}.polymorphisms.tsv.gz', sep='\t')
    f = p.daf.where(p.daf.notna(), p.maf)
    p = p[(f >= thr) & p.gene.isin(genes)]
    c = p.groupby(['gene', 'type']).size().unstack(fill_value=0).reindex(genes, fill_value=0)
    return c.get('n', 0), c.get('s', 0)


def rates(label, Pn, Ps, Dn, Ds, Ln, Ls, n_genes):
    ka, ks = Dn / Ln, Ds / Ls
    return dict(set=label, genes=n_genes, Pn=round(Pn), Ps=round(Ps), Dn=round(Dn), Ds=round(Ds),
                Ln=round(Ln), Ls=round(Ls), Ln_over_Ls=round(Ln / Ls, 3),
                PnPs_count=round(Pn / Ps, 3) if Ps else np.nan,
                DnDs_count=round(Dn / Ds, 3),
                pN=Pn / Ln, pS=Ps / Ls, pNpS=round((Pn / Ln) / (Ps / Ls), 4) if Ps else np.nan,
                Ka=ka, Ks=ks, KaKs=round(ka / ks, 4), KaKs_JC=round(jc(ka) / jc(ks), 4),
                alpha=round(1 - (Ds * Pn) / (Dn * Ps), 3) if Ps else np.nan)


rows = []
sites = {}
for focal in ['Dsub_L5_7'] + POPS:
    sp = 'Dsub' if focal.startswith('Dsub') else 'Dsuz'
    g = pd.read_csv(f'{R}/{focal}.genes.tsv', sep='\t').set_index('gene')
    g = g[g.callable > 0]
    miss = [x for x in g.index if sp not in sites or x not in sites[sp].index]
    if miss:
        sites[sp] = pd.concat([sites.get(sp), gene_sites(sp, miss)])
    s = sites[sp].reindex(g.index)
    assert s.notna().all().all(), focal
    # script 60 drops the stop codon when the alignment omits it; allow that difference
    assert ((s.ref_codons - g.codons).abs() <= 1).mean() > 0.99, focal
    frac = g.callable / g.codons
    assert (frac <= 1).all()
    g['Ls'] = s.Ls_ref * frac
    g['Ln'] = s.Ln_ref * frac
    g['Pn15'], g['Ps15'] = poly_counts(focal, list(g.index))
    t = pd.read_csv(f'{R}/{focal}.mk_tests.DAF15.tsv', sep='\t').set_index('gene')
    # tested set reproduces Table 7 sums
    tt = g.loc[t.index]
    assert int(t.Pn.sum()) == int(tt.Pn15.sum()) and int(t.Ps.sum()) == int(tt.Ps15.sum())
    assert int(t.Dn.sum()) == int(tt.Dn_focal.round().sum())
    lin = 'D. subpulchrella' if sp == 'Dsub' else f'D. suzukii ({focal[3:]})'
    rows.append(dict(focal=focal, lineage=lin, callable_frac=round(frac.mean(), 4),
                     **rates('tested (Table 7)', t.Pn.sum(), t.Ps.sum(), t.Dn.sum(), t.Ds.sum(),
                             tt.Ln.sum(), tt.Ls.sum(), len(t))))
    rows.append(dict(focal=focal, lineage=lin, callable_frac=round(frac.mean(), 4),
                     **rates('all genes', g.Pn15.sum(), g.Ps15.sum(), g.Dn_focal.sum(),
                             g.Ds_focal.sum(), g.Ln.sum(), g.Ls.sum(), len(g))))
    # concentration of Dn: drop the 1% of tested genes with the highest Dn per callable codon
    dn_rate = tt.Dn_focal / tt.callable
    top = dn_rate.rank(ascending=False, method='first') <= max(1, round(0.01 * len(tt)))
    rest = tt[~top]
    rows.append(dict(focal=focal, lineage=lin, callable_frac=round(frac.mean(), 4),
                     share_Dn_top1pct=round(tt.Dn_focal[top].sum() / tt.Dn_focal.sum(), 4),
                     **rates('tested minus top 1% Dn/codon', rest.Pn15.sum(), rest.Ps15.sum(),
                             rest.Dn_focal.sum(), rest.Ds_focal.sum(), rest.Ln.sum(),
                             rest.Ls.sum(), len(rest))))
    if focal == 'Dsub_L5_7':
        # both lineages together and the skipped (unpolarized) codons
        both = rates('both lineages, polarized codons', 0, 0, g.Dn_focal.sum() + g.Dn_sister.sum(),
                     g.Ds_focal.sum() + g.Ds_sister.sum(), g.Ln.sum(), g.Ls.sum(), len(g))
        both.update(focal=focal, lineage='Dsub + Dsuz',
                    unpolarized_codon_frac=round(g.unpolarized_codons.sum()
                                                 / (g.unpolarized_codons.sum() + g.callable.sum()), 4))
        rows.append(both)
        # D. suzukii lineage scored from the single D. suzukii reference sequence (no pool):
        # its within-species polymorphism is counted as divergence
        sis = rates('sister lineage, single reference sequence', 0, 0, g.Dn_sister.sum(),
                    g.Ds_sister.sum(), g.Ln.sum(), g.Ls.sum(), len(g))
        sis.update(focal=focal, lineage='D. suzukii (reference genome)')
        rows.append(sis)
        gsub = g

res = pd.DataFrame(rows)

# Unpolarized divergence, Dsub ingroup vs each outgroup: my recount (script 20) from Airi's
# PRANK alignments. Airi's cmktest sums are not used: 2-7 genes per table carry impossible
# counts (Ds about 4.2 million; results/L2_cmktest_impossible_counts.log) that dominate sums.
say('== count ratio vs per-site ratio across divergence (my unpolarized recount, Dsub ingroup)')
div = []
for og in ['Dsuz', 'Dbia', 'Dmel']:
    a = pd.read_csv(f'{P}/results/mk_recount/Dsub_vs_{og}.tsv', sep='\t').set_index('id')
    assert ((a.Dn + a.Ds) <= 3 * a.Ncodon).all() and ((a.Pn + a.Ps) <= 3 * a.Ncodon).all(), og
    common = a.index.intersection(gsub.index)
    sref = sites['Dsub'].loc[common]
    x = a.loc[common]
    # sites per alignment codon from the reference CDS composition
    ls = (sref.Ls_ref / sref.ref_codons * x.Ncodon).sum()
    ln = (sref.Ln_ref / sref.ref_codons * x.Ncodon).sum()
    d = rates(f'Dsub vs {og} unpolarized', x.Pn.sum(), x.Ps.sum(), x.Dn.sum(), x.Ds.sum(), ln, ls,
              len(x))
    d.update(focal='Dsub (recount)', lineage=f'Dsub-{og} total')
    div.append(d)
    if og == 'Dsuz':
        # same genes, my polarized both-lineage sums (skipped codons excluded)
        gg = gsub.loc[common]
        dd = (gg.Dn_focal.sum() + gg.Dn_sister.sum()) / (gg.Ds_focal.sum() + gg.Ds_sister.sum())
        say(f'Dsub-Dsuz genes in both: {len(gg)}; unpolarized recount Dn/Ds (count) '
            f'{x.Dn.sum() / x.Ds.sum():.3f}; polarized, both lineages summed {dd:.3f}; '
            f'Ds per kcodon unpolarized {x.Ds.sum() / x.Ncodon.sum() * 1000:.1f} vs polarized '
            f'both lineages {(gg.Ds_focal.sum() + gg.Ds_sister.sum()) / gg.callable.sum() * 1000:.1f}; '
            f'Dn per kcodon {x.Dn.sum() / x.Ncodon.sum() * 1000:.1f} vs '
            f'{(gg.Dn_focal.sum() + gg.Dn_sister.sum()) / gg.callable.sum() * 1000:.1f}')
divdf = pd.DataFrame(div)
res = pd.concat([res, divdf], ignore_index=True)
cols = ['focal', 'lineage', 'set', 'genes', 'Pn', 'Ps', 'Dn', 'Ds', 'Ln', 'Ls', 'Ln_over_Ls',
        'PnPs_count', 'pNpS', 'DnDs_count', 'Ka', 'Ks', 'KaKs', 'KaKs_JC', 'alpha',
        'share_Dn_top1pct', 'unpolarized_codon_frac', 'callable_frac']
res = res[[c for c in cols if c in res.columns]]
res.to_csv(f'{OUTD}/mk_per_site_rates.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 30)
show = res.copy()
for c in ['Ka', 'Ks']:
    show[c] = show[c].map(lambda v: f'{v:.4f}')
say(show.drop(columns=['Pn', 'Ps', 'Dn', 'Ds', 'Ln', 'Ls']).to_string(index=False))
with open(f'{OUTD}/mk_per_site_rates.log', 'w') as fh:
    fh.write('\n'.join(log) + '\n')
