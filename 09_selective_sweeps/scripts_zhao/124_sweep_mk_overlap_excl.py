#!/usr/bin/env python3
"""Overlap of sweep-peak genes with lineage-specific MK-positive genes, with the current MK sets.

Script 101 used results/L3_lineage_gene_sets.tsv, i.e. the polarized MK test
before shared polymorphisms were excluded. The manuscript now uses the EXCLUDE_SHARED=1 variant
(results/polarized_mk_excl/summary/L3_lineage_gene_sets.tsv). This script
repeats only the overlap test of script 101 (same rules: DAF15; Dsub_lineage for D. subpulchrella,
Dsuz_lineage_union for every D. suzukii pool; transcripts within 5 kb of a peak, taken from
results/F_sweep_peaks.tsv; MK-tested transcripts only; Fisher's exact test, descriptive)
for both gene-set versions. The old version must reproduce the earlier log exactly.
Output: results/F_sweep_mk_overlap.{tsv,log}
"""
import pandas as pd
from scipy.stats import fisher_exact

RES = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results'
SETS = {'original': f'{RES}/L3_lineage_gene_sets.tsv',
        'excl': f'{RES}/polarized_mk_excl/summary/L3_lineage_gene_sets.tsv'}
EXPECT_ORIGINAL = {'Dsub_L5_7': (267, 10), 'JP-Kan': (1537, 143), 'JP-Sap': (344, 32),
                   'JP-Sap2': (433, 36), 'JP-Tok': (456, 41)}   # earlier log: peak genes, MK-positive


def load_sets(path):
    gs = pd.read_csv(path, sep='\t')
    mkpos, tested = {'Dsub': set(), 'Dsuz': set()}, {'Dsub': set(), 'Dsuz': set()}
    for setname, sp in (('Dsub_lineage', 'Dsub'), ('Dsuz_lineage_union', 'Dsuz')):
        g = gs[(gs.version == 'DAF15') & (gs.set == setname)]
        assert len(g), (path, setname)
        for t, flag in zip(g.transcripts.astype(str), g.MK_positive):
            txs = [x.strip() for x in t.split(',')]
            tested[sp].update(txs)
            if flag:
                mkpos[sp].update(txs)
    return mkpos, tested


def main():
    pk = pd.read_csv(f'{RES}/F_sweep_peaks.tsv', sep='\t')
    rows = []
    with open(f'{RES}/F_sweep_mk_overlap.log', 'w') as logf:
        for ver, path in SETS.items():
            mkpos, tested = load_sets(path)
            for pool, g in pk.groupby('pool'):
                sp = 'Dsub' if pool.startswith('Dsub') else 'Dsuz'
                inpk = set(x for t in g.transcripts.fillna('').astype(str) for x in t.split(',') if x) & tested[sp]
                a = len(inpk & mkpos[sp])
                b = len(inpk) - a
                c = len(mkpos[sp] & tested[sp]) - a
                d = len(tested[sp]) - a - b - c
                orr, p = fisher_exact([[a, b], [c, d]])
                if ver == 'original':
                    assert (len(inpk), a) == EXPECT_ORIGINAL[pool], (pool, len(inpk), a)
                rows.append(dict(version=ver, pool=pool, peak_tx=len(inpk), peak_mk_pos=a,
                                 peak_frac=a / len(inpk), genome_frac=(a + c) / len(tested[sp]),
                                 tested_tx=len(tested[sp]), odds_ratio=orr, fisher_p=p))
                s = (f'{ver:8s} {pool:10s} peak transcripts {len(inpk):5d}, MK-positive {a:4d} '
                     f'({a / len(inpk):.3f}) vs genome {(a + c) / len(tested[sp]):.3f} of {len(tested[sp])}; '
                     f'odds ratio {orr:.2f}, Fisher P {p:.2g}')
                print(s)
                logf.write(s + '\n')
    pd.DataFrame(rows).to_csv(f'{RES}/F_sweep_mk_overlap.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
