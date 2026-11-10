#!/usr/bin/env python3
"""Export the Zhao & Begun (2017, PLoS Genet 13:e1007016) reference gene sets to TSV.

- pgen.1007016.s015.xlsx = Table S13, sheet "DOS": 638 D. melanogaster genes with MK P < 0.05
  and DoS > 0 (D. melanogaster-D. simulans).
- pgen.1007016.s017.xlsx = Table S15, sheet "gene": 66 genes with evidence of adaptive protein
  evolution in both the D. melanogaster/D. simulans and D. hydei/D. mojavensis comparisons;
  D. melanogaster FBgn in column Dmel_gene_name.
Needs openpyxl (run with the base anaconda python; jupyter-ml lacks it).
Writes inputs/ZhaoBegun2017_S13_DOS_638.tsv and inputs/ZhaoBegun2017_S15_both_66.tsv (FBgn, name).
"""
import openpyxl

P = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'


def rows(fname, sheet):
    ws = openpyxl.load_workbook(f'{P}/{fname}', data_only=True)[sheet]
    r = [[c.value for c in row] for row in ws.iter_rows()]
    return r[0], [x for x in r[1:] if any(v is not None for v in x)]


h, dos = rows('pgen.1007016.s015.xlsx', 'DOS')
assert h[0] == 'mel_gene' and h[2] == 'name' and h[11] == 'Direction of selection'
assert all(r[11] > 0 for r in dos)
s13 = [(r[0], r[2]) for r in dos]

h, both = rows('pgen.1007016.s017.xlsx', 'gene')
assert h[9] == 'Dmel_gene_name' and h[7] == 'DoS'
s15 = [(r[9], r[11]) for r in both]

for name, data, n in (('ZhaoBegun2017_S13_DOS_638.tsv', s13, 638),
                      ('ZhaoBegun2017_S15_both_66.tsv', s15, 66)):
    ids = [g for g, _ in data]
    assert len(ids) == n == len(set(ids)) and all(g.startswith('FBgn') for g in ids), name
    with open(f'{P}/inputs/{name}', 'w') as f:
        f.write('FBgn\tname\n')
        f.writelines(f'{g}\t{s}\n' for g, s in data)
    print(name, len(ids))
assert {g for g, _ in s15} <= {g for g, _ in s13}
