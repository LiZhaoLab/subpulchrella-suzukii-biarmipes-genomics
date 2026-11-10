"""Copy number and origin (script 144 rules) for six genes:
RpS14a, Sce, T-cp1 (=CCT1), betaTub97EF, RpL36, VhaAC39. Read-only use of project outputs."""
import importlib.util
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('do', f'{PROJ}/scripts/144_duplication_origin.py')
do = importlib.util.module_from_spec(spec)
spec.loader.exec_module(do)
ev = do.ev

LO = pd.read_csv(f'{ev.OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
E = pd.read_csv(f'{ev.OUT}/C2_per_gene.tsv.gz', sep='\t', index_col=0, low_memory=False)
LO['pgroup'] = LO.dmel.map(E.pgroup)

CG = {'RpS14a': 'CG1524', 'RpS14b': 'CG1527', 'Sce': 'CG5595', 'T-cp1/CCT1': 'CG5374',
      'betaTub97EF': 'CG4869', 'RpL36': 'CG7622', 'VhaAC39-1': 'CG2934'}
ids = {k: f'gene-Dmel_{v}' for k, v in CG.items()}
for k, v in ids.items():
    assert v in E.index, (k, v)
groups = sorted({E.loc[v, 'pgroup'] for v in ids.values()})
print('paralog groups:', {k: E.loc[v, 'pgroup'] for k, v in ids.items()})

G = do.Genomes(LO)


def introns_dmel(L):
    """Intron count of the best alignment of the locus's own D. melanogaster protein (None if absent)."""
    a = G.alns(L)
    q = 'Dmel|' + str(L.dmel)
    a = a[a.q == q]
    if not len(a):
        return None, None
    b = a.sort_values('AS').iloc[-1]
    return len(do.introns_q(b.cg, b.qs)), round((b.qe - b.qs) / b.qlen, 2)


pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 30)
for pg in groups:
    X = G.LO[G.LO.pgroup == pg].copy()
    print(f'\n=== paralog group {pg}: {sorted(X.dmel.str.replace("gene-Dmel_", "").unique())}')
    rows = []
    for _, L in X.iterrows():
        n, cov = introns_dmel(L)
        rows.append(dict(genome=L.genome, dmel=L.dmel.replace('gene-Dmel_', ''), loc=f'{L.tname}:{L.ts + 1}-{L.te}{L.strand}',
                         len_bp=L.te - L.ts, intact=L.intact, seed_cov=round(L.seed_cov, 2), fs=L.seed_fs,
                         stop=L.seed_st, ann=L.annotated_cds, introns=n, qcov=cov, w=L.w,
                         Kan=L['w_JP-Kan'], Sap=L['w_JP-Sap'], Sap2=L['w_JP-Sap2'], Tok=L['w_JP-Tok']))
    print(pd.DataFrame(rows).sort_values(['genome', 'loc']).to_string(index=False))
    for xg, rg in [('DsuzI', 'Dsub'), ('DsuzR', 'Dsub'), ('Dsub', 'Dbia'), ('DsuzI', 'Dbia'), ('DsuzR', 'Dbia')]:
        calls, note = do.classify_genome(G, pg, xg, rg)
        print(f'  {note}')
        for m, d, n in calls:
            print(f'    -> {m}: {d}')
