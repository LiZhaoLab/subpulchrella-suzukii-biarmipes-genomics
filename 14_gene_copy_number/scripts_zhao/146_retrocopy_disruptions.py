#!/usr/bin/env python3
"""Are the frameshifts and stop codons in the four D. suzukii retrocopies real and frequent in the population?

Chiu et al. (2013, G3 3:2257) reported complete ORFs for D. suzukii retrogenes of RpS14a, T-cp1, RpL36 and VhaAC39
(WT3 strain). In IsoJpt1.0 the miniprot alignments of these copies carry frameshifts or stops
(Results_Retrocopies_Six_Genes.md). This applies the event-level test of scripts/142 (artefact filters,
DNA comparison with intact copies of D. biarmipes and D. subpulchrella, frequency in the four D. suzukii pools).
The RpL36, RpS14 and CCT1 retrocopies are loci of scripts/137; the VhaAC39-1 retrocopy is hidden inside the
CG40006 locus (clustering without exon test), so a locus record is built from its D. melanogaster protein
alignment (seed = Dmel|gene-Dmel_CG2934). Run with envs/popgen/bin/python. No random numbers.
Output: results/retrocopies_six_genes/retrocopy_disruption_events.tsv
"""
import importlib.util

import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
spec = importlib.util.spec_from_file_location('ev', f'{PROJ}/scripts/142_pseudogene_event_support.py')
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)

LO = pd.read_csv(f'{ev.OUT}/C2_loci.tsv.gz', sep='\t', low_memory=False)
RETRO = [('RpL36', 'NW_027255939.1', 1918068), ('RpS14', 'NW_027255927.1', 2662679),
         ('CCT1', 'NW_027255899.1', 1547552)]
names = {}
for g, t, ts in RETRO:
    x = LO[(LO.genome == 'DsuzI') & (LO.tname == t) & (LO.ts == ts)]
    assert len(x) == 1, (g, len(x))
    names[x['name'].iloc[0]] = g
# hidden VhaAC39-1 retrocopy: locus record from the Dmel CG2934 alignment
v = dict(LO[(LO.genome == 'DsuzI') & (LO.dmel == 'gene-Dmel_CG2934')].iloc[0])
v.update(locus=-1, tname='NW_027255896.1', strand='-', ts=24702031, te=24703088, seed_q='Dmel|gene-Dmel_CG2934',
         intact=False, name='DsuzI|VhaAC39-1_retrocopy', annotated_cds=None)
LO = pd.concat([LO, pd.DataFrame([v])], ignore_index=True)
names['DsuzI|VhaAC39-1_retrocopy'] = 'VhaAC39-1'

ctx = ev.Context(LO)
rows = []
for nm, g in names.items():
    L = LO[LO['name'] == nm].iloc[0]
    for r in ev.locus_events(L, ctx):
        rows.append(dict(gene=g, **r))
E = pd.DataFrame(rows)
out = f'{PROJ}/results/retrocopies_six_genes/retrocopy_disruption_events.tsv'
E.to_csv(out, sep='\t', index=False, float_format='%.3f')
pd.set_option('display.width', 250)
pd.set_option('display.max_colwidth', 80)
cols = ['gene', 'type', 'g0', 'shift', 'qpos', 'qlen', 'intron_nt', 'dna_confirmed', 'artefact_filter'] + \
       [c for c in E.columns if c.startswith('freq_') or c.startswith('depth_')] + ['population', 'supported']
print(E[cols].to_string(index=False))
print(E[['gene', 'dna_test']].to_string(index=False))
