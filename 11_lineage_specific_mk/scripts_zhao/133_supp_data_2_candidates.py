#!/usr/bin/env python3
"""Script 133: Supplementary Data 2, candidate genes for lineage-specific adaptive protein
evolution and the category decomposition (Results_Lineage_Patterns_Candidate_Genes.md).
Inputs: results/lineage_candidates.tsv (script 130), lineage_candidates_alignment_qc.tsv
(script 131), lineage_category_decomposition.tsv (script 132). Supplementary Data 1 is not changed.
Output: Sato_Supplementary_Data_2.xlsx. Run with envs/docx_tools python (openpyxl).
"""
import pandas as pd

PROJ = '/ru-auth/local/home/lzhao/vscode_projects/suzukii'
RES = f'{PROJ}/results'
POPS = {'JP-Kan': 'Kan', 'JP-Sap': 'Sap', 'JP-Sap2': 'Sap2', 'JP-Tok': 'Tok'}
ARM = {'3R_scaffold': '3R (scaffold NW_023665577.1)'}

c = pd.read_csv(f'{RES}/lineage_candidates.tsv', sep='\t', keep_default_na=False, na_values=[''])
g = pd.read_csv(f'{RES}/lineage_genes.tsv', sep='\t').set_index('dsub_tx')
q = pd.read_csv(f'{RES}/lineage_candidates_alignment_qc.tsv', sep='\t').set_index('dsub_tx')
assert c.dsub_tx.is_unique and set(c.dsub_tx) <= set(q.index) and set(c.dsub_tx) <= set(g.index)
assert len(c) == 252

rows = []
for r in c.itertuples():
    x, a = g.loc[r.dsub_tx], q.loc[r.dsub_tx]
    flags = [f for f in str(r.qc_flags if isinstance(r.qc_flags, str) else '').split(',') if f]
    flags += [f'alignment_{f}' for f in str(a.aln_flags if isinstance(a.aln_flags, str) else '').split(',') if f]
    dn_p = r.Dn_excess_sub_P if r.candidate_class == 'Dsub_lineage' else (
        r.Dn_excess_suz_P if r.candidate_class == 'Dsuz_lineage' else float('nan'))
    if r.candidate_class != 'both' and not dn_p < 0.05:
        flags.append('no_Dn_excess_on_focal_lineage')
    passes = bool(r.robust) and not isinstance(a.aln_flags, str)
    assert passes == (len(flags) == 0), r.symbol
    row = {'class': r.candidate_class, 'passes_QC': passes, 'QC_failures': ','.join(flags),
           'Dmel_symbol': r.symbol, 'Dmel_name': r.name, 'Dmel_FBgn': r.FBgn,
           'Dsub_transcript': r.dsub_tx, 'Dsuz_transcript': r.dsuz_tx,
           'Dsub_chromosome_arm': ARM.get(r.arm, r.arm),
           'Dn_Dsub_lineage': r.Dn_sub, 'Ds_Dsub_lineage': r.Ds_sub, 'Pn_Dsub': r.Pn_sub, 'Ps_Dsub': r.Ps_sub,
           'P_Dsub_lineage': r.P_sub, 'FDR_Dsub_lineage': r.q_sub, 'alpha_Dsub_lineage': x.alpha_sub,
           'Dn_Dsuz_lineage_median': r.Dn_suz, 'Ds_Dsuz_lineage_median': r.Ds_suz,
           'Pn_Dsuz_median': r.Pn_suz, 'Ps_Dsuz_median': r.Ps_suz}
    for p, s in POPS.items():
        row[f'P_Dsuz_lineage_{s}'] = x[f'P_{p}']
    row.update({'Dsuz_populations_significant': int(r.n_pools_pos), 'FDR_Dsuz_lineage_median': r.q_suz_median,
                'lineage_test_P': r.lineage_P, 'lineage_test_FDR': r.lineage_BH,
                'lineage_test_OR': r.lineage_OR, 'Ds_balance_P': r.Ds_balance_P, 'Dn_excess_focal_P': dn_p,
                'callable_codon_ratio_Dsuz_Dsub': r.codon_ratio,
                'divergence_per_codon_Dsub_lineage': r.div_per_codon_sub,
                'divergence_per_codon_Dsuz_lineage': r.div_per_codon_suz,
                'reference_aa_changes_Dsub_lineage': int(a.ref_Dn_sub), 'max_in_30_codons_Dsub': int(a.maxwin_sub),
                'P_cluster_Dsub': a.P_cluster_sub,
                'reference_aa_changes_Dsuz_lineage': int(a.ref_Dn_suz), 'max_in_30_codons_Dsuz': int(a.maxwin_suz),
                'P_cluster_Dsuz': a.P_cluster_suz,
                'FlyAtlas2_max_tissue': r.FlyAtlas2_max, 'ZhaoBegun2017_638': bool(r.ZhaoBegun638),
                'sweep_peak_Dsub': bool(r.sweep_Dsub), 'sweep_peak_Dsuz_populations': int(r.sweep_Dsuz_pools),
                'class_at_DAF5': r.class_DAF5})
    rows.append(row)
t = pd.DataFrame(rows)
order = {'Dsub_lineage': 0, 'Dsuz_lineage': 1, 'both': 2}
t = t.sort_values(['class', 'passes_QC', 'lineage_test_P'], key=lambda s: s.map(order) if s.name == 'class' else s,
                  ascending=[True, False, True])
n = t.groupby('class').passes_QC.agg(['size', 'sum'])
assert (n.loc['Dsub_lineage'].tolist(), n.loc['Dsuz_lineage'].tolist(), n.loc['both'].tolist()) == \
    ([65, 18], [33, 11], [154, 107]), n

d = pd.read_csv(f'{RES}/lineage_category_decomposition.tsv', sep='\t')
dc = d[['category', 'genes', 'DnDs_sub', 'DnDs_suz', 'PnPs_sub', 'PnPs_suz', 'alpha_sub', 'alpha_suz',
        'RD_rel', 'RD_rel_lo', 'RD_rel_hi', 'RD_matched_P', 'DnR_rel', 'DnR_matched_P', 'DsR_rel', 'DsR_matched_P',
        'RP_rel', 'RP_rel_lo', 'RP_rel_hi', 'RP_matched_P', 'dalpha', 'dalpha_lo', 'dalpha_hi', 'dalpha_matched_P']]
dc = dc.rename(columns={'DnDs_sub': 'DnDs_Dsub_lineage', 'DnDs_suz': 'DnDs_Dsuz_lineage', 'PnPs_sub': 'PnPs_Dsub',
                        'PnPs_suz': 'PnPs_Dsuz', 'alpha_sub': 'alpha_Dsub_lineage', 'alpha_suz': 'alpha_Dsuz_lineage',
                        'dalpha': 'alpha_difference_Dsuz_minus_Dsub', 'dalpha_lo': 'alpha_difference_lo',
                        'dalpha_hi': 'alpha_difference_hi', 'dalpha_matched_P': 'alpha_difference_matched_P'})

readme = [
    ('Title', 'Supplementary Data 2 | Candidate genes for lineage-specific adaptive protein evolution, and '
              'divergence and polymorphism components of functional categories.'),
    ('Data', 'Lineage-specific McDonald–Kreitman (MK) counts of Supplementary Data 1 (shared sites removed), '
             'polymorphisms at derived-allele frequency (DAF) ≥ 15%. Genes: D. subpulchrella transcripts with a '
             'reciprocal-best-hit D. suzukii ortholog tested on both lineages (8,204 genes). D. suzukii lineage '
             'counts are medians over the four populations.'),
    ('Sheet candidates', 'The 252 genes selected by the MK criteria below, with quality checks. passes_QC = TRUE '
                         'for 18 D. subpulchrella-lineage, 11 D. suzukii-lineage and 107 both-lineage genes.'),
    ('Class Dsub_lineage', 'FDR < 0.05 and alpha > 0 on the D. subpulchrella lineage; no D. suzukii population '
                           'with P < 0.05 and alpha > 0; lineage test P < 0.05 with the higher Dn/Ds on the '
                           'D. subpulchrella lineage.'),
    ('Class Dsuz_lineage', 'Significant (P < 0.01, alpha > 0) in ≥ 3 of 4 D. suzukii populations with median '
                           'FDR < 0.05; D. subpulchrella lineage not P < 0.05 with alpha > 0; lineage test '
                           'P < 0.05 with the higher Dn/Ds on the D. suzukii lineage.'),
    ('Class both', 'FDR < 0.05 and alpha > 0 on the D. subpulchrella lineage, and significant in ≥ 3 D. suzukii '
                   'populations with median FDR < 0.05. This is stricter than the 347 both-lineage genes of '
                   'Supplementary Data 1 (P < 0.01 on both lineages).'),
    ('Lineage test', 'Two-sided Fisher’s exact test of [Dn, Ds] on the D. subpulchrella lineage against [Dn, Ds] on '
                     'the D. suzukii lineage; FDR (Benjamini–Hochberg) over 8,193 testable genes; OR = ratio of '
                     'Dn/Ds (D. subpulchrella / D. suzukii lineage), 0.5 added to each count.'),
    ('Quality checks', 'Ds_imbalance: Ds split between lineages differs from the genome share 0.503 (two-sided '
                       'binomial P < 0.05). no_Dn_excess_on_focal_lineage: Dn split not shifted toward the focal '
                       'lineage (one-sided binomial P ≥ 0.05 against the genome share 0.542; one-lineage classes '
                       'only). extreme_divergence: (Dn + Ds) per callable codon above the 99.5% quantile (0.1015). '
                       'callable_mismatch: callable-codon ratio outside the 1–99% range (0.932–1.457). '
                       'alignment_clustered: focal-lineage amino-acid changes in the three reference sequences '
                       '(PRANK codon alignments) concentrated in one 30-codon window (permutation P < 0.01, ≥ 5 '
                       'changes; 2,000 permutations, seed 20261104). alignment_near_indels: ≥ 50% of the changes '
                       'within 3 codons of a gap and ≥ 3 times the background fraction.'),
    ('Other columns', 'FlyAtlas2_max_tissue: tissue of maximal expression of the D. melanogaster ortholog (TPM > 1). '
                      'ZhaoBegun2017_638: in the 638 genes with MK P < 0.05 and DoS > 0 of Zhao and Begun (2017). '
                      'sweep_peak_*: transcript within 5 kb of a selective-sweep peak in that species. class_at_DAF5: '
                      'class with polymorphisms at DAF ≥ 5%.'),
    ('Sheet categories', 'Summed counts per category. RD = (Dn/Ds D. subpulchrella lineage) / (Dn/Ds D. suzukii '
                         'lineage); RP = (Pn/Ps D. subpulchrella) / (Pn/Ps D. suzukii); DnR, DsR = ratio of Dn or Ds '
                         'between the lineages; *_rel = category value / genome value. 95% intervals: 2,000 '
                         'bootstrap replicates over genes. matched_P: two-sided, 10,000 random gene sets of the '
                         'same size matched on chromosome arm × decile of informative counts; seed 20261104.'),
    ('Gene identifiers', 'RefSeq transcripts of RU_Dsub_v1.1 (GCF_014743375.2) and CBGP_Dsuzu_IsoJpt1.0 '
                         '(GCF_043229965.1); D. melanogaster FBgn, symbol and name from the FlyBase gene association '
                         'file (2026-05-21).'),
    ('Source', 'scripts 130, 131, 132; this file: scripts/133_supp_data_2_candidates.py.'),
]
out = f'{PROJ}/Sato_Supplementary_Data_2.xlsx'
with pd.ExcelWriter(out, engine='openpyxl') as xw:
    pd.DataFrame(readme, columns=['item', 'description']).to_excel(xw, sheet_name='README', index=False)
    t.to_excel(xw, sheet_name='candidates', index=False)
    dc.to_excel(xw, sheet_name='categories', index=False)
    for ws in xw.book.worksheets:
        ws.freeze_panes = 'A2'
    xw.book['README'].column_dimensions['A'].width = 22
    xw.book['README'].column_dimensions['B'].width = 140
print(out, len(t), 'candidate rows;', n.to_dict())
