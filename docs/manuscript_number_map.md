# Where each reported result comes from

Paragraph numbers refer to the Results and Methods of the manuscript version Sato1106 (with the tracked changes of
2026-11-06). Paths are relative to this folder. "Sato" = `scripts_sato`, "Zhao" = `scripts_zhao`.

## Results

| Manuscript statement | Script(s) | Result file(s) |
|---|---|---|
| Assembly sizes, N50, BUSCO, gene counts (assembly section) | not included (NCBI assembly reports; BUSCO) | — |
| Synteny: 99.78% / 99.85% / 99.75% on homologous Muller elements; 59, 76, 194 segments and 53, 70, 189 breakpoints; densities per 10 Mb; X densities; 33 / 19 / 1 breakpoints assigned by lineage, binomial P = 0.07 | `02_synteny`: Sato 02–08 (RepeatMasker, PROmer, one-to-one blocks); Zhao `91_synteny_breakpoints.py`, `122_synteny_paragraph_numbers.py` | `02_synteny/results_zhao/synteny_summary.log`, `synteny_segments.tsv`, `synteny_paragraph_numbers.log` |
| Six of the 33 *D. subpulchrella*-lineage breakpoints are ends of three polymorphic inversions | `08_inversions`: Zhao `111b_inversion_synteny.py` | `08_inversions/results_zhao/INV_synteny_check.tsv` |
| Centre of IsoJpt1.0 chromosome 3 (27–62 Mb): 74% TE, 15.8 genes per Mb; unmasked alignment 14% (*D. subpulchrella*) and 32% (Dsuz_RU_1.0) | Sato `11_quantify_NC092082_fusion.sh`, `11_build_synteny_TE_gene_tracks_200kb.py`; Zhao `121_chr3_centre_minimap2.sbatch`, `122_synteny_paragraph_numbers.py` | `02_synteny/results_sato/fusion/*`, `02_synteny/results_zhao/chr3_vs_*.paf`, `synteny_paragraph_numbers.log` |
| Contamination: none in the *D. subpulchrella* pool; 1.7% (1.66–1.76%) *D. subpulchrella* DNA in Tok, 1.58% on the dot chromosome, 1.15% in mtDNA; ≤ 0.04% and ≤ 0.012% in the other pools; +1% π, +2% θW in Tok | `06_contamination`: Zhao 60 (evidence for the *D. subpulchrella* pool), 65, 77, 78, 79, 87 | `06_contamination/results_zhao/L3_contamination_*`, `purity_dsuz_*` |
| Shared polymorphism: 33.0–33.8% and 23.6–23.9% of four-fold fixed differences; 6.2–6.6-fold excess; ratio 0.08–0.10 vs πN/πS 0.07–0.08 vs 0.20–0.22 | `07_shared_polymorphism_efficacy`: Zhao 97 → 112 → 113; intervals 138 | `07_.../results_zhao/alpha_shared_polymorphism.tsv`, `.log`; `shared_omega_ci.tsv` |
| Mitochondrial divergence 3.6% vs nuclear 8.4%; 105 of 15,942 fixed differences | Zhao `123_mito_vs_nuclear_divergence.py` | `07_.../results_zhao/mito_vs_nuclear_divergence.log` |
| Window means of π, θW and Tajima's D by chromosome; Fig. 2B | Sato `04_diversity_fst/scripts_sato/01_run_diversity_array.sh` (means and figure: code to be added) | `04_diversity_fst/results_sato/grenedalf_diversity/*.csv` |
| Equal depth: θW 0.0172–0.0175 (*D. suzukii*), 0.0170 (*D. subpulchrella*); 0.0210–0.0234 at full depth | Zhao 01 → 32, 33 → 40, 43 | `04_diversity_fst/results_zhao/L2_depth_standardized_windows.tsv`, `L2_diversity_means.tsv` (config P99_noRM) |
| FST means 0.0063–0.017 (autosomes), 0.010–0.027 (X) | Sato `sync_fst/01`, `02_filter_sync.py`, `04_run_fst_hudson_maf_array.sh` | `04_diversity_fst/results_sato/fst_hudson_MAF01/*/fst-*.csv` |
| Sap branch longest in 74–97% of windows; excess X/A 2.29 (2.12–2.48); drift maximum 1.78; 7,813 windows, 146 and 33 blocks | Zhao `94_fst_structure_and_outliers.py` | `04_diversity_fst/results_zhao/fst_structure.log`, `fst_windows_branches.tsv.gz` |
| Ne 3.41 × 10⁶ (3.23–3.59) and 4.51–4.92 × 10⁶; four-fold 2.69 and 4.03–4.47 × 10⁶; Fig. 2C, Table S5 | Sato `05_effective_population_size/scripts_sato` (02, 03 → 01 → 08) | `05_.../results_sato/effective_population_size_block_bootstrap.tsv` |
| 16–31% lower at short-intron sites, 27–40% at four-fold sites, no overlapping 95% intervals (also Discussion paragraph on robustness) | Zhao 31 (`SI_MIN=40`), 32, 33 (`SI_SET=_si40`, `run_si40_rerun.sh`) → `140_si40_numbers.py` | `05_.../results_zhao/si40_numbers.json`, `Ne_si40.tsv`, `depth_standardized_Ne_si40.tsv` |
| Four-fold diversity ratio by arm (0.86–0.90 on 2R … 0.48–0.49 on the 3R scaffold) | Zhao `75_per_arm_diversity_ratio.py` | `05_.../results_zhao/per_arm_diversity_ratio.log` |
| Ratio below one on every arm after excluding chromosome ends and the centre of chromosome 3 (0.52–0.78) | Zhao 76 → 86 | `05_.../results_zhao/arm_interior_ratios.tsv` |
| Divergence from *D. biarmipes* 3–5% lower in *D. subpulchrella*; divergence-normalized diversity by arm | Zhao `105_pi_over_k.py` (with `depth_common.py`, site tables of 97) | `05_.../results_zhao/B1_pi_over_k.log`, `_arm_ratios.tsv` |
| Polymorphic inversions (X 6.6 Mb, 3L 17.4 Mb, 2R 0.24 Mb) and their frequencies; Table S6 | Zhao 109 (Delly) → 111 → 111b; 110, 110b, 110c (frequency scan and sensitivity) | `08_inversions/results_zhao/INV_*` |
| Sweep-peak densities (1.49 and 2.35 on 3R and the 3R scaffold; 0.52–0.89 and 0.53–0.83 in *D. suzukii*); strongest peaks; no excess sharing; no excess overlap with MK genes; Fig. 2D–F, Fig. SX | Zhao 98 → 99 → 100 → 101 → 101b; 124; 126 | `09_selective_sweeps/results_zhao/F_sweep_*`, `F_sweep_mk_overlap.tsv`, `F_regional_pi.tsv`; `09_selective_sweeps/figures` |
| Highest πN/πS on 3R in *D. subpulchrella* | Zhao `104_pnps_windows.py` | `07_.../results_zhao/A3_pnps_arms.tsv` |
| Unpolarized MK tests: 720 genes; 393, 141, 169, 963, 59, 127, 579; 1,271–1,442 genes per population; 1,982 and 842; 319 (44.3%) reciprocal; Fig. 3A, Fig. S3, Table S8 | `10_mk_tests` (Sato pipeline; `vcf2msa_2026.py`) | `10_mk_tests/results_sato/cmktest_filtered/*.filtered.csv` (set counts: code to be added) |
| Lineage-specific MK: 0.8–0.9% of codons not polarized; shared variants 27.7% and 18.3–18.8%; 1,010 of 9,704 genes and 920–946 of 9,375–9,583; α 0.710 (0.700–0.719) and 0.703–0.714; asymptotic α; Table S9, Supplementary Data 1 | Zhao 60 (`EXCLUDE_SHARED=1`) → 64 (`run_lineage_mk_excl.sh`) | `11_lineage_specific_mk/results_zhao/per_pool/*.qc.tsv`, `*.mk_tests.DAF15.tsv`; `summary/L3_polarized_mk_counts.tsv`, `L3_polarized_mk_alpha.tsv`; `supplementary_data/Sato_Supplementary_Data_1.xlsx` |
| pN/pS 0.066 vs 0.057–0.058; dN/dS 0.229 vs 0.193–0.202; πN/πS 0.0785 (0.0754–0.0818) vs 0.0655–0.0680; Fig. 3XX | Zhao 96, 104, 138 | `11_.../results_zhao/summary/mk_per_site_rates.tsv`; `07_.../results_zhao/efficacy_per_site_ci.tsv`, `A3_pnps_arms.tsv`; `07_.../figures/Fig_shared_omega_efficacy*` |
| 347 genes positive on both lineages; 2.46-fold (P = 7.9 × 10⁻⁷¹); 1.72-fold after power matching (P = 0.0005) | Zhao `64_summarize_polarized_mk.py` | `11_.../results_zhao/summary/L3_polarized_mk_both_lineages.tsv` |
| GO enrichment (DAVID) of the unpolarized sets | DAVID web tool with `10_mk_tests/results_sato/david_mapping` (charts: to be added) | — |
| FlyAtlas2: testis 229 genes, OR 2.46; male midgut 85 genes, OR 1.77; crop and testis in common4 | Zhao `69_enrichment_trace.py` (reproduces the values from Sato's cmktest tables) | `12_enrichment/results_zhao/L1_flyatlas2.tsv` |
| Chemosensory genes: 233 genes, 23% more nonsynonymous fixed differences on the *D. suzukii* lineage (P = 0.005); α 0.76 vs 0.68; ionotropic receptor groups | Zhao 92 → 130, 132 | `11_.../results_zhao/candidates/lineage_candidates.log`, `lineage_category_decomposition.tsv` |
| Testis-maximal genes: power-matched fold 1.56 (FDR 2.1 × 10⁻¹⁶) and 1.20 (FDR 0.0047); OR 1.54 (BH P = 0.022); Dn/Ds ratio 0.98 (0.93–1.03) | Zhao 70; 132 | `11_.../results_zhao/summary/L3_lineage_contrast.tsv`; `candidates/lineage_category_decomposition.log` |
| Candidate genes; Supplementary Data 2 | Zhao 130 → 131 → 133 | `11_.../results_zhao/candidates/lineage_candidates.tsv`; `supplementary_data/Sato_Supplementary_Data_2.xlsx` |
| Zhao and Begun (2017): 601 testable, 93 genes, 2.66-fold; 124/146/133/145 genes, 1.93–2.15-fold; 66-gene set 3.11–5.37-fold; power-matched 1.52–2.08 and 1.75–2.94 (Table S12) | Zhao 68 → 69 | `12_enrichment/results_zhao/L1_zhao_begun.tsv` |
| TE content 46.4% and 50.1% of the assemblies; LTR 18.3% and 20.8% | Sato `13_transposable_elements` 01, 02 | `13_transposable_elements/results_sato/*.tbl`, `*.summary.tsv` |
| PoPoolationTE2 records (56,516; 17,063–73,299), class shares, frequency distributions | Sato 03–12 | `13_transposable_elements/results_sato/*frequency_filter*`, `*pairup*` (figures: code to be added) |
| Same-strand read pairs 12.5% (Kan) and 4.1% (Sap2) | measured from read-pair orientation during QC | `00_inputs/pool_metadata.tsv` (column `same_strand_pairs_pct`) |

## Gene copy number (section of duplication_genes_sato.docx)

| Statement | Script(s) | Result file(s) |
|---|---|---|
| Annotation counts overestimate copy-number divergence (552 annotation-based differences, most not confirmed by loci and read depth) | Zhao 135 (`135_gene_family_diamond.sbatch`) → 136 → 137 | `14_gene_copy_number/results_zhao/C2_family_copy_number.log`, `C2_validation.log`, `C2_annotation_candidates_outcome.tsv` |
| 59,235 proteins aligned with miniprot; 11,277 paralog groups; per-copy SD 0.08, 0.07 and 0.10 | Zhao 136 (`136_copy_number_validation.sbatch`), 137, 139 | `14_gene_copy_number/results_zhao/C2_validation.log`, `C2_solid_changes.log`, `C2_paralog_groups.tsv` |
| 175 high-confidence changes (172 copy-number differences, three pseudogenizations); 63, 41 + 17, 31 + 8 by lineage; 12 not polarized; 145 of 172 involve one copy; Supplementary Data 3 | Zhao 137 → 139 (with 138, 142, 144) | `14_gene_copy_number/results_zhao/C2_solid_copy_number_changes.tsv`, `C2_solid_changes.log`, `C2_copy_number_changes_all_tiers.tsv` |
| Origin of the 135 gains (55 tandem, 36 dispersed, 20 retroposition, 12 mixed, 12 ambiguous or unplaced); assemblies agree for 39 of 41 and 30 of 31 | Zhao 144 (called from 139) | `14_gene_copy_number/results_zhao/C2_solid_copy_number_changes.tsv` (columns `gain_origin`, `gain_origin_detail`), `C2_solid_changes.log` |
| Immunity enriched: 12 of 222 vs 160 of 11,055 groups; Mantel–Haenszel OR 2.77, P = 0.0017, FDR = 0.028 | Zhao 139 | `14_gene_copy_number/results_zhao/C2_solid_categories.tsv`, `C2_solid_changes.log` |
| *Or19*, *Or67a*, *Gr92a*, α-esterase, *yellow-d*, *iav* (77 kb from its parent), *Attacin*, *lectin-37Db*, *Sr-C*; *Or23a* and *Gr59a* intermediate; *Or74a* truncated in both species | Zhao 137, 139, 144 | `14_gene_copy_number/results_zhao/C2_solid_copy_number_changes.tsv`, `C2_copy_number_changes_all_tiers.tsv`, `C2_calls.tsv` |
| Pseudogenization tests; 200 control exons; thresholds 0.80 and 0.90 | Zhao 142 (called from 139) | `14_gene_copy_number/results_zhao/C2_pseudogene_events.tsv`, `C2_frameshift_estimator_calibration.tsv`, `C2_solid_readcheck_loci.tsv` |
| Retrogenes of Chiu et al. (2013): intronless copies of *RpL36*, *RpS14*, *T-cp1* (*CCT1*) and *VhaAC39-1*; none of *Sce* or *βTub97EF* | Zhao 145 (`145_six_gene_retrocopy_origin.py`), 146f | `14_gene_copy_number/results_zhao/retrogenes/six_gene_origin.out`, `cds6_vs_IsoJpt.tsv`, `parent_cds_vs_*.tsv` |
| 83–98% identity to the parent; copies in both WT3 assemblies; RU_1.0 contains only the *VhaAC39-1* copy | Zhao 146f, `145_retrocopy_flanks_junctions.sh` | `14_gene_copy_number/results_zhao/retrogenes/parent_cds_vs_*.tsv`, `rps14_cds_vs_IsoJpt.tsv`, `map_RU.paf`, `map_Dsub.paf` |
| Junction-spanning reads 0.55–1.32 of flanking reads in all four pools; flanks not Y-linked | `145_retrocopy_flanks_junctions.sh`, `145_retrocopy_flank_coverage.sh` | `14_gene_copy_number/results_zhao/retrogenes/junctions.tsv`, `point_cov.tsv`, `points.tsv` |
| Stop codon at codon 77 (*RpL36*), frameshift at codon 107 (*RpS14*), stop codon at codon 443 (*T-cp1*), frequencies 0.94–1.00; *VhaAC39-1* open reading frame 97–99% of the parental protein; same truncations in both WT3 assemblies | Zhao 146 (rules of 142), 146e (146b: *D. suzukii* copies only) | `14_gene_copy_number/results_zhao/retrogenes/retrocopy_disruption_events.tsv`, `retrocopy_orfs_3asm.tsv` |
| No exon–exon junction reads of *Sce* or *βTub97EF* (670–2,597 exonic reads per pool and gene); positive controls 41–290; no junction reads in *D. subpulchrella* (controls 30–59) | Zhao 146c, 146d | `14_gene_copy_number/results_zhao/retrogenes/junction_kmers.tsv`, `junction_kmers_dsub.tsv` |
| Gene copies merged into the locus of a larger gene (205–226 per genome; up to 35 paralog groups differ between the species) | Zhao `145_hidden_nested_copies.py` | `14_gene_copy_number/results_zhao/retrogenes/hidden_nested.tsv` |
| *iav* proteins 99.8% identical over 1,117 amino acids | separate protein alignment (code not included) | — |

## Methods

| Methods statement | Script(s) | Result file(s) |
|---|---|---|
| Depth limits (max(10, P5) to P98; P99 for the per-site analysis), indel and RepeatMasker masks | Sato `03_read_filtering_masks/scripts_sato` (03, 07, 08A; `setup_reanalysis.py`); Zhao 30, 31 | `03_read_filtering_masks/config_sato/depth_thresholds.tsv`, `results_zhao/depth_thresholds.tsv` |
| PROmer block numbers (77,307; 45,151; 26,249), breakpoints with a 200-kb minimum segment (42, 61, 141) | Zhao 120, 91 | `02_synteny/results_zhao/promer_block_audit.log`, `synteny_summary.log` |
| Depth standardization keeps 91–99% (autosomes) and 87–98% (X) of valid sites | Zhao 33 → 140 | `05_.../results_zhao/si40_numbers.json` |
| Short introns 40–65 bp; per-site estimates 2–3% lower than grenedalf | Zhao 31, 140 | `05_.../results_zhao/si40_numbers.log` |
| Sweep null: two-epoch fit within 3.6%; FDR ≤ 0.10 thresholds | Zhao 100, 101 | `09_selective_sweeps/results_zhao/F_sweep_thresholds.tsv`, `F_sweep_summary.log` |
