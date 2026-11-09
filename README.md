# Code and results: comparative population genomics of *Drosophila subpulchrella* and *D. suzukii*

Code and result tables for the comparative population genomics of the invasive fruit pest *Drosophila suzukii* and
its close relative *D. subpulchrella*, with *D. biarmipes* as outgroup: genetic diversity, effective population size,
selection, inversions, selective sweeps, gene copy number and transposable elements.

This folder contains the code and the result tables behind the manuscript "Comparative population genomics of
*Drosophila subpulchrella* and suzukii species complex" (version Sato1106, with the tracked changes of
2026-11-06) and in the gene copy-number section (duplication_genes_sato.docx with its tracked changes). It contains
only analyses whose results are reported in these texts, and only the final run of each analysis. Earlier runs, quality-control traces and analyses that are not in the manuscript are not included (see
"Not included").

`docs/manuscript_number_map.md` lists, for each Results paragraph, the script that computes each reported number
and the file that contains it.

## Layout

Folders follow the order of the Methods. Inside each folder:

| Subfolder | Content |
|---|---|
| `scripts_sato` | code written and run by Airi Sato |
| `scripts_zhao` | code run in Li Zhao's project folder |
| `scripts_peng` | `vcf2msa_2026.py` by Junhui Peng |
| `results_*`, `figures`, `supplementary_data` | result tables, logs and figure files produced by these scripts |
| `config_*` | configuration tables read by the scripts |

| Folder | Analysis (manuscript section) |
|---|---|
| `00_inputs` | pool metadata (sample sizes, read-pair orientation), Zhao and Begun (2017) gene sets, FlyAtlas2 table, chunk lists |
| `02_synteny` | PROmer synteny, collinear segments, breakpoints and their lineage (Results "Chromosome synteny"; Methods "Genome synteny analysis") |
| `03_read_filtering_masks` | duplicate removal, region-filtered BAMs, depth limits, indel and RepeatMasker masks, per-site allele counts |
| `04_diversity_fst` | grenedalf π, θW, Tajima's D and FST; depth standardization; four-population FST structure |
| `05_effective_population_size` | Ne from short-intron and four-fold sites; robustness of the Ne difference; ratios between chromosome arms; diversity relative to divergence |
| `06_contamination` | DNA of the other species in each pool (nuclear and mitochondrial) |
| `07_shared_polymorphism_efficacy` | fixed differences polymorphic in the other species; mitochondrial divergence; πN/πS; Fig. 3XX (efficacy of purifying selection) |
| `08_inversions` | polymorphic inversions in *D. subpulchrella* (Delly, allele-frequency scan, validation with synteny); data of Table S6 |
| `09_selective_sweeps` | SweepFinder2 scans, simulated null, peak densities, overlap with MK genes; Fig. 2D–F and Fig. SX |
| `10_mk_tests` | variant calling, population CDS alignments, reciprocal best hits, PRANK, cmktest (unpolarized MK tests) |
| `11_lineage_specific_mk` | polarized MK tests with *D. biarmipes* as outgroup; pN/pS and dN/dS; overlap between lineages; candidate genes; Supplementary Data 1 and 2 |
| `12_enrichment` | FlyAtlas2 tissue enrichment and overlap with Zhao and Begun (2017), with power matching (Table S12) |
| `13_transposable_elements` | RepeatMasker with the DSu-TE library; PoPoolationTE2 insertion frequencies |
| `14_gene_copy_number` | gene copy number from miniprot loci and read depth; high-confidence copy-number changes and pseudogenes; origin of gained copies (tandem, dispersed, retroposition); test of the *D. suzukii* retrogenes reported by Chiu et al. (2013) (`results_zhao/retrogenes`) |
| `environment` | package lists of the conda environments used by `scripts_zhao` |
| `MANIFEST.tsv` | every file: destination, source path on the Rockefeller HPC, size, SHA-256 |

The folder was built by `scripts/143_build_github_folder.py` of the project (rerun it to rebuild the folder after
an update; it copies files and never modifies the sources).

File names in this folder have no date prefix (in the project folder, most result files start with the date of
the run). In the scripts and logs of Li Zhao, references to file names were changed in the same way,
and dates were removed from the comments of the scripts; random seeds (for example 20261021) and release dates of
external data were not changed. Column `edited` of `MANIFEST.tsv` marks the files that differ from their source in
this way; `bytes` and `sha256` describe the files in this folder. Files of Airi Sato and Junhui Peng are unchanged.

## Running the code

- The scripts contain absolute paths of the Rockefeller University HPC (`/ru-auth/local/home/lzhao/...`,
  `/ru-auth/local/home/asato/scratch/...`). The paths are set in variables at the top of each script (for example
  `P`, `PROJ`, `QC`, `RES`, `A`). To rerun an analysis elsewhere, change these variables.
- Shell scripts with `#SBATCH` lines were run with Slurm. Most Python steps were run on many genome chunks at once
  with a job runner; the chunk lists are in `00_inputs/task_lists`.
- Script numbers are those of the project folder. Gaps in the numbering are analyses that are not included.
- Random seeds are set in the scripts (for example 20261021, 20261023, 20261030) and are recorded in their logs.

## Inputs that are not included

Large or public inputs are not copied. They are:

- Reference assemblies and annotations (NCBI RefSeq): *D. suzukii* CBGP_Dsuzu_IsoJpt1.0 (GCF_043229965.1) and
  Dsuz_RU_1.0 (GCF_037355615.1); *D. subpulchrella* RU_Dsub_v1.1 (GCF_014743375.2); *D. biarmipes* RU_DBia_V1.1
  (GCF_025231255.1); *D. melanogaster* Release 6 plus ISO1 MT (GCF_000001215.4).
- Reads: *D. subpulchrella* pool (accession to be added); *D. suzukii* pools from Olazcuaga et al. (2020) and Feng
  et al. (2024); *D. biarmipes* HiFi (SRR20215436) and RNA-seq (SRR20219023–SRR20219026); *D. subpulchrella*
  reference strain (SRR12476582).
- Intermediate files (BAM, VCF, sync, per-site allele counts, per-site tables, SweepFinder2 grids, PoPoolationTE2
  files). They are regenerated by the scripts; their locations on the HPC are in the scripts.
- Public tables: Gene Ontology (go-basic.obo) and the FlyBase GO annotation (fb.gaf.gz) used by
  `11_lineage_specific_mk/scripts_zhao/70_lineage_enrichment.py` and `14_gene_copy_number/scripts_zhao/135_gene_family_copy_number.py`;
  the DSu-TE library of Mérel et al. (2021).
- Gene copy number: DIAMOND tables, miniprot alignments (GFF/PAF, about 180 MB per genome) and read-depth tables of
  scripts 135–137 (regenerated by `135_gene_family_diamond.sbatch`, `136_copy_number_validation.sbatch` and
  `137_loci_depth.sbatch`), and the two assemblies of the *D. suzukii* WT3 strain used for the retrogene test:
  GCA_000472105.1 (Chiu et al. 2013) and GCA_013340165.1 (Paris et al. 2020), downloaded by `146f_retrocopy_blast_wt3.sh`.

## Software

Versions as stated in the Methods or recorded in the logs: BWA-MEM 0.7.15 and 0.7.19 (pools), 0.7.18 (reference
strain); SAMtools 1.9 (mpileup), 1.21; BCFtools 1.20; grenedalf 0.6.3; PRANK (codon mode); cmktest; BLASTP;
RepeatMasker 4.2.2 with Dfam 3.9 and RMBlast; MUMmer 4.0.1 (PROmer); BEDTools 2.31.1; minimap2 2.30; Delly 1.7.3;
seqtk 1.4; SweepFinder2 1.0; DIAMOND 2.1.15; miniprot 0.13-r248; gffread 0.9.9; BLAST+ 2.11.0;
Biopython 1.88; SAMtools 1.24 (copy number); moments 1.6.1; msprime 1.4.4; PoPoolationTE2 1.10.03; DAVID (web); Python 3.9.6
(NumPy 1.21.1, pandas 1.3.1) and 3.11 (NumPy 2.4.6, pandas 2.3.3). Full package lists of the environments used
by `scripts_zhao` are in `environment/`.

## Not included

These analyses were run but are not in the manuscript, or were replaced by a later run:

- Earlier runs of the diversity, FST, Ne and MK analyses (Airi Sato's 2026-08 and 2026-09 runs; direct-BAM FST;
  Ne with introns of 30–65 bp; the first cmktest run for the *D. suzukii* populations).
- Quality-control traces that recomputed manuscript numbers from Airi Sato's files (project scripts 10, 11, 41, 50,
  93), an independent MK recount used only as a library here (script 20 is included because script 60 imports it).
- The distribution of fitness effects (script 102 is included only because script 113 uses its functions and
  output tables), codon-usage and GC-biased gene conversion tests, the window regression of πN/πS, gene-family
  MK statistics, the *D. suzukii* inversion tests, an independent TE analysis with EDTA and TEMP2,
  competitive mapping (scripts 62 and 63 are included only because script 78 uses the read subsample made by 63),
  the outputs of the locus-level read check of script 138 (replaced by script 142; 138 is included because
  scripts 139 and 142 use its functions), and the scripts that wrote numbers or tables into Word documents (for
  example 71, 89, 116, 117, 125, 134, 147).

## Still to be added

These parts of the manuscript are not yet represented by code in this folder:

- Genome assemblies: hifiasm (*D. biarmipes*), Canu (Dsuz_RU_1.0) and BUSCO commands.
- Read trimming (Trimmomatic) and the initial BWA-MEM mapping of the pools; the BCFtools calling of the
  *D. subpulchrella* VCF used by the MK tests (the *D. suzukii* calling script is included).
- The CDS coordinate table used by `vcf2msa_2026.py` for *D. subpulchrella* (`Dsub.cds_map.csv`, Junhui Peng).
- Figure code for Fig. 1 (synteny, R), Fig. 2A–C, Fig. 3A, Figs. S1–S3, X1–X3 and Sx, and the code that computes the
  window means of π, θW and Tajima's D (paragraph "Genome-wide patterns of nucleotide diversity") and the MK set
  overlaps (Fig. 3A, Fig. S3, Table S8) from the included result tables.
- DAVID GO output charts (Results "Functional enrichment").

## Notes

- `05_effective_population_size/results_zhao/genes_neutral_pi.tsv` contains short-intron columns from the
  30–65-bp site set; only its four-fold columns are used (interior ratio 0.52–0.78).
- `04_diversity_fst/scripts_zhao/40_summarize_popgen.py` and `43_summarize_depth_standardized.py` read the 30–65-bp
  run by default. Their all-site window means, which the manuscript uses, are identical in the 40–65-bp run
  (checked); Ne from neutral sites comes from `140_si40_numbers.py`.
- Supplementary Data 1 was written by project script 116 from the tables in
  `11_lineage_specific_mk/results_zhao/per_pool`. Script 116 is not included because it also formats numbers
  for the Word documents.
- `14_gene_copy_number`: script 137 joins an alignment to a locus when the alignment spans overlap, without
  comparing exons, so a gene copy inside an intron of a larger gene on the same strand is merged into that gene's
  locus (for example the *VhaAC39-1* retrocopy of *D. suzukii*). `145_hidden_nested_copies.py` estimates how many
  copies this affects (`results_zhao/retrogenes/hidden_nested.tsv`); the high-confidence table has not been
  recomputed with an exon-overlap condition.
- The identity of the two *D. subpulchrella* iav proteins (99.8% over 1,117 amino acids) in the copy-number text
  was computed with a separate protein alignment whose code is not included.
