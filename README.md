# Code and results: comparative population genomics of *Drosophila subpulchrella* and *D. suzukii*

Code and result tables for the comparative population genomics of the invasive fruit pest *Drosophila suzukii* and
its close relative *D. subpulchrella*, with *D. biarmipes* as outgroup: genetic diversity, effective population size,
selection, inversions, selective sweeps, gene copy number and transposable elements.

This folder contains the code and the result tables behind the manuscript "Comparative population genomics of
*Drosophila subpulchrella* and suzukii species complex"

`docs/manuscript_number_map.md` lists, for each Results paragraph, the script that computes each reported number
and the file that contains it.

## Layout

Folders follow the order of the Methods. Inside each folder:

| Subfolder | Content |
|---|---|
| `scripts_sato` | code written and run by Airi Sato |
| `scripts_zhao` | code run in Li Zhao's project folder |
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


## Running the code

- The scripts contain absolute paths of the Rockefeller University HPC (`/ru-auth/local/home/lzhao/...`,
  `/ru-auth/local/home/asato/scratch/...`). The paths are set in variables at the top of each script (for example
  `P`, `PROJ`, `QC`, `RES`, `A`). To rerun an analysis elsewhere, change these variables.
- Shell scripts with `#SBATCH` lines were run with Slurm. Most Python steps were run on many genome chunks at once
  with a job runner; the chunk lists are in `00_inputs/task_lists`.
- Script numbers are those of the project folder. Gaps in the numbering are analyses that are not included.
- Random seeds are set in the scripts (for example 20261021, 20261023, 20261030) and are recorded in their logs.

## Software

Versions as stated in the Methods or recorded in the logs: BWA-MEM 0.7.15 and 0.7.19 (pools), 0.7.18 (reference
strain); SAMtools 1.9 (mpileup), 1.21; BCFtools 1.20; grenedalf 0.6.3; PRANK (codon mode); cmktest; BLASTP;
RepeatMasker 4.2.2 with Dfam 3.9 and RMBlast; MUMmer 4.0.1 (PROmer); BEDTools 2.31.1; minimap2 2.30; Delly 1.7.3;
seqtk 1.4; SweepFinder2 1.0; DIAMOND 2.1.15; miniprot 0.13-r248; gffread 0.9.9; BLAST+ 2.11.0;
Biopython 1.88; SAMtools 1.24 (copy number); moments 1.6.1; msprime 1.4.4; PoPoolationTE2 1.10.03; DAVID (web); Python 3.9.6
(NumPy 1.21.1, pandas 1.3.1) and 3.11 (NumPy 2.4.6, pandas 2.3.3). 
