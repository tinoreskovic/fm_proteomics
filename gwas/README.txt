Scripts for de novo UK Biobank GWAS of fat mass and fat-free mass

These scripts document the REGENIE GWAS workflow used to produce the fat-mass and fat-free-mass GWAS summary statistics that are used as input by 03_proteins_fat_mass_MR_and_plots.R and 04_proteins_fat_free_mass_MR_and_plots.R. 

The cleaned scripts assume these files already exist in controlled_access/gwas/ or the equivalent RAP project folder:

n_fm_ppp_complement_4regenie.phe
Phenotype/covariate file for the fat-mass GWAS sample. This file should contain the non-overlapping UKB complement described in the manuscript: participants with measured fat mass, excluding heterozygosity/missingness outliers, sex chromosome aneuploidy, discordant genetic and registry sex, non-White-British ancestries as defined by Bycroft et al., and the 32,594 participants selected for proteomic profiling. Expected final n = 367,620.

n_ffm_ppp_complement_4regenie.phe
Phenotype/covariate file for the fat-free-mass GWAS sample after the same exclusions, plus exclusion of participants missing fat-free mass. Expected final n = 367,586.

Both files should contain FID and IID columns and the relevant phenotype column:
fat_mass or fat_free_mass.

Both files should also contain the manuscript covariates:
sex, age, age2, age_sex, age2_sex, pc1-pc20, genotype_measurement_batch, centre, and genetic_array.

Execution order

Run from the project root:

1. bash gwas/01_prepare_genotyped_variants_for_step1.sh
2. bash gwas/02_run_regenie_step1.sh fat_mass
3. bash gwas/02_run_regenie_step1.sh fat_free_mass
4. bash gwas/03_prepare_imputed_variants_for_step2.sh fat_mass
5. bash gwas/03_prepare_imputed_variants_for_step2.sh fat_free_mass
6. bash gwas/04_run_regenie_step2.sh fat_mass
7. bash gwas/04_run_regenie_step2.sh fat_free_mass
8. bash gwas/05_merge_regenie_results.sh fat_mass
9. bash gwas/05_merge_regenie_results.sh fat_free_mass

The merge script writes:
controlled_access/fat_mass_gwas_regenie.txt
controlled_access/fat_free_mass_gwas_regenie.txt
