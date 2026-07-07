Code accompanying "Analyses of protein levels in relation to fat mass and fat-free mass in the UK Biobank implicate several proteins as potential drug targets" by Oreskovic et al (2026).

These scripts are intended to document the analysis workflow.

1. Folder layout and execution order

project_root/
  input/                 lookup/index files described below
  controlled_access/     user-supplied UK Biobank and GWAS inputs/outputs; not shared here
  summary_stats/         user-downloaded pQTL summary statistics; not shared here
  ld_reference/          user-supplied 1000 Genomes EUR PLINK reference; not shared here
  gwas/                  scripts for de novo fat-mass and fat-free-mass GWAS
  selection_output/      created by 01_stability_selection.R
  positive_control_output/ created by 04b_positive_control_GLP1R_MR.R
  ShareProColoc/         created by the colocalisation scripts
  plots/                 figures produced by the analysis scripts
  *.R / *.py / *.ipynb   scripts listed below

Run the scripts from project_root/ in filename order.

The gwas/ subfolder contains the cleaned REGENIE scripts used to generate controlled_access/fat_mass_gwas_regenie.txt and controlled_access/fat_free_mass_gwas_regenie.txt. These GWAS outputs are not included in this repository.

01_stability_selection.R
Loads protein and fat-mass data, merges Olink plate metadata, imputes missing NPX values, inverse-rank-normalises them, and runs randomised LASSO stability selection with cross-validated hyperparameters. Outputs selected proteins, cross-validation summaries, and the stability-selection plot.

02_synapse_download_summary_stats.ipynb
Retrieves EUR-ancestry UKB-PPP pQTL summary-statistic files from Synapse for the selected proteins and stores them under summary_stats/pqtls_as_exposures/ for the MR scripts.

03_proteins_fat_mass_MR_and_plots.R
Performs cis MR of selected proteins on fat mass using IVW estimates for proteins with multiple genome-wide-significant cis-pQTLs, Wald-ratio estimates for proteins with one such cis-pQTL, and PCA-GMM only as a fallback where no primary cis instrument is available. MR-Egger slope and the MR-Egger intercept estimates are saved for pleiotropy assessment. 

04_proteins_fat_free_mass_MR_and_plots.R
Runs the analogous MR workflow for fat-free mass, with the same IVW/Wald/PCA-GMM fallback structure and pleiotropy checks.

04b_positive_control_GLP1R_MR.R
Runs the GLP1R positive-control PCA-GMM analysis for fat mass and fat-free mass. Place the GLP1R pQTL archive under summary_stats/pqtls_as_exposures_positive_control/ or summary_stats/pqtls_as_exposures/. Outputs are written to positive_control_output/.

05_coloc_preprocessing.R
For proteins with MR evidence for fat mass, prepares SharePro inputs: matched protein and fat-mass summary-statistic files and LD matrices.

05_colocalisation_sharepro.py
Runs SharePro for each prepared protein/fat-mass pair. Set SHAREPRO_SCRIPT if the SharePro entry-point script is not available at SharePro_coloc/src/SharePro/sharepro_coloc.py, and set SHAREPRO_PYTHON if it should be run with a specific Python executable.

05_postprocessing_coloc.R
Reads SharePro outputs and summarises proteins with posterior probability of colocalisation with fat mass greater than 0.8.

06_coloc_plots.R
Generates two-panel locus plots for proteins with colocalisation evidence.

07_drug_db.R
Queries the OpenTargets API for tractability and existing-drug information for MR-supported proteins.

07_go_kegg_enrichment.R
Runs GO biological-process and KEGG enrichment analyses for MR-supported proteins using the proteins screened in the randomised LASSO as the background universe. It also retrieves UniProt FUNCTION text.

2. Files in input/

all_synapse_ukb_ppp_filenames.csv
Index of Synapse IDs and filenames for UKB-PPP protein GWAS archives, used by script 02.

coding143.tsv
UK Biobank coding 143 table mapping Olink assay IDs to protein names, used for labels and output tables.

olink_assay.txt
Olink protein-panel categorisation, used by scripts 01, 03, and 04.

olink_batch_number.txt
Plate ID to batch mapping, used by script 01 — not provided here but available in the .dat format from the UKB Showcase.

olink_protein_map_3k_v1.tsv
Lookup table containing assay IDs, gene coordinates, Ensembl IDs, UniProt IDs, and related protein metadata, used by scripts 03-07.

3. Controlled and downloaded files not shared

UK Biobank individual-level data, UK Biobank GWAS outputs, 1000 Genomes PLINK reference files, and the full UKB-PPP/Sun et al. pQTL summary statistics are not bundled here because access is controlled by their respective platforms. The scripts refer to these with relative paths only: controlled_access/, summary_stats/, and ld_reference/.

Researchers can apply for access through:
a) UK Biobank Research Analysis Platform (RAP), for individual-level genetic, proteomic, and body-composition data.
b) Synapse / UKB-PPP, for pQTL summary statistics, also browsable at http://ukb-ppp.gwas.eu.

Please see the manuscript for details.
