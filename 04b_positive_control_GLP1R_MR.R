library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(MendelianRandomization)

options(scipen = 999)

input_dir <- "input"
controlled_access_dir <- "controlled_access"
summary_stats_dir <- "summary_stats"
output_dir <- "positive_control_output"
dir.create(output_dir, showWarnings = FALSE)

positive_control_id <- "GLP1R"
exposure_n <- 34557
pc_thresh <- 0.99
bfile <- file.path("ld_reference", "1kg_eur_hg38")
rsid_pos_map_dir_path <- file.path(summary_stats_dir, "rsid_pos_map")
extract_dir_path <- file.path(summary_stats_dir, "extracted_pqtls_positive_control_glp1r")
pqtl_tar_dirs <- c(
  file.path(summary_stats_dir, "pqtls_as_exposures_positive_control"),
  file.path(summary_stats_dir, "pqtls_as_exposures")
)

read_cis_pqtl_gz_from_tar <- function(tar_dirs, chr_value, id, extract_dir) {
  tar_files <- unlist(lapply(tar_dirs, function(tar_dir) {
    if (!dir.exists(tar_dir)) {
      return(character())
    }
    list.files(tar_dir, pattern = paste0(id, ".*\\.tar$"), full.names = TRUE)
  }))
  
  if (length(tar_files) == 0) {
    stop("No GLP1R pQTL .tar file found. Place it in summary_stats/pqtls_as_exposures_positive_control/ or summary_stats/pqtls_as_exposures/.")
  }
  
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE)
  }
  
  untar(tar_files[1], exdir = extract_dir)
  
  gz_file_pattern <- paste0(".*", chr_value, ".*", id, ".*\\.gz$")
  gz_file <- list.files(extract_dir, pattern = gz_file_pattern, recursive = TRUE, full.names = TRUE)
  
  if (length(gz_file) == 0) {
    stop("No GLP1R pQTL .gz file found after extracting the .tar archive.")
  }
  
  data.frame(fread(gz_file[1], sep = " "))
}

read_rsid_map_gz_from_directory <- function(rsid_pos_map_dir, chr_value) {
  gz_file_pattern <- paste0(".*", chr_value, ".*\\.gz$")
  gz_file <- list.files(rsid_pos_map_dir, pattern = gz_file_pattern, full.names = TRUE)
  
  if (length(gz_file) == 0) {
    stop("No matching rsid position-map .gz file found.")
  }
  
  data.frame(fread(gz_file[1], sep = "\t"))
}

load_outcome_gwas <- function(file_name) {
  outcome <- fread(file.path(controlled_access_dir, file_name))
  
  filtered_variant_file <- file.path(controlled_access_dir, "filtered_imputation_variants.csv")
  if (file.exists(filtered_variant_file)) {
    imputation_variants <- fread(filtered_variant_file, header = FALSE)
    outcome <- subset(outcome, ID %in% unique(imputation_variants$V3))
    rm(imputation_variants)
  } else {
    warning("controlled_access/filtered_imputation_variants.csv not found; using all variants in ", file_name)
  }
  
  data.frame(outcome)
}

prepare_glp1r_exposure <- function() {
  pqtl_map <- read.delim(file.path(input_dir, "olink_protein_map_3k_v1.tsv"))
  pqtl_map <- subset(pqtl_map, Assay == positive_control_id)
  
  chr_value <- paste0("chr", pqtl_map$chr[1])
  cis_pqtls <- read_cis_pqtl_gz_from_tar(pqtl_tar_dirs, chr_value, positive_control_id, extract_dir_path)
  rsid_pos <- read_rsid_map_gz_from_directory(rsid_pos_map_dir_path, chr_value)
  
  cis_pqtls <- merge(
    cis_pqtls,
    rsid_pos,
    by.x = c("ID", "ALLELE0", "ALLELE1"),
    by.y = c("ID", "REF", "ALT"),
    all.x = TRUE
  )
  
  cis_pqtls$chr <- as.numeric(cis_pqtls$CHROM)
  cis_pqtls$position <- as.numeric(cis_pqtls$GENPOS)
  cis_pqtls$beta <- cis_pqtls$BETA
  cis_pqtls$se <- cis_pqtls$SE
  cis_pqtls$id.exposure <- positive_control_id
  cis_pqtls$SNP <- cis_pqtls$rsid
  cis_pqtls$pval.exposure <- 10^(-cis_pqtls$LOG10P)
  
  start_pos <- pqtl_map$gene_start[1]
  end_pos <- pqtl_map$gene_end[1]
  cis_pqtls <- cis_pqtls[
    cis_pqtls$position >= start_pos - 1000000 &
      cis_pqtls$position <= end_pos + 1000000,
  ]
  
  cis_pqtls <- cis_pqtls %>%
    filter(!is.na(rsid), !is.na(pval.exposure)) %>%
    arrange(pval.exposure) %>%
    distinct(rsid, .keep_all = TRUE)
  
  keep <- ieugwasr::ld_clump_local(
    dplyr::tibble(rsid = cis_pqtls$rsid, pval = cis_pqtls$pval.exposure, id = cis_pqtls$id.exposure),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = bfile,
    clump_kb = 10000,
    clump_r2 = 0.95,
    clump_p = 1
  )
  
  if (nrow(keep) <= 1) {
    stop("Not enough GLP1R cis variants retained for PCA-GMM.")
  }
  
  cis_pqtls <- subset(cis_pqtls, rsid %in% unique(keep$rsid))
  cis_pqtls$exposure <- positive_control_id
  cis_pqtls$beta.exposure <- cis_pqtls$beta
  cis_pqtls$se.exposure <- cis_pqtls$se
  cis_pqtls$effect_allele.exposure <- cis_pqtls$ALLELE1
  cis_pqtls$other_allele.exposure <- cis_pqtls$ALLELE0
  cis_pqtls$eaf.exposure <- cis_pqtls$A1FREQ
  cis_pqtls$pos <- cis_pqtls$POS19
  
  list(
    cis_pqtls = cis_pqtls,
    diagnostics = data.frame(
      exposure = positive_control_id,
      min_cis_p = min(cis_pqtls$pval.exposure, na.rm = TRUE),
      any_genome_wide_cis_pqtl = any(cis_pqtls$pval.exposure < 5e-08, na.rm = TRUE),
      cis_variants_after_clumping = nrow(cis_pqtls)
    )
  )
}

run_glp1r_pca_gmm <- function(outcome_gwas, outcome_name, outcome_label, outcome_n) {
  exposure <- prepare_glp1r_exposure()
  cis_pqtls <- exposure$cis_pqtls
  
  chr_value <- unique(cis_pqtls$CHROM)[1]
  if (chr_value == "X") {
    outcome <- subset(outcome_gwas, CHROM == 23)
  } else {
    outcome <- subset(outcome_gwas, CHROM == as.numeric(chr_value))
  }
  
  outcome$P <- 10^(-outcome$LOG10P)
  outcome$pval.outcome <- outcome$P
  outcome$id.outcome <- outcome_label
  outcome$pos <- outcome$GENPOS
  
  outcome <- format_data(
    outcome,
    type = "outcome",
    snps = cis_pqtls$SNP,
    header = TRUE,
    snp_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "A1FREQ",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    pval_col = "pval.outcome",
    min_pval = 1e-200,
    chr_col = "CHROM",
    pos_col = "pos",
    log_pval = FALSE
  )
  
  dat <- harmonise_data(cis_pqtls, outcome)
  dat <- subset(dat, mr_keep == TRUE)
  dat <- dat[!duplicated(dat$SNP), ]
  
  x_corr <- ieugwasr::ld_matrix(
    variants = unique(dat$SNP),
    with_alleles = TRUE,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = bfile
  )
  
  snpnames <- sapply(strsplit(rownames(x_corr), split = "_"), `[`, 1)
  x_corr <- x_corr[order(snpnames), order(snpnames)]
  dat <- dat[order(dat$SNP), ]
  
  ld_matrix <- harmonise_ld_dat(dat, x_corr)
  x_corr <- ld_matrix[["ld"]]
  found_snps <- ld_matrix[["x"]][["SNP"]]
  dat <- dat[match(found_snps, dat$SNP), ]
  
  mr_input <- MendelianRandomization::mr_input(
    bx = dat$beta.exposure,
    bxse = dat$se.exposure,
    by = dat$beta.outcome,
    byse = dat$se.outcome,
    correlation = matrix(),
    exposure = positive_control_id,
    outcome = outcome_label,
    snps = found_snps,
    effect_allele = dat$effect_allele.exposure,
    other_allele = dat$other_allele.exposure,
    eaf = dat$eaf.exposure
  )
  mr_input@correlation <- x_corr
  
  pca_gmm <- MendelianRandomization::mr_pcgmm(
    mr_input,
    nx = exposure_n,
    ny = outcome_n,
    thres = pc_thresh,
    robust = TRUE,
    alpha = 0.05
  )
  
  data.frame(
    exposure = positive_control_id,
    outcome = outcome_name,
    method = "PCA-GMM",
    PCs = pca_gmm@PCs,
    nsnp = nrow(pca_gmm@Correlation),
    b = pca_gmm@Estimate,
    se = pca_gmm@StdError,
    pval = pca_gmm@Pvalue,
    ci_lwr = pca_gmm@CILower,
    ci_upr = pca_gmm@CIUpper,
    fstat = pca_gmm@Fstat,
    hetr_stat = pca_gmm@Heter.Stat,
    min_cis_p = exposure$diagnostics$min_cis_p,
    any_genome_wide_cis_pqtl = exposure$diagnostics$any_genome_wide_cis_pqtl,
    cis_variants_after_clumping = exposure$diagnostics$cis_variants_after_clumping,
    stringsAsFactors = FALSE
  )
}

fat_mass_gwas <- load_outcome_gwas("fat_mass_gwas_regenie.txt")
fat_free_mass_gwas <- load_outcome_gwas("fat_free_mass_gwas_regenie.txt")

positive_control_results <- rbind(
  run_glp1r_pca_gmm(
    outcome_gwas = fat_mass_gwas,
    outcome_name = "fat_mass",
    outcome_label = "Whole body fat mass",
    outcome_n = 367620
  ),
  run_glp1r_pca_gmm(
    outcome_gwas = fat_free_mass_gwas,
    outcome_name = "fat_free_mass",
    outcome_label = "Whole body fat-free mass",
    outcome_n = 367586
  )
)

positive_control_results$`95% CI` <- paste0(
  round(positive_control_results$ci_lwr, 2),
  " to ",
  round(positive_control_results$ci_upr, 2)
)

write.csv(
  positive_control_results,
  file.path(output_dir, "glp1r_positive_control_pca_gmm_results.csv"),
  row.names = FALSE
)
