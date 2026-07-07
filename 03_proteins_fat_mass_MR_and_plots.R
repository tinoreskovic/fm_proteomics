library(tidyverse)
library(broom)
library(dplyr)
library(readxl)
library(TwoSampleMR)
library(ggplot2)
library(stringr)
library(R.utils)
library(remotes)
library(ckbplotr)
library(ieugwasr)
library(flextable)
library(grid)
library(gridGraphics)
library(ggtext)
library(ggrepel)
library(curl)
library(readr)
library(corrplot)
library(irlba)
library(glmnetUtils)
library(Rmpfr)
library(data.table)
library(R.utils)
library(susieR)
library(coloc)
library(Rfast)
library(MendelianRandomization)
library(dplyr)
library(genetics.binaRies)
library(ieugwasr)
#remotes::install_github("RfastOfficial/Rfast")
library(Rfast)

options(scipen = 999)

input_dir <- "input"
controlled_access_dir <- "controlled_access"
summary_stats_dir <- "summary_stats"

#settings
restarted = T #T F
run = T #T F
proteins = "selected"
clump_or_pc = "lax" #lax, strict

tar_dir_path <- file.path(summary_stats_dir, "pqtls_as_exposures")
extract_dir_path <- file.path(summary_stats_dir, "extracted_pqtls_protein_to_fat_mass")

filtered_selection_df <- read.csv(file.path("selection_output", "selected_proteins.csv"))

filtered_selection_df$feature <- toupper(filtered_selection_df$feature)
pqtl_map <- read.delim(file.path(input_dir, "olink_protein_map_3k_v1.tsv"))

pqtl_map <- subset(pqtl_map, pqtl_map$Assay %in% filtered_selection_df$feature)

ids <- unique(pqtl_map$Assay)

output_folder <- ""
prefix <- ""

rsid_pos_map_dir_path <- file.path(summary_stats_dir, "rsid_pos_map")

if(restarted){
  fat_mass_regenie <- fread(file.path(controlled_access_dir, "fat_mass_gwas_regenie.txt"))
  
  gc()
  
  #only if sub-setting by info score for the first time
  if(F){
    library(data.table)
    
    # Function to read and filter MFI files
    read_and_filter_mfi <- function(chr) {
      file_path <- file.path(controlled_access_dir, "imputation_info", paste0("imputation_info_chr", chr, ".txt"))
      filtered_data <- fread(file_path, header=FALSE)[V8 > 0.7]
      return(filtered_data)
    }
    

    filtered_mfi_list <- list()
    
    for (chr in c(1:22, "X")) {
      filtered_mfi_list[[as.character(chr)]] <- read_and_filter_mfi(chr)
    }
    
    # Combine all filtered data tables
    ukb_mfi <- rbindlist(filtered_mfi_list)
    
    write.csv(ukb_mfi, file.path(controlled_access_dir, "filtered_imputation_variants.csv"))
    
  }
  
  ukb_mfi <- fread(file.path(controlled_access_dir, "filtered_imputation_variants.csv"), header=FALSE)
  info_07 <- unique(ukb_mfi$V3)
  remove(ukb_mfi)
  fat_mass_regenie <- subset(fat_mass_regenie, fat_mass_regenie$ID %in% info_07)
  remove(info_07)
  gc()
  
}

if(clump_or_pc=="lax"){
    clump_r2 = 0.4
    clump_r2_name = "4"
    pc_thresh = 0.99      
    pc_thresh_name = "99"  
    }else{
      clump_r2 = 0.01
      clump_r2_name = "01"
      pc_thresh = 0.99      
      pc_thresh_name = "99"
    }



bfile <- file.path("ld_reference", "1kg_eur_hg38")

mr_results_file <- paste0(output_folder, prefix, "fat_mass_mr_results_", clump_r2_name, "_", pc_thresh_name, ".csv")
mr_results_by_clump_file <- paste0(output_folder, prefix, "fat_mass_mr_results_", clump_r2_name, ".csv")
pca_gmm_results_file <- paste0(output_folder, prefix, "fat_mass_pca_gmm_results_", pc_thresh_name, ".csv")
fm_layout_file <- paste0(output_folder, prefix, "fat_mass_forest_plot_layout.rds")
lax_mr_results_file <- paste0(output_folder, prefix, "fat_mass_mr_results_4_", pc_thresh_name, ".csv")

make_forest_axis_limits <- function(plot_data, step = 0.5) {
  ci_limits <- range(c(plot_data$b - 1.96 * plot_data$se,
                      plot_data$b + 1.96 * plot_data$se),
                    na.rm = TRUE)
  c(floor(ci_limits[1] / step) * step, ceiling(ci_limits[2] / step) * step)
}

read_cis_pqtl_gz_from_tar <- function(tar_dir, chr_value, id, extract_dir) {
  # Ensure the full path without tilde (~)
  tar_dir <- normalizePath(tar_dir)
  extract_dir <- normalizePath(extract_dir)
  
  # list all .tar files in the directory
  tar_files <- list.files(tar_dir, pattern = paste0(id, ".*\\.tar$"), full.names = TRUE)
  
  if (length(tar_files) == 0) {
    stop("No matching .tar file found.")
  }
  
  tar_file <- tar_files[1]  # assuming there's only one matching .tar file
  
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir)
  }
  

  untar(tar_file, exdir = extract_dir)
  
  gz_file_pattern <- paste0(".*", chr_value, ".*", id, ".*\\.gz$")
  gz_file <- list.files(extract_dir, pattern = gz_file_pattern, recursive = TRUE, full.names = TRUE)
  
  if (length(gz_file) == 0) {
    stop("No matching .gz file found.")
  }
  
  data <- fread(gz_file[1], sep=" ")
  data <- data.frame(data)
  
  return(data)
}

#reading rsid=position map function
read_rsid_map_gz_from_directory <- function(rsid_pos_map_dir, chr_value) {

  rsid_pos_map_dir <- normalizePath(rsid_pos_map_dir)
  
  gz_file_pattern <- paste0(".*", chr_value, ".*\\.gz$")
  gz_file <- list.files(rsid_pos_map_dir, pattern = gz_file_pattern, full.names = TRUE)
  
  if (length(gz_file) == 0) {
    stop("No matching .gz file found.")
  }
  
  #data <- read.csv(gz_file[1], sep="")
  data <- fread(gz_file[1], sep="\t")
  data <- data.frame(data)
  
  return(data)
}

#IVW/Wald ratio MR protein to FM run 
if(run==T){
  
  gc()
  mro_fm <- data.frame()
  no_cis_pqtls <- c()
  two_or_fewer_pqtls <- c()
  rsqrs_df <- data.frame()
  
  for(id in ids){try({
    iteration_start_time <- Sys.time()
    message("Starting processing for id: ", id, " at ", iteration_start_time)
    
    chr_value <- pqtl_map$chr[pqtl_map$Assay == id]
    if(length(chr_value)>1){
      chr_value <- chr_value[1]
    }
    chr <- chr_value
    chr_value <- paste0("chr", chr_value)
    
    tar_dir <- tar_dir_path 
  
    extract_dir <- extract_dir_path
    
    cis_pqtls <- read_cis_pqtl_gz_from_tar(tar_dir, chr_value, id, extract_dir)
    
    rsid_pos_map_dir <- rsid_pos_map_dir_path
    
    rsid_pos <- read_rsid_map_gz_from_directory(rsid_pos_map_dir, chr_value)
    
    cis_pqtls <- merge(cis_pqtls, rsid_pos, by.x=c("ID", "ALLELE0", "ALLELE1"),
                       by.y=c("ID", "REF", "ALT"), all.x=TRUE)
    
    cis_pqtls$chr <- cis_pqtls$CHROM
    cis_pqtls$chr <- as.numeric(cis_pqtls$chr)
    cis_pqtls$position <- as.numeric(cis_pqtls$GENPOS)
    cis_pqtls$beta <- cis_pqtls$BETA
    cis_pqtls$se <- cis_pqtls$SE
    cis_pqtls$id.exposure <- id
    cis_pqtls$SNP <- cis_pqtls$rsid
    
    cis_pqtls <- 
      cis_pqtls %>% mutate(pval.exposure = (10^(-LOG10P)))
    
    start_pos <- pqtl_map$gene_start[pqtl_map$Assay == id]
    end_pos <- pqtl_map$gene_end[pqtl_map$Assay == id]
    
    if(clump_or_pc=="strict"){
      cis_pqtls <- cis_pqtls[cis_pqtls$position >= start_pos-250000 & 
                               cis_pqtls$position <= end_pos+250000, ] 
    }
    if(clump_or_pc=="lax"){
      cis_pqtls <- cis_pqtls[cis_pqtls$position >= start_pos-1000000 & 
                               cis_pqtls$position <= end_pos+1000000, ] 
    }

    
    clump_kb <- 10000
    clump_r2 <- clump_r2
    clump_p <- 5e-08
    #remotes::install_github("MRCIEU/genetics.binaRies")
    library(genetics.binaRies)
    
    keep <- data.frame()
    tryCatch({
      keep <- ieugwasr::ld_clump_local(
        dplyr::tibble(rsid = cis_pqtls$rsid, pval = cis_pqtls$pval, id = cis_pqtls$id.exposure),
        plink_bin = genetics.binaRies::get_plink_binary(),
        bfile = bfile,
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = clump_p
      )
    }, error = function(e) {
      message("Error during clumping")
      keep <- data.frame()
    })
    
    if(nrow(keep) == 0){
      min_pval_snp <- cis_pqtls[which.min(cis_pqtls$pval.exposure), ]
      if(min_pval_snp$pval.exposure < 5e-08){
        keep <- data.frame(rsid = min_pval_snp$rsid)
      }else{
        print("no cis-pqtls for")
        print(id)
        no_cis_pqtls <- c(no_cis_pqtls, id)
        next
      }
    }
    
    cis_pqtls <- subset(cis_pqtls, cis_pqtls$rsid %in% unique(keep$rsid))
    
    pqtl_map_bit <- pqtl_map[pqtl_map$Assay == id,]
    pqtl_map_bit <- pqtl_map_bit[, c("Assay", "ensembl_id")]
    
    cis_pqtls <- merge(cis_pqtls, pqtl_map_bit, by.x=c("id.exposure"),
                       by.y=c("Assay"), all.x=TRUE)
    
    cis_pqtls$exposure <- cis_pqtls$id.exposure
    cis_pqtls$beta.exposure <- cis_pqtls$beta
    cis_pqtls$se.exposure <- cis_pqtls$se
    cis_pqtls$effect_allele.exposure <- cis_pqtls$ALLELE1
    cis_pqtls$other_allele.exposure <- cis_pqtls$ALLELE0
    cis_pqtls$eaf.exposure <- cis_pqtls$A1FREQ
    
	    cis_pqtls <- cis_pqtls %>%
	      mutate(F_statistic = (beta.exposure^2) / (se.exposure^2))
	    # F-statistics are retained for diagnostics/reporting only; instruments
	    # are not filtered by F-statistic in the revised primary analyses.
    
    gc()
    
    if(chr=="X"){
      outcome <- subset(fat_mass_regenie, fat_mass_regenie$CHROM==23)
    }else{
      outcome <- subset(fat_mass_regenie, fat_mass_regenie$CHROM==as.numeric(chr))
    }
    
    outcome <- 
      outcome %>% mutate(P = (10^(-LOG10P)))
    
    outcome$pval.outcome <- outcome$P
    outcome$id.outcome <- "Whole body fat mass"
    
    outcome$chr_name <- outcome$CHROM
    outcome$chrom_start <- outcome$GENPOS
    outcome$chr <- ifelse(outcome$chr_name==23,
                          "chrX", paste0("chr", outcome$chr_name))

    gc()
  
    outcome <- data.frame(outcome)
    
    cis_pqtls$pos <- cis_pqtls$POS19
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
      log_pval = FALSE)
    
    dat <- harmonise_data(cis_pqtls, outcome)
    dat <- subset(dat, dat$mr_keep==TRUE)
    dat <- dat[!duplicated(dat$SNP),]
    
    dat2 <- dat
    
    dat2$samplesize.exposure <- 34557
    dat2$samplesize.outcome <- 367620 
    
    dat2$r_exp <- get_r_from_bsen(
      dat2$beta.exposure,
      dat2$se.exposure,
      dat2$samplesize.exposure)
    
    dat2$r_exp2 <- get_r_from_pn(
      dat2$pval.exposure,
      dat2$samplesize.exposure
    )
    
    dat2$r_out <- get_r_from_bsen(
      dat2$beta.outcome,
      dat2$se.outcome,
      dat2$samplesize.outcome)
    
    dat2$r_out2 <- get_r_from_pn(
      dat2$pval.outcome,
      dat2$samplesize.outcome
    )
 
    direction <- mr_steiger2(dat2$r_exp, dat2$r_out, dat2$samplesize.exposure,
                             dat2$samplesize.outcome)
  
    rsqrs <- data.frame(exposure = id, rsq.exposure = direction$r2_exp, rsq.outcome = direction$r2_out, direction = direction$steiger_test)
    
    rsqrs_df <- rbind(rsqrs_df, rsqrs)
    
    rsqrs_df <- rsqrs_df[!duplicated(rsqrs_df$exposure), ]
  
	    write.csv(rsqrs_df, paste0(output_folder, prefix, "fat_mass_steiger_directionality_", clump_or_pc,  ".csv"), row.names = FALSE)
    
    if(length(dat$SNP) == 1){
      two_or_fewer_pqtls <- c(two_or_fewer_pqtls, id)
      mro <- mr_wald_ratio(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome)
      mro <- data.frame(mro)
      mro$id.exposure <- id
      mro$exposure <- id
      mro$id.outcome <- "Whole body fat mass"
      mro$outcome <- "outcome"
      mro$method <- "Wald ratio"
      
      mro$direction <- direction$correct_causal_direction
      mro$direction_p <- direction$steiger_test
      
      mro_fm <- rbind(mro_fm, mro)
      mro_fm <- subset(mro_fm, !is.na(mro_fm$b))
      
	      write.csv(mro_fm, mr_results_by_clump_file)
      
    } else if(length(dat$SNP) > 1){
      
      if(length(dat$SNP)==2){
        two_or_fewer_pqtls <- c(two_or_fewer_pqtls, id)
      }
      
      mro <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
      
      tryCatch({
        intercept <- data.frame(mr_pleiotropy_test(dat))
        intercept$method <- "MR Egger intercept"
        intercept$nsnp <- NA
        names(intercept)[names(intercept) == "egger_intercept"] <- "b"
        mro <- rbind(mro, intercept)
        
        
      }, error = function(e) {
        message("Error in MR Egger intercept calculation: ", e$message)
      })
      
      mro$direction <- direction$correct_causal_direction
      mro$direction_p <- direction$steiger_test
      
      mro_fm <- rbind(mro_fm, mro)
      mro_fm <- subset(mro_fm, !is.na(mro_fm$b))
      
	      write.csv(mro_fm, mr_results_by_clump_file)
      
    } else {
      mro <- data.frame()
    }
  })
    iteration_end_time <- Sys.time()  # End time of the iteration
    duration <- iteration_end_time - iteration_start_time
    message("Finished processing for id: ", id, " at ", iteration_end_time, " (Duration: ", duration, " seconds)")
  }
  
	  write.csv(no_cis_pqtls, paste0(output_folder, prefix, "fat_mass_no_genome_wide_cis_pqtls_", clump_r2_name, ".csv"))
	  write.csv(two_or_fewer_pqtls, paste0(output_folder, prefix, "fat_mass_two_or_fewer_cis_pqtls_", clump_r2_name, ".csv"))
  
}


gc()

# PCA-GMM fallback for proteins without genome-wide-significant cis instruments.
if(run==T){
  
  mro_fm_pc <- data.frame()
  
  for(id in unique(no_cis_pqtls)){
    try({
      iteration_start_time <- Sys.time()
      message("Starting processing for id: ", id, " at ", iteration_start_time)
      
      chr_value <- pqtl_map$chr[pqtl_map$Assay == id]
      if(length(chr_value)>1){
        chr_value <- chr_value[1]
      }
      chr <- chr_value
      chr_value <- paste0("chr", chr_value)
      
      tar_dir <- tar_dir_path 
      
      # Define the extraction directory
      extract_dir <- extract_dir_path
      
      # Call the function
      cis_pqtls <- read_cis_pqtl_gz_from_tar(tar_dir, chr_value, id, extract_dir)
      
      rsid_pos_map_dir <- rsid_pos_map_dir_path
      
      # Call the function
      rsid_pos <- read_rsid_map_gz_from_directory(rsid_pos_map_dir, chr_value)
      
      cis_pqtls <- merge(cis_pqtls, rsid_pos, by.x = c("ID", "ALLELE0", "ALLELE1"), by.y = c("ID", "REF", "ALT"), all.x = TRUE)
      
      cis_pqtls$chr <- as.numeric(cis_pqtls$CHROM)
      cis_pqtls$position <- as.numeric(cis_pqtls$GENPOS)
      cis_pqtls$beta <- cis_pqtls$BETA
      cis_pqtls$se <- cis_pqtls$SE
      cis_pqtls$id.exposure <- id
	      cis_pqtls$SNP <- cis_pqtls$rsid
	      
	      cis_pqtls <- cis_pqtls %>% mutate(pval.exposure = (10^(-LOG10P)))
	      cis_pqtls <- cis_pqtls %>%
	        filter(!is.na(rsid), !is.na(pval.exposure)) %>%
	        arrange(pval.exposure) %>%
	        distinct(rsid, .keep_all = TRUE)
      
      start_pos <- pqtl_map$gene_start[pqtl_map$Assay == id]
      end_pos <- pqtl_map$gene_end[pqtl_map$Assay == id]
      
      
      if(clump_or_pc=="strict"){
        cis_pqtls <- cis_pqtls[cis_pqtls$position >= start_pos - 250000 & cis_pqtls$position <= end_pos + 250000, ] #1000000
      }
      if(clump_or_pc=="lax"){
        cis_pqtls <- cis_pqtls[cis_pqtls$position >= start_pos - 1000000 & cis_pqtls$position <= end_pos + 1000000, ] #1000000
      }
      
      
      clump_kb <- 10000
      clump_r2 <- 0.95
      # The LD threshold
      clump_p <- 1
      #remotes::install_github("MRCIEU/genetics.binaRies")
      library(genetics.binaRies)
      #remotes::#install_github("explodecomputer/genetics.binaRies")
      
      keep <- data.frame()
	      tryCatch({
	        keep <- ieugwasr::ld_clump_local(
	          dplyr::tibble(rsid = cis_pqtls$rsid, pval = cis_pqtls$pval.exposure, id = cis_pqtls$id.exposure),
	          plink_bin = genetics.binaRies::get_plink_binary(),
	          bfile = bfile,
          clump_kb = clump_kb,
          clump_r2 = clump_r2,
          clump_p = clump_p
        )
      }, error = function(e) {
        message("Error during clumping")
        keep <- data.frame()
      })
      
	      if(nrow(keep) <= 1){
	        message("Not enough cis variants for PCA-GMM for id: ", id)
	        next
	      }
      
      cis_pqtls <- subset(cis_pqtls, cis_pqtls$rsid %in% unique(keep$rsid))
      
      pqtl_map_bit <- pqtl_map[pqtl_map$Assay == id,]
      pqtl_map_bit <- pqtl_map_bit[, c("Assay", "ensembl_id")]
      
      cis_pqtls <- merge(cis_pqtls, pqtl_map_bit, by.x=c("id.exposure"),
                         by.y=c("Assay"), all.x=TRUE)
      
      cis_pqtls$exposure <- cis_pqtls$id.exposure
      cis_pqtls$beta.exposure <- cis_pqtls$beta
      cis_pqtls$se.exposure <- cis_pqtls$se
      cis_pqtls$effect_allele.exposure <- cis_pqtls$ALLELE1
      cis_pqtls$other_allele.exposure <- cis_pqtls$ALLELE0
      cis_pqtls$eaf.exposure <- cis_pqtls$A1FREQ
      
      gc()
      
      if(chr=="X"){
        outcome <- subset(fat_mass_regenie, fat_mass_regenie$CHROM==23)
      }else{
        outcome <- subset(fat_mass_regenie, fat_mass_regenie$CHROM==as.numeric(chr))
      }
      
      outcome <- 
        outcome %>% mutate(P = (10^(-LOG10P)))
      
      outcome$pval.outcome <- outcome$P
      outcome$id.outcome <- "Whole body fat mass"
      
      outcome$chr_name <- outcome$CHROM
      outcome$chrom_start <- outcome$GENPOS
      outcome$chr <- ifelse(outcome$chr_name==23,
                            "chrX", paste0("chr", outcome$chr_name))
      
      gc()
      
      outcome <- data.frame(outcome)
      
      cis_pqtls$pos <- cis_pqtls$POS19
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
        log_pval = FALSE)
      
      dat <- harmonise_data(cis_pqtls, outcome)
      dat <- subset(dat, dat$mr_keep==TRUE)
      dat <- dat[!duplicated(dat$SNP),]
      
      snps <- unique(dat$SNP)
      
      if (length(dat$SNP) >= 1) {
        tryCatch({
          x.corr <- ieugwasr::ld_matrix(
            variants = snps,
            with_alleles = TRUE,
            plink_bin = genetics.binaRies::get_plink_binary(),
            bfile = bfile
          )
          
          snpnames <- sapply(strsplit(rownames(x.corr), split = "_"), `[`, 1)
          order_indices <- order(snpnames)
          x.corr <- x.corr[order_indices, order_indices]
          
          dat <- dat[order(dat$SNP), ]
          
          ld_matrix <- harmonise_ld_dat(dat, x.corr)
          x.corr <- ld_matrix[["ld"]]
          
          found_snps_modified <- ld_matrix[["x"]][["SNP"]]
          dat <- dat[dat$SNP %in% found_snps_modified, ]
          
          dat <- dat[match(found_snps_modified, dat$SNP), ]
          
          outcome <- outcome[outcome$SNP %in% found_snps_modified, ]
          
          dat <- MendelianRandomization::mr_input(
            bx = dat$beta.exposure,
            bxse = dat$se.exposure,
            by = dat$beta.outcome,
            byse = dat$se.outcome,
            correlation = matrix(),
            exposure = id,
            outcome = "Whole body fat mass",
            snps = found_snps_modified,
            effect_allele = dat$effect_allele.exposure,
            other_allele = dat$other_allele.exposure,
            eaf = dat$eaf.exposure
          )
          
          dat@correlation <- x.corr
          
          mro <- MendelianRandomization::mr_pcgmm(
            dat,
            nx = 34557,
            ny = 367620,
            thres = pc_thresh,
            robust = TRUE,
            alpha = 0.05
          )
          
          mro_dat <- data.frame(
            exposure = mro@Exposure,
            outcome = mro@Outcome,
            PCs = mro@PCs,
            b = mro@Estimate,
            se = mro@StdError,
            pval = mro@Pvalue,
            ci_lwr = mro@CILower,
            ci_upr = mro@CIUpper,
            fstat = mro@Fstat,
            hetr_stat = mro@Heter.Stat, 
            nsnp = nrow(mro@Correlation),
            stringsAsFactors = FALSE
          )
          
	          mro_fm_pc <- rbind(mro_fm_pc, mro_dat)
	          write.csv(mro_fm_pc, pca_gmm_results_file)
          
        }, error = function(e) {
          message("Error during MR PCGMM: ", e$message)
          # Proceed with the rest of the loop
        }, warning = function(w) {
          message("Warning during MR PCGMM: ", w$message)
          # Proceed with the rest of the loop
        })
        
      } else {
        message("Not enough SNPs to perform MR PCGMM for id: ", id)
      }
    })
    iteration_end_time <- Sys.time()  # End time of the iteration
    duration <- iteration_end_time - iteration_start_time
    message("Finished processing for id: ", id, " at ", iteration_end_time, " (Duration: ", duration, " seconds)")
  }
  
  if(nrow(mro_fm_pc) == 0){
    mro_fm_pc <- data.frame(
      exposure = character(), outcome = character(), PCs = numeric(),
      b = numeric(), se = numeric(), pval = numeric(), ci_lwr = numeric(),
      ci_upr = numeric(), fstat = numeric(), hetr_stat = numeric(),
      nsnp = numeric(), stringsAsFactors = FALSE
    )
  }
  write.csv(mro_fm_pc, pca_gmm_results_file)

  mro_fm$code <- tolower(mro_fm$exposure)
  mro_fm$PCs <- NA
  mro_fm$fstat <- NA
  mro_fm$hetr_stat <- NA
  mro_fm$ci_lwr <- mro_fm$b - (1.96*mro_fm$se)
  mro_fm$ci_upr <- mro_fm$b + (1.96*mro_fm$se)
  
  mro_fm_pc$code <- tolower(mro_fm_pc$exposure)
  mro_fm_pc$id.exposure <- mro_fm_pc$exposure
  mro_fm_pc$id.outcome <- mro_fm_pc$outcome
  mro_fm_pc$method <- "PCA-GMM"
  mro_fm_pc$se <- (mro_fm_pc$ci_upr - mro_fm_pc$b)/1.96
  
  mro_fm_pc$direction <- NA
  mro_fm_pc$direction_p <- NA

  
  mro_fm <- rbind(mro_fm, mro_fm_pc)
  
  mro_fm <- subset(mro_fm, !is.na(mro_fm$b))
  
  
  mro_fm$method_group <- ifelse(mro_fm$method == "Wald ratio" |
                                  mro_fm$method == "PCA-GMM" |
                                  mro_fm$method == "Inverse variance weighted", "main", mro_fm$method)
  
  if(clump_or_pc=="strict"){
  
    mro_fm <- subset(mro_fm, mro_fm$method!="PCA-GMM")
    
  }
  
  mro_fm <- mro_fm %>%
    group_by(method_group) %>%
    dplyr::mutate(fdr = p.adjust(pval, method = "fdr"))
  
  mro_fm <- subset(mro_fm, !is.na(mro_fm$b))
  
  lines <- readLines(file.path(input_dir, "coding143.tsv"))
  
  data <- lapply(lines, function(line) {
    parts <- unlist(strsplit(line, ";", fixed = TRUE))
    code <- parts[1]
    description <- paste(parts[-1], collapse = ";")
    return(c(code, description))
  })
  
  coding143_edited <- do.call(rbind, data)
  colnames(coding143_edited) <- c("code", "Protein")
  coding143_edited <- data.frame(coding143_edited)
  
  coding143_edited$code <- gsub('^"|"$', '', coding143_edited$code)
  coding143_edited$Protein <- gsub('^"|"$', '', coding143_edited$Protein)
  
  coding143_edited$code <- tolower(coding143_edited$code)
  coding143_edited$code <- gsub('-', '_', coding143_edited$code)
  coding143_edited <- coding143_edited[-1,]
  
  mro_fm$code <- tolower(mro_fm$exposure)
  mro_fm <- merge(mro_fm, coding143_edited, by.x = "code", by.y = "code", all.x = TRUE)
  
  write.csv(mro_fm, mr_results_file)
  
}


#reading in MR
if(run==F){
  mro_fm <- read.csv(mr_results_file)
}

mro_fm$method <- ifelse(mro_fm$method == "PC-GMM", "PCA-GMM", mro_fm$method)

if(clump_or_pc=="strict"){
  
  mro_fm <- subset(mro_fm, mro_fm$method!="PCA-GMM")
  
}

mro_fm <- mro_fm %>%
  group_by(method_group) %>%
  dplyr::mutate(fdr = p.adjust(pval, method = "fdr"))


#filtering significant results selecting significant results

sig_proteins <- subset(mro_fm, (mro_fm$fdr < 0.05 & mro_fm$method=="Inverse variance weighted") |
                         (mro_fm$fdr < 0.05 & mro_fm$method=="Wald ratio") |
                         (mro_fm$fdr < 0.05 & mro_fm$method=="PCA-GMM")
                       )

sig_proteins <- unique(sig_proteins$exposure)

pleiotropic_total <- unique(subset(mro_fm, mro_fm$fdr < 0.05 & mro_fm$method == "MR Egger intercept")$exposure)

egger_sig <- unique(subset(mro_fm, mro_fm$fdr < 0.05 & mro_fm$method=="MR Egger")$exposure)
  
pleiotropic <- setdiff(pleiotropic_total, egger_sig)

mro_fm_concordance <- subset(mro_fm, (method!="MR Egger intercept") &
                               (method!="MR Egger" |
                                  mro_fm$exposure %in% pleiotropic_total))

mro_fm_concordance$positive <- ifelse(mro_fm_concordance$b > 0, 1, 0)
mro_fm_concordance <- mro_fm_concordance %>% group_by(exposure) %>% 
  dplyr::mutate(positive_sum = sum(positive))
mro_fm_concordance <- mro_fm_concordance %>% group_by(exposure) %>% 
  dplyr::mutate(num_methods = n())
mro_fm_concordance$positive_ratio <- mro_fm_concordance$positive_sum/mro_fm_concordance$num_methods
concordant <- mro_fm_concordance %>% filter(positive_ratio==1 | positive_ratio==0)
concordant <- unique(concordant$exposure)

sig <- mro_fm %>% filter(!(exposure %in% pleiotropic) & (exposure %in% sig_proteins) & (exposure %in% concordant))

sig_main <- mro_fm %>% filter(!(exposure %in% pleiotropic) & (exposure %in% sig_proteins) & (exposure %in% concordant))


if (proteins == "selected"){
  sig_ids <- unique(sig$exposure)
  sig_main_ids <- unique(sig_main$exposure)
}else{
  sig_ids <- sig_proteins
}
  

gc()


#volcano plot 
if(T){
  try({
    
    volcano_data <- subset(mro_fm, mro_fm$method=="Inverse variance weighted" | mro_fm$method=="Wald ratio" | mro_fm$method=="PCA-GMM")
    
    volcano_data <- volcano_data %>%
      mutate(neg_log_p = -log(pval, base=exp(1)),
             significant = exposure %in% sig_ids,
             point_color = ifelse(significant, ifelse(b < 0, "blue", "red"), "grey"))
    
    volcano_data <- volcano_data %>%
      arrange(desc(significant), desc(point_color == "red"), desc(point_color == "blue"))
    
    
    if(proteins=="selected"){
      closest_fdr_value <- max(volcano_data$fdr[volcano_data$fdr < 0.05])
      corresponding_neg_log_p <- volcano_data$neg_log_p[volcano_data$fdr == closest_fdr_value]
      cutoff <- 0.99 * corresponding_neg_log_p
      
      y_breaks <- pretty(volcano_data$neg_log_p)
      y_breaks <- sort(unique(c(y_breaks, cutoff)))
      
      y_labels <- ifelse(y_breaks == cutoff, "FDR-P=0.05", as.character(y_breaks))
      
    }else{
      
      y_breaks <- pretty(volcano_data$neg_log_p)
      cutoff <- -log(0.05, base=exp(1))
      y_breaks <- sort(unique(c(y_breaks, cutoff)))
      
      y_labels <- ifelse(y_breaks == cutoff, "p=0.05", as.character(y_breaks))
      
    }
    
    plt <- ggplot(volcano_data, aes(x = b, y = neg_log_p)) +
      geom_point(size = ifelse(volcano_data$point_color == "grey", 3, 4),  # Increase size for significant points
                 aes(color = volcano_data$point_color, fill = volcano_data$point_color), 
                 stroke = 0.5, shape = 21, alpha = ifelse(volcano_data$point_color == "grey", 0.6, 0.9)) +
      scale_color_manual(
        values = c("blue" = "#1F77B1", "red" = "#D62744", "grey" = "grey"),
        labels = c("blue" = "Negative coefficient", "red" = "Positive coefficient", "grey" = "Not statistically significant"),
        guide = guide_legend(override.aes = list(alpha = 1, size=4, fill = c("#1F77B1", "grey", "#D62744" )))
      ) +
      scale_fill_manual(
        values = c("blue" = "#1F77B1", "red" = "#D62744", "grey" = "grey"),
        guide = "none"
      ) +
      theme_minimal() +
      labs(
        title = element_blank(),
        x = "Regression coefficient (kg per SD difference in NPX)",
        y = "-log(P-value)"
      ) +
      geom_text_repel(
        aes(label = ifelse(fdr < 0.05, as.character(code), '')),
        box.padding = 0.2,
        size = 4,  # Increase text size
        point.padding = 0.3,
        segment.color = 'grey50',
        max.overlaps = 20
      ) +
      geom_hline(yintercept = cutoff, linetype = "dashed", color = "black") +  # Dashed line at cutoff
      scale_y_continuous(
        breaks = y_breaks,
        labels = y_labels
      ) +
      scale_x_continuous(
        limits = c(-1, 2),
        breaks = seq(-1, 2, by = 0.5)
      ) 
    
    library(ckbplotr)
    plt + theme_classic() + theme(axis.line = element_line(),
                                  legend.position = c(0.83, 0.90),
                                  legend.title = element_blank(),
                                  legend.background = element_rect(color = "black", fill = "white", size = 0.25, linetype = "solid"),
                                  text = element_text(size = 12),  # Increase overall text size
                                  axis.title = element_text(size = 12),  # Increase axis title size
                                  axis.text = element_text(size = 12),  # Increase axis labels size
                                  axis.title.y = element_text(vjust=-19.5))  # Adjust margin for y-axis title to bring it closer to the axis
    
    
    ggsave(paste0(output_folder, "plots/", prefix, 'volcano_mr_plot_', clump_r2_name, "_",  pc_thresh_name, '.pdf'), width=8, height=6)
    
    ggsave(paste0(output_folder, "plots/", prefix, 'wide_volcano_mr_plot_', clump_r2_name, "_",   pc_thresh_name, '.pdf'), width=10, height=6)
    
  })
}

#tables
if(proteins=="selected"){
  
  mro_fm$"95% CI" <- paste0(round(mro_fm$ci_lwr, 2), " - ", round(mro_fm$ci_upr, 2))
  mro_fm$"SNPs/PCs" <- ifelse(mro_fm$method == "PCA-GMM", mro_fm$PCs, mro_fm$nsnp)
  
  if("X" %in% colnames(mro_fm)){
    mro_fm <- subset(mro_fm, select=-c(X))
  }
  
  mro_fm <- mro_fm %>%
    group_by(method_group) %>%
    dplyr::mutate(fdr = p.adjust(pval, method = "fdr"))
  
  primary_mr_summary <- subset(mro_fm, mro_fm$method=="Inverse variance weighted" | mro_fm$method=="Wald ratio" | mro_fm$method=="PCA-GMM") 
  primary_mr_summary <- primary_mr_summary[, c(5, 28, 6, 30, 8, 29, 27, 12)]
  colnames(primary_mr_summary) <- c("Olink Assay Code", "Protein name", "Primary MR method", "SNPs/PCs", "Beta", "95% CI", "FDR P", "MR Steiger P")
  primary_mr_summary$`MR Steiger P` <- round(primary_mr_summary$`MR Steiger P`, 2)
  primary_mr_summary$Beta <- round(primary_mr_summary$Beta, 2)
  primary_mr_summary$"FDR P" <- round(primary_mr_summary$"FDR P", 2)
  
  egger_intercept_summary <- subset(mro_fm, mro_fm$method=="MR Egger intercept") 
  egger_intercept_summary <- egger_intercept_summary[, c(5, 27)]
  colnames(egger_intercept_summary) <- c("Olink Assay Code", "MR-Egger intercept FDR P")
  egger_intercept_summary$"MR-Egger intercept FDR P" <- round(egger_intercept_summary$"MR-Egger intercept FDR P", 2)
  
  primary_mr_summary <- merge(primary_mr_summary, egger_intercept_summary, by.x="Olink Assay Code", all.x=T)
  
  write.csv(primary_mr_summary, paste0(output_folder, prefix, "fat_mass_primary_mr_summary_", clump_r2_name, "_", pc_thresh_name, ".csv"))
  
  protein_annotation <- pqtl_map[, c(5, 6, 4, 10, 11, 12, 13)]
  colnames(protein_annotation) <- c("Olink Assay Code", "Olink Panel", "UniProt ID", "Coding gene ENSEMBL ID", "Chromosome", "Coding gene start", "Coding gene end")
  protein_annotation_names <- coding143_edited
  colnames(protein_annotation_names) <- c("code", "UniProt name")
  protein_annotation_names$code <- toupper(protein_annotation_names$code)
  protein_annotation <- merge(protein_annotation, protein_annotation_names, by.x="Olink Assay Code", by.y="code")
  protein_annotation <- protein_annotation[, c(1, 8, 2, 3, 4, 5, 6, 7)]
  write.csv(protein_annotation, paste0(output_folder, prefix, "selected_protein_annotation_", clump_r2_name, "_", pc_thresh_name, ".csv"))
  
}

#plotting/summary forest plot
if(T){
  lax_primary_results <- data.frame()
  if(file.exists(lax_mr_results_file)){
    lax_primary_results <- read.csv(lax_mr_results_file)
    lax_primary_results <- subset(
      lax_primary_results,
      method %in% c("Inverse variance weighted", "Wald ratio") & fdr < 0.05
    )
  }
  
  if(proteins=="selected" & clump_or_pc=="lax"){
    plot_dta <- subset(mro_fm, (mro_fm$method=="Inverse variance weighted" |
                                  mro_fm$method=="Wald ratio" |
                                  mro_fm$method=="PCA-GMM") &
                         (mro_fm$exposure %in% sig_ids))}
  
  if(proteins=="selected" & clump_or_pc=="strict"){
    lax_sig_ids <- if(nrow(lax_primary_results) > 0) unique(lax_primary_results$exposure) else unique(mro_fm$exposure)
    mro_fm$fdr <- round(mro_fm$fdr, 2)
    plot_dta <- subset(mro_fm, (mro_fm$method=="Inverse variance weighted" |
                                  mro_fm$method=="Wald ratio")
                       & (mro_fm$exposure %in% lax_sig_ids)
    )}
  
  plot_dta$mthd <- ifelse(plot_dta$method=="Inverse variance weighted", "IVW",
                          ifelse(plot_dta$method=="Wald ratio", "WR", "PCA-GMM"))
  
  plot_dta$ivs <- ifelse(plot_dta$mthd=="IVW", plot_dta$nsnp,
                         ifelse(plot_dta$mthd=="WR", plot_dta$nsnp, plot_dta$PCs))
  
  
  by_panel <- read.csv(file.path("selection_output", "selected_proteins_by_panel.csv"))
  by_panel <- by_panel[, -c(1)]
  
  plot_dta <- merge(plot_dta, by_panel[, c(1, 7)], by.x="exposure", by.y="Assay")
  
  
  plot_dta$fill <- "white"
  
  plot_dta$shape <- 22
  
  plot_dta <- plot_dta %>% arrange(desc(b))
  
  protein_order_FM <- plot_dta$id.exposure  
  
  plot_dta$key <- plot_dta$id.exposure
  
  if(proteins=="selected"){
    xlim_vals <- make_forest_axis_limits(plot_dta)
    xticks_vals <- seq(from = xlim_vals[1], to = xlim_vals[2], by = 0.5)
  }else{
    xlim_vals <- c(-1.5, 1)
    xticks_vals <- seq(from = xlim_vals[1], to = xlim_vals[2], by = 0.5)
  }
  saveRDS(list(order = protein_order_FM, xlim = xlim_vals, xticks = xticks_vals), file = fm_layout_file)
  
  if(clump_or_pc=="lax"){
    plot_dta$fill <- ifelse(plot_dta$Panel=="Inflammation", "#FF9999",
                            ifelse(plot_dta$Panel=="Cardiometabolic", "#07CADF",
                                   ifelse(plot_dta$Panel=="Oncology", "#FFFF99", "#99ABFF")))
  }else{
    plot_dta$fill <- ifelse((plot_dta$fdr < 0.05), ifelse(
      plot_dta$Panel=="Inflammation", "#FF9999",
      ifelse(plot_dta$Panel=="Cardiometabolic", "#07CADF",
             ifelse(plot_dta$Panel=="Oncology", "#FFFF99", "#99ABFF"))), "white")
    
    plot_dta$Panel <- "            "
    
  }
  
  if (length(unique(plot_dta$method)) > 1) {
    plot_cols <- c("ivs", "mthd")
    plot_col_labels <- c("SNPs/PCs", "Method")
    
  } else {
    
    plot_cols <- list("nsnp")
    plot_col_labels <- list("SNPs")
  }
  
  
  ckbplotr::forest_plot(panels = list(plot_dta),
                        panel.headings = "FAT MASS",
                        col.key = "key",
                        panel.names = c("FAT MASS"),
                        row.labels.heading = "Protein",
                        col.estimate = "b",
                        col.stderr = "se",
                        showcode = FALSE,
                        exponentiate = FALSE,
                        col.right.heading = c("Beta (95% CI)"),
                        xlab = "Regression coefficient (kg per SD difference in NPX)",
                        col.left = plot_cols,
                        col.left.heading = plot_col_labels,
                        pointsize = 5,
                        scalepoints = FALSE,
                        shape = "shape", 
                        stroke = 0.4, 
                        base_line_size = 0.4,
                        ciunder = FALSE, 
                        fill = "fill",
                        nullval = 0,
                        xlim = xlim_vals,
                        xticks = xticks_vals,
                        col.left.hjust = 0.2,
                        base_size = 10) + theme(family = "serif")
  
  ggsave(paste0(output_folder, "plots/", prefix, 'panel_protein_to_fat_mass_mr_forest_plot_', clump_r2_name, "_", pc_thresh_name, '.pdf'), 
         width = ifelse(nrow(plot_dta)>3,
                        8, 5), height = ifelse(nrow(plot_dta)>6,
                                               nrow(plot_dta)*0.6, nrow(plot_dta)),
         unit = 'in', dpi = 500)
  
  ggsave(paste0(output_folder, "plots/", prefix, 'wide_panel_protein_to_fat_mass_mr_forest_plot_', clump_r2_name, "_", pc_thresh_name, '.pdf'), 
         width = ifelse(nrow(plot_dta)>3,
                        9, 5), height = ifelse(nrow(plot_dta)>6,
                                               nrow(plot_dta)*0.4, nrow(plot_dta)),
         unit = 'in', dpi = 500)
  
  
  
}

  
