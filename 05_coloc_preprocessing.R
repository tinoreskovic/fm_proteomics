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

run = T #T F

#writing files for SharePro
if(run==T){
  for (id in sig_ids){
    tryCatch({
      message("Writing ShareProColoc files for protein: ", id)
    
      if (file.exists(paste0(output_folder, "ShareProColoc/", prefix, id, ".txt"))){
        message("Files already exists. Skipping.")
        next
      }
      
      gc()
      chr_value <- pqtl_map$chr[pqtl_map$Assay == id]
      if(length(chr_value)>1){
        chr_value <- chr_value[1]
      }
      chr <- chr_value
      chr_value <- paste0("chr", chr_value)

      extract_dir <- extract_dir_path
      tar_dir <- tar_dir_path 
      
      cis_pqtls <- read_cis_pqtl_gz_from_tar(tar_dir, chr_value, id, extract_dir)
      
      rsid_pos_map_dir <- rsid_pos_map_dir_path

      rsid_pos <- read_rsid_map_gz_from_directory(rsid_pos_map_dir, chr_value)
      
      # Merge the pQTL data with the rsid position map
      cis_pqtls <- merge(cis_pqtls, rsid_pos, by.x = c("ID", "ALLELE0", "ALLELE1"), by.y = c("ID", "REF", "ALT"), all.x = TRUE)
      
      cis_pqtls$chr <- as.numeric(cis_pqtls$CHROM)
      cis_pqtls$position <- as.numeric(cis_pqtls$GENPOS)
      cis_pqtls$beta <- cis_pqtls$BETA
      cis_pqtls$se <- cis_pqtls$SE
      cis_pqtls$id.exposure <- id
      cis_pqtls$SNP <- cis_pqtls$rsid

      cis_pqtls <- cis_pqtls %>% mutate(pval.exposure = (10^(-LOG10P)))
      
      start_pos <- pqtl_map$gene_start[pqtl_map$Assay == id]
      end_pos <- pqtl_map$gene_end[pqtl_map$Assay == id]
      
      cis_pqtls <- cis_pqtls[cis_pqtls$position >= start_pos - 1000000 & cis_pqtls$position <= end_pos + 1000000, ] 
      
      
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
      
      common_snps <- dat$SNP
      
      ld_matrix <- ieugwasr::ld_matrix(
        common_snps,
        with_alleles = TRUE,
        plink_bin = genetics.binaRies::get_plink_binary(),
        bfile = bfile
      )
      
      snpnames <- sapply(strsplit(rownames(ld_matrix), split = "_"), `[`, 1)
      order_indices <- order(snpnames)
      ld_matrix <- ld_matrix[order_indices, order_indices]
      
      dat <- dat[order(dat$SNP), ]
      
      ld_matrix <- harmonise_ld_dat(dat, ld_matrix)
      common_snps <- ld_matrix[["x"]][["SNP"]]
      common_snps <- ld_matrix[["x"]][["SNP"]]
      
      cis_pqtls <- subset(cis_pqtls, cis_pqtls$rsid %in% common_snps)
      fat_mass_snps <- data.frame(outcome)
      fat_mass_snps <- subset(fat_mass_snps, fat_mass_snps$SNP %in% common_snps)
      fat_mass_snps <- fat_mass_snps[!duplicated(fat_mass_snps$SNP),]
      cis_pqtls <- cis_pqtls[!duplicated(cis_pqtls$SNP),]
      
      gc()
      
      # the order of SNPs has to be the same in both datasets
      cis_pqtls <- cis_pqtls[order(match(cis_pqtls$rsid, common_snps)), ]
      fat_mass_snps <- fat_mass_snps[order(match(fat_mass_snps$SNP, common_snps)), ]
      
      cis_pqtls$z <- cis_pqtls$beta / cis_pqtls$se
      
      ld <- ld_matrix[["ld"]]
      
      
      cis_pqtls <- cis_pqtls[, c(16, 3, 4, 7, 11, 12, 24)]
      names(cis_pqtls) <- c("SNP", "A1", "A2", "EAF", "BETA", "SE", "P")
      cis_pqtls$N <- 34557
      
      write.table(cis_pqtls, paste0(output_folder, "ShareProColoc/", prefix, id, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
      
      fat_mass_snps <- fat_mass_snps[, c(2, 3, 4, 5, 6, 7, 8)]
      names(fat_mass_snps) <- c("SNP", "A1", "A2", "EAF", "BETA", "SE", "P")
      fat_mass_snps$N <- 367620
      
      write.table(fat_mass_snps, paste0(output_folder, "ShareProColoc/", prefix, id, "_fat_mass", ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
      
      write.table(ld, paste0(output_folder, "ShareProColoc/", prefix, id, "_ld", ".ld"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    })
  }
}
