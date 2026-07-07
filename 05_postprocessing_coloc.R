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


#reading in coloc sharepro files
coloc_confirmed <- list()
combined_coloc_df <- data.frame()

for(id in sig_ids) {
  # Construct the file path
  file_path <- file.path(paste0(output_folder, "ShareProColoc/", prefix, id, ".sharepro.txt"))

  if (file.exists(file_path)) {
    print(file_path)
    df <- read.table(file_path, header = TRUE, sep = "\t")
    if (any(df$share > 0.8)) {
      try({
        coloc_confirmed <- append(coloc_confirmed, id)
      })
      
    }
    
    try({
      df$exposure <- id
      combined_coloc_df <- rbind(combined_coloc_df, df)
    })
    
  } else {
    warning(paste("File not found:", file_path))
  }
}

coloc_confirmed <- unique(unlist(coloc_confirmed))

length(unique(combined_coloc_df$exposure))

write.csv(combined_coloc_df, paste0(output_folder, prefix, "sharepro_colocalisation_results_", clump_r2_name, "_", pc_thresh_name, ".csv"), row.names = FALSE)
