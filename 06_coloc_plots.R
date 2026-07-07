library(gridExtra)
library(ggplot2)

rsid_pos_map_dir_path <- file.path("summary_stats", "rsid_pos_map")

for (i in coloc_confirmed){
  coloc_pqtl_map <- subset(pqtl_map, pqtl_map$Assay == i)

  chr_value <- coloc_pqtl_map$chr[coloc_pqtl_map$Assay == i]
  if(length(chr_value)>1){
    chr_value <- chr_value[1]
  }
  chr <- chr_value
  chr_value <- paste0("chr", chr_value)
  
  rsid_pos_map_dir <- rsid_pos_map_dir_path
  
  rsid_pos <- read_rsid_map_gz_from_directory(rsid_pos_map_dir, chr_value)
  
  fm_g <- read.delim(paste0("ShareProColoc/", i, "_fat_mass.txt"), header=TRUE)
  fm_g <- merge(fm_g, rsid_pos, by.x = c("SNP", "A1", "A2"), by.y = c("rsid", "REF", "ALT"), all.x=TRUE)
  
  p_g <- read.delim(paste0("ShareProColoc/", i, ".txt"), header=TRUE)
  p_g <- merge(p_g, rsid_pos, by.x = c("SNP", "A1", "A2"), by.y = c("rsid", "REF", "ALT"), all.x=TRUE)
  
  coloc_df <- subset(combined_coloc_df, combined_coloc_df$exposure == i &
                       combined_coloc_df$share > 0.8)
  
  rsid_list <- unlist(strsplit(as.character(coloc_df$cs), split = "/"))
  group_lengths <- sapply(strsplit(as.character(coloc_df$cs), split = "/"), length)
  rsid_df <- data.frame(rsid = rsid_list, group = rep(1:nrow(coloc_df), times = group_lengths))
  
  fm_g <- fm_g %>%
    left_join(rsid_df, by = c("SNP" = "rsid")) %>%
    mutate(group = as.factor(group))
  
  p_g <- p_g %>%
    left_join(rsid_df, by = c("SNP" = "rsid")) %>%
    mutate(group = as.factor(group))
  
  fm_g$POS38_MB <- fm_g$POS38 / 1e6
  p_g$POS38_MB <- p_g$POS38 / 1e6
  
  pastel_colors <- c("#E69F00", "#E57373", "#009E73", "#F0E442", "#FFC107", "#D55E00", "#CC79A7", "black", "#F4A582", "#8E44AD")
  
  pastel_colors <- c(
    "#F0E442",
    "#E69F99",
    "#009E73",
    "#F28E2F",
    "#F0E442",
    "#D55E00",
    "#999999",
    "#CC79A7"
  )
  
  names(pastel_colors) <- levels(fm_g$group)
  
  plot_p <- ggplot(p_g, aes(x = POS38_MB, y = -log10(P))) +
    geom_point(color = '#1F77B1', alpha=0.65) + 
    geom_point(data = p_g %>% filter(!is.na(group)), aes(fill = group), alpha=0.90, size = 4, shape = 21, stroke = 0.3, color = "black") +  # Highlight group points
    scale_fill_manual(values = pastel_colors, na.translate = FALSE) +
    guides(color = "none", alpha="none") +
    labs(title=paste0("Associations with ", i), fill = "Effect group") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(0, 1),
          legend.justification = c(0, 1))
  
  plot_fm <- ggplot(fm_g, aes(x = POS38_MB, y = -log10(P))) +
    geom_point(color = '#1F77B1', alpha=0.65) +
    geom_point(data = fm_g %>% filter(!is.na(group)), aes(fill = group), alpha=0.9, size = 3.5, shape = 21, stroke = 0.3, color = "black") +  # Highlight group points
    scale_fill_manual(values = pastel_colors, na.translate = FALSE) +
    guides(color = "none", alpha="none") +
    labs(x = paste0("Position in cis region, on chromosome ", chr, " (MB)"), y = "-log10(P)",
         title="Associations with fat mass", fill = "Effect group") +
    theme_minimal() +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1))
  
  pdf(paste0("plots/coloc_plot_", i, ".pdf"), width = 10, height = 8)
  grid.arrange(plot_p, plot_fm, ncol = 1, nrow = 2, heights = c(1, 1))
  dev.off()
}
