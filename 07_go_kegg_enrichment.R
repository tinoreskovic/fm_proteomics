library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(dplyr)
library(tidyr)
library(stringr)
library(httr)

dir.create("plots", showWarnings = FALSE)

protein_map <- read.delim(file.path("input", "olink_protein_map_3k_v1.tsv"))

# These proteins were not included in the randomised LASSO screening set.
excluded_from_screen <- c("AMY2B", "CST1", "CTSS", "NPM1", "PCOLCE", "TACSTD2")
screened_proteins <- subset(protein_map, !(Assay %in% excluded_from_screen))

background_uniprot <- unique(screened_proteins$UniProt)
background_entrez <- bitr(
  background_uniprot,
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

mr_supported_proteins <- c(
  "ADAMTS15", "ADM", "AGER", "APOD", "CALB1", "CCL27", "CCN5",
  "CD300LG", "CLMP", "CRIM1", "DNER", "GHRL", "NEFL", "PRRT3",
  "PSPN", "ROBO1", "RTN4R", "SCGB3A1", "SCGB3A2", "SHBG", "WFIKKN2"
)

supported_map <- subset(screened_proteins, Assay %in% mr_supported_proteins)
supported_uniprot <- unique(supported_map$UniProt)
supported_entrez <- bitr(
  supported_uniprot,
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

write.csv(
  background_entrez,
  "enrichment_background_screened_proteins.csv",
  row.names = FALSE
)

go_enrichment <- enrichGO(
  gene = supported_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_results <- go_enrichment@result
go_results$p.adjust <- round(go_results$p.adjust, 3)
write.csv(
  go_results,
  "go_biological_process_enrichment_screened_background.csv",
  row.names = FALSE
)

kegg_enrichment <- enrichKEGG(
  gene = supported_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05
)

kegg_results <- kegg_enrichment@result
kegg_results$p.adjust <- round(kegg_results$p.adjust, 3)
write.csv(
  kegg_results,
  "kegg_enrichment_screened_background.csv",
  row.names = FALSE
)

significant_go <- subset(go_results, p.adjust < 0.05)

if(nrow(significant_go) > 1){
  chord_data <- significant_go %>%
    mutate(geneID = strsplit(as.character(geneID), "/")) %>%
    unnest(geneID) %>%
    dplyr::select(geneID, Description)
  
  names(chord_data) <- c("Protein", "Biological_Process")
  
  connection_matrix <- table(chord_data$Protein, chord_data$Biological_Process)
  connection_matrix <- as.matrix(connection_matrix)
  connection_matrix <- connection_matrix[rev(order(rownames(connection_matrix))), ]
  
  wrap_labels <- function(labels, width = 20) {
    sapply(labels, function(x) paste(strwrap(x, width = width), collapse = "\n"))
  }
  
  rownames(connection_matrix) <- wrap_labels(rownames(connection_matrix), width = 20)
  colnames(connection_matrix) <- wrap_labels(colnames(connection_matrix), width = 20)
  
  pdf("plots/go_enrichment_chord_diagram_screened_background.pdf", width = 16, height = 16)
  circos.clear()
  par(mar = c(1, 1, 1, 1))
  chordDiagram(
    connection_matrix,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.2)
  )
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 1.35
    )
  }, bg.border = NA)
  dev.off()
}else{
  message("One or fewer GO biological-process terms pass FDR < 0.05; no chord diagram was generated.")
}

extract_function_text <- function(text) {
  match <- str_match(text, "-!- FUNCTION:((?:.|\\n)*?)-!-")
  if(!is.na(match[1, 2])){
    return(str_trim(str_remove(match[1, 2], "-!-$")))
  }
  NA
}

uniprot_function <- data.frame(
  UniProt = character(),
  Function = character(),
  stringsAsFactors = FALSE
)

for(id in supported_uniprot){
  response <- GET(paste0("https://rest.uniprot.org/uniprotkb/", id, ".txt"))
  if(status_code(response) == 200){
    content_text <- content(response, as = "text", encoding = "UTF-8")
    function_text <- extract_function_text(content_text)
  }else{
    function_text <- NA
  }
  uniprot_function <- rbind(
    uniprot_function,
    data.frame(UniProt = id, Function = function_text, stringsAsFactors = FALSE)
  )
}

uniprot_function <- merge(
  uniprot_function,
  supported_map[, c("UniProt", "Assay")],
  by = "UniProt",
  all.x = TRUE
)

write.csv(
  uniprot_function,
  "uniprot_function_annotations_mr_supported_proteins.csv",
  row.names = FALSE
)
