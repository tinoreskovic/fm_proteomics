library(dbparser)
library(otargen)
library(httr)
library(jsonlite)
library(dplyr)

pqtl_map <- read.delim(file.path("input", "olink_protein_map_3k_v1.tsv"))


sig_ids <-  c("ADAMTS15", "ADM", "AGER", "APOD",
              "CALB1", "CCL27", "CCN5", "CD300LG",
              "CLMP", "CRIM1", "DNER", "GHRL", "NEFL", 
              "PRRT3", "PSPN", "ROBO1", "RTN4R",
              "SCGB3A1", "SCGB3A2", "SHBG", "WFIKKN2")

pqtl_map <- subset(pqtl_map, pqtl_map$Assay %in% sig_ids)

ensemblIDs <- unique(pqtl_map$ensembl_id)

results_list <- list()

# API request for information on known drugs
get_known_drugs <- function(ensembl_id) {
  query_string <- '
  query KnownDrugsQuery($ensgId: String!, $cursor: String, $freeTextQuery: String, $size: Int = 10) {
    target(ensemblId: $ensgId) {
      id
      knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
        count
        cursor
        rows {
          phase
          status
          urls {
            name
            url
          }
          disease {
            id
            name
          }
          drug {
            id
            name
            mechanismsOfAction {
              rows {
                actionType
                targets {
                  id
                }
              }
            }
          }
          drugType
          mechanismOfAction
        }
      }
    }
  }
  '
  
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  variables <- list("ensgId" = ensembl_id)
  post_body <- list(query = query_string, variables = variables)
  
  response <- httr::POST(url=base_url, body=post_body, encode='json')
  
  if (response$status_code == 200) {
    content_data <- content(response, "parsed", simplifyDataFrame = TRUE)
    return(content_data$data$target$knownDrugs$rows)
  } else {
    warning(paste("Failed to retrieve data for Ensembl ID:", ensembl_id))
    return(NULL)
  }
}

# Loop through each Ensembl ID and collect the results
for (ensembl_id in ensemblIDs) {
  drug_data <- get_known_drugs(ensembl_id)
  if (!is.null(drug_data) && length(drug_data) > 0) {
    results_list[[ensembl_id]] <- drug_data
  }
}

# Combine the results into a dataframe
clinical_precedence <- bind_rows(
  lapply(names(results_list), function(ensembl_id) {
    df <- as.data.frame(results_list[[ensembl_id]])
    if (nrow(df) > 0) {
      df$ensembl_id <- ensembl_id
      df$olink_assay_id <- pqtl_map$Assay[pqtl_map$ensembl_id == ensembl_id]
    }
    return(df)
  })
)

# Transform the urls column
clinical_precedence$urls <- sapply(clinical_precedence$urls, function(urls_df) {
  if (is.data.frame(urls_df) && nrow(urls_df) > 0) {
    urls <- sapply(urls_df$url, function(url) {
      sub(".*\\/([^\\/]+)$", "\\1", url)
    })
    return(paste(urls, collapse = ", "))
  } else {
    return(NA)
  }
})


results_list <- list()



# API request for additional information on target
get_target_annotation <- function(ensembl_id) {
  query_string <- "
    query target($ensemblId: String!){
      target(ensemblId: $ensemblId){
        id
        approvedSymbol
        biotype
        geneticConstraint {
          constraintType
          exp
          obs
          score
          oe
          oeLower
          oeUpper
        }
        tractability {
          label
          modality
          value
        }
      }
    }
  "
  
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  variables <- list("ensemblId" = ensembl_id)
  post_body <- list(query = query_string, variables = variables)
  
  response <- httr::POST(url=base_url, body=post_body, encode='json')
  
  if (response$status_code == 200) {
    content_data <- content(response, "parsed", simplifyDataFrame = TRUE)
    return(content_data$data$target)
  } else {
    warning(paste("Failed to retrieve data for Ensembl ID:", ensembl_id))
    return(NULL)
  }
}

# Loop through each Ensembl ID and collect the results
for (ensembl_id in ensemblIDs) {
  target_data <- get_target_annotation(ensembl_id)
  if (!is.null(target_data)) {
    results_list[[ensembl_id]] <- target_data
  }
}

open_targets <- data.frame()
for (i in 1:length(ensemblIDs)){
  ensembl_id <- ensemblIDs[i]
  tractability <- results_list[[ensembl_id]][["tractability"]]
  tractability$ensembl_id <- ensembl_id
  tractability$approvedSymbol <- results_list[[ensembl_id]][["approvedSymbol"]]
  open_targets <- rbind(open_targets, tractability)
}

clinical_open_targets <- subset(open_targets, 
                                open_targets$label=="Approved Drug" |
                                open_targets$label=="Advanced Clinical" |
                                open_targets$label=="Phase 1 Clinical")


sm_open_targets <- subset(open_targets, open_targets$modality=="SM") 
