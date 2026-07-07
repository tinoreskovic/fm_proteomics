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
library(monaLisa)
library(impute)
library(caret)
library(glmnet)
library(ggplot2)
library(scales)  
library(dplyr)
library(ggrepel)
library(parallel)
library(stabs)
library(monaLisa)

options(scipen = 999)

input_dir <- "input"
controlled_access_dir <- "controlled_access"
selection_dir <- "selection_output"
dir.create(selection_dir, showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

fm_ppp_complement <- read.delim(file.path(controlled_access_dir, "fat_mass_proteomics_phenotypes.tsv"))


#https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1016
olink_batch_number <- read_delim(file.path(input_dir, "olink_batch_number.txt"), 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

olink_batch_number$PlateID <- as.character(olink_batch_number$PlateID)

olink_batch_number$PlateID <- sub("^0+", "", olink_batch_number$PlateID)


fm_ppp_complement$protein_plate <- as.character(fm_ppp_complement$protein_plate)

fm_ppp_complement <- merge(fm_ppp_complement, olink_batch_number, how="left", by.x = "protein_plate",
                           by.y = "PlateID", all.x = TRUE)

fm_ppp_complement <- merge(fm_ppp_complement, olink_batch_number, how="left", by.x = "protein_plate_2",
                           by.y = "PlateID", all.x = TRUE)

fm_ppp_complement <- merge(fm_ppp_complement, olink_batch_number, how="left", by.x = "protein_plate_3",
                           by.y = "PlateID", all.x = TRUE)

#32594
prot_obs <- subset(fm_ppp_complement, fm_ppp_complement$Batch.x<7 &
                     is.na(fm_ppp_complement$consortium_selected) &
                     fm_ppp_complement$ethnic_group==1 &
                     !is.na(fm_ppp_complement$fat_mass) &
                     is.na(fm_ppp_complement$sex_chromosome_aneuploidy) &
                     is.na(fm_ppp_complement$outliers_for_heterozygosity_or_missing) &
                     (fm_ppp_complement$genetic_sex==fm_ppp_complement$sex))



keep_cols <- c("IID", "fat_mass")

prot_obs <- subset(prot_obs, select=keep_cols)

protein_npx <- read.csv(file.path(controlled_access_dir, "protein_npx.csv"))


prot_obs <- merge(prot_obs, protein_npx, how="left", by.x = "IID", by.y = "eid",
                  all.x = TRUE)


missing_percentage <- sapply(prot_obs[, 3:2925], function(x) sum(is.na(x)) / length(x))


#amy2b cst1 ctss (glipr1) npm1 pcolce tacstd2 have more than 10% missing values
columns_to_exclude <- names(which(missing_percentage > 0.1))
# Exclude these columns from the dataframe
prot_obs <- prot_obs[, !(names(prot_obs) %in% columns_to_exclude)]

predictors <- prot_obs[, 3:2918]
outcome <- prot_obs$fat_mass

predictors_matrix <- as.matrix(predictors)


# Impute missing values with the mean of each column
imputed_data <- impute::impute.knn(predictors_matrix, rng.seed=362436069)
imputed_data <- data.frame(imputed_data$data)

inverse_rank_normalize <- function(x) {
  n <- sum(!is.na(x))
  ranks <- rank(x, na.last = "keep")
  percentiles <- (ranks - 0.5) / n
  qnorm(percentiles)
}

scaled_data <- apply(imputed_data, 2, inverse_rank_normalize)

x_matrix <- as.matrix(scaled_data)
y_vector <- as.vector(outcome)


nrow(x_matrix)
ncol(x_matrix)

x_matrix <- x_matrix[, 2:2917]


set.seed(1234)
num_cores <- detectCores() - 2

# Stability selection parameters grid
cutoff_grid <- seq(0.6, 0.9, by = 0.05)
weakness_grid <- seq(0.2, 0.8, by = 0.1)
PFER <- 2  

library(monaLisa)
# Cross-validation setup
k_folds <- 10
cv_folds <- sample(1:k_folds, nrow(x_matrix), replace = TRUE)
cv_results <- array(NA, dim = c(length(cutoff_grid), length(weakness_grid), k_folds))

cv_r2_results <- array(NA, dim = c(length(cutoff_grid), length(weakness_grid), k_folds))

# Stability selection with specific cutoff and weakness
perform_stability_selection <- function(cutoff, weakness, X_train, y_train) {
  stab_results <- randLassoStabSel(x = X_train, y = y_train, weakness = weakness, cutoff = cutoff, PFER = PFER, mc.cores = num_cores)
  max_selection_frequencies <- stab_results@metadata[["stabsel.params.max"]]
  selection_df <- data.frame(
    feature = colnames(X_train),
    frequency = max_selection_frequencies
  )
  selected_features <- subset(selection_df, selection_df$frequency > cutoff)$feature
  return(selected_features)
}

# Evaluate stability selection across folds
for (i in seq_along(cutoff_grid)) {
  for (j in seq_along(weakness_grid)) {
    cutoff <- cutoff_grid[i]
    weakness <- weakness_grid[j]
    
    for (fold in 1:k_folds) {
      X_train <- x_matrix[cv_folds != fold, ]
      y_train <- y_vector[cv_folds != fold]
      X_test <- x_matrix[cv_folds == fold, ]
      y_test <- y_vector[cv_folds == fold]
      
      selected_features <- perform_stability_selection(cutoff, weakness, X_train, y_train)
      
      if (length(selected_features) > 0) {
        # Fit OLS regression on selected features
        fit <- lm(y_train ~ ., data = as.data.frame(X_train[, selected_features, drop = FALSE]))
        preds <- predict(fit, as.data.frame(X_test[, selected_features, drop = FALSE]))
        cv_results[i, j, fold] <- mean((y_test - preds)^2)  # MSE
        
        # R-squared
        ss_res <- sum((y_test - preds)^2)
        ss_tot <- sum((y_test - mean(y_test))^2)
        cv_r2_results[i, j, fold] <- 1 - (ss_res / ss_tot)
        
      } else {
        cv_results[i, j, fold] <- NA
        
        cv_r2_results[i, j, fold] <- NA
      }
    }
  }
}

# Mean CV MSE for each cutoff and weakness
mean_cv_errors <- apply(cv_results, c(1, 2), mean, na.rm = TRUE)

mean_cv_r2 <- apply(cv_r2_results, c(1, 2), mean, na.rm = TRUE)

optimal_indices <- which(mean_cv_errors == min(mean_cv_errors), arr.ind = TRUE)
optimal_cutoff <- cutoff_grid[optimal_indices[1]]
optimal_weakness <- weakness_grid[optimal_indices[2]]

rownames(mean_cv_errors) <- paste("cutoff", cutoff_grid, sep = "_")
colnames(mean_cv_errors) <- paste("weakness", weakness_grid, sep = "_")

write.csv(mean_cv_errors, file = file.path(selection_dir, "stability_selection_cv_mse.csv"))

rownames(mean_cv_r2) <- paste("cutoff", cutoff_grid, sep = "_")
colnames(mean_cv_r2) <- paste("weakness", weakness_grid, sep = "_")

write.csv(mean_cv_r2, file = file.path(selection_dir, "stability_selection_cv_r2.csv"))

# Perform final stability selection with the optimal cutoff and weakness

final_stab_results <- randLassoStabSel(x = x_matrix, y = y_vector, weakness = optimal_weakness, cutoff = optimal_cutoff, PFER = PFER, mc.cores = num_cores)

save(final_stab_results, file = file.path(selection_dir, "stability_selection_fit.rda"))

#load(file.path(selection_dir, "stability_selection_fit.rda"))
optimal_cutoff <- 0.9
optimal_weakness <- 0.8


mat <- as.matrix(SummarizedExperiment::colData(final_stab_results))
#Maximum value in each row for columns 4 to 40
max_values <- apply(mat[, 4:43], 1, max)
# Combine the first three columns with the max values
mat <- cbind(mat[, 1:3], max_values)
colnames(mat)[4] <- "max_value"
mat <- data.frame(mat)

max_selection_frequencies <- max_values


# A dataframe of features and their maximum selection frequencies

selection_df <- data.frame(
  feature = colnames(x_matrix),
  frequency = max_selection_frequencies
)


write.csv(selection_df, file = file.path(selection_dir, "stability_selection_frequencies.csv"))

# Select features based on the optimal cutoff
final_selected_features <- subset(selection_df, selection_df$frequency > optimal_cutoff)$feature


#https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=143
lines <- readLines(file.path(input_dir, "coding143.tsv"))

# Parse lines manually
data <- lapply(lines, function(line) {
  parts <- unlist(strsplit(line, ";", fixed = TRUE))
  # Assuming the first part is the code and the rest is the description
  code <- parts[1]
  description <- paste(parts[-1], collapse = ";")  # Rejoin remaining parts if any
  return(c(code, description))
})

# Convert list to data frame
coding143 <- do.call(rbind, data)
colnames(coding143) <- c("code", "Protein")
coding143 <- data.frame(coding143)

coding143$code <- gsub('^"|"$', '', coding143$code)
coding143$Protein <- gsub('^"|"$', '', coding143$Protein)

coding143$code <- tolower(coding143$code)
coding143$code <- gsub('-', '_', coding143$code)
coding143 <- coding143[-1,]

selection_df <- read.csv(file.path(selection_dir, "stability_selection_frequencies.csv"))

# Merge the selection_df with coding143 to get the Protein names
merged_selection_df <- merge(selection_df, coding143, by.x = "feature", by.y = "code", all.x = TRUE)

# Filter for plotting
filtered_selection_df <- merged_selection_df %>% filter(frequency > 0.5)

filtered_selection_df <- filtered_selection_df %>% arrange(desc(frequency))

optimal_cutoff=0.9 
label_df <- filtered_selection_df %>% filter(frequency > optimal_cutoff)

# Plot the selection frequencies 
ggplot(filtered_selection_df, aes(x = frequency, y = reorder(feature, frequency))) +
  geom_point(size = 1) +  
  labs(title = "",
       x = "Selection frequency",
       y = "") +  
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid.minor.y = element_blank(),  
    panel.grid.major.y = element_blank()
  ) +
  geom_text(data = label_df,  
            aes(x = 0.3, y = reorder(feature, frequency), label = feature),
            hjust = 1, 
            size = 2.5,
            inherit.aes = FALSE,
            family = "serif") +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red") +  
  scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, by = 0.1)) +  
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major.y = element_line(color = "gray90", size = 0.1)  
  ) +
  coord_cartesian(clip = "off") 

ggsave("plots/stability_protein_selection_rlasso_cv.pdf", width = 9, height = 8)

filtered_selection_df <- subset(filtered_selection_df, frequency > optimal_cutoff)

write.csv(filtered_selection_df, file.path(selection_dir, "selected_proteins.csv"))

olink_assay <- read.delim(file.path(input_dir, "olink_assay.txt"))
rmvd <- c("AMY2B", "CST1", "CTSS", "GLIPR1", "NPM1", "PCOLCE", "TACSTD2")
olink_assay <- subset(olink_assay, !(olink_assay$Assay %in% rmvd))
table(olink_assay$Panel)


by_panel <- filtered_selection_df
by_panel$Assay <- toupper(by_panel$feature)
by_panel <- merge(by_panel, olink_assay,
                  by.x = "Assay", by.y = "Assay", all.x = TRUE)

by_panel$Panel <- sapply(str_split(by_panel$Panel, " "),
                         function(x) x[1])

table(by_panel$Panel)

write.csv(by_panel, file.path(selection_dir, "selected_proteins_by_panel.csv"))
