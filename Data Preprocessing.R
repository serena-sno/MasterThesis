
### -------- DATA PRE-PROCESSING -------------- ###

# Library load
load_libraries <- function() {
  library(tidyverse) # for data manipulation
  library(missForest) # for data imputation for UKB data
  library(caret) # for correlation filteration
}

# Set working directory
set_working_directory <- function(directory_path) {
  setwd(dir = directory_path)
}

# Function to read data for Cohort 1 (UK Biobank)
read_cohort1_data <- function(features_file, diagnosis_file) {
  BBdata_metadata <- read.csv(features_file, header = TRUE, row.names = NULL)
  BBdiagnosis_metadata <- read.csv(diagnosis_file, header = TRUE, row.names = NULL)
  
  return(list(BBdata_metadata = BBdata_metadata, BBdiagnosis_metadata = BBdiagnosis_metadata))
}

# Function to preprocess Cohort 1 (UK Biobank) data
preprocess_cohort1_data <- function(BBdata_metadata, BBdiagnosis_metadata) {
  # Cleaning and Preprocessing
  # ---------------------------------
  
  # Subsetting the metabolite section
  my_string <- "Metabolite"
  my_cols <- grep(my_string, colnames(BBdata_metadata), value = TRUE)
  BBmetabolomics <- BBdata_metadata[, my_cols]
  rownames(BBmetabolomics) <- BBdata_metadata$f.eid
  
  # Formatting the names of the metabolites
  my_string <- ".0.0.Metabolites"
  new_colnames <- gsub(my_string, "", colnames(BBmetabolomics))
  colnames(BBmetabolomics) <- new_colnames
  
  # Cleaning unwanted objects
  rm(list = c("my_string", "my_cols", "new_colnames"))
  
  # Subset to keep only metabolite samples for UK Biobank
  max_null_percentage <- 50  # The metabolite samples are over 90% filled for all metabolites
  null_percentages <- rowMeans(is.na(BBmetabolomics)) * 100
  BBmetabolomics <- BBmetabolomics[null_percentages < max_null_percentage, ]
  
  # Formatting the Metadata
  BBmetadata <- BBdata_metadata[, c((ncol(BBdata_metadata) - 5):ncol(BBdata_metadata))]
  BBmetadata <- BBmetadata[BBmetadata$f.eid %in% rownames(BBmetabolomics), ]
  BBdiagnosis_metadata <- BBdiagnosis_metadata[BBdiagnosis_metadata$f.eid %in% rownames(BBmetabolomics), ]
  
  head(BBdiagnosis_metadata)
  count(distinct(BBdiagnosis_metadata, f.eid))  # How many are distinct?
  
  # Function to convert date format to a unified one of YYYY-MM-DD
  convert_date_format <- function(date) {
    if (grepl("^\\d{4}\\.\\d$", date)) {
      date <- str_replace(date, "\\.", "-0")
      date <- paste(date, "15", sep = "-")
    } else if (grepl("^\\d{4}\\.\\d{2}$", date)) {
      date <- str_replace(date, "\\.", "-")
      date <- paste(date, "15", sep = "-")
    } else if (grepl("^\\d{4}$", date)) {
      date <- paste(date, "06-15", sep = "-")
    }
    return(date)
  }
  
  # Applying the conversion function to the 'Date_value' column
  BBdiagnosis_metadata$Date_value <- sapply(BBdiagnosis_metadata$Date_value, convert_date_format)
  
  # Subset the dataframe to keep only the sample duplicates with the latest dates
  BBdiagnosis_metadata$Date_value <- as.Date(BBdiagnosis_metadata$Date_value)  # Convert to Date
  
  # Filtering to keep only the rows with the maximum Date_value within each f.eid group
  BBdiagnosis_metadata <- BBdiagnosis_metadata %>%
    group_by(f.eid) %>%
    arrange(desc(Date_value)) %>%
    slice(if (max(Date_value) %in% Date_value) 1 else 0)
  
  # Combining the diagnosis data with the metadata to form the new BBmetadata
  BBmetadata <- merge(BBmetadata, BBdiagnosis_metadata, by = 'f.eid', all = TRUE)
  rownames(BBmetadata) <- BBmetadata$f.eid
  BBmetadata <- select(BBmetadata, -f.eid)
  
  # Formatting BBmetadata columns
  BBmetadata$Ethnicity <- as.character(BBmetadata$Ethnicity)
  BBmetadata$Sex <- as.character(BBmetadata$Sex)
  BBmetadata$casecontrol <- as.character(BBmetadata$casecontrol)
  # change name of "casecontrol" to "Label" for MetaboAnalyst format
  names(BBmetadata)[names(BBmetadata) == "casecontrol"] <- "Label"
  
  # Formatting the BBdata_metadata
  rownames(BBdata_metadata) <- BBdata_metadata$f.eid
  BBdata_metadata <- BBdata_metadata[rownames(BBdata_metadata) %in% rownames(BBmetabolomics), ]
  
  # Subset to the Complete Cases of Diet Samples
  # --------------------------------------------
  
  # Dealing with the diet data
  my_string <- "Diet"
  my_cols <- grep(my_string, colnames(BBdata_metadata), value = TRUE)
  BBdiet <- BBdata_metadata[, my_cols]
  
  # Formatting the names of the metabolites
  my_string <- ".Diet...0.0" 
  new_colnames <- gsub(my_string, "", colnames(BBdiet))
  colnames(BBdiet) <- new_colnames
  
  # filter out columns not to be used (have numerous NAs {type} or hard to quantify)
  columns_to_remove <- grepl("type|Hot|Variation|Major|questionnaires|Never", colnames(BBdiet))
  BBdiet.feat <- BBdiet[, !columns_to_remove]
  
  # filter out incomplete samples across the remaining columns
  BBdiet.feat.samp <- BBdiet.feat[complete.cases(BBdiet.feat), ]
  rows_to_remove <- apply(BBdiet.feat.samp == -3 | BBdiet.feat.samp == -1, 1, any)
  BBdiet.feat.samp <- BBdiet.feat.samp[!rows_to_remove, ]
  
  # subset the metadata and metabolite to the common samples
  BBmetabolomics <- BBmetabolomics[rownames(BBmetabolomics) %in% rownames(BBdiet.feat.samp), ]
  BBmetadata <- BBmetadata[rownames(BBmetadata) %in% rownames(BBdiet.feat.samp), ]
  
  # Removing non-baseline IBD individuals
  # ------------------------------------------
  
  # Adding a column of difference between age at diagnosis and age at measurement
  BBmetadata <- BBmetadata %>%
    mutate(difference = ifelse(!is.na(Age_value) & !is.na(Age.at.measurement), Age.at.measurement - Age_value, NA))
  
  # Storing the nonIBD sample names
  nonIBD_sample_names <- rownames(BBmetadata[which(BBmetadata$Label == "0"), ])
  
  # Keeping IBD samples that had a history of IBD at baseline
  IBD_sample_names_to_keep <- rownames(BBmetadata[which(!is.na(BBmetadata$difference) & BBmetadata$difference >= 0), ])
  
  # Combining the two sets of sample names
  all_sample_names_to_keep <- union(nonIBD_sample_names, IBD_sample_names_to_keep)
  
  # Subsetting BBmetabolomics and BBmetadata and BBdiet based on the combined sample names
  BBmetabolomics <- BBmetabolomics[rownames(BBmetabolomics) %in% all_sample_names_to_keep, ]
  BBmetadata <- BBmetadata[rownames(BBmetadata) %in% all_sample_names_to_keep, ]
  BBdiet.feat.samp <- BBdiet.feat.samp[rownames(BBdiet.feat.samp) %in% all_sample_names_to_keep, ]
  
  # Continue editting the Diet Data
  # -------------------------------
  
  BBdiet_numerical_features_vector <- c(
    'Cooked_vegetable_intake', 
    'Salad_._raw_vegetable_intake',  
    'Fresh_fruit_intake', 
    'Dried_fruit_intake', 
    'Bread_intake', 
    'Cereal_intake', 
    'Tea_intake', 
    'Coffee_intake', 
    'Water_intake' 
  )
  BBdiet_categorical_features_vector <- colnames(BBdiet.feat.samp)[!(colnames(BBdiet.feat.samp) %in% 
                                                                       BBdiet_numerical_features_vector)]
  BBdiet_numerical <- BBdiet.feat.samp[, BBdiet_numerical_features_vector]
  BBdiet_categorical <- BBdiet.feat.samp[, BBdiet_categorical_features_vector]
  # replacing -10 with a numerical average equivalence of "less than one"
  BBdiet_numerical <- replace(BBdiet_numerical, (BBdiet_numerical == "-10"), 0.5)
  
  # Robust Scaling function
  robust_scale <- function(x) {
    (x - median(x)) / IQR(x)
  }
  
  BBdiet_numerical_robust_scaled <- as.data.frame(apply(BBdiet_numerical, 2, robust_scale))
  
  # to replace the the values of some categories
  value_mapping <- list(
    `0` = 0, 
    `1` = 0.07, 
    `2` = 0.14, 
    `3` = 0.43, 
    `4` = 0.79, 
    `5` = 1
  )
  # Replace values in all columns except "salt added to food"
  BBdiet_categorical_replaced <- BBdiet_categorical %>%
    mutate_at(vars(-Salt_added_to_food), ~ as.numeric(value_mapping[as.character(.)]))
  
  # Replace values in "salt added to food" column
  salt_mapping <- list(
    `1` = 0, 
    `2` = 0.28, 
    `3` = 0.71, 
    `4` = 1
  )
  BBdiet_categorical_replaced$Salt_added_to_food <- as.numeric(salt_mapping[as.character(BBdiet_categorical$Salt_added_to_food)])
  BBdiet_categorical_replaced
  
  BBdiet_final <- cbind(BBdiet_numerical_robust_scaled, BBdiet_categorical_replaced)
  
  return(list(BBmetabolomics = BBmetabolomics, BBmetadata = BBmetadata, BBdiet = BBdiet_final))
  
}

# Function to impute, standardize, and transform Cohort 1 (UK Biobank) data
impute_standardize_transform_cohort1 <- function(BBmetabolomics) {
  # Calculate the percentage of missing values for each column in BBmetabolomics
  BBmetabolomics_missing_perc <- mean(is.na(BBmetabolomics)) * 100
  
  # Impute missing values using missForest
  BBmetabolomics.imp <- missForest(BBmetabolomics, maxiter = 10, ntree = 100)
  BBmetabolomics.imp <- BBmetabolomics.imp$ximp

  # Perform Pareto scaling (column-wise)
  BBmetabolomics.imp.sc <- scale(BBmetabolomics.imp, center = FALSE, scale = sqrt(apply(BBmetabolomics.imp, 2, var)))
  
  # Apply log transformation after shifting by 1
  BBmetabolomics.imp.sc.trans <- as.data.frame(log(BBmetabolomics.imp.sc + abs(min(BBmetabolomics.imp.sc) + 1)))
  
  # Feature Selection
  # Remove columns that are a sum of other columns to avoid collinearity and 
  # to help the correlation remover make better decisions
  # i am keeping the individuals to help in pathway analysis
  
  # Find column indices containing the string "Total"
  columns_to_remove <- grep("Total", colnames(BBmetabolomics.imp.sc.trans))
  
  # Remove columns with "Total" in their names
  BBmetabolomics.imp.sc.trans.filtered <- BBmetabolomics.imp.sc.trans[, -columns_to_remove]
  
  # Remove highly correlated features
  # use spearman because of the possibly non linear nature of the data
  cor_mat <- cor(BBmetabolomics.imp.sc.trans.filtered, method = "spearman")
  
  highly_corr_feat <- findCorrelation(cor_mat, cutoff=0.9)
  highly_corr_feat = sort(highly_corr_feat)
  
  removed_highly_cor_feat_data <- BBmetabolomics.imp.sc.trans.filtered[, -c(highly_corr_feat)]
  return(list(BBprocessed = BBmetabolomics.imp.sc.trans,
              BBprocessed_filtered = removed_highly_cor_feat_data,
              BBmissingno = BBmetabolomics_missing_perc))
}

# Function to read Cohort 2 (HMP2 data)
read_cohort2_data <- function(hmp2_metabolomics_file, hmp2_metadata_file) {
  HMPmetabolomics <- read.csv(hmp2_metabolomics_file, header = TRUE, row.names = NULL)
  HMPmetadata <- read.csv(hmp2_metadata_file, header = TRUE, row.names = NULL)
  
  return(list(HMPmetabolomics = HMPmetabolomics, 
              HMPmetadata = HMPmetadata))
}

# Function to preprocess Cohort 2 (HMP2 data)
preprocess_cohort2_data <- function(HMPmetabolomics, HMPmetadata) {
  
  # Preprocess the Metadata
  # -------------------------------
  
  # Change name of "diagnosis" to "Label" for MetaboAnalyst format
  names(HMPmetadata)[names(HMPmetadata) == "diagnosis"] <- "Label"
  
  # Get the complete cases of the HMPdiet
  HMPdiet <- HMPmetadata[, c(which(colnames(HMPmetadata) == "External.ID"), 72:81, 83:92, 101)]
  HMPdiet[HMPdiet == ""] <- NA
  HMPdiet <- HMPdiet[complete.cases(HMPdiet), ]
  HMPmetadata <- HMPmetadata[(HMPmetadata$External.ID %in% HMPdiet$External.ID), ]
  
  # Remove non-baseline IBD individuals
  
  HMPmetadata <- HMPmetadata %>% 
    mutate(difference = ifelse(!is.na(consent_age) & !is.na(Age.at.diagnosis), consent_age - Age.at.diagnosis, NA))
  
  # Storing the nonIBD sample names
  nonIBD_sample_names <- HMPmetadata[which(HMPmetadata$Label == "nonIBD"), "External.ID"]
  
  # Keeping IBD samples that had a history of IBD at baseline
  IBD_sample_names_to_keep <- HMPmetadata[which(!is.na(HMPmetadata$difference) & HMPmetadata$difference >= 0), "External.ID"]
  
  # Combining the two sets of sample names
  all_sample_names_to_keep <- union(nonIBD_sample_names, IBD_sample_names_to_keep)
  
  # Subset HMPmetadata based on the combined sample names
  HMPmetadata <- HMPmetadata[HMPmetadata$External.ID %in% all_sample_names_to_keep, ]
  
  # Sub HMP metadata to the ones with the metabolites -----
  HMPmetadata <- HMPmetadata %>% filter(data_type %in% "metabolomics")
  rownames(HMPmetadata) <- HMPmetadata$External.ID
  HMPmetadata$External.ID <- NULL
 
  # Edit the metadata to fill in the NA sections
  HMPmetadata <- replace(HMPmetadata, is.na(HMPmetadata), "NA")
  HMPmetadata <- replace(HMPmetadata, HMPmetadata == "", "NA")
  
  # Pre-process the Metabolomics Data
  # ----------------------------------
  
  # Select one metabolomic method (HILIC-pos) in HMPmetabolomics
  HMPmetabolomics <- subset(HMPmetabolomics, Method == "HILIC-pos")
  
  # Subset HMPmetabolomics to only the metabolites that exist
  HMPmetabolomics <- subset(HMPmetabolomics, Metabolite != "")
  
  # Keep only the samples and the metabolite columns in HMPmetabolomics
  HMPmetabolomics <- HMPmetabolomics[, -c(1:5,7), drop = FALSE]
  
  # Remove duplicated rows (metabolites), keeping the one with fewer NA values
  na_counts <- rowSums(is.na(HMPmetabolomics))
  duplicated_rows <- duplicated(HMPmetabolomics$Metabolite) | duplicated(HMPmetabolomics$Metabolite, fromLast = TRUE)
  keep_rows <- !duplicated_rows | (duplicated_rows & na_counts == min(na_counts[duplicated_rows]))
  HMPmetabolomics <- HMPmetabolomics[keep_rows, ]
  
  # Set row names and drop the metabolite column
  row.names(HMPmetabolomics) <- HMPmetabolomics$Metabolite
  HMPmetabolomics <- HMPmetabolomics %>% select(-Metabolite)
  
  # Transpose HMPmetabolomics to fit the form of BBmetabolomics (samples as rows, metabolites as columns)
  HMPmetabolomics <- as.data.frame(t(HMPmetabolomics))
  
  ## ---> INSERT FILTERING STEP HERE (E.G. SAMPLES THAT HAVE >50% MISSING)
  ## This metabolomic dataset has about <5% missing values, so no need for filtering
  
  # Sub HMPmetabolomics to available samples in HMPmetadata
  HMPmetabolomics <- HMPmetabolomics[rownames(HMPmetabolomics) %in% rownames(HMPmetadata), ]
  
  # In case there was some filtering step in metabolomics that took out more samples
  HMPmetadata <- HMPmetadata[rownames(HMPmetadata) %in% rownames(HMPmetabolomics), ]

  
  # Continue editing the Diet Data
  # ------------------------------
  
  # remove duplicated rows
  HMPdiet <- HMPdiet[!duplicated(HMPdiet), ]
  
  # Sub the diet data to the final samples in the metadata
  HMPdiet <- HMPdiet[HMPdiet$External.ID %in% rownames(HMPmetadata), ]
  rownames(HMPdiet) <- HMPdiet$External.ID
  HMPdiet$External.ID <- NULL
  
  value_mapping <- list(
    `No, I did not consume these products in the last 7 days` = 0,
    `Within the past 4 to 7 days` = 0.2,
    `Within the past 2 to 3 days` = 0.58,
    `Yesterday, 1 to 2 times` = 0.9,
    `Yesterday, 3 or more times` = 1
  )
  
  HMPdiet_replaced <- HMPdiet %>%
    mutate_at(vars(everything()), ~ as.numeric(value_mapping[as.character(.)]))
  
  # Edit column names to be shorter
  # Extract the name before the ".."
  new_col_names <- sub("\\.\\..*", "", colnames(HMPdiet))
  
  # Define patterns and replacements
  patterns <- c("Yogurt.or.other.foods.containing.active.bacterial.cultures",
                "Tea.or.coffee.no.sugar.and.no.sugar.replacement")
  replacements <- c("Yoghurt.or.other", "Tea.or.coffee.no.sugar")
  
  # Loop through patterns and replacements
  short_col_names <- new_col_names  # Initialize with new_col_names
  for (i in seq_along(patterns)) {
    short_col_names <- str_replace(short_col_names, patterns[i], replacements[i])
  }
  
  # Replace original column names with the short ones
  colnames(HMPdiet_replaced) <- short_col_names
  colnames(HMPdiet) <- short_col_names
  
  # Trim down the metadata
  # ---------------------------
  
  # Define the threshold for removing columns with one unique value
  unique_threshold <- 1  # All values are the same
  # Define the threshold for removing columns with more than 50% "NA" values
  NA_threshold <- 0.5  # More than 50% "NA"
  # Combine all removal operations in a single pipeline
  HMPmetadata_selected_names <- HMPmetadata %>%
    select(-where(~ n_distinct(.) <= unique_threshold)) %>%
    select(-where(~ mean(. == "NA", na.rm = TRUE) > NA_threshold))
  
  # DIET IDs to remove
  dietIDs_to_remove <- c(17:26, 28:37, 46)
  
  # Unnecessary columns to remove
  unnecessary_columns_to_remove <- c(1:10, 55, 57:63)
  
  HMPmetadata_selected_names <- HMPmetadata_selected_names %>%
    select(-all_of(union(dietIDs_to_remove, unnecessary_columns_to_remove)))
  
  HMPmetadata_names <- names(HMPmetadata_selected_names)
  
  return(list(HMPdiet = HMPdiet_replaced,
              HMPmetabolomics = HMPmetabolomics, 
              HMPmetadata = HMPmetadata_selected_names))
}

# Function to impute, standardize, and transform Cohort 2 (HMP2 data) data
impute_standardize_transform_cohort2 <- function(HMPmetabolomics) {
  # Calculate missing and zero percentages
  HMPmetabolomics_missing_perc <- mean(is.na(HMPmetabolomics)) * 100
  HMPmetabolomics_zero_perc <- mean(HMPmetabolomics == "0", na.rm = TRUE) * 100
  
  # Change zero values with half of the minimum positive value
  replacezero <- function(x) "[<-"(x, !x | is.na(x), min(x[x > 0], na.rm = TRUE) / 2)
  HMPmetabolomics.imp <- as.data.frame(apply(HMPmetabolomics, 2, replacezero))
  
  # Perform scaling and transform
  HMPmetabolomics.imp.sc <- scale(HMPmetabolomics.imp, center = FALSE, scale = sqrt(apply(HMPmetabolomics.imp, 2, var)))
  HMPmetabolomics.imp.sc.trans <- as.data.frame(log(HMPmetabolomics.imp.sc + abs(min(HMPmetabolomics.imp.sc)) + 1))
  
  # Remove highly collinear features
  cor_mat <- cor(HMPmetabolomics.imp.sc.trans, method = "spearman")
  
  highly_corr_feat <- findCorrelation(cor_mat, cutoff=0.9)
  highly_corr_feat = sort(highly_corr_feat)
  
  removed_highly_cor_feat_data <- HMPmetabolomics.imp.sc.trans[, -c(highly_corr_feat)]
  return(list(HMPprocessed = HMPmetabolomics.imp.sc.trans,
              HMPprocessed_filtered = removed_highly_cor_feat_data,
              HMPmissingno = HMPmetabolomics_missing_perc,
              HMPzerono = HMPmetabolomics_zero_perc))
}

# Function to create datasets from intersecting metabolites across both Cohorts
create_common_datasets <- function(HMPmetabolomics, BBmetabolomics) {
  # Clean: Change the metabolites to same case
  colnames(HMPmetabolomics) <- tolower(colnames(HMPmetabolomics))
  colnames(BBmetabolomics) <- tolower(colnames(BBmetabolomics))
  
  # HMP subset to only the common metabolites
  HMPmetabolomics_common <- HMPmetabolomics[, colnames(HMPmetabolomics)[which(colnames(HMPmetabolomics) %in% colnames(BBmetabolomics))]]
  
  # Biobank subset to only the common metabolites
  BBmetabolomics_common <- BBmetabolomics[, colnames(BBmetabolomics)[which(colnames(BBmetabolomics) %in% colnames(HMPmetabolomics))]]
  
  return(list(BBmetabolomics_common, HMPmetabolomics_common))
}

# Function to impute, standardize, and transform the common datasets
impute_standardize_transform_common <- function(BBmetabolomics_common, HMPmetabolomics_common) {
  
  # UK BB common with HMP
  BBmetabolomics_common_imp <- missForest(BBmetabolomics_common, maxiter = 10, ntree = 100)
  BBmetabolomics_common_imp <- BBmetabolomics_common_imp$ximp
  BBmetabolomics_common_imp_sc <- scale(BBmetabolomics_common_imp, center = FALSE, scale = sqrt(apply(BBmetabolomics_common_imp, 2, var)))
  BBmetabolomics_common_imp_sc_trans <- as.data.frame(log(BBmetabolomics_common_imp_sc + abs(min(BBmetabolomics_common_imp_sc) + 1)))
  
  
  # HMP common with UK BB
  replacezero <- function(x) "[<-"(x, !x | is.na(x), min(x[x > 0], na.rm = TRUE) / 2)
  HMPmetabolomics_common_imp <- as.data.frame(apply(HMPmetabolomics_common, 2, replacezero))
  HMPmetabolomics_common_imp_sc <- scale(HMPmetabolomics_common_imp, center = FALSE, scale = sqrt(apply(HMPmetabolomics_common_imp, 2, var)))
  HMPmetabolomics_common_imp_sc_trans <- as.data.frame(log(HMPmetabolomics_common_imp_sc + abs(min(HMPmetabolomics_common_imp_sc)) + 1))
  
  return(list(BBmetabolomics_common_imp_sc_trans, 
              HMPmetabolomics_common_imp_sc_trans))
}

# Function to format labels for both Cohort 1 and Cohort 2 data
format_labels <- function(cohort_metadata) {
  # Define a function for find and replace
  find_replace <- function(x) {
    x <- gsub("CD", "IBD", x)
    x <- gsub("UC", "IBD", x)
    x <- gsub("1", "IBD", x)
    x <- gsub("0", "nonIBD", x)
    return(x)
  }
  
  # Change the IBD and nonIBD labels to a unified form
  cohort_metadata$Label <- find_replace(cohort_metadata$Label)
  
  # Retrieve the label column from the Metadata
  label_data <- cohort_metadata %>% select(Label)
  
  return(list(label_data, cohort_metadata))
}

# Function to merge the data with the labels
merge_data_with_label <- function(BBmetabolomics_label, BBmetabolomics_imp_sc_trans_filtered, BBmetabolomics_imp_sc_trans,
                                  HMPmetabolomics_label, HMPmetabolomics_imp_sc_trans_filtered, HMPmetabolomics_imp_sc_trans,
                                  BBmetabolomics_HMPcommon_imp_sc_trans, HMPmetabolomics_BBcommon_imp_sc_trans) {
  
  # For Pathway Analysis
  BBmetabolomics.imp.sc.trans_with_label <- merge(BBmetabolomics_label, BBmetabolomics_imp_sc_trans, by = "row.names")
  HMPmetabolomics.imp.sc.trans_with_label <- merge(HMPmetabolomics_label, HMPmetabolomics_imp_sc_trans, by = "row.names")
  
  # For Machine Learning and Exploration
  BBmetabolomics.imp.sc.trans.filt_with_label <- merge(BBmetabolomics_label, BBmetabolomics_imp_sc_trans_filtered, by = "row.names")
  HMPmetabolomics.imp.sc.trans.filt_with_label <- merge(HMPmetabolomics_label, HMPmetabolomics_imp_sc_trans_filtered, by = "row.names")
  BBmetabolomics.HMPcommon.imp.sc.trans_with_label <- merge(BBmetabolomics_label, BBmetabolomics_HMPcommon_imp_sc_trans, by = "row.names")
  HMPmetabolomics.BBcommon.imp.sc.trans_with_label <- merge(HMPmetabolomics_label, HMPmetabolomics_BBcommon_imp_sc_trans, by = "row.names")
 
  # Name the row names to "Sample" to accommodate MetaboAnalyst or MicrobiomeAnalyst format
  names(BBmetabolomics.imp.sc.trans_with_label)[names(BBmetabolomics.imp.sc.trans_with_label) == "Row.names"] <- "Sample"
  names(HMPmetabolomics.imp.sc.trans_with_label)[names(HMPmetabolomics.imp.sc.trans_with_label) == "Row.names"] <- "Sample"
  names(BBmetabolomics.imp.sc.trans.filt_with_label)[names(BBmetabolomics.imp.sc.trans.filt_with_label) == "Row.names"] <- "Sample"
  names(HMPmetabolomics.imp.sc.trans.filt_with_label)[names(HMPmetabolomics.imp.sc.trans.filt_with_label) == "Row.names"] <- "Sample"
  names(BBmetabolomics.HMPcommon.imp.sc.trans_with_label)[names(BBmetabolomics.HMPcommon.imp.sc.trans_with_label) == "Row.names"] <- "Sample"
  names(HMPmetabolomics.BBcommon.imp.sc.trans_with_label)[names(HMPmetabolomics.BBcommon.imp.sc.trans_with_label) == "Row.names"] <- "Sample"
  
  # Return the merged datasets with labels
  return(list(BBmetabolomics_imp_sc_trans_with_label = BBmetabolomics.imp.sc.trans_with_label,
              HMPmetabolomics_imp_sc_trans_with_label = HMPmetabolomics.imp.sc.trans_with_label,
              BBmetabolomics_imp_sc_trans_filt_with_label = BBmetabolomics.imp.sc.trans.filt_with_label,
              HMPmetabolomics_imp_sc_trans_filt_with_label = HMPmetabolomics.imp.sc.trans.filt_with_label,
              BBmetabolomics_HMPcommon_imp_sc_trans_with_label = BBmetabolomics.HMPcommon.imp.sc.trans_with_label,
              HMPmetabolomics_BBcommon_imp_sc_trans_with_label = HMPmetabolomics.BBcommon.imp.sc.trans_with_label))
}

# Function to export processed data
format_and_export_data <- function(data_list, BBmetadata, HMPmetadata, BBdiet, HMPdiet, output_directory) {
  
  BBmetabolomics.imp.sc.trans_with_label_sample <- data_list$BBmetabolomics_imp_sc_trans_with_label
  HMPmetabolomics.imp.sc.trans_with_label_sample <- data_list$HMPmetabolomics_imp_sc_trans_with_label
  BBmetabolomics.imp.sc.trans.filt_with_label_sample <- data_list$BBmetabolomics_imp_sc_trans_filt_with_label
  HMPmetabolomics.imp.sc.trans.filt_with_label_sample <- data_list$HMPmetabolomics_imp_sc_trans_filt_with_label
  BBmetabolomics.HMPcommon.imp.sc.trans_with_label_sample <- data_list$BBmetabolomics_HMPcommon_imp_sc_trans_with_label
  HMPmetabolomics.BBcommon.imp.sc.trans_with_label_sample <- data_list$HMPmetabolomics_BBcommon_imp_sc_trans_with_label
 
  # Transpose to accommodate MetaboAnalyst format
  BBmetabolomics.imp.sc.trans_with_label.T <- as.data.frame(t(BBmetabolomics.imp.sc.trans_with_label_sample))
  HMPmetabolomics.imp.sc.trans_with_label.T <- as.data.frame(t(HMPmetabolomics.imp.sc.trans_with_label_sample))
  BBmetabolomics.imp.sc.trans.filt_with_label.T <- as.data.frame(t(BBmetabolomics.imp.sc.trans.filt_with_label_sample))
  HMPmetabolomics.imp.sc.trans.filt_with_label.T <- as.data.frame(t(HMPmetabolomics.imp.sc.trans.filt_with_label_sample))
  BBmetabolomics.HMPcommon.imp.sc.trans_with_label.T <- as.data.frame(t(BBmetabolomics.HMPcommon.imp.sc.trans_with_label_sample))
  HMPmetabolomics.BBcommon.imp.sc.trans_with_label.T <- as.data.frame(t(HMPmetabolomics.BBcommon.imp.sc.trans_with_label_sample))
 
  # Export the processed data for MetaboAnalyst to CSV files
  write.table(BBmetabolomics.imp.sc.trans_with_label.T, file = file.path(output_directory, "For MetaboAnalyst/BB_imp_sc_trans_label_T.csv"), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
  write.table(HMPmetabolomics.imp.sc.trans_with_label.T, file = file.path(output_directory, "For MetaboAnalyst/HMP_imp_sc_trans_label_T.csv"), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
  write.table(BBmetabolomics.imp.sc.trans.filt_with_label.T, file = file.path(output_directory, "For MetaboAnalyst/BB_imp_sc_trans_filt_label_T.csv"), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
  write.table(HMPmetabolomics.imp.sc.trans.filt_with_label.T, file = file.path(output_directory, "For MetaboAnalyst/HMP_imp_sc_trans_filt_label_T.csv"), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
  write.table(BBmetabolomics.HMPcommon.imp.sc.trans_with_label.T, file = file.path(output_directory, "For MetaboAnalyst/BB_HMP_common_imp_sc_trans_label_T.csv"), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
  write.table(HMPmetabolomics.BBcommon.imp.sc.trans_with_label.T, file = file.path(output_directory, "For MetaboAnalyst/HMP_BB_common_imp_sc_trans_label_T.csv"), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
 
  # Set the rownames as the Sample IDs
  rownames(BBmetabolomics.imp.sc.trans_with_label_sample) <- BBmetabolomics.imp.sc.trans_with_label_sample$Sample
  BBmetabolomics.imp.sc.trans_with_label_sample$Sample <- NULL
  rownames(HMPmetabolomics.imp.sc.trans_with_label_sample) <- HMPmetabolomics.imp.sc.trans_with_label_sample$Sample
  HMPmetabolomics.imp.sc.trans_with_label_sample$Sample <- NULL
  rownames(BBmetabolomics.imp.sc.trans.filt_with_label_sample) <- BBmetabolomics.imp.sc.trans.filt_with_label_sample$Sample
  BBmetabolomics.imp.sc.trans.filt_with_label_sample$Sample <- NULL
  rownames(HMPmetabolomics.imp.sc.trans.filt_with_label_sample) <- HMPmetabolomics.imp.sc.trans.filt_with_label_sample$Sample
  HMPmetabolomics.imp.sc.trans.filt_with_label_sample$Sample <- NULL
  rownames(BBmetabolomics.HMPcommon.imp.sc.trans_with_label_sample) <- BBmetabolomics.HMPcommon.imp.sc.trans_with_label_sample$Sample
  BBmetabolomics.HMPcommon.imp.sc.trans_with_label_sample$Sample <- NULL
  rownames(HMPmetabolomics.BBcommon.imp.sc.trans_with_label_sample) <- HMPmetabolomics.BBcommon.imp.sc.trans_with_label_sample$Sample
  HMPmetabolomics.BBcommon.imp.sc.trans_with_label_sample$Sample <- NULL 
  
  # Export the normal data
  write.csv(BBmetabolomics.imp.sc.trans_with_label_sample, file = file.path(output_directory, "BB_imp_sc_trans_label.csv"), row.names = TRUE)
  write.csv(HMPmetabolomics.imp.sc.trans_with_label_sample, file = file.path(output_directory, "HMP_imp_sc_trans_label.csv"), row.names = TRUE)
  write.csv(BBmetabolomics.imp.sc.trans.filt_with_label_sample, file = file.path(output_directory, "BB_imp_sc_trans_filt_label.csv"), row.names = TRUE)
  write.csv(HMPmetabolomics.imp.sc.trans.filt_with_label_sample, file = file.path(output_directory, "HMP_imp_sc_trans_filt_label.csv"), row.names = TRUE)
  write.csv(BBmetabolomics.HMPcommon.imp.sc.trans_with_label_sample, file = file.path(output_directory, "BB_HMPcommon_imp_sc_trans_label.csv"), row.names = TRUE)
  write.csv(HMPmetabolomics.BBcommon.imp.sc.trans_with_label_sample, file = file.path(output_directory, "HMP_BBcommon_imp_sc_trans_label.csv"), row.names = TRUE)
 
  # Export the metadata
  write.csv(BBmetadata[, c(1,2,4,5)], file = file.path(output_directory, "BBmetadata.csv"), row.names = TRUE)
  write.csv(HMPmetadata, file = file.path(output_directory, "HMPmetadata.csv"), row.names = TRUE)
  
  # Export the Diet Data
  write.csv(BBdiet, file = file.path(output_directory, "BBdiet.csv"), row.names = TRUE)
  write.csv(HMPdiet, file = file.path(output_directory, "HMPdiet.csv"), row.names = TRUE)
}



#### ------------ BEGIN RUNNING -----------------
#### --------------------------------------------


load_libraries()
# >> input working directory path <<
set_working_directory("path/to/Masters Thesis/")

# Step 1: Read data for Cohort 1 (UK Biobank)
cohort1_data <- read_cohort1_data("Cohort Data/UKBiobank Data/FeaturesData.csv", "Cohort Data/UKBiobank Data/DiagnosisData.csv")
BBdata_metadata <- cohort1_data$BBdata_metadata
BBdiagnosis_metadata <- cohort1_data$BBdiagnosis_metadata

# Step 2: Preprocess data for Cohort 1 (UK Biobank)
cohort1_processed_data <- preprocess_cohort1_data(BBdata_metadata, BBdiagnosis_metadata)
BBmetabolomics <- cohort1_processed_data$BBmetabolomics
BBdiet <- cohort1_processed_data$BBdiet
BBmetadata <- cohort1_processed_data$BBmetadata

# Step 3: Impute, standardize, and transform data for Cohort 1 (UK Biobank)
cohort1_transformed_data <- impute_standardize_transform_cohort1(BBmetabolomics)
BBmetabolomics_missing_perc <- cohort1_transformed_data$BBmissingno
BBmetabolomics_imp_sc_trans <- cohort1_transformed_data$BBprocessed
columns_to_remove <- grep("Total|Average|Degree", colnames(BBmetabolomics_imp_sc_trans))
BBmetabolomics_imp_sc_trans <- BBmetabolomics_imp_sc_trans[, -columns_to_remove]
BBmetabolomics_imp_sc_trans_filtered <- cohort1_transformed_data$BBprocessed_filtered
columns_to_remove <- grep("Total|Average|Degree", colnames(BBmetabolomics_imp_sc_trans_filtered))
BBmetabolomics_imp_sc_trans_filtered <- BBmetabolomics_imp_sc_trans_filtered[, -columns_to_remove]


# Step 4: Read in data for Cohort 2 (HMP2 data)
cohort2_data <- read_cohort2_data("Cohort Data/HMP2DB Data/iHMP_metabolomics.csv", "Cohort Data/HMP2DB Data/hmp2_metadata.csv")
HMPmetabolomics <- cohort2_data$HMPmetabolomics
HMPmetadata <- cohort2_data$HMPmetadata

# Step 5: Preprocess data for Cohort 2 (HMP2 data)
cohort2_processed_data <- preprocess_cohort2_data(HMPmetabolomics, HMPmetadata)
HMPmetabolomics <- cohort2_processed_data$HMPmetabolomics
HMPdiet <- cohort2_processed_data$HMPdiet
HMPdiet <- HMPdiet[, -which(colnames(HMPdiet) == "Probiotic")]
HMPmetadata <- cohort2_processed_data$HMPmetadata 

# Step 6: Impute, standardize, and transform metabolomics data for Cohort 2 (HMP2 data)
cohort2_transformed_data <- impute_standardize_transform_cohort2(HMPmetabolomics)
HMPmetabolomics_missing_perc <- cohort2_transformed_data$HMPmissingno
HMPmetabolomics_zero_perc <- cohort2_transformed_data$HMPzerono
HMPmetabolomics_imp_sc_trans <- cohort2_transformed_data$HMPprocessed
HMPmetabolomics_imp_sc_trans_filtered <- cohort2_transformed_data$HMPprocessed_filtered

# Step 7: Subset datasets based on common metabolites
common_datasets <- create_common_datasets(HMPmetabolomics_imp_sc_trans, BBmetabolomics_imp_sc_trans)
BBmetabolomics_common <- common_datasets[[1]]
HMPmetabolomics_common <- common_datasets[[2]]

# Step 8: Impute, standardize, and transform combined datasets
processed_common_data <- impute_standardize_transform_common(BBmetabolomics_common, HMPmetabolomics_common)
BBmetabolomics_common_imp_sc_trans <- processed_common_data[[1]]
HMPmetabolomics_common_imp_sc_trans <- processed_common_data[[2]]

# Step 9: Format labels for both Cohort 1 and Cohort 2 data
BBmetabolomics_label <- format_labels(BBmetadata)[[1]]
BBmetadata <- format_labels(BBmetadata)[[2]]
HMPmetabolomics_label <- format_labels(HMPmetadata)[[1]]
HMPmetadata <- format_labels(HMPmetadata)[[2]]

# Step 10: Merge processed data with labels
data_with_labels <- merge_data_with_label(BBmetabolomics_label, BBmetabolomics_imp_sc_trans_filtered, BBmetabolomics_imp_sc_trans,
                                          HMPmetabolomics_label, HMPmetabolomics_imp_sc_trans_filtered, HMPmetabolomics_imp_sc_trans,
                                          BBmetabolomics_common_imp_sc_trans, HMPmetabolomics_common_imp_sc_trans)

# Step 11: Format into Metabo/MicrobiomeAnalyst Format and export processed data
format_and_export_data(data_with_labels, BBmetadata, HMPmetadata, BBdiet, HMPdiet, "R Processed Data/")
