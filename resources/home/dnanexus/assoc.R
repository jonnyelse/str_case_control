#!/usr/bin/env Rscript
library(purrr)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
AC_min_threshold <- as.numeric(args[1])

samples_list <- read.table("samples_list.txt", stringsAsFactors = F)$V1
covar_pheno <- read.table("covar_pheno.tsv", header=T, stringsAsFactors = F)
covar_pheno$IID <- as.factor(covar_pheno$IID)
covar_names <- setdiff(colnames(covar_pheno), c("disease_status", "FID","IID"))
cases_covar <- covar_pheno %>% filter(disease_status == 1)
cases_ids <- cases_covar$IID
controls_covar <- covar_pheno %>% filter(disease_status == 0)
controls_ids <- controls_covar$IID

custom_round <- function(x, significant_digits = 4) {
  if (length(x) > 1) {
    return(sapply(x, custom_round, significant_digits = significant_digits))  # Apply recursively for vectors
  }
  if (x == 0 || is.na(x)) return(x)  # Handle zeros or NA
  if (abs(x) >= 1) {
    return(round(x, digits = significant_digits))  # Fixed decimal places for large numbers
  }
  digits <- -floor(log10(abs(x))) + significant_digits  # Dynamically adjust for small numbers
  round(x, digits = digits)
}

low_cols <- function(locus_data, alt_colnames){
  cases_locus_data <- locus_data %>% filter(IID %in% cases_ids)
  cases_locus_data$all_alleles <- rowSums(cases_locus_data[alt_colnames])
  all_allele_count_cases <- sum(cases_locus_data$all_alleles, na.rm = T)
  no_refs_cases <- 2 * sum(!is.na(cases_locus_data$all_alleles)) - all_allele_count_cases
  
  controls_locus_data <- locus_data %>% filter(IID %in% controls_ids)
  controls_locus_data$all_alleles <- rowSums(controls_locus_data[alt_colnames])
  all_allele_count_controls <- sum(controls_locus_data$all_alleles, na.rm = T)
  no_refs_controls <- 2 * sum(!is.na(controls_locus_data$all_alleles)) - all_allele_count_controls
  if(no_refs_cases < AC_min_threshold | no_refs_controls < AC_min_threshold) {
    return("low_ref_count")
  }
  cols_to_combine <- c()
  datasets <- list(cases_locus_data, controls_locus_data)
  cols_AC <- data.frame(casesAC=rep(NA, length(alt_colnames)), controlsAC=rep(NA, length(alt_colnames)))
  min_to_keep <- 1e20
  for(i in seq(length(alt_colnames))){
    for(j in seq(2)){
      allele_num <- 2 * sum(!is.na(datasets[[j]][[alt_colnames[i]]]))
      lowAC <- sum(datasets[[j]][[alt_colnames[i]]], na.rm = T)
      highAC <- allele_num - lowAC
      minAC <- min(c(lowAC, highAC))
      cols_AC[i,j] <- minAC
      if(minAC < AC_min_threshold){
        cols_to_combine <- c(cols_to_combine,alt_colnames[i])
      }
    }
  }
  if(length(cols_to_combine) == 0){
    return(NULL)
  }
  else if (length(unique(cols_to_combine)) == length(alt_colnames)){
    cols_to_combine <- unique(cols_to_combine)
    sorted_cols <- cols_to_combine[order(as.numeric(sub(".*_", "", cols_to_combine)))]
    return(sorted_cols)
  }
  else{
    cols_AC_filtered <- cols_AC %>%
      mutate(hiAC = if_else(casesAC >= AC_min_threshold & controlsAC >= AC_min_threshold, "Y", "N"))
    min_row <- which.min(ifelse(cols_AC_filtered$hiAC == "Y", 
                                pmin(cols_AC_filtered$casesAC, cols_AC_filtered$controlsAC), 
                                Inf))
    hiAC_col_to_merge <- paste0("allele_counts_", min_row)
    cols_to_combine <- unique(c(cols_to_combine, hiAC_col_to_merge))
    
    # Sort based on the numeric part at the end of each string
    sorted_cols <- cols_to_combine[order(as.numeric(sub(".*_", "", cols_to_combine)))]
    return(sorted_cols) 
  }
}

process_and_regress <- function(genotype_chunk, covar_pheno, no_alts) {
  print("running progress and regress")
  
  results_df <- data.frame()
  long_gt <- reshape2::melt(genotype_chunk, id.vars=c("locus","alt_lengths"), variable.name="IID", value.name="allele_counts")
  merged_data <- long_gt %>%
    left_join(covar_pheno, by=c("IID" = "IID"))
  if(no_alts == 1){
    for(locus_id in unique(merged_data$locus)) {
      print(locus_id)
      current_data <- merged_data[merged_data$locus == locus_id, ]
      current_data <- current_data %>% filter(!is.na(allele_counts))
      allele_col_names <- paste0(locus_id, "_", current_data$alt_lengths[1])
      cases_current_data <- current_data %>% filter(IID %in% cases_ids)
      AN_cases <- 2* sum(!is.na(cases_current_data$allele_counts))
      controls_current_data <- current_data %>% filter(IID %in% controls_ids)
      AN_controls <- 2* sum(!is.na(controls_current_data$allele_counts))
      if(sum(cases_current_data$allele_counts) < AC_min_threshold | sum(controls_current_data$allele_counts) < AC_min_threshold){
        result_row <- data.frame(
          locus = locus_id,
          Alleles = paste(allele_col_names, collapse = ", "),
          OR = NA,
          Std_Errors = NA,
          Z_Values = NA,
          P_Values = NA,
          Smallest_P_Value = "low_AC",
          N_Samples = NA)
      }
      else if((AN_cases - sum(cases_current_data$allele_counts)) < AC_min_threshold | (AN_controls - sum(controls_current_data$allele_counts)) < AC_min_threshold){
        result_row <- data.frame(
          locus = locus_id,
          Alleles = paste(allele_col_names, collapse = ", "),
          OR = NA,
          Std_Errors = NA,
          Z_Values = NA,
          P_Values = NA,
          Smallest_P_Value = "low_ref_AC",
          N_Samples = NA)
      } 
      else{
        covariate_cols <- covar_names
        formula_str <- paste0("disease_status ~ allele_counts + ", paste(covariate_cols, collapse=" + "))
        formula <- as.formula(formula_str)
        model <- glm(formula, data=current_data, family=binomial)
        n_samples <- nobs(model)
        model_summary <- summary(model)$coefficients
        allele_summary <- model_summary["allele_counts", ]
        betas <- custom_round(unname(allele_summary["Estimate"]),6)
        odds_ratios <- custom_round(exp(betas),6)
        std_errors <-  custom_round(unname(allele_summary["Std. Error"]),6)
        z_values <-  custom_round(unname(allele_summary["z value"]),4)
        p_values <- custom_round(unname(allele_summary["Pr(>|z|)"]),100)
        result_row <- data.frame(
          locus = locus_id,
          Alleles = paste(allele_col_names, collapse = ", "),
          OR = paste(odds_ratios, collapse = ", "),
          Std_Errors = paste(std_errors, collapse = ", "),
          Z_Values = paste(z_values, collapse = ", "),
          P_Values = paste(p_values, collapse = ", "),
          Smallest_P_Value = p_values,
          N_Samples = n_samples)
        
      }
      
      results_df <- rbind(results_df, result_row)
      
    }
    
  }
  else{
    alt_colnames <- paste0("allele_counts_", 1:no_alts)
    merged_data <- merged_data %>%
      # Use 'separate' to split the 'allele_counts' column into multiple new columns
      separate(col = "allele_counts", 
               into = alt_colnames, 
               sep = ",", 
               remove = TRUE, # Set to FALSE if you want to keep the original column
               convert = TRUE)
    
    for(locus_id in unique(merged_data$locus)) {
      print(locus_id)
      current_data <- merged_data[merged_data$locus == locus_id, ]
      cols_to_merge <- low_cols(current_data, alt_colnames)
      #check if cols_to_merge returns "low_ref_count
      if ("low_ref_count" %in% cols_to_merge){
        allele_lengths <- current_data$alt_lengths[1]
        split_alleles <- str_split(allele_lengths, pattern=",")[[1]]
        alleles <- paste(locus_id, split_alleles, sep = "_")
        result_row <- data.frame(
          locus = locus_id,
          Alleles = paste(alleles, collapse = ", "),
          OR = NA,
          Std_Errors = NA,
          Z_Values = NA,
          P_Values = NA,
          Smallest_P_Value = "low_ref_AC",
          N_Samples = NA)
      }
      else{
        if(is.null(cols_to_merge)){
          allele_lengths <- current_data$alt_lengths[1]
          split_alleles <- str_split(allele_lengths, pattern=",")[[1]]
          alleles <- paste(locus_id, split_alleles, sep = "_")
          covariate_cols <- covar_names
          formula_str <- paste0("disease_status ~ ", paste(alt_colnames, collapse = " + "), " + ", paste(covariate_cols, collapse=" + "))
          formula <- as.formula(formula_str)
          model <- glm(formula, data = current_data, family = binomial)
          n_samples <- nobs(model)
          model_summary <- summary(model)$coefficients
          existing_colnames <- alt_colnames
          allele_summary <- model_summary[existing_colnames, , drop=F]
          odds_ratios <- paste(custom_round(exp(as.numeric(allele_summary[, "Estimate"])), 6), collapse = ',')
          std_errors <- paste(custom_round(as.numeric(allele_summary[, "Std. Error"]), 6), collapse = ',')
          z_values <- paste(custom_round(as.numeric(allele_summary[, "z value"]), 6), collapse = ',')
          p_values <- paste(sapply(as.numeric(allele_summary[, "Pr(>|z|)"]), custom_round, significant_digits = 6), collapse = ',')
          numeric_p_values <- as.numeric(strsplit(p_values, ",")[[1]])
          smallest_p <- min(numeric_p_values, na.rm=T)
          
          result_row <- data.frame(
            locus = locus_id,
            Alleles = paste(alleles, collapse = ','),
            OR = odds_ratios,
            Std_Errors = std_errors,
            Z_Values = z_values,
            P_Values = p_values,
            Smallest_P_Value = smallest_p,
            N_Samples = n_samples)
        }
        else{
          numbers <- sub(".*_", "", cols_to_merge)
          new_col_name <- paste("allele_counts", paste(numbers, collapse = "_"), sep = "_")
          current_data[[new_col_name]] <- rowSums(current_data[cols_to_merge])
          current_data <- current_data %>% select(-cols_to_merge)
          indices <- as.numeric(numbers)
          allele_lengths <- current_data$alt_lengths[1]
          split_alleles <- str_split(allele_lengths, pattern = ",")[[1]]
          alleles_to_combine <- split_alleles[indices]
          ordered_alleles <- sort(as.numeric(alleles_to_combine))
          combined_allele_lengths <- paste(ordered_alleles, collapse = "_")
          allele_identifier <- paste(locus_id, combined_allele_lengths, sep = "_")
          alleles_to_leave <- split_alleles[-indices]
          if(length(cols_to_merge) == length(alt_colnames)){
            alleles <- allele_identifier
          }
          else{
            alleles <- paste(locus_id, alleles_to_leave, sep = "_")
            alleles <- c(alleles, allele_identifier)
          }
          cases_current_data <- current_data %>% filter(IID %in% cases_ids)
          controls_current_data <- current_data %>% filter(IID %in% controls_ids)
          minAC_combined_cases <- sum(cases_current_data[new_col_name], na.rm = T)
          minAC_combined_controls <- sum(controls_current_data[new_col_name], na.rm = T)
          min_all <- min(c(minAC_combined_cases,minAC_combined_controls))
          if(length(cols_to_merge) == length(alt_colnames) & min_all < AC_min_threshold){
            result_row <- data.frame(
              locus = locus_id,
              Alleles = paste(alleles, collapse = ','),
              OR = NA,
              Std_Errors = NA,
              Z_Values = NA,
              P_Values = NA,
              Smallest_P_Value = "lowAC",
              N_Samples = NA)
          }
          
          else{
            
            covariate_cols <- covar_names
            if(length(cols_to_merge) == length(alt_colnames)){
              formula_str <- paste0("disease_status ~ ", new_col_name, " + ", paste(covariate_cols, collapse=" + "))
            }
            else{
              formula_str <- paste0("disease_status ~ ", paste(alt_colnames[-indices], collapse = " + ")," + ", new_col_name, " + ", paste(covariate_cols, collapse=" + "))
            }
            
            formula <- as.formula(formula_str)
            model <- glm(formula, data=current_data, family=binomial)
            n_samples <- nobs(model)
            model_summary <- summary(model)$coefficients
            existing_colnames <- c(paste(alt_colnames[-indices]), new_col_name)
            allele_summary <- model_summary[existing_colnames, , drop=F]
            odds_ratios <- paste(custom_round(exp(as.numeric(allele_summary[, "Estimate"])), 6), collapse = ',')
            std_errors <- paste(custom_round(as.numeric(allele_summary[, "Std. Error"]), 6), collapse = ',')
            z_values <- paste(custom_round(as.numeric(allele_summary[, "z value"]), 6), collapse = ',')
            p_values <- paste(sapply(as.numeric(allele_summary[, "Pr(>|z|)"]), custom_round, significant_digits = 6), collapse = ',')
            numeric_p_values <- as.numeric(strsplit(p_values, ",")[[1]])
            smallest_p <- min(numeric_p_values, na.rm=T)
            
            
            result_row <- data.frame(
              locus = locus_id,
              Alleles = paste(alleles, collapse = ','),
              OR = odds_ratios,
              Std_Errors = std_errors,
              Z_Values = z_values,
              P_Values = p_values,
              Smallest_P_Value = smallest_p,
              N_Samples = n_samples)
            
          }
        }
      }
      results_df <- rbind(results_df, result_row)}
    
    
  }
  return(results_df)
  
}

directory <- "QC"
pattern = "_alt_alleles\\.txt$"

process_grouped_files <- function(directory, pattern, covar_pheno) {
  files <- list.files(path = directory, pattern = pattern, full.names = TRUE)
  all_results <- list()
  column_names <- c("locus", "alt_lengths", samples_list)
  for (file in files) {
    print(file)
    no_alts <- as.numeric(gsub("[^0-9]", "", basename(file)))
    total_rows <- nrow(fread(file, select = 1))  # Get total number of rows in the file
    chunk_size <- 100
    chunks <- ceiling(total_rows / chunk_size)  # Calculate the number of chunks
    
    # Initialize an empty list to store results of each chunk
    chunk_results_list <- list()
    for (chunk_idx in 1:chunks) {
      print(chunk_idx)
      start_row <- (chunk_idx - 1) * chunk_size + 1
      nrows <- ifelse(chunk_idx * chunk_size > total_rows, total_rows - (chunk_idx - 1) * chunk_size, chunk_size)
      genotype_chunk <- fread(file, skip = start_row - 1, nrows = nrows, sep='\t',header = FALSE, dec = '.')  # Read chunk_size lines
      colnames(genotype_chunk) <- column_names
      # Process each chunk
      chunk_results <- process_and_regress(genotype_chunk, covar_pheno, no_alts)
      chunk_results_list[[chunk_idx]] <- chunk_results
    }
    
    # Combine results from all chunks for the current file
    results <- rbindlist(chunk_results_list, use.names = TRUE)
    fwrite(results, paste0("results/chunks/", no_alts, "_no_alts_results.tsv"), sep = "\t")
    all_results[[no_alts]] <- results
  }
  
  # Combine results from all files and write to a single file
  final_results_df <- rbindlist(all_results, use.names = TRUE)
  fwrite(final_results_df, paste0("results", "/final_combined_regression_results.tsv"), sep = "\t")
}
process_grouped_files(directory = "QC", pattern = "_alt_alleles\\.txt$", covar_pheno)


