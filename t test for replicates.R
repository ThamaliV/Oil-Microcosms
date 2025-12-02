library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Set working directory (adjust path as needed)
setwd("E:/Thami-Uni/Oil Microcosm Data")

# Read the data
cat("Loading data from working_data.csv...\n")
data <- read.csv("working_data.csv", stringsAsFactors = FALSE)

# Examine data structure
cat("Data dimensions:", dim(data), "\n")
cat("Column names:\n")
print(colnames(data))

# Check unique values in key columns
cat("\nUnique replicates:", unique(data$Replicate), "\n")
cat("Unique corrected sample names:", length(unique(data$Corrected_Sample_name)), "\n")

# Display first few rows
cat("\nFirst few rows of data:\n")
print(head(data[, c("Corrected_Sample_name", "Replicate", "Abundance", "ASV", "Phylum")]))

# Create a function to extract sample base name (remove REP1/REP2)
extract_sample_base <- function(sample_name) {
  gsub("_REP[12]$", "", sample_name)
}

# Add sample base name column
data$Sample_base <- extract_sample_base(data$Corrected_Sample_name)

# Check which samples have both replicates
replicate_summary <- data %>%
  group_by(Sample_base) %>%
  summarise(
    n_replicates = n_distinct(Replicate),
    replicates = paste(sort(unique(Replicate)), collapse = ", "),
    total_abundance = sum(Abundance, na.rm = TRUE),
    n_asvs = n_distinct(ASV)
  ) %>%
  arrange(desc(n_replicates))

cat("\nReplicate summary:\n")
print(replicate_summary)

# Filter for samples with both REP1 and REP2
samples_with_both_reps <- replicate_summary %>%
  filter(n_replicates == 2) %>%
  pull(Sample_base)

cat("\nSamples with both replicates:", length(samples_with_both_reps), "\n")
print(samples_with_both_reps)

# Function to perform t-test for a specific ASV between replicates

perform_asv_ttest <- function(asv_data) {
  # Initialize default values
  

  
  p_value <- NA
  t_stat <- NA
  significant <- FALSE
  
  # Check if we have enough data
  if (nrow(asv_data) < 2) {
    return(data.frame(
      Sample_base = sample_base,
      p_value = NA,
      t_stat = NA,
      significant = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract replicate data
  rep1_data <- asv_data$Abundance[asv_data$Replicate == "REP1"]
  rep2_data <- asv_data$Abundance[asv_data$Replicate == "REP2"]
  
  # Perform t-test if we have data for both replicates
  if (length(rep1_data) > 0 && length(rep2_data) > 0) {
    # Use log transformation for abundance data (add 1 to avoid log(0))
    rep1_log <- log(rep1_data + 1)
    rep2_log <- log(rep2_data + 1)
    
    # Perform t-test
    tryCatch({
      ttest_result <- t.test(rep1_log, rep2_log, var.equal = FALSE)
      p_value <- ttest_result$p.value
      t_stat <- as.numeric(ttest_result$statistic)
      significant <- p_value < 0.05
    }, error = function(e) {
      # If t-test fails, set default values
      p_value <<- 1.0
      t_stat <<- 0.0
      significant <<- FALSE
    })
  }
  
  # Return simplified results
  return(data.frame(
    
    p_value = p_value,
    t_stat = t_stat,
    significant = significant,
    stringsAsFactors = FALSE
  ))
}

filtered_data <- data %>%
  filter(Sample_base %in% samples_with_both_reps)

# Perform t-tests
ttest_results <- filtered_data %>%
  group_by( Sample_base) %>%
  group_modify(~ perform_asv_ttest(.x))
print(ttest_results)
ttest_results%>%
count(Sample_base)
# Remove rows with NA p-values
ttest_results <- ttest_results %>%
  filter(!is.na(p_value))


# Summary of significant differences
significant_results <- ttest_results %>%
  filter(significant == TRUE)

cat("Number of significant differences (p < 0.05):", nrow(significant_results), "\n")

# Summary by sample
sample_summary <- ttest_results %>%
  group_by(Sample_base) %>%
  summarise(
    total_asvs_tested = n(),
    significant_asvs = sum(significant, na.rm = TRUE),
    percent_significant = round(100 * significant_asvs / total_asvs_tested, 2),
    mean_p_value = round(mean(p_value, na.rm = TRUE), 4),
    median_p_value = round(median(p_value, na.rm = TRUE), 4)
  ) %>%
  arrange(desc(percent_significant))

library(officer)
library(flextable)
ft <- flextable(ttest_results)

# Create a Word document and add the table
doc <- read_docx() %>%
  body_add_par("T-test Results", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal")

# Save to file
print(doc, target = "ttest_results.docx")
