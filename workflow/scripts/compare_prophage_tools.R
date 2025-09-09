library(readr)
library(tidyverse)
library(Biostrings)
library(VennDiagram)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Function to check if 2 regions overlap
check_overlap <- function(region1, region2, tolerance = 0) {
  if (region1[2] >= (region2[1] - tolerance) && region1[1] <= (region2[2] + tolerance)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Function to calculate overlap percentage
calculate_overlap_percent <- function(region1, region2) {
  overlap_start <- max(region1[1], region2[1])
  overlap_end <- min(region1[2], region2[2])
  if (overlap_start <= overlap_end) {
    overlap_length <- overlap_end - overlap_start + 1
    total_length <- max(region1[2], region2[2]) - min(region1[1], region2[1]) + 1
    return(overlap_length / total_length * 100)
  }
  return(0)
}

# Load and process geNomad predictions
print("Loading geNomad predictions...")
genomad_path <- file.path(snakemake@input[["genomad"]], 
                          "final_filtered_contigs_find_proviruses", 
                          "final_filtered_contigs_provirus.tsv")
genomad <- read_tsv(genomad_path) %>%
  extract(source_seq, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
  as.data.frame() %>%
  select(contig, start, end) %>%
  mutate(tool = "genomad", length = end - start + 1)

# Load and process PhiSpy predictions
print("Loading PhiSpy predictions...")
phispy_path <- file.path(snakemake@input[["phispy"]], "prophage.tsv") 
phispy <- read_tsv(phispy_path) %>%
  separate_wider_delim('Prophage number', '_', names=c('pp', 'pp_number')) %>%
  extract(Contig, into = "contig", regex = "contig_(\\d+)", remove = FALSE) %>%
  as.data.frame()

# Handle PhiSpy column names
if (!all(c('Start', 'Stop') %in% colnames(phispy))) {
  stop("Columns 'Start' and 'Stop' not found in phispy data.")
}

phispy <- phispy %>%
  dplyr::rename(start = Start, end = Stop) %>%
  select(contig, start, end, pp_number) %>%
  mutate(tool = "phispy", length = end - start + 1)

# Load and process PIDE predictions
print("Loading PIDE predictions...")
# PIDE outputs cluster.csv with prophage island predictions
pide_path <- file.path(snakemake@input[["pide"]], "cluster.csv")
if (file.exists(pide_path)) {
  pide <- read_csv(pide_path) %>%
    extract(Contig, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
    as.data.frame() %>%
    select(contig, Start, End) %>%
    dplyr::rename(start = Start, end = End) %>%
    mutate(tool = "pide", length = end - start + 1)
} else {
  # Create empty dataframe if PIDE output doesn't exist
  print("Warning: PIDE predictions file not found. Creating empty dataset.")
  pide <- data.frame(contig = character(), start = numeric(), end = numeric(), 
                     tool = character(), length = numeric(), stringsAsFactors = FALSE)
}

# Combine all predictions
all_predictions <- bind_rows(
  genomad %>% select(contig, start, end, tool, length),
  phispy %>% select(contig, start, end, tool, length),
  pide %>% select(contig, start, end, tool, length)
)

# Convert contig to numeric for consistency
all_predictions$contig <- as.numeric(all_predictions$contig)

# Summary statistics
print("Generating summary statistics...")
tool_stats <- all_predictions %>%
  group_by(tool) %>%
  summarise(
    total_predictions = n(),
    mean_length = mean(length, na.rm = TRUE),
    median_length = median(length, na.rm = TRUE),
    sd_length = sd(length, na.rm = TRUE),
    min_length = min(length, na.rm = TRUE),
    max_length = max(length, na.rm = TRUE),
    total_bp_predicted = sum(length, na.rm = TRUE),
    unique_contigs = n_distinct(contig)
  )

# Overlap analysis
print("Performing overlap analysis...")
overlap_results <- data.frame()
contigs_with_predictions <- unique(all_predictions$contig)

for (contig_id in contigs_with_predictions) {
  contig_preds <- all_predictions %>% filter(contig == contig_id)
  
  if (nrow(contig_preds) > 1) {
    # Check all pairwise overlaps within this contig
    tools_present <- unique(contig_preds$tool)
    
    for (i in 1:(nrow(contig_preds)-1)) {
      for (j in (i+1):nrow(contig_preds)) {
        pred1 <- contig_preds[i,]
        pred2 <- contig_preds[j,]
        
        overlap_exists <- check_overlap(c(pred1$start, pred1$end), c(pred2$start, pred2$end))
        overlap_pct <- calculate_overlap_percent(c(pred1$start, pred1$end), c(pred2$start, pred2$end))
        
        overlap_results <- rbind(overlap_results, data.frame(
          contig = contig_id,
          tool1 = pred1$tool,
          tool2 = pred2$tool,
          start1 = pred1$start, end1 = pred1$end, length1 = pred1$length,
          start2 = pred2$start, end2 = pred2$end, length2 = pred2$length,
          overlap = overlap_exists,
          overlap_percent = overlap_pct
        ))
      }
    }
  }
}

# Tool agreement analysis
agreement_summary <- overlap_results %>%
  filter(overlap == TRUE) %>%
  group_by(tool1, tool2) %>%
  summarise(
    overlapping_pairs = n(),
    mean_overlap_percent = mean(overlap_percent),
    .groups = 'drop'
  )

# Unique predictions by tool
print("Identifying unique predictions...")
unique_predictions <- data.frame()

for (tool_name in c("genomad", "phispy", "pide")) {
  tool_preds <- all_predictions %>% filter(tool == tool_name)
  
  for (i in 1:nrow(tool_preds)) {
    current_pred <- tool_preds[i,]
    other_tools_preds <- all_predictions %>% 
      filter(tool != tool_name, contig == current_pred$contig)
    
    is_unique <- TRUE
    if (nrow(other_tools_preds) > 0) {
      for (j in 1:nrow(other_tools_preds)) {
        other_pred <- other_tools_preds[j,]
        if (check_overlap(c(current_pred$start, current_pred$end), 
                         c(other_pred$start, other_pred$end))) {
          is_unique <- FALSE
          break
        }
      }
    }
    
    if (is_unique) {
      unique_predictions <- rbind(unique_predictions, current_pred)
    }
  }
}

unique_summary <- unique_predictions %>%
  group_by(tool) %>%
  summarise(
    unique_predictions = n(),
    mean_unique_length = mean(length, na.rm = TRUE),
    total_unique_bp = sum(length, na.rm = TRUE)
  )

# Generate plots
print("Creating visualizations...")

# Length distribution plot
length_plot <- ggplot(all_predictions, aes(x = length, fill = tool)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_x_log10() +
  facet_wrap(~tool, scales = "free_y") +
  labs(title = "Prophage Length Distribution by Tool",
       x = "Prophage Length (bp, log scale)",
       y = "Count") +
  theme_minimal()

# Prediction counts by tool
count_plot <- ggplot(tool_stats, aes(x = tool, y = total_predictions, fill = tool)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Prophage Predictions by Tool",
       x = "Tool", y = "Number of Predictions") +
  theme_minimal() +
  theme(legend.position = "none")

# Write outputs
print("Writing output files...")

# Raw predictions with tool info
write_tsv(all_predictions, snakemake@output[["raw_predictions"]])

# Tool statistics
write_tsv(tool_stats, snakemake@output[["tool_stats"]])

# Overlap analysis
write_tsv(overlap_results, snakemake@output[["overlap_analysis"]])

# Agreement summary
write_tsv(agreement_summary, snakemake@output[["agreement_summary"]])

# Unique predictions
write_tsv(unique_predictions, snakemake@output[["unique_predictions"]])

# Unique summary
write_tsv(unique_summary, snakemake@output[["unique_summary"]])

# Save plots
pdf(snakemake@output[["plots"]], width = 12, height = 8)
print(length_plot)
print(count_plot)
dev.off()

print("Comparison analysis complete!")
print(paste("Total predictions: geNomad =", nrow(genomad), 
            ", PhiSpy =", nrow(phispy), 
            ", PIDE =", nrow(pide)))