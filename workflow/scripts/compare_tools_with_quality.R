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

print("Loading tool predictions...")

# Load and process geNomad predictions
print("Loading geNomad predictions...")
genomad_path <- file.path(snakemake@input[["genomad"]], 
                          "final_filtered_contigs_find_proviruses", 
                          "final_filtered_contigs_provirus.tsv")
genomad <- read_tsv(genomad_path) %>%
  extract(source_seq, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
  mutate(
    tool = "genomad",
    length = end - start + 1,
    sequence_id = paste0("genomad_", row_number()),
    prediction_id = paste0(source_seq, "_", start, "_", end)
  ) %>%
  select(contig, start, end, tool, length, sequence_id, prediction_id, source_seq)

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
  mutate(
    tool = "phispy",
    length = end - start + 1,
    sequence_id = paste0("phispy_", pp_number),
    prediction_id = paste0(Contig, "_", start, "_", end)
  ) %>%
  select(contig, start, end, tool, length, sequence_id, prediction_id, pp_number)

# Load and process PIDE predictions
print("Loading PIDE predictions...")
pide_path <- file.path(snakemake@input[["pide"]], "cluster.csv")
if (file.exists(pide_path)) {
  pide <- read_csv(pide_path) %>%
    extract(Contig, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
    as.data.frame() %>%
    select(contig, Start, End, Contig) %>%
    dplyr::rename(start = Start, end = End, contig_id = Contig) %>%
    mutate(
      tool = "pide",
      length = end - start + 1,
      sequence_id = paste0("pide_", row_number()),
      prediction_id = paste0(contig_id, "_", start, "_", end)
    ) %>%
    select(contig, start, end, tool, length, sequence_id, prediction_id, contig_id)
} else {
  print("Warning: PIDE predictions file not found. Creating empty dataset.")
  pide <- data.frame(
    contig = character(), start = numeric(), end = numeric(), 
    tool = character(), length = numeric(), sequence_id = character(),
    prediction_id = character(), contig_id = character(),
    stringsAsFactors = FALSE
  )
}

# Load CheckV quality results for each tool
print("Loading CheckV quality results...")

load_checkv_results <- function(checkv_dir, tool_name) {
  quality_file <- file.path(checkv_dir, "quality_summary.tsv")
  if (file.exists(quality_file)) {
    checkv_data <- read_tsv(quality_file) %>%
      mutate(tool = tool_name) %>%
      select(contig_id, contig_length, provirus, proviral_length, gene_count, 
             viral_genes, host_genes, checkv_quality, miuvig_quality,
             completeness, contamination, tool)
    return(checkv_data)
  } else {
    print(paste("Warning: CheckV quality file not found for", tool_name))
    return(data.frame())
  }
}

checkv_genomad <- load_checkv_results(snakemake@input[["checkv_genomad"]], "genomad")
checkv_phispy <- load_checkv_results(snakemake@input[["checkv_phispy"]], "phispy")
checkv_pide <- load_checkv_results(snakemake@input[["checkv_pide"]], "pide")

# Combine all CheckV results
all_checkv <- bind_rows(checkv_genomad, checkv_phispy, checkv_pide)

# Create sequence mapping for joining with CheckV results
# This is tricky - we need to map tool predictions to CheckV sequence IDs
# For now, we'll do a simplified mapping based on sequence names/coordinates

# Combine all predictions
all_predictions <- bind_rows(
  genomad %>% select(contig, start, end, tool, length, sequence_id, prediction_id),
  phispy %>% select(contig, start, end, tool, length, sequence_id, prediction_id),
  pide %>% select(contig, start, end, tool, length, sequence_id, prediction_id)
)

# Convert contig to numeric for consistency
all_predictions$contig <- as.numeric(all_predictions$contig)

# Enhanced analysis with quality metrics
print("Performing enhanced analysis with quality metrics...")

# Basic tool statistics (preserve original functionality)
basic_stats <- all_predictions %>%
  group_by(tool) %>%
  summarise(
    total_predictions = n(),
    mean_length = mean(length, na.rm = TRUE),
    median_length = median(length, na.rm = TRUE),
    sd_length = sd(length, na.rm = TRUE),
    min_length = min(length, na.rm = TRUE),
    max_length = max(length, na.rm = TRUE),
    total_bp_predicted = sum(length, na.rm = TRUE),
    unique_contigs = n_distinct(contig),
    .groups = 'drop'
  )

# Quality-enhanced statistics
if (nrow(all_checkv) > 0) {
  quality_stats <- all_checkv %>%
    group_by(tool) %>%
    summarise(
      high_quality = sum(checkv_quality == "High-quality", na.rm = TRUE),
      medium_quality = sum(checkv_quality == "Medium-quality", na.rm = TRUE),
      low_quality = sum(checkv_quality == "Low-quality", na.rm = TRUE),
      not_determined = sum(checkv_quality == "Not-determined", na.rm = TRUE),
      mean_completeness = mean(completeness, na.rm = TRUE),
      mean_contamination = mean(contamination, na.rm = TRUE),
      mean_viral_genes = mean(viral_genes, na.rm = TRUE),
      mean_host_genes = mean(host_genes, na.rm = TRUE),
      viral_gene_ratio = mean(viral_genes / (viral_genes + host_genes), na.rm = TRUE),
      .groups = 'drop'
    )
} else {
  quality_stats <- data.frame(tool = character())
}

# Biological plausibility filters
bio_filtered <- all_predictions %>%
  mutate(
    reasonable_length = length >= 5000 & length <= 200000,
    very_short = length < 5000,
    very_long = length > 200000
  )

bio_stats <- bio_filtered %>%
  group_by(tool) %>%
  summarise(
    reasonable_length_count = sum(reasonable_length),
    very_short_count = sum(very_short),
    very_long_count = sum(very_long),
    reasonable_length_pct = mean(reasonable_length) * 100,
    .groups = 'drop'
  )

# Tool performance comparison
tool_performance <- basic_stats %>%
  left_join(bio_stats, by = "tool") %>%
  left_join(quality_stats, by = "tool") %>%
  mutate(
    quality_ratio = ifelse(is.na(high_quality + medium_quality), NA, 
                          (high_quality + medium_quality) / (high_quality + medium_quality + low_quality + not_determined)),
    recommendations = case_when(
      reasonable_length_pct > 80 & (is.na(quality_ratio) | quality_ratio > 0.6) ~ "High confidence tool",
      reasonable_length_pct > 60 & (is.na(quality_ratio) | quality_ratio > 0.4) ~ "Medium confidence tool", 
      TRUE ~ "Needs validation"
    )
  )

# Overlap analysis (preserve from original)
print("Performing overlap analysis...")
overlap_results <- data.frame()
contigs_with_predictions <- unique(all_predictions$contig)

for (contig_id in contigs_with_predictions) {
  contig_preds <- all_predictions %>% filter(contig == contig_id)
  
  if (nrow(contig_preds) > 1) {
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

# High confidence predictions (multi-tool agreement + quality)
high_confidence <- overlap_results %>%
  filter(overlap == TRUE, overlap_percent > 50) %>%
  select(contig, tool1, tool2) %>%
  distinct() %>%
  mutate(agreement_type = "multi_tool")

# Enhanced visualizations
print("Creating enhanced visualizations...")

# Quality distribution plot
if (nrow(all_checkv) > 0) {
  quality_plot <- ggplot(all_checkv, aes(x = tool, fill = checkv_quality)) +
    geom_bar(position = "stack") +
    labs(title = "CheckV Quality Distribution by Tool",
         x = "Tool", y = "Number of Predictions", fill = "CheckV Quality") +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set2")
} else {
  quality_plot <- ggplot() + 
    annotate("text", x = 1, y = 1, label = "No CheckV quality data available") +
    theme_minimal()
}

# Length vs quality scatter (if quality data available)
if (nrow(all_checkv) > 0) {
  length_quality_plot <- all_checkv %>%
    filter(!is.na(checkv_quality)) %>%
    ggplot(aes(x = contig_length, y = completeness, color = tool, shape = checkv_quality)) +
    geom_point(alpha = 0.7) +
    scale_x_log10() +
    labs(title = "Sequence Length vs Completeness",
         x = "Sequence Length (bp, log scale)", 
         y = "Completeness (%)",
         color = "Tool", shape = "Quality") +
    theme_minimal()
} else {
  length_quality_plot <- ggplot() + 
    annotate("text", x = 1, y = 1, label = "No quality data for length analysis") +
    theme_minimal()
}

# Enhanced length distribution
length_plot <- ggplot(all_predictions, aes(x = length, fill = tool)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(5000, 200000), linetype = "dashed", color = "red") +
  annotate("text", x = 5000, y = Inf, label = "5kb", vjust = 2, hjust = -0.1) +
  annotate("text", x = 200000, y = Inf, label = "200kb", vjust = 2, hjust = 1.1) +
  scale_x_log10() +
  facet_wrap(~tool, scales = "free_y") +
  labs(title = "Prophage Length Distribution by Tool (with biological filters)",
       x = "Prophage Length (bp, log scale)",
       y = "Count") +
  theme_minimal()

# Tool performance comparison plot
performance_plot <- tool_performance %>%
  select(tool, total_predictions, reasonable_length_count, high_quality, medium_quality) %>%
  pivot_longer(cols = -tool, names_to = "metric", values_to = "count") %>%
  filter(!is.na(count)) %>%
  ggplot(aes(x = tool, y = count, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Tool Performance Comparison",
       x = "Tool", y = "Count", fill = "Metric") +
  theme_minimal()

# Write outputs
print("Writing output files...")

write_tsv(all_predictions, snakemake@output[["raw_predictions"]])
write_tsv(quality_stats, snakemake@output[["quality_stats"]])
write_tsv(tool_performance, snakemake@output[["tool_performance"]])
write_tsv(high_confidence, snakemake@output[["high_confidence"]])

# Save plots
pdf(snakemake@output[["plots"]], width = 12, height = 8)
print(quality_plot)
print(length_quality_plot)  
print(length_plot)
print(performance_plot)
dev.off()

print("Enhanced comparison analysis complete!")
print(paste("Results saved to:", dirname(snakemake@output[["plots"]])))