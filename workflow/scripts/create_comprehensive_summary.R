#!/usr/bin/env Rscript

# Comprehensive Prophage Pipeline Summary Generator
# Creates data warehouse tables and basic visualizations

library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(knitr)
library(scales)

# Helper function to create directory structure
create_output_structure <- function(outdir) {
  subdirs <- c("tables", "plots", "per_sample_summaries", "raw")
  for (subdir in subdirs) {
    dir.create(file.path(outdir, "summary_results", subdir), 
               showWarnings = FALSE, recursive = TRUE)
  }
  return(file.path(outdir, "summary_results"))
}

# Function to collect all prophage data with detailed information
collect_master_prophage_catalog <- function(outdir) {
  message("Creating master prophage catalog...")
  
  # Find all prophage table files
  prophage_files <- list.files(
    path = file.path(outdir, "phage_analysis"),
    pattern = "final_prophage_table_with_host_taxonomy.tsv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  master_catalog <- map_dfr(prophage_files, function(file) {
    sample_name <- str_extract(file, "phage_analysis/([^/]+)/", group = 1)
    
    if (file.exists(file)) {
      df <- read_tsv(file, show_col_types = FALSE) %>%
        mutate(
          sample_id = sample_name,
          prophage_id = paste(sample_name, contig, start, end, sep = "_"),
          length = end - start + 1
        ) %>%
        select(prophage_id, sample_id, contig, start, end, length, tool, everything())
    } else {
      tibble()
    }
  })
  
  return(master_catalog)
}

# Function to collect CheckV quality data and merge
collect_checkv_data <- function(outdir) {
  message("Collecting CheckV quality data...")
  
  checkv_files <- list.files(
    path = file.path(outdir, "phage_analysis"),
    pattern = "quality_summary.tsv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  checkv_data <- map_dfr(checkv_files, function(file) {
    sample_name <- str_extract(file, "phage_analysis/([^/]+)/", group = 1)
    
    if (file.exists(file)) {
      read_tsv(file, show_col_types = FALSE) %>%
        mutate(sample_id = sample_name)
    } else {
      tibble()
    }
  })
  
  return(checkv_data)
}

# Function to create host-prophage relationships table
create_host_prophage_table <- function(master_catalog) {
  message("Creating host-prophage relationships table...")
  
  # Extract taxonomic columns
  tax_cols <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  available_tax_cols <- tax_cols[tax_cols %in% colnames(master_catalog)]
  
  if (length(available_tax_cols) > 0) {
    host_prophage <- master_catalog %>%
      select(prophage_id, sample_id, contig, all_of(available_tax_cols)) %>%
      # Create simplified taxonomy string
      rowwise() %>%
      mutate(
        simplified_taxonomy = paste(
          na.omit(c_across(all_of(available_tax_cols))), 
          collapse = ";"
        )
      ) %>%
      ungroup()
  } else {
    host_prophage <- master_catalog %>%
      select(prophage_id, sample_id, contig) %>%
      mutate(simplified_taxonomy = "Unknown")
  }
  
  return(host_prophage)
}

# Function to create sample-level summary statistics
create_sample_summary <- function(master_catalog, checkv_data) {
  message("Creating sample-level summary statistics...")
  
  sample_summary <- master_catalog %>%
    group_by(sample_id) %>%
    summarise(
      total_prophages = n(),
      unique_contigs_with_prophages = n_distinct(contig),
      genomad_prophages = sum(tool == "genomad", na.rm = TRUE),
      phispy_prophages = sum(tool == "phispy", na.rm = TRUE),
      mean_prophage_length = round(mean(length, na.rm = TRUE), 1),
      median_prophage_length = round(median(length, na.rm = TRUE), 1),
      min_prophage_length = min(length, na.rm = TRUE),
      max_prophage_length = max(length, na.rm = TRUE),
      unique_phyla = n_distinct(phylum, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add CheckV quality metrics if available
  if (nrow(checkv_data) > 0) {
    checkv_summary <- checkv_data %>%
      group_by(sample_id) %>%
      summarise(
        high_quality_prophages = sum(checkv_quality == "High-quality", na.rm = TRUE),
        medium_quality_prophages = sum(checkv_quality == "Medium-quality", na.rm = TRUE),
        complete_prophages = sum(checkv_quality == "Complete", na.rm = TRUE),
        .groups = "drop"
      )
    
    sample_summary <- sample_summary %>%
      left_join(checkv_summary, by = "sample_id")
  }
  
  return(sample_summary)
}

# Function to create tool detection comparison
create_tool_comparison <- function(master_catalog) {
  message("Creating tool detection comparison...")
  
  # Overall tool detection statistics
  tool_summary <- master_catalog %>%
    count(tool, name = "total_detections") %>%
    mutate(percentage = round(total_detections / sum(total_detections) * 100, 1))
  
  # Tool overlap analysis
  tool_overlap <- master_catalog %>%
    group_by(sample_id, contig) %>%
    summarise(
      tools_detected = paste(sort(unique(tool)), collapse = ","),
      num_tools = n_distinct(tool),
      .groups = "drop"
    ) %>%
    count(tools_detected, name = "contig_count") %>%
    mutate(percentage = round(contig_count / sum(contig_count) * 100, 1))
  
  # Per-sample tool detection counts
  tool_by_sample <- master_catalog %>%
    group_by(sample_id, tool) %>%
    summarise(detections = n(), .groups = "drop") %>%
    pivot_wider(names_from = tool, values_from = detections, values_fill = 0)
  
  return(list(
    overall = tool_summary,
    overlap = tool_overlap,
    by_sample = tool_by_sample
  ))
}

# Function to create basic visualizations
create_basic_visualizations <- function(master_catalog, sample_summary, tool_comparison, plots_dir) {
  message("Creating basic visualizations...")
  
  # Set theme
  theme_set(theme_bw(base_size = 12))
  
  # 1. Prophage counts per sample
  p1 <- ggplot(sample_summary, aes(x = reorder(sample_id, total_prophages), y = total_prophages)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = "Total Prophages per Sample", x = "Sample", y = "Number of Prophages") +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(plots_dir, "prophages_per_sample.png"), p1, width = 10, height = 6, dpi = 300)
  
  # 2. Tool detection comparison by sample (stacked bar)
  tool_data <- master_catalog %>%
    count(sample_id, tool) %>%
    complete(sample_id, tool, fill = list(n = 0))
  
  p2 <- ggplot(tool_data, aes(x = sample_id, y = n, fill = tool)) +
    geom_col(position = "stack", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Tool Detection Comparison by Sample", x = "Sample", y = "Number of Prophages", fill = "Tool") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  
  ggsave(file.path(plots_dir, "tool_detection_by_sample.png"), p2, width = 12, height = 6, dpi = 300)
  
  # 3. Overall tool detection comparison
  p3 <- ggplot(tool_comparison$overall, aes(x = tool, y = total_detections, fill = tool)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0(total_detections, "\n(", percentage, "%)")), vjust = -0.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Overall Tool Detection Comparison", x = "Detection Tool", y = "Total Detections") +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  ggsave(file.path(plots_dir, "overall_tool_detection_comparison.png"), p3, width = 8, height = 6, dpi = 300)
  
  # 4. Prophage length distribution
  p4 <- ggplot(master_catalog, aes(x = length)) +
    geom_histogram(bins = 30, fill = "darkorange", alpha = 0.8, color = "black") +
    scale_x_log10(labels = comma) +
    labs(title = "Prophage Length Distribution", x = "Prophage Length (bp)", y = "Count") +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(plots_dir, "length_distribution.png"), p4, width = 10, height = 6, dpi = 300)
  
  # 5. Host taxonomy distribution (top phyla)
  if ("phylum" %in% colnames(master_catalog)) {
    phylum_data <- master_catalog %>%
      filter(!is.na(phylum)) %>%
      count(phylum, sort = TRUE) %>%
      slice_head(n = 10)
    
    p5 <- ggplot(phylum_data, aes(x = reorder(phylum, n), y = n)) +
      geom_col(fill = "forestgreen", alpha = 0.8) +
      coord_flip() +
      labs(title = "Top 10 Host Phyla", x = "Phylum", y = "Number of Prophages") +
      theme(panel.grid.minor = element_blank())
    
    ggsave(file.path(plots_dir, "host_phyla_distribution.png"), p5, width = 10, height = 6, dpi = 300)
  }
  
  message("Basic visualizations saved to:", plots_dir)
}

# Function to create per-sample summary files
create_per_sample_summaries <- function(master_catalog, per_sample_dir) {
  message("Creating per-sample summary files...")
  
  samples <- unique(master_catalog$sample_id)
  
  for (sample in samples) {
    sample_data <- master_catalog %>%
      filter(sample_id == sample)
    
    # Create sample-specific summary
    sample_file <- file.path(per_sample_dir, paste0(sample, "_prophage_summary.tsv"))
    write_tsv(sample_data, sample_file)
  }
  
  message("Per-sample summaries saved to:", per_sample_dir)
}

# Main function
main <- function(outdir) {
  message("Starting comprehensive prophage summary analysis...")
  
  # Create output structure
  summary_dir <- create_output_structure(outdir)
  tables_dir <- file.path(summary_dir, "tables")
  plots_dir <- file.path(summary_dir, "plots")
  per_sample_dir <- file.path(summary_dir, "per_sample_summaries")
  
  # Collect all data
  master_catalog <- collect_master_prophage_catalog(outdir)
  checkv_data <- collect_checkv_data(outdir)
  
  if (nrow(master_catalog) == 0) {
    stop("No prophage data found. Please check that the pipeline has completed successfully.")
  }
  
  # Create summary tables
  host_prophage_table <- create_host_prophage_table(master_catalog)
  sample_summary <- create_sample_summary(master_catalog, checkv_data)
  tool_comparison <- create_tool_comparison(master_catalog)
  
  # Save main tables
  write_tsv(master_catalog, file.path(tables_dir, "master_prophage_catalog.tsv"))
  write_tsv(host_prophage_table, file.path(tables_dir, "host_prophage_relationships.tsv"))
  write_tsv(sample_summary, file.path(tables_dir, "sample_level_summary.tsv"))
  write_tsv(tool_comparison$overall, file.path(tables_dir, "tool_detection_overall.tsv"))
  write_tsv(tool_comparison$overlap, file.path(tables_dir, "tool_overlap_analysis.tsv"))
  write_tsv(tool_comparison$by_sample, file.path(tables_dir, "tool_detection_by_sample.tsv"))
  
  # Create visualizations
  create_basic_visualizations(master_catalog, sample_summary, tool_comparison, plots_dir)
  
  # Create per-sample summaries
  create_per_sample_summaries(master_catalog, per_sample_dir)
  
  # Generate markdown report
  # Read the main tables
  markdown_content <- paste0('
# Prophage Pipeline Summary Report

**Generated on:** ', Sys.Date(), '  
**Total samples:** ', length(unique(master_catalog$sample_id)), '  
**Total prophages detected:** ', nrow(master_catalog), '  

## Overview

This report summarizes the results of prophage detection across all samples in your dataset. The pipeline detected prophages using multiple tools and provided comprehensive quality assessment.

### Key Statistics

- **Average prophages per sample:** ', round(nrow(master_catalog) / length(unique(master_catalog$sample_id)), 1), '
- **Total unique contigs with prophages:** ', length(unique(master_catalog$contig)), '
- **Prophage length range:** ', min(master_catalog$length, na.rm = TRUE), ' - ', max(master_catalog$length, na.rm = TRUE), ' bp
- **Median prophage length:** ', median(master_catalog$length, na.rm = TRUE), ' bp

## Sample-Level Results

', kable(sample_summary, format = "markdown"), '

## Tool Detection Comparison

', kable(tool_comparison$overall, format = "markdown"), '

## Output Files Description

### Main Data Tables (`tables/` directory)

- **`master_prophage_catalog.tsv`** - Complete catalog of all prophages detected across all samples
- **`host_prophage_relationships.tsv`** - Host taxonomy information for each prophage
- **`sample_level_summary.tsv`** - Per-sample summary statistics
- **`tool_detection_*.tsv`** - Tool detection statistics and overlap analyses

### Visualizations (`plots/` directory)
- `prophages_per_sample.png` - Sample comparison chart
- `tool_detection_by_sample.png` - Tool detection comparison by sample
- `overall_tool_detection_comparison.png` - Overall tool detection comparison
- `length_distribution.png` - Prophage length distribution
- `host_phyla_distribution.png` - Host taxonomy overview

## Usage Examples

### Load main dataset
```r
prophages <- read_tsv("tables/master_prophage_catalog.tsv")
```

### Filter by quality
```r
high_quality <- prophages %>% filter(checkv_quality %in% c("High-quality", "Complete"))
```

### Sample comparisons
```r
sample_stats <- prophages %>% group_by(sample_id) %>% summarise(n = n(), mean_length = mean(length))
```

---
*Generated by the Prophage Detection Pipeline*
')

  # Write the markdown file
  report_file <- file.path(summary_dir, "SUMMARY_REPORT.md")
  writeLines(markdown_content, report_file)
  
  # Print summary statistics
  message("\n=== COMPREHENSIVE SUMMARY STATISTICS ===")
  message(paste("Total samples analyzed:", length(unique(master_catalog$sample_id))))
  message(paste("Total prophages detected:", nrow(master_catalog)))
  message(paste("Average prophages per sample:", round(nrow(master_catalog) / length(unique(master_catalog$sample_id)), 1)))
  message(paste("Unique contigs with prophages:", length(unique(master_catalog$contig))))
  
  if ("phylum" %in% colnames(master_catalog)) {
    message(paste("Unique host phyla:", length(unique(master_catalog$phylum[!is.na(master_catalog$phylum)]))))
  }
  
  message("\n=== FILES CREATED ===")
  message("Main tables:")
  message("  - master_prophage_catalog.tsv")
  message("  - host_prophage_relationships.tsv") 
  message("  - sample_level_summary.tsv")
  message("  - tool_detection_*.tsv")
  message("\nVisualization plots:")
  message("  - prophages_per_sample.png")
  message("  - tool_detection_by_sample.png")
  message("  - overall_tool_detection_comparison.png")
  message("  - length_distribution.png")
  message("  - host_phyla_distribution.png")
  
  message(paste("\nAll results saved to:", summary_dir))
  
  return(summary_dir)
}

# Run if called as script or from Snakemake
if (!interactive()) {
  if (exists("snakemake")) {
    # Extract output directory from snakemake
    outdir <- dirname(snakemake@output[["summary_dir"]])
    main(outdir)
  } else {
    # Command line usage
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
      stop("Usage: Rscript create_comprehensive_summary.R <output_directory>")
    }
    main(args[1])
  }
}