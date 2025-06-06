#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(plotly)
library(DT)
library(htmlwidgets)
library(knitr)
library(kableExtra)

# Function to collect prophage data from all samples
collect_prophage_data <- function(outdir) {
  message("Collecting prophage data...")
  
  # Find all prophage table files
  prophage_files <- list.files(
    path = file.path(outdir, "phage_analysis"),
    pattern = "final_prophage_table_with_host_taxonomy.tsv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Extract sample names from file paths
  sample_names <- str_extract(prophage_files, "phage_analysis/([^/]+)/", group = 1)
  
  # Read and combine all prophage data
  prophage_data <- map2_dfr(prophage_files, sample_names, function(file, sample) {
    if (file.exists(file)) {
      read_tsv(file, show_col_types = FALSE) %>%
        mutate(sample = sample)
    } else {
      tibble()
    }
  })
  
  return(prophage_data)
}

# Function to collect CheckV quality data
collect_checkv_data <- function(outdir) {
  message("Collecting CheckV quality data...")
  
  checkv_files <- list.files(
    path = file.path(outdir, "phage_analysis"),
    pattern = "quality_summary.tsv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  sample_names <- str_extract(checkv_files, "phage_analysis/([^/]+)/", group = 1)
  
  checkv_data <- map2_dfr(checkv_files, sample_names, function(file, sample) {
    if (file.exists(file)) {
      read_tsv(file, show_col_types = FALSE) %>%
        mutate(sample = sample)
    } else {
      tibble()
    }
  })
  
  return(checkv_data)
}

# Function to collect CheckM bin quality data
collect_checkm_data <- function(outdir) {
  message("Collecting CheckM bin quality data...")
  
  checkm_files <- list.files(
    path = file.path(outdir, "binning"),
    pattern = "checkm_out.tsv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  sample_names <- str_extract(checkm_files, "binning/([^/]+)/", group = 1)
  
  checkm_data <- map2_dfr(checkm_files, sample_names, function(file, sample) {
    if (file.exists(file)) {
      read_tsv(file, show_col_types = FALSE) %>%
        mutate(sample = sample)
    } else {
      tibble()
    }
  })
  
  return(checkm_data)
}

# Function to collect coverage data
collect_coverage_data <- function(outdir) {
  message("Collecting coverage data...")
  
  coverage_files <- list.files(
    path = file.path(outdir, "coverm"),
    pattern = "_stats.txt",
    recursive = TRUE,
    full.names = TRUE
  )
  
  sample_names <- str_extract(coverage_files, "coverm/([^/]+)/", group = 1)
  
  coverage_data <- map2_dfr(coverage_files, sample_names, function(file, sample) {
    if (file.exists(file)) {
      read_tsv(file, show_col_types = FALSE) %>%
        mutate(sample = sample)
    } else {
      tibble()
    }
  })
  
  return(coverage_data)
}

# Function to create summary statistics
create_summary_stats <- function(prophage_data, checkv_data, checkm_data, coverage_data) {
  message("Creating summary statistics...")
  
  # Prophage summary by sample
  prophage_summary <- prophage_data %>%
    group_by(sample) %>%
    summarise(
      total_prophages = n(),
      genomad_prophages = sum(tool == "genomad", na.rm = TRUE),
      phispy_prophages = sum(tool == "phispy", na.rm = TRUE),
      unique_contigs = n_distinct(contig),
      .groups = "drop"
    )
  
  # Tool detection summary
  tool_summary <- prophage_data %>%
    count(tool, name = "prophages_detected") %>%
    mutate(percentage = round(prophages_detected / sum(prophages_detected) * 100, 1))
  
  # Taxonomy summary at multiple levels
  taxonomy_summary <- list()
  
  # Define taxonomic levels to summarize
  tax_levels <- c("phylum", "class", "order", "family", "genus")
  
  for (level in tax_levels) {
    if (level %in% colnames(prophage_data)) {
      taxonomy_summary[[level]] <- prophage_data %>%
        filter(!is.na(.data[[level]])) %>%
        count(.data[[level]], sort = TRUE, name = "count") %>%
        rename(!!level := 1) %>%
        slice_head(n = 15) %>%
        mutate(
          percentage = round(count / sum(prophage_data %>% filter(!is.na(.data[[level]])) %>% nrow()) * 100, 1),
          level = level
        )
    }
  }
  
  # Keep the original phylum summary for backward compatibility
  if ("phylum" %in% colnames(prophage_data)) {
    taxonomy_summary_phylum <- prophage_data %>%
      filter(!is.na(phylum)) %>%
      count(phylum, sort = TRUE) %>%
      slice_head(n = 10)
  } else {
    taxonomy_summary_phylum <- tibble(phylum = character(), n = integer())
  }
  
  # CheckV quality summary
  if (nrow(checkv_data) > 0 && "checkv_quality" %in% colnames(checkv_data)) {
    quality_summary <- checkv_data %>%
      count(checkv_quality, name = "count") %>%
      mutate(percentage = round(count / sum(count) * 100, 1))
  } else {
    quality_summary <- tibble(checkv_quality = character(), count = integer(), percentage = numeric())
  }
  
  return(list(
    prophage_summary = prophage_summary,
    tool_summary = tool_summary,
    taxonomy_summary = taxonomy_summary_phylum,  # For plots
    taxonomy_summary_all = taxonomy_summary,      # All levels
    quality_summary = quality_summary
  ))
}

# Function to create visualizations
create_visualizations <- function(prophage_data, checkv_data, summary_stats) {
  plots <- list()
  
  # Set a clean theme with white background
  theme_set(theme_bw(base_size = 12))
  
  # 1. Prophages per sample
  plots$prophages_per_sample <- ggplot(summary_stats$prophage_summary, aes(x = reorder(sample, total_prophages), y = total_prophages)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = "Total Prophages Detected per Sample",
         x = "Sample", y = "Number of Prophages") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  # 2. Tool comparison
  plots$tool_comparison <- ggplot(summary_stats$tool_summary, aes(x = tool, y = prophages_detected, fill = tool)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0(prophages_detected, "\n(", percentage, "%)")), 
              vjust = -0.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Prophage Detection by Tool",
         x = "Detection Tool", y = "Number of Prophages") +
    theme_bw() +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  # 3. Tool comparison by sample (stacked)
  tool_by_sample <- prophage_data %>%
    group_by(sample, tool) %>%
    summarise(count = n(), .groups = "drop")
  
  plots$tool_by_sample <- ggplot(tool_by_sample, aes(x = sample, y = count, fill = tool)) +
    geom_col(position = "stack", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Prophage Detection Tools by Sample",
         x = "Sample", y = "Number of Prophages", fill = "Tool") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  
  # 4. Host taxonomy distribution (top 10)
  if (nrow(summary_stats$taxonomy_summary) > 0) {
    plots$taxonomy_distribution <- ggplot(summary_stats$taxonomy_summary, aes(x = reorder(phylum, n), y = n)) +
      geom_col(fill = "forestgreen", alpha = 0.8) +
      coord_flip() +
      labs(title = "Top 10 Host Phyla for Prophages",
           x = "Phylum", y = "Number of Prophages") +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
  }
  
  # 5. CheckV quality distribution
  if (nrow(summary_stats$quality_summary) > 0) {
    plots$quality_distribution <- ggplot(summary_stats$quality_summary, aes(x = checkv_quality, y = count, fill = checkv_quality)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = paste0(count, "\n(", percentage, "%)")), 
                vjust = -0.5) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      labs(title = "Prophage Quality Distribution (CheckV)",
           x = "Quality Category", y = "Number of Prophages") +
      theme_bw() +
      theme(legend.position = "none", 
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank())
  }
  
  # 6. Prophage length distribution
  if ("provirus_length" %in% colnames(checkv_data)) {
    plots$length_distribution <- ggplot(checkv_data, aes(x = provirus_length)) +
      geom_histogram(bins = 30, fill = "darkorange", alpha = 0.8, color = "black") +
      labs(title = "Prophage Length Distribution",
           x = "Prophage Length (bp)", y = "Count") +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      scale_x_log10(labels = scales::comma)
  }
  
  return(plots)
}

# Function to generate HTML report
generate_html_report <- function(outdir, prophage_data, checkv_data, summary_stats, plots) {
  message("Generating HTML report...")
  
  report_file <- file.path(outdir, "prophage_summary_report.html")
  
  # Create HTML content
  html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
<title>Prophage Pipeline Summary Report</title>
<style>
body { font-family: Arial, sans-serif; margin: 40px; }
.summary-box { background-color: #f0f0f0; padding: 15px; margin: 20px 0; border-radius: 5px; }
.plot-container { margin: 30px 0; text-align: center; }
table { border-collapse: collapse; width: 100%; margin: 20px 0; }
th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
th { background-color: #f2f2f2; }
</style>
</head>
<body>

<h1>Prophage Detection Pipeline Summary Report</h1>
<p>Generated on: ', Sys.Date(), '</p>

<div class="summary-box">
<h2>Overview</h2>
<ul>
<li><strong>Total Samples:</strong> ', length(unique(prophage_data$sample)), '</li>
<li><strong>Total Prophages Detected:</strong> ', nrow(prophage_data), '</li>
<li><strong>Total Unique Contigs with Prophages:</strong> ', n_distinct(prophage_data$contig), '</li>
<li><strong>Average Prophages per Sample:</strong> ', round(nrow(prophage_data) / length(unique(prophage_data$sample)), 1), '</li>
</ul>
</div>

<h2>Sample Summary</h2>
', kable(summary_stats$prophage_summary, format = "html", table.attr = 'class="summary-table"') %>%
  kable_styling(bootstrap_options = c("striped", "hover")), '

<h2>Detection Tool Performance</h2>
', kable(summary_stats$tool_summary, format = "html", table.attr = 'class="summary-table"') %>%
  kable_styling(bootstrap_options = c("striped", "hover")), '

</body>
</html>')
  
  writeLines(html_content, report_file)
  
  message(paste("HTML report saved to:", report_file))
  return(report_file)
}

# Main function
main <- function(outdir) {
  message("Starting prophage pipeline summary analysis...")
  message(paste("Output directory:", outdir))
  
  # Collect all data
  prophage_data <- collect_prophage_data(outdir)
  checkv_data <- collect_checkv_data(outdir)
  checkm_data <- collect_checkm_data(outdir)
  coverage_data <- collect_coverage_data(outdir)
  
  if (nrow(prophage_data) == 0) {
    stop("No prophage data found. Please check that the pipeline has completed successfully.")
  }
  
  # Create summary statistics
  summary_stats <- create_summary_stats(prophage_data, checkv_data, checkm_data, coverage_data)
  
  # Create visualizations
  plots <- create_visualizations(prophage_data, checkv_data, summary_stats)
  
  # Save plots
  plot_dir <- file.path(outdir, "summary_plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (plot_name in names(plots)) {
    if (!is.null(plots[[plot_name]])) {
      ggsave(
        filename = file.path(plot_dir, paste0(plot_name, ".png")),
        plot = plots[[plot_name]],
        width = 10, height = 6, dpi = 300
      )
      message(paste("Saved plot:", plot_name))
    }
  }
  
  # Generate HTML report
  report_file <- generate_html_report(outdir, prophage_data, checkv_data, summary_stats, plots)
  
  # Save data summaries as TSV
  write_tsv(summary_stats$prophage_summary, file.path(outdir, "prophage_summary_by_sample.tsv"))
  write_tsv(summary_stats$tool_summary, file.path(outdir, "tool_detection_summary.tsv"))
  
  # Save original phylum summary for backward compatibility
  if (nrow(summary_stats$taxonomy_summary) > 0) {
    write_tsv(summary_stats$taxonomy_summary, file.path(outdir, "host_taxonomy_summary.tsv"))
  }
  
  # Save all taxonomic levels in a single file
  if (length(summary_stats$taxonomy_summary_all) > 0) {
    # Combine all levels into one dataframe
    all_taxonomy <- bind_rows(summary_stats$taxonomy_summary_all)
    write_tsv(all_taxonomy, file.path(outdir, "host_taxonomy_all_levels.tsv"))
    
    # Also save each level separately for convenience
    for (level_name in names(summary_stats$taxonomy_summary_all)) {
      if (nrow(summary_stats$taxonomy_summary_all[[level_name]]) > 0) {
        write_tsv(
          summary_stats$taxonomy_summary_all[[level_name]], 
          file.path(outdir, paste0("host_taxonomy_", level_name, ".tsv"))
        )
      }
    }
  }
  
  message("Summary analysis completed!")
  message(paste("Results saved to:", outdir))
  message(paste("HTML report:", report_file))
  message(paste("Plots directory:", plot_dir))
  
  return(list(
    prophage_data = prophage_data,
    summary_stats = summary_stats,
    plots = plots,
    report_file = report_file
  ))
}

# Run if called as script or from Snakemake
if (!interactive()) {
  # Check if running from Snakemake
  if (exists("snakemake")) {
    # Extract output directory from snakemake output paths
    outdir <- dirname(snakemake@output[["report"]])
    
    # Run main analysis
    result <- main(outdir)
    
    # Ensure outputs match Snakemake expectations
    if (!file.exists(snakemake@output[["report"]])) {
      stop("HTML report was not created successfully")
    }
    
  } else {
    # Command line usage
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
      stop("Usage: Rscript create_summary_report.R <output_directory>")
    }
    
    main(args[1])
  }
}