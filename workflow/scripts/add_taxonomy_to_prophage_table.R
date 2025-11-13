library(readr)
library(tidyverse)
library(taxonomizr)

cat("=== HYBRID TAXONOMY INTEGRATION (GTDB-Tk + MMseqs) ===\n")

# Read the basic prophage table (now includes bin column)
cat("1. Reading basic prophage table...\n")
final_prophage_table <- read_tsv(snakemake@input[["basic_table"]])
cat("Basic prophage table:", nrow(final_prophage_table), "prophages\n")
cat("Columns:", paste(colnames(final_prophage_table), collapse=", "), "\n")
cat("Binned prophages:", sum(final_prophage_table$bin != "none"), "of", nrow(final_prophage_table), "\n")
cat("Unbinned prophages:", sum(final_prophage_table$bin == "none"), "of", nrow(final_prophage_table), "\n")

# Initialize taxonomy data structure
taxonomy_columns <- c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
empty_taxonomy <- data.frame(
  contig = character(),
  superkingdom = character(), phylum = character(), class = character(), 
  order = character(), family = character(), genus = character(), species = character(),
  taxonomy_source = character()
)

# Read GTDB-Tk results for binned contigs
cat("\n2. Reading GTDB-Tk taxonomy for binned contigs...\n")
gtdbtk_taxonomy <- empty_taxonomy

gtdbtk_path <- snakemake@input[["gtdbtk_taxonomy"]]
bac_summary <- file.path(gtdbtk_path, "gtdbtk.bac120.summary.tsv")
ar_summary <- file.path(gtdbtk_path, "gtdbtk.ar53.summary.tsv")

gtdbtk_data <- data.frame()

# Read bacterial results
if (file.exists(bac_summary) && file.size(bac_summary) > 0) {
  bac_data <- read_tsv(bac_summary) %>%
    filter(!is.na(classification)) %>%
    select(user_genome, classification) %>%
    rename(bin = user_genome, gtdbtk_classification = classification)
  gtdbtk_data <- bind_rows(gtdbtk_data, bac_data)
  cat("GTDB-Tk bacterial results:", nrow(bac_data), "bins\n")
}

# Read archaeal results  
if (file.exists(ar_summary) && file.size(ar_summary) > 0) {
  ar_data <- read_tsv(ar_summary) %>%
    filter(!is.na(classification)) %>%
    select(user_genome, classification) %>%
    rename(bin = user_genome, gtdbtk_classification = classification)
  gtdbtk_data <- bind_rows(gtdbtk_data, ar_data)
  cat("GTDB-Tk archaeal results:", nrow(ar_data), "bins\n")
}

# Parse GTDB-Tk classifications
if (nrow(gtdbtk_data) > 0) {
  gtdbtk_taxonomy <- gtdbtk_data %>%
    # Parse GTDB-Tk classification string (d__Bacteria;p__Firmicutes;...)
    separate(gtdbtk_classification, 
             into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), 
             sep = ";", fill = "right") %>%
    # Remove GTDB prefixes (d__, p__, c__, etc.)
    mutate(
      superkingdom = str_remove(domain, "^[a-z]__"),
      phylum = str_remove(phylum, "^[a-z]__"),
      class = str_remove(class, "^[a-z]__"), 
      order = str_remove(order, "^[a-z]__"),
      family = str_remove(family, "^[a-z]__"),
      genus = str_remove(genus, "^[a-z]__"),
      species = str_remove(species, "^[a-z]__"),
      taxonomy_source = "GTDB-Tk"
    ) %>%
    select(bin, superkingdom, phylum, class, order, family, genus, species, taxonomy_source)
  
  cat("GTDB-Tk taxonomy processed for", nrow(gtdbtk_taxonomy), "bins\n")
  cat("Sample GTDB-Tk results:\n")
  print(head(gtdbtk_taxonomy, 3))
} else {
  cat("No GTDB-Tk results found\n")
}

# Read MMseqs results for unbinned contigs  
cat("\n3. Reading MMseqs taxonomy for unbinned contigs...\n")
mmseqs_taxonomy <- empty_taxonomy

mmseqs_path <- file.path(snakemake@input[["mmseqs_taxonomy"]], "contig.taxonomy")
taxonomizr_db <- snakemake@params[["taxonomizr_db"]]

if (file.exists(mmseqs_path)) {
  # Read MMseqs output (9 columns: contig, taxid, rank, name, 4 numeric values, lineage_ids)
  mmseqs_raw <- read_tsv(mmseqs_path,
                         col_names = c("contig_full", "taxid", "rank", "name", "v1", "v2", "v3", "v4", "lineage_ids"),
                         show_col_types = FALSE) %>%
    # Extract contig number
    extract(contig_full, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
    filter(!is.na(contig))

  cat("MMseqs taxonomy loaded for", nrow(mmseqs_raw), "contigs\n")

  # Convert lineage IDs to lineage names using taxonomizr
  if (file.exists(taxonomizr_db)) {
    cat("Converting taxonomy IDs to names using taxonomizr...\n")

    # Split the lineage IDs and convert each to names
    mmseqs_taxonomy <- mmseqs_raw %>%
      rowwise() %>%
      mutate(
        # Split lineage IDs
        taxid_list = list(as.numeric(strsplit(lineage_ids, ";")[[1]])),
        # Get taxonomy for all IDs in the lineage
        lineage_tax = list(tryCatch({
          tax_info <- getTaxonomy(taxid_list, taxonomizr_db)
          # getTaxonomy returns a matrix/dataframe with one row per taxid
          # We want the most specific (last) entry for each rank
          if (nrow(tax_info) > 0) {
            # For each rank, get the last non-NA value
            sapply(taxonomy_columns, function(rank) {
              vals <- tax_info[, rank]
              vals <- vals[!is.na(vals)]
              if (length(vals) > 0) tail(vals, 1) else NA_character_
            })
          } else {
            setNames(rep(NA_character_, length(taxonomy_columns)), taxonomy_columns)
          }
        }, error = function(e) {
          setNames(rep(NA_character_, length(taxonomy_columns)), taxonomy_columns)
        }))
      ) %>%
      ungroup() %>%
      # Extract taxonomy columns
      mutate(
        superkingdom = sapply(lineage_tax, function(x) x["superkingdom"]),
        phylum = sapply(lineage_tax, function(x) x["phylum"]),
        class = sapply(lineage_tax, function(x) x["class"]),
        order = sapply(lineage_tax, function(x) x["order"]),
        family = sapply(lineage_tax, function(x) x["family"]),
        genus = sapply(lineage_tax, function(x) x["genus"]),
        species = sapply(lineage_tax, function(x) x["species"]),
        taxonomy_source = "MMseqs"
      ) %>%
      select(contig, all_of(taxonomy_columns), taxonomy_source)

    cat("Taxonomy ID conversion completed\n")
  } else {
    cat("Warning: Taxonomizr database not found at:", taxonomizr_db, "\n")
    cat("Creating empty taxonomy for MMseqs results\n")
    mmseqs_taxonomy <- mmseqs_raw %>%
      mutate(
        superkingdom = NA_character_,
        phylum = NA_character_,
        class = NA_character_,
        order = NA_character_,
        family = NA_character_,
        genus = NA_character_,
        species = NA_character_,
        taxonomy_source = NA_character_
      ) %>%
      select(contig, all_of(taxonomy_columns), taxonomy_source)
  }

  cat("MMseqs taxonomy processed for", nrow(mmseqs_taxonomy), "contigs\n")
  
  # Show taxonomy summary
  mmseqs_summary <- mmseqs_taxonomy %>%
    summarise(
      total_contigs = n(),
      with_superkingdom = sum(!is.na(superkingdom) & superkingdom != "unknown"),
      with_species = sum(!is.na(species) & species != "unknown" & !str_detect(species, "^uc_"))
    )
  cat("MMseqs summary:", mmseqs_summary$total_contigs, "contigs,", 
      mmseqs_summary$with_superkingdom, "with superkingdom,", 
      mmseqs_summary$with_species, "with species\n")
} else {
  cat("Warning: MMseqs taxonomy file not found\n")
}

# Create hybrid taxonomy integration
cat("\n4. Integrating hybrid taxonomy...\n")

# First, create bin-to-taxonomy mapping from GTDB-Tk
bin_taxonomy <- gtdbtk_taxonomy %>%
  select(bin, all_of(taxonomy_columns), taxonomy_source)

# Join with prophage table
final_prophage_table_tax <- final_prophage_table %>%
  mutate(contig = as.character(contig)) %>%
  # First join with GTDB-Tk taxonomy by bin (for binned contigs)
  left_join(bin_taxonomy, by = 'bin') %>%
  # Then join with MMseqs taxonomy by contig (for unbinned contigs where GTDB-Tk is NA)
  left_join(
    mmseqs_taxonomy %>% mutate(contig = as.character(contig)) %>%
      rename_with(~paste0(.x, "_mmseqs"), all_of(taxonomy_columns)) %>%
      select(contig, all_of(paste0(taxonomy_columns, "_mmseqs")), taxonomy_source_mmseqs = taxonomy_source),
    by = 'contig'
  ) %>%
  # Use GTDB-Tk for binned contigs, MMseqs for unbinned contigs
  mutate(
    superkingdom = ifelse(is.na(superkingdom), superkingdom_mmseqs, superkingdom),
    phylum = ifelse(is.na(phylum), phylum_mmseqs, phylum),
    class = ifelse(is.na(class), class_mmseqs, class),
    order = ifelse(is.na(order), order_mmseqs, order),
    family = ifelse(is.na(family), family_mmseqs, family),
    genus = ifelse(is.na(genus), genus_mmseqs, genus),
    species = ifelse(is.na(species), species_mmseqs, species),
    taxonomy_source = ifelse(is.na(taxonomy_source), taxonomy_source_mmseqs, taxonomy_source)
  ) %>%
  # Clean up temporary columns
  select(-ends_with("_mmseqs"))

# Final statistics
cat("\n5. Hybrid taxonomy integration results...\n")
cat("Final prophage table with hybrid taxonomy:", nrow(final_prophage_table_tax), "prophages\n")

tax_stats <- final_prophage_table_tax %>%
  summarise(
    total = n(),
    with_gtdbtk = sum(!is.na(taxonomy_source) & taxonomy_source == "GTDB-Tk"),
    with_mmseqs = sum(!is.na(taxonomy_source) & taxonomy_source == "MMseqs"),
    with_any_taxonomy = sum(!is.na(taxonomy_source)),
    binned_contigs = sum(bin != "none"),
    unbinned_contigs = sum(bin == "none")
  )

cat("Taxonomy sources:\n")
cat("  GTDB-Tk (binned):", tax_stats$with_gtdbtk, "/", tax_stats$binned_contigs, "binned prophages\n")
cat("  MMseqs (unbinned):", tax_stats$with_mmseqs, "/", tax_stats$unbinned_contigs, "unbinned prophages\n") 
cat("  Total with taxonomy:", tax_stats$with_any_taxonomy, "/", tax_stats$total, "prophages\n")

# Show breakdown by tool and taxonomy source
tool_tax_summary <- final_prophage_table_tax %>%
  group_by(tool, taxonomy_source) %>%
  summarise(count = n(), .groups = 'drop')
cat("\nTaxonomy by detection tool and source:\n")
print(tool_tax_summary)

# Write output file
write.table(final_prophage_table_tax, snakemake@output[["table_with_taxonomy"]], 
            row.names=FALSE, sep="\t", quote=FALSE)

cat("\n=== HYBRID TAXONOMY INTEGRATION COMPLETE ===\n")
cat("- Output table with hybrid taxonomy:", snakemake@output[["table_with_taxonomy"]], "\n")
cat("- Schema: contig | start | end | tool | bin | superkingdom | phylum | class | order | family | genus | species | taxonomy_source\n")