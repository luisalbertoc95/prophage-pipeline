library(readr)
library(tidyverse)

cat("=== ADDING TAXONOMY TO PROPHAGE TABLE ===\n")

# Read the basic prophage table
cat("1. Reading basic prophage table...\n")
final_prophage_table <- read_tsv(snakemake@input[["basic_table"]])
cat("Basic prophage table:", nrow(final_prophage_table), "prophages\n")
print(final_prophage_table)

# Get taxonomy data
cat("\n2. Processing taxonomy data...\n")
if ("taxonomy" %in% names(snakemake@input)) {
  cat("Using MMseqs taxonomy...\n")
  # MMseqs taxonomy parsing
  mmseqs_path <- file.path(snakemake@input[["taxonomy"]], "contig.taxonomy")
  
  if (file.exists(mmseqs_path)) {
    taxonomy_data <- read_tsv(mmseqs_path, col_names = c("contig_full", "taxid", "rank", "name", "retained", "assigned", "agreement", "confidence", "lineage", "lineage_names")) %>%
      # Extract contig number using extract() - revert to working method
      extract(contig_full, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
      filter(!is.na(contig)) %>%
      select(contig, lineage) %>%
      separate_wider_delim(lineage, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species)
    
    cat("MMseqs taxonomy loaded:", nrow(taxonomy_data), "records\n")
    cat("Unique contigs in taxonomy:", length(unique(taxonomy_data$contig)), "\n")
    
    # Show taxonomy summary
    tax_summary <- taxonomy_data %>%
      summarise(
        total_contigs = n(),
        with_superkingdom = sum(!is.na(superkingdom) & superkingdom != "unknown"),
        with_species = sum(!is.na(species) & species != "unknown" & !str_detect(species, "^uc_"))
      )
    cat("Taxonomy summary: ", tax_summary$total_contigs, " contigs,", 
        tax_summary$with_superkingdom, " with superkingdom,", 
        tax_summary$with_species, " with species\n")
  } else {
    cat("Warning: MMseqs taxonomy file not found, creating empty taxonomy\n")
    taxonomy_data <- data.frame(contig = character(), superkingdom = character(), phylum = character(), 
                               class = character(), order = character(), family = character(), 
                               genus = character(), species = character())
  }
} else {
  cat("Warning: No taxonomy input specified, creating empty taxonomy\n")
  taxonomy_data <- data.frame(contig = character(), superkingdom = character(), phylum = character(), 
                             class = character(), order = character(), family = character(), 
                             genus = character(), species = character())
}

# Join taxonomy with prophage table
cat("\n3. Joining taxonomy with prophage data...\n")
final_prophage_table_tax <- final_prophage_table %>%
  mutate(contig = as.character(contig)) %>%
  left_join(taxonomy_data %>% mutate(contig = as.character(contig)), by = 'contig')

cat("Final prophage table with taxonomy:", nrow(final_prophage_table_tax), "prophages\n")
cat("Prophages with taxonomy:", sum(!is.na(final_prophage_table_tax$superkingdom) & final_prophage_table_tax$superkingdom != "unknown"), "of", nrow(final_prophage_table_tax), "\n")

# Show breakdown by tool
tool_summary <- final_prophage_table_tax %>%
  group_by(tool) %>%
  summarise(
    total = n(),
    with_taxonomy = sum(!is.na(superkingdom) & superkingdom != "unknown"),
    .groups = 'drop'
  )
cat("Taxonomy by tool:\n")
for(i in 1:nrow(tool_summary)) {
  cat("  ", tool_summary$tool[i], ":", tool_summary$with_taxonomy[i], "/", tool_summary$total[i], "\n")
}

# Write output file
write.table(final_prophage_table_tax, snakemake@output[["table_with_taxonomy"]], row.names=FALSE, sep="\t", quote=FALSE)

cat("\n=== TAXONOMY ADDITION COMPLETE ===\n")
cat("- Output table with taxonomy:", snakemake@output[["table_with_taxonomy"]], "\n")