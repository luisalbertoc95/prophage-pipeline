library(readr)
library(tidyverse)
library(Biostrings)

cat("=== BEDTOOLS-BASED PROPHAGE MERGE ===\n")

# Read the merged BED file (already contains overlap detection results)
cat("1. Reading merged prophage predictions from BEDtools...\n")
merged_bed <- read_tsv(snakemake@input[["merged_bed"]], 
                      col_names = c("contig", "start", "end", "tool", "id"),
                      col_types = cols(
                        contig = col_character(),
                        start = col_integer(),
                        end = col_integer(), 
                        tool = col_character(),
                        id = col_character()
                      ))

cat("Merged prophage predictions:")
print(merged_bed)
cat("\n")

# Extract unique PhiSpy FASTA sequences
cat("2. Extracting unique PhiSpy prophage sequences...\n")
if (file.exists(snakemake@input[["phispy_unique_ids"]]) && 
    file.size(snakemake@input[["phispy_unique_ids"]]) > 0) {
  
  # Read the unique PhiSpy IDs from BEDtools output
  unique_ids <- readLines(snakemake@input[["phispy_unique_ids"]])
  cat("Unique PhiSpy IDs to extract:", paste(unique_ids, collapse=", "), "\n")
  
  # Read PhiSpy FASTA file
  fasta_path <- file.path(snakemake@input[["phispy"]], "phage.fasta")
  if (file.exists(fasta_path)) {
    fasta <- readDNAStringSet(fasta_path)
    cat("Total PhiSpy sequences available:", length(fasta), "\n")
    cat("FASTA sequence names (first 3):", paste(names(fasta)[1:min(3, length(fasta))], collapse=", "), "...\n")
    
    # Get unique PhiSpy bed entries to match by contig and coordinates
    unique_bed <- merged_bed %>% filter(tool == "phispy")
    cat("Unique PhiSpy predictions to extract:", nrow(unique_bed), "\n")
    
    # Match FASTA sequences by contig number (more robust than pp numbers)
    matched_sequences <- DNAStringSet()
    for (i in 1:nrow(unique_bed)) {
      contig_num <- unique_bed$contig[i]
      start_coord <- unique_bed$start[i]
      end_coord <- unique_bed$end[i]
      
      # Look for FASTA headers containing this contig number
      pattern <- paste0("NODE_", contig_num, "_")
      matching_idx <- grep(pattern, names(fasta))
      
      if (length(matching_idx) > 0) {
        # If multiple matches, try to find one with matching coordinates
        coord_pattern <- paste0("_", start_coord, "_", end_coord)
        coord_matches <- grep(coord_pattern, names(fasta)[matching_idx])
        
        if (length(coord_matches) > 0) {
          # Found exact coordinate match
          final_idx <- matching_idx[coord_matches[1]]
          matched_sequences <- c(matched_sequences, fasta[final_idx])
          cat("Matched contig", contig_num, "with coordinates", start_coord, "-", end_coord, "\n")
        } else {
          # Take first contig match (coordinates might be slightly different)
          final_idx <- matching_idx[1]
          matched_sequences <- c(matched_sequences, fasta[final_idx])
          cat("Matched contig", contig_num, "(coordinates may differ)\n")
        }
      } else {
        cat("Warning: No FASTA sequence found for contig", contig_num, "\n")
      }
    }
    
    if (length(matched_sequences) > 0) {
      writeXStringSet(matched_sequences, filepath=snakemake@output[["fasta"]])
      cat("Successfully wrote", length(matched_sequences), "unique PhiSpy sequences to FASTA\n")
    } else {
      # Create empty FASTA if no sequences found
      writeXStringSet(DNAStringSet(), filepath=snakemake@output[["fasta"]])
      cat("No matching PhiSpy sequences found - wrote empty FASTA\n")
    }
  } else {
    stop("PhiSpy FASTA file not found: ", fasta_path)
  }
} else {
  # No unique PhiSpy sequences - create empty FASTA
  writeXStringSet(DNAStringSet(), filepath=snakemake@output[["fasta"]])
  cat("No unique PhiSpy sequences - wrote empty FASTA\n")
}

# Get taxonomy data (keep existing complex logic)
cat("\n3. Processing taxonomy data...\n")
if ("mmseqs" %in% names(snakemake@input)) {
  cat("Using MMseqs2 taxonomy...\n")
  # MMseqs2 taxonomy parsing
  mmseqs_path <- file.path(snakemake@input[["taxonomy"]], "contig.taxonomy")
  
  if (file.exists(mmseqs_path)) {
    taxonomy_data <- read_tsv(mmseqs_path, col_names = c("contig_full", "taxid", "rank", "name", "retained", "assigned", "agreement", "confidence", "lineage", "lineage_names")) %>%
      # Extract contig number using str_extract for consistency
      mutate(contig = str_extract(contig_full, "(?<=NODE_)\\d+(?=_)")) %>%
      filter(!is.na(contig)) %>%
      select(contig, lineage) %>%
      separate_wider_delim(lineage, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species)
    
    cat("MMseqs2 taxonomy loaded:", nrow(taxonomy_data), "records\n")
  } else {
    cat("Warning: MMseqs2 taxonomy file not found, creating empty taxonomy\n")
    taxonomy_data <- data.frame(contig = character(), superkingdom = character(), phylum = character(), 
                               class = character(), order = character(), family = character(), 
                               genus = character(), species = character())
  }
  
} else if ("gtdbtk" %in% names(snakemake@input)) {
  cat("Using GTDB-Tk taxonomy...\n")
  # GTDB-Tk bin-level taxonomy parsing (keep existing complex logic)
  gtdbtk_path <- file.path(snakemake@input[["taxonomy"]], "gtdbtk.bac120.summary.tsv")
  gtdbtk_ar_path <- file.path(snakemake@input[["taxonomy"]], "gtdbtk.ar53.summary.tsv")
  gtdbtk_genomes_dir <- file.path(snakemake@input[["taxonomy"]], "genomes")
  
  # Function to create bin-to-contig mapping from FASTA files
  create_bin_contig_mapping <- function(genomes_dir) {
    bin_contig_map <- data.frame(bin_file = character(), contig = character(), stringsAsFactors = FALSE)
    
    if (dir.exists(genomes_dir)) {
      bin_files <- list.files(genomes_dir, pattern = "\\.fa$", full.names = TRUE)
      
      for (bin_file in bin_files) {
        bin_name <- basename(bin_file)
        if (file.exists(bin_file) && file.size(bin_file) > 0) {
          fasta_lines <- readLines(bin_file)
          header_lines <- fasta_lines[grepl("^>", fasta_lines)]
          
          for (header in header_lines) {
            # Extract contig number from header (NODE format)
            contig_match <- str_extract(header, "(?<=NODE_)\\d+(?=_)")
            if (!is.na(contig_match)) {
              bin_contig_map <- rbind(bin_contig_map, data.frame(bin_file = bin_name, contig = contig_match, stringsAsFactors = FALSE))
            }
          }
        }
      }
    }
    return(bin_contig_map)
  }
  
  # Create mapping from bins to contigs
  bin_contig_mapping <- create_bin_contig_mapping(gtdbtk_genomes_dir)
  cat("Created bin-to-contig mapping with", nrow(bin_contig_mapping), "entries\n")
  
  taxonomy_data <- data.frame()
  
  # Read bacterial taxonomy if file exists
  if (file.exists(gtdbtk_path) && file.size(gtdbtk_path) > 0) {
    gtdbtk_bac_joined <- read_tsv(gtdbtk_path, show_col_types = FALSE) %>%
      select(user_genome, classification) %>%
      mutate(user_genome_fa = paste0(user_genome, ".fa")) %>%
      left_join(bin_contig_mapping, by = c("user_genome_fa" = "bin_file")) %>%
      select(contig, classification) %>%
      filter(!is.na(contig) & !is.na(classification) & classification != "")
    
    cat("Found", nrow(gtdbtk_bac_joined), "bacterial classifications with matching contigs\n")
    
    if (nrow(gtdbtk_bac_joined) > 0) {
      gtdbtk_bac <- gtdbtk_bac_joined %>%
        separate_wider_delim(classification, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
        # Remove GTDB prefixes
        mutate(
          superkingdom = gsub("^d__", "", superkingdom),
          phylum = gsub("^p__", "", phylum),
          class = gsub("^c__", "", class),
          order = gsub("^o__", "", order),
          family = gsub("^f__", "", family),
          genus = gsub("^g__", "", genus),
          species = gsub("^s__", "", species)
        ) %>%
        select(contig, superkingdom, phylum, class, order, family, genus, species)
      taxonomy_data <- rbind(taxonomy_data, gtdbtk_bac)
    }
  }
  
  # Read archaeal taxonomy if file exists
  if (file.exists(gtdbtk_ar_path) && file.size(gtdbtk_ar_path) > 0) {
    gtdbtk_ar_joined <- read_tsv(gtdbtk_ar_path, show_col_types = FALSE) %>%
      select(user_genome, classification) %>%
      mutate(user_genome_fa = paste0(user_genome, ".fa")) %>%
      left_join(bin_contig_mapping, by = c("user_genome_fa" = "bin_file")) %>%
      select(contig, classification) %>%
      filter(!is.na(contig) & !is.na(classification) & classification != "")
    
    cat("Found", nrow(gtdbtk_ar_joined), "archaeal classifications with matching contigs\n")
    
    if (nrow(gtdbtk_ar_joined) > 0) {
      gtdbtk_ar <- gtdbtk_ar_joined %>%
        separate_wider_delim(classification, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
        # Remove GTDB prefixes
        mutate(
          superkingdom = gsub("^d__", "", superkingdom),
          phylum = gsub("^p__", "", phylum),
          class = gsub("^c__", "", class),
          order = gsub("^o__", "", order),
          family = gsub("^f__", "", family),
          genus = gsub("^g__", "", genus),
          species = gsub("^s__", "", species)
        ) %>%
        select(contig, superkingdom, phylum, class, order, family, genus, species)
      taxonomy_data <- rbind(taxonomy_data, gtdbtk_ar)
    }
  }
  
  cat("GTDB-Tk taxonomy loaded:", nrow(taxonomy_data), "records\n")
} else {
  cat("Warning: No taxonomy method specified, creating empty taxonomy\n")
  taxonomy_data <- data.frame(contig = character(), superkingdom = character(), phylum = character(), 
                             class = character(), order = character(), family = character(), 
                             genus = character(), species = character())
}

# Create output tables
cat("\n4. Creating output tables...\n")

# Basic table (coordinates only)
final_prophage_table <- merged_bed %>%
  select(contig, start, end, tool)

cat("Final prophage table:", nrow(final_prophage_table), "prophages\n")
print(final_prophage_table)

# Table with taxonomy
final_prophage_table_tax <- final_prophage_table %>%
  mutate(contig = as.character(contig)) %>%
  left_join(taxonomy_data %>% mutate(contig = as.character(contig)), by = 'contig')

cat("Final prophage table with taxonomy:", nrow(final_prophage_table_tax), "prophages\n")
cat("Prophages with taxonomy:", sum(!is.na(final_prophage_table_tax$superkingdom)), "of", nrow(final_prophage_table_tax), "\n")

# Write output files
write.table(final_prophage_table, snakemake@output[["table"]], row.names=FALSE, sep="\t", quote=FALSE)
write.table(final_prophage_table_tax, snakemake@output[["table_with_taxonomy"]], row.names=FALSE, sep="\t", quote=FALSE)

cat("\n=== BEDTOOLS MERGE COMPLETE ===\n")
cat("- Output table:", snakemake@output[["table"]], "\n")
cat("- Output table with taxonomy:", snakemake@output[["table_with_taxonomy"]], "\n") 
cat("- Output FASTA:", snakemake@output[["fasta"]], "\n")