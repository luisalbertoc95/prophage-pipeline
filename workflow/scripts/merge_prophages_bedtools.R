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

# Get taxonomy data using hybrid approach (GTDB-Tk + MMseqs)
cat("\n3. Processing taxonomy data with hybrid approach...\n")

# Helper function to create bin-to-contig mapping from FASTA files
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

# Load GTDB-Tk taxonomy (for binned contigs)
gtdbtk_taxonomy <- data.frame()
if ("gtdbtk" %in% names(snakemake@input)) {
  cat("Loading GTDB-Tk taxonomy for binned contigs...\n")
  gtdbtk_path <- file.path(snakemake@input[["gtdbtk"]], "gtdbtk.bac120.summary.tsv")
  gtdbtk_ar_path <- file.path(snakemake@input[["gtdbtk"]], "gtdbtk.ar53.summary.tsv")
  gtdbtk_genomes_dir <- file.path(snakemake@input[["gtdbtk"]], "genomes")
  
  # Create bin-to-contig mapping
  bin_contig_mapping <- create_bin_contig_mapping(gtdbtk_genomes_dir)
  cat("Created bin-to-contig mapping with", nrow(bin_contig_mapping), "entries\n")
  
  # Read bacterial taxonomy
  if (file.exists(gtdbtk_path) && file.size(gtdbtk_path) > 0) {
    gtdbtk_bac <- read_tsv(gtdbtk_path, show_col_types = FALSE) %>%
      select(user_genome, classification) %>%
      mutate(user_genome_fa = paste0(user_genome, ".fa")) %>%
      left_join(bin_contig_mapping, by = c("user_genome_fa" = "bin_file")) %>%
      select(contig, classification) %>%
      filter(!is.na(contig) & !is.na(classification) & classification != "") %>%
      separate_wider_delim(classification, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      mutate(
        superkingdom = gsub("^d__", "", superkingdom),
        phylum = gsub("^p__", "", phylum),
        class = gsub("^c__", "", class),
        order = gsub("^o__", "", order),
        family = gsub("^f__", "", family),
        genus = gsub("^g__", "", genus),
        species = gsub("^s__", "", species),
        source = "gtdbtk"
      ) %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species, source)
    
    gtdbtk_taxonomy <- rbind(gtdbtk_taxonomy, gtdbtk_bac)
    cat("Loaded", nrow(gtdbtk_bac), "bacterial classifications\n")
  }
  
  # Read archaeal taxonomy
  if (file.exists(gtdbtk_ar_path) && file.size(gtdbtk_ar_path) > 0) {
    gtdbtk_ar <- read_tsv(gtdbtk_ar_path, show_col_types = FALSE) %>%
      select(user_genome, classification) %>%
      mutate(user_genome_fa = paste0(user_genome, ".fa")) %>%
      left_join(bin_contig_mapping, by = c("user_genome_fa" = "bin_file")) %>%
      select(contig, classification) %>%
      filter(!is.na(contig) & !is.na(classification) & classification != "") %>%
      separate_wider_delim(classification, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      mutate(
        superkingdom = gsub("^d__", "", superkingdom),
        phylum = gsub("^p__", "", phylum),
        class = gsub("^c__", "", class),
        order = gsub("^o__", "", order),
        family = gsub("^f__", "", family),
        genus = gsub("^g__", "", genus),
        species = gsub("^s__", "", species),
        source = "gtdbtk"
      ) %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species, source)
    
    gtdbtk_taxonomy <- rbind(gtdbtk_taxonomy, gtdbtk_ar)
    cat("Loaded", nrow(gtdbtk_ar), "archaeal classifications\n")
  }
}

cat("Total GTDB-Tk taxonomy records:", nrow(gtdbtk_taxonomy), "\n")

# Load MMseqs taxonomy (for unbinned contigs)
mmseqs_taxonomy <- data.frame()
if ("mmseqs" %in% names(snakemake@input)) {
  cat("Loading MMseqs taxonomy for unbinned contigs...\n")
  mmseqs_path <- file.path(snakemake@input[["mmseqs"]], "contig.taxonomy")
  
  if (file.exists(mmseqs_path) && file.size(mmseqs_path) > 0) {
    mmseqs_taxonomy <- read_tsv(mmseqs_path, col_names = c("contig_full", "taxid", "rank", "name", "retained", "assigned", "agreement", "confidence", "lineage", "lineage_names")) %>%
      mutate(contig = str_extract(contig_full, "(?<=NODE_)\\d+(?=_)")) %>%
      filter(!is.na(contig)) %>%
      select(contig, lineage) %>%
      separate_wider_delim(lineage, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      mutate(source = "mmseqs") %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species, source)
    
    cat("Loaded", nrow(mmseqs_taxonomy), "MMseqs taxonomy records\n")
  } else {
    cat("MMseqs taxonomy file not found or empty\n")
  }
}

# Merge taxonomies with GTDB-Tk priority
cat("Merging taxonomies (GTDB-Tk priority for binned contigs)...\n")
taxonomy_data <- gtdbtk_taxonomy %>%
  full_join(mmseqs_taxonomy, by = "contig", suffix = c("_gtdb", "_mmseqs")) %>%
  mutate(
    # GTDB-Tk takes priority, MMseqs fills gaps
    superkingdom = coalesce(superkingdom_gtdb, superkingdom_mmseqs),
    phylum = coalesce(phylum_gtdb, phylum_mmseqs),
    class = coalesce(class_gtdb, class_mmseqs),
    order = coalesce(order_gtdb, order_mmseqs),
    family = coalesce(family_gtdb, family_mmseqs),
    genus = coalesce(genus_gtdb, genus_mmseqs),
    species = coalesce(species_gtdb, species_mmseqs),
    source = coalesce(source_gtdb, source_mmseqs)
  ) %>%
  select(contig, superkingdom, phylum, class, order, family, genus, species, source)

cat("Final hybrid taxonomy:", nrow(taxonomy_data), "records\n")
cat("From GTDB-Tk:", sum(taxonomy_data$source == "gtdbtk", na.rm = TRUE), "contigs\n")
cat("From MMseqs:", sum(taxonomy_data$source == "mmseqs", na.rm = TRUE), "contigs\n")

# Remove source column for final output (keep compatibility)
taxonomy_data <- taxonomy_data %>% select(-source)

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