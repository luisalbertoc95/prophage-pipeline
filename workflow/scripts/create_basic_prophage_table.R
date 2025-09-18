library(readr)
library(tidyverse)
library(Biostrings)

cat("=== CREATING BASIC PROPHAGE TABLE ===\n")

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

# Read binning information from final DAS Tool bins
cat("\n3. Reading binning information from final DAS Tool bins...\n")
sample_name <- basename(dirname(dirname(snakemake@output[["table"]])))

bin_mapping <- data.frame(contig = character(), bin = character())

# Look for final DAS Tool bin files (copied to GTDB-Tk genomes directory)
dastool_bins_dir <- file.path(dirname(dirname(dirname(snakemake@output[["table"]]))), "taxonomy", "gtdbtk", "genomes")

if (dir.exists(dastool_bins_dir)) {
  # Find all bin FASTA files
  bin_files <- list.files(dastool_bins_dir, pattern = "bin\\..*\\.fa$", full.names = TRUE)
  
  if (length(bin_files) > 0) {
    cat("Found", length(bin_files), "final DAS Tool bins\n")
    
    # Extract contig information from each bin file
    for (bin_file in bin_files) {
      bin_name <- str_remove(basename(bin_file), "\\.fa$")  # e.g., "bin.1"
      
      # Read FASTA headers to get contig names
      fasta_lines <- readLines(bin_file)
      header_lines <- fasta_lines[grepl("^>", fasta_lines)]
      
      # Extract contig numbers from headers
      for (header in header_lines) {
        # Remove the ">" and extract NODE number
        contig_match <- str_extract(header, "NODE_(\\d+)_")
        if (!is.na(contig_match)) {
          contig_num <- str_extract(contig_match, "\\d+")
          if (!is.na(contig_num)) {
            bin_mapping <- rbind(bin_mapping, data.frame(contig = contig_num, bin = bin_name))
          }
        }
      }
    }
    
    cat("Final binning information:\n")
    cat("- Total binned contigs:", nrow(bin_mapping), "\n")
    cat("- Number of final bins:", length(unique(bin_mapping$bin)), "\n")
    cat("- Bin names:", paste(unique(bin_mapping$bin), collapse=", "), "\n")
    
    # Write mapping file for reference
    mapping_file <- file.path(dirname(snakemake@output[["table"]]), "final_dastool_contig2bin.tsv")
    write_tsv(bin_mapping, mapping_file)
    cat("- Mapping file written to:", mapping_file, "\n")
    
  } else {
    cat("No final DAS Tool bin files found in:", dastool_bins_dir, "\n")
  }
} else {
  cat("DAS Tool bins directory not found:", dastool_bins_dir, "\n")
  cat("Proceeding without binning information - all contigs will be marked as 'none'\n")
}

# Create basic prophage table with bin information
cat("\n4. Creating basic prophage table with bin information...\n")
final_prophage_table <- merged_bed %>%
  select(contig, start, end, tool) %>%
  mutate(contig = as.character(contig)) %>%
  left_join(bin_mapping %>% mutate(contig = as.character(contig)), by = 'contig') %>%
  mutate(bin = ifelse(is.na(bin), "none", bin))

cat("Final prophage table:", nrow(final_prophage_table), "prophages\n")
cat("Prophages in bins:", sum(final_prophage_table$bin != "none"), "of", nrow(final_prophage_table), "\n")
cat("Prophages unbinned:", sum(final_prophage_table$bin == "none"), "of", nrow(final_prophage_table), "\n")

# Show breakdown by tool and binning status
bin_tool_summary <- final_prophage_table %>%
  mutate(is_binned = bin != "none") %>%
  group_by(tool, is_binned) %>%
  summarise(count = n(), .groups = 'drop')
cat("Prophage distribution by tool and binning status:\n")
print(bin_tool_summary)

# Write output file
write.table(final_prophage_table, snakemake@output[["table"]], row.names=FALSE, sep="\t", quote=FALSE)

cat("\n=== BASIC PROPHAGE TABLE WITH BINNING COMPLETE ===\n")
cat("- Output table (with bin column):", snakemake@output[["table"]], "\n")
cat("- Output FASTA:", snakemake@output[["fasta"]], "\n")