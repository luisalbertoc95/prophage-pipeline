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

# Create basic prophage table (coordinates only)
cat("\n3. Creating basic prophage table...\n")
final_prophage_table <- merged_bed %>%
  select(contig, start, end, tool)

cat("Final prophage table:", nrow(final_prophage_table), "prophages\n")
print(final_prophage_table)

# Write output file
write.table(final_prophage_table, snakemake@output[["table"]], row.names=FALSE, sep="\t", quote=FALSE)

cat("\n=== BASIC PROPHAGE TABLE COMPLETE ===\n")
cat("- Output table:", snakemake@output[["table"]], "\n")
cat("- Output FASTA:", snakemake@output[["fasta"]], "\n")