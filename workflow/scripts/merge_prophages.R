library(readr)
library(tidyverse)
library(Biostrings)

# Function to check if 2 regions overlap
check_overlap <- function(region1, region2) {
  if (region1[2] >= region2[1] && region1[1] <= region2[2]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

genomad_path <- file.path(snakemake@input[["genomad"]], 
                                          "final_filtered_contigs_find_proviruses", 
                                          "final_filtered_contigs_provirus.tsv")
genomad <- read_tsv(genomad_path) %>%
  # Extract contig number from format: NovaSeq_N1028_metagenomics_I14076_FGT_Metagenomic_FRESH_41460_NODE_2333_length_43637_cov_146.0889
  extract(source_seq, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
  as.data.frame() %>%
  select(contig, start, end) %>%
  mutate(tool = "genomad")


phispy_path <- file.path(snakemake@input[["phispy"]],
                         "prophage.tsv") 
phispy <- read_tsv(phispy_path) %>%
  separate_wider_delim('Prophage number', '_', names=c('pp', 'pp_number')) %>%
  # Extract contig number from format: contig_27
  extract(Contig, into = "contig", regex = "contig_(\\d+)", remove = FALSE) %>%
  as.data.frame()

# Debugging: print column names after separation
print("Column names in phispy after separation:")
print(colnames(phispy))

# Ensure columns 'Start' and 'Stop' exist before renaming
if (!all(c('Start', 'Stop') %in% colnames(phispy))) {
  stop("Columns 'Start' and 'Stop' not found in phispy data after separation.")
  }

phispy <- phispy %>%
  dplyr::rename(start = Start, end = Stop) %>%
  mutate(tool = "phispy")

# Debugging: Print column names after renaming
print("Column names in phispy after renaming:")
print(colnames(phispy))

# Set up empty data frame to hold unique-to-phispy prophages
phispy_unique <- data.frame(matrix(ncol=4, nrow=0))
cols <- c("contig", "start", "end", "tool")
colnames(phispy_unique) <- cols
shared_contigs <- intersect(genomad$contig, phispy$contig)

# Get all phispy prophages on unique contigs
for (row in 1:nrow(phispy)) {
  if (!(phispy$contig[row] %in% genomad$contig)) {
    phispy_unique <- rbind(phispy_unique, phispy[row, ])
  }
}

# Get phispy prophages on shared prophages but with no overlap with any genomad prophage
for (c in shared_contigs){
  genomad_shared <- genomad %>%
    filter(contig == c)
  phispy_shared <- phispy %>%
    filter(contig == c)
  for (i in 1:nrow(phispy_shared)){
    overlap_found <- FALSE
    for (j in 1:nrow(genomad_shared)){
      r1 <- c(phispy_shared$start[i], phispy_shared$end[i])
      r2 <- c(genomad_shared$start[j], genomad_shared$end[j])
      if (check_overlap(r1, r2) == TRUE){
        overlap_found <- TRUE
        break
      }
    }
    if (overlap_found == FALSE){
      phispy_unique <- bind_rows(phispy_unique, phispy_shared[i, ])
    }
  }
}


# Get fasta file of unique-to-phispy contigs
pp_num <- as.numeric(phispy_unique$pp_number)
fasta <- readDNAStringSet(file.path(snakemake@input[["phispy"]], 
                                          "phage.fasta"))
fasta_unique <- fasta[pp_num]
writeXStringSet(fasta_unique, filepath=snakemake@output[["fasta"]])


phispy_unique <- phispy_unique %>%
  select(contig, start, end, tool)

# Get taxonomy output - handle both MMseqs2 and GTDB-Tk formats
if ("mmseqs" %in% names(snakemake@input)) {
  # MMseqs2 taxonomy parsing
  mmseqs_path <- file.path(snakemake@input[["mmseqs"]], "contig.taxonomy")
  
  # Debugging: print the first few lines of the MMseqs2 file
  print("First few lines of MMseqs2 taxonomy file:")
  print(readLines(mmseqs_path, n = 3))
  
  taxonomy_data <- read_tsv(mmseqs_path, col_names = c("contig_full", "taxid", "rank", "name", "retained", "assigned", "agreement", "confidence", "lineage", "lineage_names")) %>%
    # Extract contig number from format: NovaSeq_N1028_metagenomics_I14076_FGT_Metagenomic_FRESH_41460_NODE_22_length_1060_cov_5.0000
    extract(contig_full, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%  
    as.data.frame() %>%
    select(contig, lineage) %>%
    separate_wider_delim(lineage, ';', names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
    select(contig, superkingdom, phylum, class, order, family, genus, species)
    
  print("Using MMseqs2 taxonomy data")
  
} else if ("gtdbtk" %in% names(snakemake@input)) {
  # GTDB-Tk taxonomy parsing
  gtdbtk_path <- file.path(snakemake@input[["gtdbtk"]], "gtdbtk.bac120.summary.tsv")
  
  # Check if archaeal file exists and combine with bacterial
  gtdbtk_ar_path <- file.path(snakemake@input[["gtdbtk"]], "gtdbtk.ar53.summary.tsv")
  
  print("Reading GTDB-Tk taxonomy file")
  
  taxonomy_data <- data.frame()
  
  # Read bacterial taxonomy if file exists
  if (file.exists(gtdbtk_path)) {
    gtdbtk_bac <- read_tsv(gtdbtk_path) %>%
      # Extract contig number from user_genome column
      extract(user_genome, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
      as.data.frame() %>%
      select(contig, classification) %>%
      separate_wider_delim(classification, ';', names=c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      rename(superkingdom = domain) %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species)
    taxonomy_data <- rbind(taxonomy_data, gtdbtk_bac)
  }
  
  # Read archaeal taxonomy if file exists
  if (file.exists(gtdbtk_ar_path)) {
    gtdbtk_ar <- read_tsv(gtdbtk_ar_path) %>%
      extract(user_genome, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
      as.data.frame() %>%
      select(contig, classification) %>%
      separate_wider_delim(classification, ';', names=c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      rename(superkingdom = domain) %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species)
    taxonomy_data <- rbind(taxonomy_data, gtdbtk_ar)
  }
  
  print("Using GTDB-Tk taxonomy data")
} else {
  stop("Neither mmseqs nor gtdbtk input found in snakemake inputs")
}

# Debugging: print the structure of the parsed taxonomy
print("Structure of parsed taxonomy data:")
print(str(taxonomy_data))
print("First few rows of taxonomy_data:")
print(head(taxonomy_data))

# Merge the tables
final_prophage_table <- rbind(genomad, phispy_unique)

final_prophage_table_tax <- rbind(genomad, phispy_unique) %>%
  mutate(contig = as.numeric(contig)) %>%
  merge(taxonomy_data, by='contig', all.x = TRUE)
  
write.table(final_prophage_table, snakemake@output[["table"]], row.names=FALSE, sep="\t", quote=FALSE)
write.table(final_prophage_table_tax, snakemake@output[["table_with_taxonomy"]], row.names=FALSE, sep="\t", quote=FALSE)
