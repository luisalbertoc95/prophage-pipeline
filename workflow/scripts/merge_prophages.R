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
    separate_wider_delim(lineage, ";", names=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
    select(contig, superkingdom, phylum, class, order, family, genus, species)
    
  print("Using MMseqs2 taxonomy data")
  
} else if ("gtdbtk" %in% names(snakemake@input)) {
  # GTDB-Tk bin-level taxonomy parsing
  gtdbtk_path <- file.path(snakemake@input[["gtdbtk"]], "gtdbtk.bac120.summary.tsv")
  gtdbtk_ar_path <- file.path(snakemake@input[["gtdbtk"]], "gtdbtk.ar53.summary.tsv")
  gtdbtk_genomes_dir <- file.path(snakemake@input[["gtdbtk"]], "genomes")
  
  print("Reading GTDB-Tk bin-level taxonomy files")
  
  # Function to create NODE->contig mapping from Bakta input/output
  create_node_contig_mapping <- function() {
    node_contig_mapping <- data.frame(node_num = character(), contig_num = character(), stringsAsFactors = FALSE)
    
    # Path to original contigs file (Bakta input)
    orig_input_path <- file.path(dirname(dirname(snakemake@input[["gtdbtk"]])), "binning", "final_filt_contigs_5000.fasta")
    # Path to Bakta output with renumbered contigs
    bakta_output_path <- file.path(dirname(dirname(snakemake@input[["gtdbtk"]])), "phage_analysis", "bakta", "final_filt_contigs_5000.fna")
    
    if (file.exists(orig_input_path) && file.exists(bakta_output_path)) {
      # Read original headers to get NODE names in order
      orig_lines <- readLines(orig_input_path)
      orig_headers <- orig_lines[grepl("^>", orig_lines)]
      
      # Read Bakta headers to get contig_X names in same order
      bakta_lines <- readLines(bakta_output_path)
      bakta_headers <- bakta_lines[grepl("^>", bakta_lines)]
      
      # Create mapping: NODE number -> contig number
      for (i in seq_along(orig_headers)) {
        if (i <= length(bakta_headers)) {
          # Extract NODE number from original header
          node_match <- regmatches(orig_headers[i], regexpr("NODE_(\\d+)_", orig_headers[i], perl = TRUE))
          if (length(node_match) > 0) {
            node_num <- gsub("NODE_(\\d+)_", "\\1", node_match, perl = TRUE)
            contig_num <- as.character(i)  # Bakta uses contig_1, contig_2, etc. (1-indexed)
            node_contig_mapping <- rbind(node_contig_mapping, data.frame(node_num = node_num, contig_num = contig_num, stringsAsFactors = FALSE))
          }
        }
      }
      
      print(paste("Created NODE->contig mapping with", nrow(node_contig_mapping), "entries"))
    } else {
      print("Warning: Could not find files for NODE->contig mapping")
    }
    
    return(node_contig_mapping)
  }

  # Function to create bin-to-contig mapping using NODE->contig translation
  create_bin_contig_mapping <- function(genomes_dir, node_contig_mapping) {
    bin_contig_map <- data.frame(bin_file = character(), contig = character(), stringsAsFactors = FALSE)
    
    if (dir.exists(genomes_dir)) {
      bin_files <- list.files(genomes_dir, pattern = "\\.fa$", full.names = TRUE)
      
      for (bin_file in bin_files) {
        bin_name <- basename(bin_file)
        # Read FASTA headers to get contig names
        fasta_lines <- readLines(bin_file)
        header_lines <- fasta_lines[grepl("^>", fasta_lines)]
        
        for (header in header_lines) {
          # Extract NODE number from bin file header
          contig_match <- regmatches(header, regexpr("NODE_(\\d+)_", header, perl = TRUE))
          if (length(contig_match) > 0) {
            node_num <- gsub("NODE_(\\d+)_", "\\1", contig_match, perl = TRUE)
            # Look up corresponding contig number using NODE->contig mapping
            contig_mapping <- node_contig_mapping[node_contig_mapping$node_num == node_num, ]
            if (nrow(contig_mapping) > 0) {
              contig_num <- contig_mapping$contig_num[1]
              bin_contig_map <- rbind(bin_contig_map, data.frame(bin_file = bin_name, contig = contig_num, stringsAsFactors = FALSE))
            }
          }
        }
      }
    }
    return(bin_contig_map)
  }
  
  # Create NODE->contig mapping first
  node_contig_mapping <- create_node_contig_mapping()
  
  # Create mapping from bins to contigs using NODE->contig translation
  bin_contig_mapping <- create_bin_contig_mapping(gtdbtk_genomes_dir, node_contig_mapping)
  
  taxonomy_data <- data.frame()
  
  # Read bacterial taxonomy if file exists
  if (file.exists(gtdbtk_path) && file.size(gtdbtk_path) > 0) {
    gtdbtk_bac <- read_tsv(gtdbtk_path, show_col_types = FALSE) %>%
      as.data.frame() %>%
      select(user_genome, classification) %>%
      # Join with bin-to-contig mapping (add .fa extension to match bin files)
      mutate(user_genome_fa = paste0(user_genome, ".fa")) %>%
      left_join(bin_contig_mapping, by = c("user_genome_fa" = "bin_file")) %>%
      select(contig, classification) %>%
      filter(!is.na(contig)) %>%
      separate_wider_delim(classification, ";", names=c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      # Remove GTDB prefixes (d__, p__, c__, o__, f__, g__, s__)
      mutate(
        domain = gsub("^d__", "", domain),
        phylum = gsub("^p__", "", phylum),
        class = gsub("^c__", "", class),
        order = gsub("^o__", "", order),
        family = gsub("^f__", "", family),
        genus = gsub("^g__", "", genus),
        species = gsub("^s__", "", species)
      ) %>%
      rename(superkingdom = domain) %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species)
    taxonomy_data <- rbind(taxonomy_data, gtdbtk_bac)
  }
  
  # Read archaeal taxonomy if file exists
  if (file.exists(gtdbtk_ar_path) && file.size(gtdbtk_ar_path) > 0) {
    gtdbtk_ar <- read_tsv(gtdbtk_ar_path, show_col_types = FALSE) %>%
      as.data.frame() %>%
      select(user_genome, classification) %>%
      # Join with bin-to-contig mapping (add .fa extension to match bin files)
      mutate(user_genome_fa = paste0(user_genome, ".fa")) %>%
      left_join(bin_contig_mapping, by = c("user_genome_fa" = "bin_file")) %>%
      select(contig, classification) %>%
      filter(!is.na(contig)) %>%
      separate_wider_delim(classification, ";", names=c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), too_few = "align_start", too_many = "drop") %>%
      # Remove GTDB prefixes (d__, p__, c__, o__, f__, g__, s__)
      mutate(
        domain = gsub("^d__", "", domain),
        phylum = gsub("^p__", "", phylum),
        class = gsub("^c__", "", class),
        order = gsub("^o__", "", order),
        family = gsub("^f__", "", family),
        genus = gsub("^g__", "", genus),
        species = gsub("^s__", "", species)
      ) %>%
      rename(superkingdom = domain) %>%
      select(contig, superkingdom, phylum, class, order, family, genus, species)
    taxonomy_data <- rbind(taxonomy_data, gtdbtk_ar)
  }
  
  print(paste("Using GTDB-Tk bin-level taxonomy data with", nrow(taxonomy_data), "classified contigs"))
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
