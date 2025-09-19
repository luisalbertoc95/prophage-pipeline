#!/usr/bin/env Rscript

library(tidyverse)

# Test the regex extraction on actual data
test_string <- "NovaSeq_N1028_metagenomics_I14076_FGT_Metagenomic_FRESH_41460_NODE_1662_length_43637_cov_146.0889"

cat("Testing regex extraction methods:\n")
cat("Input string:", test_string, "\n\n")

# Method 1: str_extract with lookbehind/lookahead (current failing method)
method1 <- str_extract(test_string, "(?<=NODE_)\\d+(?=_)")
cat("Method 1 (str_extract with lookbehind/lookahead):", method1, "\n")

# Method 2: extract with capturing group (what we reverted to)
df <- data.frame(contig_full = test_string)
method2 <- df %>%
  extract(contig_full, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE)
cat("Method 2 (extract with capturing group):", method2$contig, "\n")

# Test on actual taxonomy data
cat("\n--- Testing on actual taxonomy file ---\n")
tax_data <- read_tsv("contig.taxonomy", col_names = c("contig_full", "taxid", "rank", "name", "retained", "assigned", "agreement", "confidence", "lineage", "lineage_names"))

cat("First few taxonomy entries:\n")
print(head(tax_data$contig_full, 3))

# Test extraction
tax_extracted <- tax_data %>%
  extract(contig_full, into = "contig", regex = "NODE_(\\d+)_", remove = FALSE) %>%
  filter(!is.na(contig))

cat("Extracted contig IDs from taxonomy:\n")
print(head(tax_extracted$contig, 10))
cat("Number of successfully extracted contigs:", nrow(tax_extracted), "out of", nrow(tax_data), "\n")

# Test join with prophage data
cat("\n--- Testing join with prophage data ---\n")
prophage_data <- read_tsv("final_prophage_table_with_host_taxonomy.tsv")
cat("Prophage contig IDs:\n")
print(unique(prophage_data$contig))

# Perform the join
joined <- prophage_data %>%
  select(contig, start, end, tool) %>%
  mutate(contig = as.character(contig)) %>%
  left_join(tax_extracted %>% 
            select(contig, lineage) %>%
            mutate(contig = as.character(contig)), 
            by = 'contig')

cat("Join results:\n")
cat("Rows with taxonomy data:", sum(!is.na(joined$lineage)), "out of", nrow(joined), "\n")
print(joined)