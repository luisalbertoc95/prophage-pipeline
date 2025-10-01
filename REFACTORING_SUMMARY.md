# Prophage Pipeline Refactoring Summary

## Overview
Split the post-binning workflow into two separate paths (MAGs and unbinned contigs) to handle prophage prediction and taxonomy assignment appropriately for each data type.

## Key Changes

### New Files Created
1. **workflow/rules/prophage_mags.smk** - Prophage prediction for MAG bins
   - `bakta_per_mag` - Annotate each MAG individually
   - `phispy_per_mag` - Predict prophages in each MAG
   - `genomad_per_mag` - Predict prophages in each MAG with genomad
   - `get_mag_bins` - Checkpoint to dynamically discover bins
   - `merge_mag_prophages` - Merge genomad + unique phispy predictions (same bedtools logic)
   - `extract_mag_prophage_sequences` - Extract prophage FASTAs and create table

2. **workflow/rules/prophage_unbinned.smk** - Prophage and free phage prediction for unbinned contigs
   - `genomad_unbinned` - Predict prophages AND free phages
   - `extract_free_phages` - Extract free phage sequences from virus calls
   - `extract_unbinned_prophages` - Extract prophage sequences
   - `mask_prophage_regions` - Mask prophage coordinates with N's using bbtools
   - `checkv_unbinned` - Quality check viral sequences

3. **workflow/envs/bbtools_env.yaml** - Conda environment for bbtools/bbmask

### Modified Files

#### workflow/Snakefile
- Added includes for `prophage_mags.smk` and `prophage_unbinned.smk`

#### workflow/rules/phages.smk
- **Removed**: All prophage prediction rules (moved to new files)
- **Kept**: Database downloads (genomad_db, bakta_db, checkv_db), checkm rule
- **Added**: `merge_prophage_tables` - Combine MAG and unbinned tables
- **Modified**: `add_taxonomy_to_prophage_table` - Now uses merged table
- **Added**: `final_prophage_outputs` - Creates three separate outputs:
  - `prophages_from_mags.fasta` - Prophages from MAG bins
  - `prophages_from_unbinned.fasta` - Prophages from unbinned contigs
  - `free_phages.fasta` - Free phages from unbinned contigs
  - `all_prophages_combined.fasta` - All prophages (MAG + unbinned, no free phages)
- **Modified**: `run_everything` - Updated inputs for new outputs

#### workflow/rules/taxonomy.smk
- **Modified**: `gtdbtk_classify_bins` - Now depends on MAG prophage table (ensures prophage prediction completes first)
- **Modified**: `mmseqs_taxonomy_all_contigs` - Now uses masked unbinned contigs
- **Modified**: `mmseqs_taxonomy_prophage_only` - Now uses masked unbinned contigs and filters for unbinned prophages only (bin=="none")

#### config/config.yaml
- **Added**: `bbtools: "../envs/bbtools_env.yaml"` to conda_envs

## Workflow Logic

### PATH A - MAGs (Binned Contigs)
```
binning → per-MAG bakta → per-MAG phispy
                       ↘
                         merge prophages (all genomad + unique phispy)
                       ↗
          per-MAG genomad → prophage table & sequences
                         ↓
                    gtdb-tk taxonomy (on ORIGINAL unmasked bins)
```

**Why no masking for MAGs?**
- GTDB-Tk uses phylogenetic markers that are robust to prophage contamination
- MAGs have sufficient genomic context

### PATH B - Unbinned Contigs
```
unbinned contigs → genomad → prophages → mask regions with N's → mmseqs taxonomy
                           ↘
                             free phages → checkv
                           ↗
                      prophages → checkv
```

**Why mask unbinned contigs?**
- Unbinned contigs may be predominantly prophage sequence
- MMseqs similarity search could misclassify them without masking

## Dependency Flow

### Critical Dependencies
1. **gtdb-tk waits for MAG prophage prediction** to complete before running (ensures proper execution order)
2. **mmseqs waits for prophage masking** to complete (uses masked contigs as input)

### Data Sources in Final Table
Each prophage is tagged with a `source` column:
- `MAG` - Prophages from binned contigs (has bin assignment like "bin.1")
- `unbinned` - Prophages from unbinned contigs (bin == "none")

Additionally, free phages are tracked separately (not in prophage table).

## Testing Recommendations

1. **Dry-run test**:
   ```bash
   snakemake --dry-run
   ```

2. **Check for**:
   - Proper wildcard expansion for bin.{bin_num}
   - Checkpoint resolution for dynamic bin discovery
   - Dependency resolution (especially taxonomy after prophage prediction)
   - Proper file paths in all rules

3. **Verify outputs**:
   - `mags/prophage_table.tsv` - Should have columns: contig, start, end, tool, bin, source
   - `unbinned/prophage_table.tsv` - Same structure, bin should be "none", source should be "unbinned"
   - `final_prophage_table.tsv` - Merged version
   - `final_prophage_table_with_host_taxonomy.tsv` - With taxonomy added

## Potential Issues to Watch

1. **Empty bins**: Handles case where no bins exist (creates empty outputs)
2. **No prophages found**: Handles case where genomad finds nothing (creates empty BED/tables)
3. **Wildcard constraints**: May need to add bin_num constraints if issues arise
4. **BBtools masking format**: Verify bbmask.sh accepts the coordinate format we're using

## Next Steps

1. Test on actual data with `snakemake --dry-run`
2. Check logs for any path/dependency issues
3. May need to adjust bbmask.sh coordinate format based on testing
4. Verify taxonomy integration still works with new table structure
