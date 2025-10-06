# GeNomad Comparison Branch

## Purpose
This branch implements a **per-MAG + per-unbinned** GeNomad approach to compare against the **complete assembly** approach on the main branch.

## Changes Made

### 1. Standardized GeNomad Parameters
All GeNomad runs now use **identical parameters** for fair comparison:
- `--cleanup`
- `--threads 24`
- No `--splits` parameter (let GeNomad decide automatically)

### 2. New Rules in `prophage_mags.smk`
- **`genomad_per_mag`**: Runs GeNomad on each individual MAG bin
- **`collect_genomad_per_mag`**: Collects prophage predictions from all per-MAG GeNomad runs

### 3. New Rules in `prophage_unbinned.smk`
- **`extract_unbinned_contigs`**: Creates FASTA file containing only unbinned contigs
- **`genomad_unbinned`**: Runs GeNomad only on unbinned contigs
- Updated `identify_unbinned_genomad` and `extract_free_phages` to use the new unbinned-only GeNomad output

### 4. Modified Rules in `phages.smk`
- **`genomad_complete_assembly_ORIGINAL`**: Renamed from `genomad_complete_assembly`, kept for comparison
- **`compare_genomad_approaches`**: New optional rule that compares both approaches

### 5. Comparison Rule
The `compare_genomad_approaches` rule:
- Collects prophage predictions from per-MAG approach (MAGs + unbinned separately)
- Collects prophage predictions from complete assembly approach
- Generates a comparison report with counts and differences
- Creates TSV files for detailed manual comparison

## How to Use

### Run the per-MAG approach (this branch):
```bash
cd workflow
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output
```

### Run the comparison (requires running ORIGINAL rule first):
```bash
# First, trigger the original complete assembly run
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output \
  {outdir}/{sample}/phage_analysis/genomad_complete_ORIGINAL

# Then run the comparison
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output \
  {outdir}/{sample}/phage_analysis/genomad_comparison_report.txt
```

### View comparison results:
```bash
cat {outdir}/{sample}/phage_analysis/genomad_comparison_report.txt
```

## Output Files for Comparison

### Per-MAG Approach:
- `{sample}/phage_analysis/mags/genomad/bin.{bin_num}/` - Per-MAG GeNomad outputs
- `{sample}/phage_analysis/unbinned/genomad/` - Unbinned GeNomad output
- `{sample}/phage_analysis/per_mag_all_prophages.tsv` - Combined prophage list

### Complete Assembly Approach:
- `{sample}/phage_analysis/genomad_complete_ORIGINAL/` - Complete assembly GeNomad output
- `{sample}/phage_analysis/original_all_prophages.tsv` - Prophage list

### Comparison Report:
- `{sample}/phage_analysis/genomad_comparison_report.txt` - Summary statistics

## Expected Observations

Based on your initial observations, we expect:
- **Complete assembly approach**: More prophage predictions
- **Per-MAG approach**: Fewer prophage predictions

The comparison will help determine:
1. How many prophages are "lost" in the per-MAG approach
2. Which specific prophages differ between approaches
3. Whether the difference is consistent across samples

## Next Steps

1. Run both approaches on the same sample
2. Review the comparison report
3. Investigate specific prophages that differ
4. Determine if the difference is due to:
   - Edge effects (prophages at contig boundaries)
   - Context loss (ML model needs surrounding sequence)
   - Binning artifacts
   - Other technical factors
