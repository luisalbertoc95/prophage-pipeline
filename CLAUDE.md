# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Running the Pipeline
```bash
# Basic pipeline execution from the workflow directory
cd workflow
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output

# Dry run to see planned jobs
snakemake -n --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output

# Run with custom parameters
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output fastp_min_sequence_length=100

# Generate workflow visualization
snakemake --dag | dot -Tpng > workflow.png
```

### Development Commands
```bash
# Validate Snakefile syntax
snakemake --lint

# Run specific rule for testing
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output -R rule_name

# Check which files will be created
snakemake --summary --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output
```

## Architecture

This is a Snakemake-based bioinformatics pipeline for prophage detection in metagenomic samples. The workflow processes paired-end sequencing reads through seven main stages:

1. **Preprocessing**: Quality control and host decontamination
2. **Assembly**: Metagenomic assembly into contigs
3. **Binning**: Grouping contigs into genomic bins
4. **Bin Refinement**: Optimizing bins and quality assessment
5. **Coverage Analysis**: Calculating coverage statistics across samples
6. **Taxonomy**: Taxonomic classification of contigs
7. **Prophage Analysis**: Multi-tool prophage detection and merging

### Key Design Patterns

- **Modular Rules**: Each processing stage is in a separate `.smk` file under `workflow/rules/`
- **Conda Integration**: Each tool has its own environment specification in `workflow/envs/`
- **Resource Management**: SLURM profile defines memory/CPU requirements per rule
- **Error Tolerance**: Some rules use `|| true` to continue on non-critical failures
- **Result Integration**: Custom R script merges prophage predictions from multiple tools

### Critical Files

- `workflow/Snakefile`: Main workflow orchestrator that includes all rule modules
- `config/config.yaml`: Default parameters and database paths
- `workflow/scripts/merge_prophages.R`: Custom logic for combining geNomad and PhiSpy results
- `profile/slurm/config.v8+.yaml`: HPC resource allocations per rule

### Output Structure

Each sample creates five subdirectories in the output directory:
- `assembly/`: SPAdes assembly results
- `binning/`: Bins from CONCOCT, MaxBin2, MetaBAT2, and DAS Tool
- `coverm/`: Coverage statistics
- `taxonomy/`: CAT/BAT taxonomic classifications
- `phage_analysis/`: Prophage predictions (main results)

The primary output is `phage_analysis/final_prophage_table_with_host_taxonomy.tsv` containing all identified prophages with their genomic coordinates and host taxonomy.