
# Prophage Pipeline

<img width="930" height="580" alt="updated_workflow_072525 drawio" src="https://github.com/user-attachments/assets/4da9d036-4329-403c-a740-8934d393b127" />

A Snakemake pipeline for comprehensive prophage detection and host taxonomy assignment from metagenomic data. The pipeline combines multiple prophage detection tools, performs binning-based host assignment, and integrates hybrid taxonomy (GTDB-Tk + MMseqs) for accurate host identification.

## Quick Start

### Prerequisites
- Snakemake 8+
- [mamba](https://anaconda.org/conda-forge/mamba)
- [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)

### Running the Pipeline
```bash
cd workflow
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output
```

**Note:** Paired-end reads must have identical names in R1 and R2 files.

## Configuration Options

### Required Parameters
- `reads`: Directory containing paired-end FASTQ files (suffixes: `_1.fastq.gz`, `_2.fastq.gz`)
- `outdir`: Output directory for all results

### Optional Parameters
- `fastq_names_1`: R1 file pattern (default: `{sample}_1.fastq.gz`)
- `fastq_names_2`: R2 file pattern (default: `{sample}_2.fastq.gz`)
- `fastp_min_sequence_length`: Minimum read length after trimming (default: 120)
- `taxonomy_scope`: Taxonomy analysis scope - `"prophage_only"` or `"all_contigs"` (default: `"prophage_only"`)

### Database Paths
- `genomad_database`: GeNomad viral database
- `bakta_database`: Bakta annotation database  
- `gtdbtk_database`: GTDB-Tk reference database
- `mmseqs_database`: MMseqs NR database
- `checkv_database`: CheckV quality database
- `human_ref`: Human reference genome for decontamination

## Key Features

### Prophage Detection
- **GeNomad**: ML-based prophage detection
- **PhiSpy**: HMM-based prophage detection  
- **Overlap Resolution**: Combines predictions, prioritizing GeNomad with unique PhiSpy additions

### Host Taxonomy Assignment
- **Hybrid Approach**: GTDB-Tk for high-quality binned contigs, MMseqs for unbinned contigs
- **Binning Integration**: Uses DAS Tool consensus bins from CONCOCT, MaxBin, and MetaBAT
- **Configurable Scope**: Taxonomy on prophage contigs only or all contigs

### Quality Control
- **CheckV**: Prophage completeness and contamination assessment
- **CheckM**: Bin quality evaluation
- **Coverage Analysis**: CoverM mapping statistics

## Outputs

Each sample generates the following directory structure:
```
{sample}/
├── assembly/           # Megahit assembly results
├── binning/           # DAS Tool consensus bins and CheckM results  
├── coverm/            # Read mapping and coverage statistics
├── taxonomy/          # GTDB-Tk and MMseqs taxonomy results
└── phage_analysis/    # Prophage detection and analysis
```

### Key Output Files

#### Prophage Tables
- `final_prophage_table.tsv`: Prophage coordinates with bin assignments
- `final_prophage_table_with_host_taxonomy.tsv`: Prophages with host taxonomy (hybrid GTDB-Tk/MMseqs)

#### Sequences  
- `final_prophage.fasta`: Extracted prophage sequences
- `contigs_with_prophages.fasta`: Full contigs containing prophages

#### Quality Assessment
- `checkv/`: Prophage quality and completeness analysis
- `checkm/`: Bin quality statistics

## Example Commands

### Basic Run
```bash
snakemake --profile ../profile/slurm/ --config \
  reads=/scratch/data/fastq_files \
  outdir=/scratch/results/prophage_analysis
```

### All-Contigs Taxonomy Mode  
```bash
snakemake --profile ../profile/slurm/ --config \
  reads=/scratch/data/fastq_files \
  outdir=/scratch/results/prophage_analysis \
  taxonomy_scope=all_contigs
```

