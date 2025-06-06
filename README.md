
![updated_workflow_091024 drawio](https://github.com/user-attachments/assets/4f59c6d9-a453-4985-a11e-f8eed6714539)


# Prophage Detection Pipeline

A Snakemake-based bioinformatics pipeline for comprehensive prophage detection in metagenomic samples. The pipeline processes paired-end sequencing reads through seven main stages: preprocessing, assembly, binning, bin refinement, coverage analysis, taxonomy assignment, and multi-tool prophage detection.

## Requirements

- Snakemake version 8+
- [mamba](https://anaconda.org/conda-forge/mamba) 
- [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)

## How to Run

```bash
# Navigate to workflow directory
cd workflow

# Run the pipeline (without summaries)
snakemake --profile ../profile/slurm/ --config [options]

# Run the pipeline with automatic summary generation
snakemake all_with_summary --profile ../profile/slurm/ --config [options]

# Run only the summary generation (after pipeline completion)
snakemake all_summaries --profile ../profile/slurm/ --config [options]

# Dry run to see planned jobs
snakemake -n --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output

# Generate workflow visualization
snakemake --dag | dot -Tpng > workflow.png
```

## Configuration Options

### Required Parameters

- `reads`: Path to directory containing paired-end fastq reads (with suffixes _1.fastq.gz and _2.fastq.gz)
- `outdir`: Path to directory where all outputs will be created

### Optional Parameters

- `fastq_names_1`: Default is `{sample}_1.fastq.gz`
- `fastq_names_2`: Default is `{sample}_2.fastq.gz`
- `fastp_min_sequence_length`: Length threshold (in bp) for fastp step (default: 120)

### Database Paths (with defaults)

- `human_ref`: `/ref/sahlab/data/GRCh38.fna.gz`
- `genomad_database`: `/ref/sahlab/data/viral_analysis_DBs/genomad_DBs/genomad_db`
- `bakta_database`: `/ref/sahlab/data/bakta_db`
- `cat_database`: `/ref/sahlab/data/CAT_prepare_20210107`
- `checkv_database`: `/ref/sahlab/data/viral_analysis_DBs/checkV_DB/checkv-db-v1.4`

### Example Command

```bash
snakemake --profile ../profile/slurm/ --config reads=/scratch/sahlab/Megan/test_reads outdir=/scratch/sahlab/Megan/pipeline_test_out
```

## Important Notes

- Paired-end reads belonging to the same pair must have identical names in the r1 and r2 fastq files
- The pipeline uses centralized conda environment management for reproducibility
- Conda environments are stored persistently at: `/ref/sahlab/software/miniforge3/envs/smk_envs_prophage_pipeline`

## Output Directory Structure

The pipeline uses a process-first directory structure, organizing outputs by analysis type rather than by sample:

```
outputs/
├── preprocessing/          # Quality control and host decontamination
│   ├── sample1/
│   └── sample2/
├── assembly/              # SPAdes metagenomic assembly
│   ├── sample1/
│   └── sample2/
├── binning/               # Genomic binning (CONCOCT, MaxBin2, MetaBAT2, DAS Tool)
│   ├── sample1/
│   └── sample2/
├── coverm/                # Coverage statistics
│   ├── sample1/
│   └── sample2/
├── taxonomy/              # CAT/BAT taxonomic classification
│   ├── sample1/
│   └── sample2/
└── phage_analysis/        # Prophage detection results (main output)
    ├── sample1/
    └── sample2/
```

## Key Output Files

### In the phage_analysis/{sample}/ directory:

- `final_prophage_table.tsv`: Prophage predictions with genomic coordinates
- `final_prophage_table_with_host_taxonomy.tsv`: **Primary output** - prophages with host taxonomy information
- `final_prophage.fasta`: Sequences of all identified prophage regions

## Pipeline Stages

1. **Preprocessing**: Quality control (fastp) and host decontamination (KneadData)
2. **Assembly**: Metagenomic assembly using SPAdes
3. **Binning**: Multiple binning tools (CONCOCT, MaxBin2, MetaBAT2) refined with DAS Tool
4. **Bin Refinement**: Quality assessment with CheckM
5. **Coverage Analysis**: Coverage calculation using CoverM
6. **Taxonomy**: Taxonomic classification using CAT/BAT
7. **Prophage Analysis**: Multi-tool prophage detection using geNomad and PhiSpy, with results merged using custom R script

## Automatic Summary Generation

The pipeline can automatically generate comprehensive summaries and visualizations after all samples are processed:

### Option 1: Run Pipeline with Automatic Summaries
```bash
snakemake all_with_summary --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output
```

### Option 2: Generate Summaries After Pipeline Completion
```bash
snakemake all_summaries --profile ../profile/slurm/ --config outdir=/path/to/output
```

### Manual Summary Generation
You can also run the summary scripts independently:

```bash
# Generate Python visualizations
python workflow/scripts/create_summary_plots.py /path/to/output

# Generate R HTML report
Rscript workflow/scripts/create_summary_report.R /path/to/output
```

### Summary Outputs

The summary generation creates:
- `summary_plots/` - Python-generated publication-quality plots (PNG)
- `summary_plots_R/` - R-generated plots
- `prophage_summary_report.html` - Interactive HTML report
- `prophage_summary_by_sample.tsv` - Per-sample statistics
- `tool_detection_summary.tsv` - Tool performance metrics

