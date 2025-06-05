# Prophage Pipeline Changes Summary

This document summarizes the major improvements implemented in the prophage detection pipeline on the `restructure-directories` branch.

## Overview

The pipeline underwent significant organizational improvements to enhance usability, maintainability, and HPC performance. These changes make cross-sample comparisons easier and provide better control over computational environments.

## 1. Process-First Directory Structure

### Previous Structure (Sample-First)
```
outputs/
├── sample1/
│   ├── preprocessing/
│   ├── assembly/
│   ├── binning/
│   ├── taxonomy/
│   └── phage_analysis/
└── sample2/
    ├── preprocessing/
    ├── assembly/
    ├── binning/
    ├── taxonomy/
    └── phage_analysis/
```

### New Structure (Process-First)
```
outputs/
├── preprocessing/
│   ├── sample1/
│   └── sample2/
├── assembly/
│   ├── sample1/
│   └── sample2/
├── binning/
│   ├── sample1/
│   └── sample2/
├── taxonomy/
│   ├── sample1/
│   └── sample2/
└── phage_analysis/
    ├── sample1/
    └── sample2/
```

### Benefits
- **Easier comparative analysis**: All results for a specific analysis step are grouped together
- **Better for downstream tools**: MultiQC and custom analysis scripts can more easily aggregate results
- **Cleaner workflow organization**: Researchers can focus on specific analysis types across all samples
- **Improved collaboration**: Different team members can work on different analysis steps independently

### Implementation
- Updated all output paths in 7 rule files (`preprocessing.smk`, `assembly.smk`, `binning.smk`, `refine_bins.smk`, `coverm.smk`, `taxonomy.smk`, `phages.smk`)
- Modified main `Snakefile` target rule
- Updated ~50+ file path references throughout the pipeline
- Validated with successful dry-run testing

## 2. Centralized Conda Environment Management

### Previous Approach
Each rule had hardcoded conda environment paths:
```yaml
rule genomad:
    conda: "../envs/genomad_env.yaml"
```

### New Approach
Centralized configuration in `config.yaml`:
```yaml
# Conda environments
conda_envs:
  bakta: "../envs/bakta_env.yaml"
  cat: "../envs/cat_env.yaml"
  checkm: "../envs/checkm_env.yaml"
  checkv: "../envs/checkv_env.yaml"
  concoct: "../envs/concoct_env.yaml"
  coverm: "../envs/coverm_env.yaml"
  dastool: "../envs/dastool_env.yaml"
  fastp: "../envs/fastp_test.yaml"
  genomad: "../envs/genomad_env.yaml"
  kneaddata: "../envs/kneaddata.yaml"
  maxbin: "../envs/maxbin_env.yaml"
  metabat: "../envs/metabat_env.yaml"
  minimap: "../envs/minimap_env.yaml"
  phage_all: "../envs/phage_all_env.yaml"
  phispy: "../envs/phispy_env.yaml"
  spades: "../envs/spades_env.yml"
```

Rules now reference the config:
```yaml
rule genomad:
    conda: config["conda_envs"]["genomad"]
```

### Benefits
- **Centralized management**: Single location to modify all environment paths
- **Flexibility**: Easy to switch between YAML files and pre-built environments
- **HPC optimization**: Can easily point to shared environments (e.g., `/ref/sahlab/software/miniforge3/envs/env_name`)
- **Maintenance**: Simplified debugging and updates
- **Environment consistency**: Ensures all rules use the intended environment versions

### Implementation
- Added `conda_envs` section to `config/config.yaml`
- Updated 16 conda environment references across 7 rule files
- All rule files now use `config["conda_envs"]["env_name"]` syntax

## 3. Persistent Conda Environment Configuration

### Configuration Added
Added to SLURM profile (`profile/slurm/config.v8+.yaml`):
```yaml
conda-prefix: /ref/sahlab/software/miniforge3/envs/smk_envs_prophage_pipeline
```

### Benefits
- **Persistent environments**: Environments are stored in a permanent location instead of temporary scratch space
- **Faster execution**: Environments don't need to be recreated for each run
- **Shared environments**: Multiple pipeline runs can share the same environment installations
- **Resource efficiency**: Reduces conda environment setup time and storage usage
- **HPC optimization**: Environments persist beyond job completion

## 4. Development Documentation

### CLAUDE.md Added
Created comprehensive documentation for future development:
- **Commands section**: Key commands for running, testing, and developing the pipeline
- **Architecture section**: High-level overview of the seven-stage workflow
- **Critical files**: Important components and their purposes
- **Output structure**: Explanation of the new process-first organization

### Benefits
- **Faster onboarding**: New developers can quickly understand the pipeline structure
- **Consistent development**: Standardized commands and practices
- **Maintenance guidance**: Clear documentation of key components and design patterns

## Files Modified

### Configuration Files
- `config/config.yaml` - Added centralized conda environment configuration
- `profile/slurm/config.v8+.yaml` - Added persistent conda environment location

### Workflow Files
- `workflow/Snakefile` - Updated target rule for process-first structure
- `workflow/rules/preprocessing.smk` - Updated paths and conda references
- `workflow/rules/assembly.smk` - Updated paths and conda references
- `workflow/rules/binning.smk` - Updated paths and conda references
- `workflow/rules/refine_bins.smk` - Updated paths and conda references
- `workflow/rules/coverm.smk` - Updated paths and conda references
- `workflow/rules/taxonomy.smk` - Updated paths and conda references
- `workflow/rules/phages.smk` - Updated paths and conda references

### Documentation
- `CLAUDE.md` - Added development guidance (new file)
- `PIPELINE_CHANGES.md` - This summary document (new file)

## Validation

### Testing Performed
- **Syntax validation**: All Snakemake files parse correctly
- **Dry-run testing**: Pipeline generates expected workflow with new directory structure
- **Configuration testing**: Centralized conda environments work as expected
- **Path verification**: All output paths follow the new process-first structure

### Results
✅ Pipeline successfully restructured
✅ All conda environments properly configured  
✅ Process-first directories implemented
✅ SLURM profile configured for persistent environments
✅ Documentation created for future development

## Usage

### Running the Pipeline
The pipeline usage remains the same:
```bash
cd workflow
snakemake --profile ../profile/slurm/ --config reads=/path/to/reads outdir=/path/to/output
```

### Expected Output Structure
Results will now be organized by analysis type:
- `outdir/preprocessing/` - Quality control and host removal results
- `outdir/assembly/` - SPAdes assembly outputs
- `outdir/binning/` - Genomic bins and quality metrics
- `outdir/taxonomy/` - CAT taxonomic classifications
- `outdir/phage_analysis/` - **Primary prophage detection results**

### Key Output Files
- `outdir/phage_analysis/{sample}/final_prophage_table_with_host_taxonomy.tsv` - Main results table
- `outdir/phage_analysis/{sample}/final_prophage.fasta` - All identified prophage sequences

## Future Considerations

### Potential Enhancements
1. **Environment flexibility**: Easy switching between YAML files and pre-built environments for different HPC systems
2. **Cross-sample analysis**: The new structure facilitates development of cross-sample comparison rules
3. **Modular execution**: Process-first structure enables running specific analysis steps across all samples
4. **Resource optimization**: Persistent environments reduce computational overhead

### Migration Notes
- Existing users should update their output directory expectations
- Analysis scripts may need updates to reflect new directory structure
- The changes are backward-compatible with existing configuration files

## Conclusion

These improvements significantly enhance the pipeline's usability and maintainability while optimizing it for HPC environments. The process-first directory structure makes comparative analysis much more intuitive, while centralized conda management provides better control and flexibility for different computational environments.

---

*Changes implemented on the `restructure-directories` branch*  
*Generated with Claude Code assistance*