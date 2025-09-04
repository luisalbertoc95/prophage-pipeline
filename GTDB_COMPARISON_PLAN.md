# GTDB Taxonomy Comparison Plan

## Context
Currently using MMseqs2 with NR database for taxonomy assignment. Coworker recommended switching to GTDB-Tk or MMseqs2 with GTDB database for better performance and more up-to-date taxonomic structure.

## Current Implementation
- **File**: `workflow/rules/taxonomy.smk` (lines 1-34)
- **Database**: `/ref/sahlab/data/nr/nr_mmseqs15_DB` (NR database)
- **Parsing**: `workflow/scripts/merge_prophages.R` (lines 98-127)
- **Environment**: `workflow/envs/mmseqs_env.yaml`

## Plan: Benchmark GTDB-Tk vs MMseqs2+GTDB

### 1. Branch Setup
- Create new branch `gtdb-comparison` from current `development`
- This allows testing both approaches without affecting main pipeline

### 2. Implement Both Approaches

#### Option A: GTDB-Tk
- Create `workflow/envs/gtdbtk_env.yaml`
- Add new rule `gtdbtk_taxonomy` in `workflow/rules/taxonomy.smk`
- Download/configure GTDB-Tk database
- Update config for GTDB-Tk paths

#### Option B: MMseqs2 + GTDB
- Keep existing MMseqs2 setup
- Download GTDB sequences and convert to MMseqs2 format
- Create alternative rule `mmseqs_gtdb_taxonomy`
- Update config for GTDB database path

### 3. Make Both Approaches Configurable
- Add config option to switch between:
  - `mmseqs_nr` (current)
  - `gtdbtk` (new)
  - `mmseqs_gtdb` (new)
- Conditional rules based on config selection

### 4. Benchmark Framework
- Use Snakemake's built-in benchmarking (already configured)
- Compare:
  - Runtime (wall clock time)
  - Memory usage
  - CPU utilization
  - Result quality/consistency

### 5. Test & Compare
- Run same sample through both approaches
- Analyze benchmark files and taxonomy outputs
- Document differences in taxonomic assignments

## GTDB-Tk Benefits (from research)
- Uses phylogenetically consistent GTDB database
- 50%+ performance improvement with ANI screening
- Specifically designed for bacterial/archaeal genome classification
- More up-to-date taxonomic structure than NR
- Works well with MAGs and contigs

## Implementation Notes
- GTDB-Tk output format will be different from MMseqs2
- Will need to update parsing in `merge_prophages.R`
- Both approaches should produce comparable taxonomic lineages
- Snakemake benchmarks already configured - will automatically capture performance metrics

## Next Steps When Resuming
1. `git checkout -b gtdb-comparison`
2. Start with creating GTDB-Tk environment file
3. Add configurable taxonomy method to config.yaml
4. Implement both taxonomy rules
5. Test with small dataset first

## Files to Modify
- `config/config.yaml` - add taxonomy method selection
- `workflow/rules/taxonomy.smk` - add new rules
- `workflow/envs/gtdbtk_env.yaml` - new file
- `workflow/scripts/merge_prophages.R` - update parsing logic