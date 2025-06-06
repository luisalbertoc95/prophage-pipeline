# SLURM Job Naming in Snakemake v8

## The Issue
Snakemake v8 changed how SLURM job naming works, making it harder to get meaningful job names in `squeue`. This is a known limitation that affects many users.

## Current Status
With the default configuration, jobs will appear with generic names like:
```
smk-generic-0
smk-generic-1
```

## Workarounds

### Option 1: Add to your ~/.bashrc
```bash
alias snakemake-prophage='snakemake --executor slurm --default-resources slurm_partition=normal'
```

### Option 2: Use job groups
Add to your Snakemake command:
```bash
--groups spades=group1 genomad=group2
```

### Option 3: Use squeue with custom format
Create an alias to see more job details:
```bash
alias squeue-detail='squeue -u $USER -o "%.18i %.9P %.50j %.8u %.2t %.10M %.6D %R %Z"'
```

### Option 4: Check job details after submission
The SLURM output files will contain the rule and sample information:
```bash
ls slurm-*.out
```

## Future Solutions
- Snakemake developers are aware of this issue
- Future versions may restore better job naming functionality
- Consider using Snakemake v7 for production pipelines if job names are critical

## Finding Your Jobs
Even with generic names, you can identify jobs by:
1. Submission time
2. Resource allocation (memory, CPUs)
3. SLURM output files which contain rule information
4. The order of submission (first jobs are usually preprocessing)