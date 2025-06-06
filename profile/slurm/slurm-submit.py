#!/usr/bin/env python3
"""
Custom SLURM submission script for Snakemake with meaningful job names
"""
import sys
import re
import subprocess
from snakemake.utils import read_job_properties

# Read job properties from Snakemake
jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# Extract rule name and wildcards
rule = job_properties.get("rule", "unknown")
wildcards = job_properties.get("wildcards", {})

# Create job name
if wildcards:
    # Get the sample wildcard if it exists
    sample = wildcards.get("sample", "")
    if sample:
        job_name = f"{rule}_{sample}"
    else:
        # Use first wildcard value if no sample
        first_wildcard = list(wildcards.values())[0] if wildcards else ""
        job_name = f"{rule}_{first_wildcard}" if first_wildcard else rule
else:
    job_name = rule

# Truncate job name to SLURM's limit (typically 256 chars, but we'll use 50 for readability)
job_name = job_name[:50]

# Build sbatch command
sbatch_cmd = ["sbatch"]

# Add job name
sbatch_cmd.extend(["--job-name", job_name])

# Add resources from job properties
resources = job_properties.get("resources", {})
if "mem_mb" in resources:
    sbatch_cmd.extend(["--mem", f"{resources['mem_mb']}M"])
if "runtime" in resources:
    sbatch_cmd.extend(["--time", f"{resources['runtime']}"])
if "cores" in resources:
    sbatch_cmd.extend(["--cpus-per-task", str(resources["cores"])])
if "slurm_account" in resources:
    sbatch_cmd.extend(["--account", resources["slurm_account"]])

# Add any cluster parameters
cluster_params = job_properties.get("cluster", {})
for key, value in cluster_params.items():
    sbatch_cmd.extend([f"--{key}", str(value)])

# Add the job script
sbatch_cmd.append(jobscript)

# Submit the job
subprocess.run(sbatch_cmd)