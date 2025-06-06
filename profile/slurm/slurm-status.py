#!/usr/bin/env python3
"""
Check SLURM job status for Snakemake
"""
import sys
import subprocess

jobid = sys.argv[1]

try:
    # Query job status
    result = subprocess.run(
        ["sacct", "-j", jobid, "-o", "State", "-n"],
        capture_output=True,
        text=True,
        check=True
    )
    
    # Get the status (first line, first word)
    status = result.stdout.strip().split()[0]
    
    # Map SLURM states to Snakemake states
    status_map = {
        "RUNNING": "running",
        "PENDING": "running",
        "COMPLETING": "running",
        "COMPLETED": "success",
        "FAILED": "failed",
        "TIMEOUT": "failed",
        "CANCELLED": "failed",
        "OUT_OF_MEMORY": "failed",
    }
    
    print(status_map.get(status, "running"))
    
except Exception:
    # If we can't get status, assume it's still running
    print("running")