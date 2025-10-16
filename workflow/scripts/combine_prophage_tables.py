#!/usr/bin/env python3
"""
Combine prophage tables from multiple samples into a single table.
Adds a 'sample' column to identify the source sample for each prophage.
"""

import pandas as pd
import sys
from pathlib import Path


def combine_prophage_tables(sample_files, output_file):
    """
    Combine prophage tables from multiple samples.

    Args:
        sample_files: List of tuples (sample_name, file_path)
        output_file: Path to output combined table
    """
    combined_data = []

    for sample_name, file_path in sample_files:
        # Check if file exists
        if not Path(file_path).exists():
            print(f"Warning: File not found for sample {sample_name}: {file_path}", file=sys.stderr)
            continue

        # Read the prophage table
        try:
            df = pd.read_csv(file_path, sep='\t')

            # Check if table is empty
            if df.empty:
                print(f"Warning: Empty table for sample {sample_name}", file=sys.stderr)
                continue

            # Add sample column as the first column
            df.insert(0, 'sample', sample_name)

            combined_data.append(df)
            print(f"Added {len(df)} prophages from sample {sample_name}", file=sys.stderr)

        except Exception as e:
            print(f"Error reading file for sample {sample_name}: {e}", file=sys.stderr)
            continue

    # Combine all dataframes
    if combined_data:
        combined_df = pd.concat(combined_data, ignore_index=True)

        # Write to output file
        combined_df.to_csv(output_file, sep='\t', index=False)
        print(f"Combined table written to {output_file}", file=sys.stderr)
        print(f"Total prophages: {len(combined_df)} from {len(combined_data)} samples", file=sys.stderr)
    else:
        print("Warning: No data to combine", file=sys.stderr)
        # Create empty file with header
        empty_df = pd.DataFrame(columns=['sample', 'contig', 'start', 'end', 'tool', 'bin', 'source',
                                         'superkingdom', 'phylum', 'class', 'order', 'family',
                                         'genus', 'species', 'taxonomy_source'])
        empty_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    # Parse command line arguments
    # Format: python combine_prophage_tables.py output_file sample1:file1 sample2:file2 ...

    if len(sys.argv) < 3:
        print("Usage: combine_prophage_tables.py <output_file> <sample1:file1> <sample2:file2> ...", file=sys.stderr)
        sys.exit(1)

    output_file = sys.argv[1]

    # Parse sample:file pairs
    sample_files = []
    for arg in sys.argv[2:]:
        if ':' not in arg:
            print(f"Warning: Invalid format for argument '{arg}', expected 'sample:file'", file=sys.stderr)
            continue
        sample_name, file_path = arg.split(':', 1)
        sample_files.append((sample_name, file_path))

    if not sample_files:
        print("Error: No valid sample files provided", file=sys.stderr)
        sys.exit(1)

    # Combine tables
    combine_prophage_tables(sample_files, output_file)
