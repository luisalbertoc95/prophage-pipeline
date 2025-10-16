#!/usr/bin/env python3
"""
Combine prophage tables from multiple samples into a single table.
Adds a 'sample' column to identify the source sample for each prophage.
"""

import sys
import os


def combine_prophage_tables(sample_files, output_file):
    """
    Combine prophage tables from multiple samples.

    Args:
        sample_files: List of tuples (sample_name, file_path)
        output_file: Path to output combined table
    """
    header_written = False
    total_prophages = 0
    samples_processed = 0

    with open(output_file, 'w') as outfile:
        for sample_name, file_path in sample_files:
            # Check if file exists
            if not os.path.exists(file_path):
                print(f"Warning: File not found for sample {sample_name}: {file_path}", file=sys.stderr)
                continue

            try:
                with open(file_path, 'r') as infile:
                    lines = infile.readlines()

                    if not lines:
                        print(f"Warning: Empty file for sample {sample_name}", file=sys.stderr)
                        continue

                    # Process header
                    if not header_written:
                        # Add 'sample' column to header
                        header = lines[0].strip()
                        outfile.write(f"sample\t{header}\n")
                        header_written = True

                    # Process data lines
                    sample_count = 0
                    for line in lines[1:]:
                        line = line.strip()
                        if line:  # Skip empty lines
                            outfile.write(f"{sample_name}\t{line}\n")
                            sample_count += 1
                            total_prophages += 1

                    if sample_count > 0:
                        samples_processed += 1
                        print(f"Added {sample_count} prophages from sample {sample_name}", file=sys.stderr)
                    else:
                        print(f"Warning: No data rows in file for sample {sample_name}", file=sys.stderr)

            except Exception as e:
                print(f"Error reading file for sample {sample_name}: {e}", file=sys.stderr)
                continue

        # If no data was processed, write an empty file with header
        if not header_written:
            print("Warning: No data to combine, creating empty file with header", file=sys.stderr)
            outfile.write("sample\tcontig\tstart\tend\ttool\tbin\tsource\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttaxonomy_source\n")

    if samples_processed > 0:
        print(f"Combined table written to {output_file}", file=sys.stderr)
        print(f"Total prophages: {total_prophages} from {samples_processed} samples", file=sys.stderr)


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
