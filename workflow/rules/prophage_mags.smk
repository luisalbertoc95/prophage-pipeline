# Prophage prediction workflow for MAG bins
# Runs bakta, phispy, and genomad on individual bins

import os

# Helper function to get all bin files for a sample
def get_bin_files(wildcards):
    bins_dir = os.path.join(config["outdir"], wildcards.sample, "binning", "dastool", f"{wildcards.sample}_DASTool_bins")
    if os.path.exists(bins_dir):
        bin_files = [f for f in os.listdir(bins_dir) if f.startswith("bin.") and f.endswith(".fa")]
        return [os.path.join(bins_dir, f) for f in bin_files]
    return []

# Helper function to extract bin number from filename
def get_bin_number(bin_file):
    """Extract bin number from bin.X.fa format"""
    import re
    match = re.search(r'bin\.(\d+)\.fa', os.path.basename(bin_file))
    return match.group(1) if match else None

rule bakta_per_mag:
    input:
        bin_file = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}_DASTool_bins", "bin.{bin_num}.fa"),
        db = config["bakta_database"]
    threads: 24
    conda: config["conda_envs"]["bakta"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bakta", "bin.{bin_num}"))
    log:
        os.path.join(config["outdir"], "logs", "bakta_per_mag", "{sample}_bin.{bin_num}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "bakta_per_mag", "{sample}_bin.{bin_num}_bmrk.txt")
    shell:
        """
        bakta --db {input.db}/db --force --skip-plot --keep-contig-headers --output {output} \
        --threads {threads} {input.bin_file} 2> {log}
        """

rule phispy_per_mag:
    input:
        bakta_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bakta", "bin.{bin_num}")
    threads: 8
    conda: config["conda_envs"]["phispy"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "phispy", "bin.{bin_num}"))
    log:
        os.path.join(config["outdir"], "logs", "phispy_per_mag", "{sample}_bin.{bin_num}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "phispy_per_mag", "{sample}_bin.{bin_num}_bmrk.txt")
    shell:
        """
        gbff_file=$(find {input.bakta_dir} -name "*.gbff" | head -1)
        if [ -f "$gbff_file" ]; then
            PhiSpy.py "$gbff_file" -o {output} --output_choice 63 2> {log} || true
        else
            echo "No .gbff file found in {input.bakta_dir}" > {log}
            mkdir -p {output}
            touch {output}/.empty
        fi
        """

# Checkpoint to dynamically determine which bins exist
checkpoint get_mag_bins:
    input:
        bins_done = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins")
    output:
        bin_list = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bin_list.txt")
    shell:
        """
        # Find all bin.*.fa files and extract just the bin identifier
        find {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins \
        -name "bin.*.fa" -type f -exec basename {{}} \\; | \
        sed 's/^bin\\.//;s/\\.fa$//' > {output.bin_list}
        """

def aggregate_mag_prophage_inputs(wildcards):
    """Aggregate all per-MAG prophage prediction outputs"""
    checkpoint_output = checkpoints.get_mag_bins.get(**wildcards).output.bin_list

    with open(checkpoint_output) as f:
        bin_nums = [line.strip() for line in f if line.strip()]

    return {
        "phispy": expand(os.path.join(config["outdir"], wildcards.sample, "phage_analysis", "mags", "phispy", "bin.{bin_num}"), bin_num=bin_nums),
        "bin_list": checkpoint_output
    }

rule map_genomad_to_bins:
    input:
        genomad_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad_complete"),
        bin_list = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bin_list.txt")
    conda: config["conda_envs"]["phage_all"]
    output:
        genomad_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "genomad_all.bed"),
        genomad_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "genomad_prophages.fasta")
    log:
        os.path.join(config["outdir"], "logs", "map_genomad_to_bins", "{sample}.log")
    shell:
        """
        # Find GeNomad provirus predictions
        provirus_tsv=$(find {input.genomad_dir} -name "*_provirus.tsv" | head -1)
        provirus_fna=$(find {input.genomad_dir} -name "*_provirus.fna" | head -1)

        # Initialize output files
        touch {output.genomad_bed}
        touch {output.genomad_fasta}

        if [ -f "$provirus_tsv" ] && [ -s "$provirus_tsv" ]; then
            # Build contig-to-bin mapping from bin FASTA files
            bin_dir={config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins

            # Create temporary mapping file
            temp_map=$(mktemp)

            # For each bin, extract contig IDs
            while IFS= read -r bin_num; do
                bin_file="$bin_dir/bin.$bin_num.fa"
                if [ -f "$bin_file" ]; then
                    grep "^>" "$bin_file" | sed 's/^>//' | while read contig_name; do
                        # Extract NODE number from contig name
                        if [[ "$contig_name" =~ NODE_([0-9]+)_ ]]; then
                            echo "${{BASH_REMATCH[1]}}\tbin.$bin_num"
                        fi
                    done >> "$temp_map"
                fi
            done < {input.bin_list}

            # Parse GeNomad provirus predictions and assign to bins
            awk 'NR>1 {{
                if (match($2, /NODE_([0-9]+)_/, arr)) {{
                    print arr[1] "\t" $3 "\t" $4 "\tgenomad\t" NR-1
                }}
            }}' "$provirus_tsv" | while IFS=$'\t' read -r contig start end tool pred_id; do
                # Look up bin for this contig
                bin=$(awk -v c="$contig" '$1 == c {{print $2; exit}}' "$temp_map")
                if [ -n "$bin" ]; then
                    echo -e "$contig\t$start\t$end\t$tool\t$pred_id\t$bin" >> {output.genomad_bed}
                fi
            done

            # Extract only binned prophage sequences from FASTA
            if [ -f "$provirus_fna" ] && [ -s {output.genomad_bed} ]; then
                # Get list of prediction IDs that are in bins
                temp_binned_ids=$(mktemp)
                awk '{{print $5}}' {output.genomad_bed} > "$temp_binned_ids"

                # Extract those sequences from FASTA
                # GeNomad names sequences like: final_filtered_contigs|provirus_1
                while read pred_id; do
                    grep -A1 "|provirus_$pred_id" "$provirus_fna" || true
                done < "$temp_binned_ids" >> {output.genomad_fasta}

                rm -f "$temp_binned_ids"
            fi

            rm -f "$temp_map"
        fi 2> {log}
        """

rule merge_mag_prophages:
    input:
        unpack(aggregate_mag_prophage_inputs),
        genomad_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "genomad_all.bed"),
        bin_list = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bin_list.txt")
    params:
        phispy_dirs = lambda wildcards: aggregate_mag_prophage_inputs(wildcards)["phispy"]
    conda: config["conda_envs"]["phage_all"]
    output:
        merged_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "merged_prophages.bed"),
        phispy_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "phispy_all.bed"),
        phispy_unique_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "phispy_unique.bed"),
        phispy_unique_ids = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "phispy_unique_ids.txt")
    log:
        os.path.join(config["outdir"], "logs", "merge_mag_prophages", "{sample}.log")
    shell:
        """
        # Initialize empty BED file
        touch {output.phispy_bed}

        # Extract PhiSpy predictions from all bins
        printf '%s\n' {params.phispy_dirs} | while IFS= read -r phispy_dir; do
            bin_num=$(basename "$phispy_dir" | sed 's/bin\\.//;s/\\..*$//')
            tsv_file="$phispy_dir/prophage.tsv"

            if [ -f "$tsv_file" ]; then
                awk -v bin="bin.$bin_num" 'NR>1 {{
                    if (match($2, /NODE_([0-9]+)_/, arr)) {{
                        split($1, pp_parts, "_")
                        print arr[1] "\t" $3 "\t" $4 "\tphispy\t" pp_parts[2] "\t" bin
                    }}
                }}' "$tsv_file" >> {output.phispy_bed}
            fi
        done 2>> {log}

        # Find PhiSpy predictions that DON'T overlap with GeNomad
        if [ -s {input.genomad_bed} ] && [ -s {output.phispy_bed} ]; then
            bedtools intersect -a {output.phispy_bed} -b {input.genomad_bed} -v > {output.phispy_unique_bed} 2>> {log}
        elif [ -s {output.phispy_bed} ]; then
            # No genomad predictions, all phispy are unique
            cp {output.phispy_bed} {output.phispy_unique_bed}
        else
            # No phispy predictions
            touch {output.phispy_unique_bed}
        fi

        # Create merged BED: all GeNomad + unique PhiSpy
        cat {input.genomad_bed} {output.phispy_unique_bed} > {output.merged_bed} 2>> {log}

        # Extract PhiSpy unique IDs for FASTA extraction
        if [ -s {output.phispy_unique_bed} ]; then
            awk '{{print $5}}' {output.phispy_unique_bed} > {output.phispy_unique_ids} 2>> {log}
        else
            touch {output.phispy_unique_ids}
        fi
        """

rule extract_mag_prophage_sequences:
    input:
        merged_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "merged_prophages.bed"),
        genomad_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "genomad_prophages.fasta"),
        phispy_unique_ids = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "phispy_unique_ids.txt"),
        bin_list = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bin_list.txt")
    params:
        phispy_dirs = lambda wildcards: aggregate_mag_prophage_inputs(wildcards)["phispy"]
    conda: config["conda_envs"]["phage_all"]
    output:
        prophage_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "all_prophages.fasta"),
        prophage_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "prophage_table.tsv")
    log:
        os.path.join(config["outdir"], "logs", "extract_mag_prophages", "{sample}.log")
    shell:
        """
        # Initialize output FASTA
        > {output.prophage_fasta}

        # Copy GeNomad prophage sequences (already filtered to binned contigs)
        if [ -s {input.genomad_fasta} ]; then
            cat {input.genomad_fasta} >> {output.prophage_fasta} 2> {log}
        fi

        # Collect unique PhiSpy prophage sequences
        if [ -s {input.phispy_unique_ids} ]; then
            printf '%s\n' {params.phispy_dirs} | while IFS= read -r phispy_dir; do
                fasta_file="$phispy_dir/phage.fasta"
                if [ -f "$fasta_file" ]; then
                    cat "$fasta_file" >> {output.prophage_fasta}
                fi
            done 2>> {log}
        fi

        # Create prophage table from merged BED with source column
        echo -e "contig\tstart\tend\ttool\tbin\tsource" > {output.prophage_table}
        awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $4, $6, "MAG"}}' {input.merged_bed} >> {output.prophage_table}
        """
