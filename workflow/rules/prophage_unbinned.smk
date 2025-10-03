# Prophage and free phage prediction workflow for unbinned contigs
# Runs genomad, extracts free phages, masks prophage regions

import os

rule identify_unbinned_genomad:
    input:
        genomad_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad_complete"),
        bin_list = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bin_list.txt")
    conda: config["conda_envs"]["phage_all"]
    output:
        prophage_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.bed"),
        prophage_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.fasta"),
        prophage_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophage_table.tsv")
    log:
        os.path.join(config["outdir"], "logs", "identify_unbinned_genomad", "{sample}.log")
    shell:
        """
        # Find GeNomad provirus predictions
        provirus_tsv=$(find {input.genomad_dir} -name "*_provirus.tsv" | head -1)
        provirus_fna=$(find {input.genomad_dir} -name "*_provirus.fna" | head -1)

        # Initialize output files
        touch {output.prophage_bed}
        touch {output.prophage_fasta}
        echo -e "contig\tstart\tend\ttool\tbin\tsource" > {output.prophage_table}

        if [ -f "$provirus_tsv" ] && [ -s "$provirus_tsv" ]; then
            # Build set of binned contigs from bin FASTA files
            bin_dir={config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins
            temp_binned_contigs=$(mktemp)

            # Extract all contig IDs from bins
            while IFS= read -r bin_num; do
                bin_file="$bin_dir/bin.$bin_num.fa"
                if [ -f "$bin_file" ]; then
                    grep "^>" "$bin_file" | sed 's/^>//' | while read contig_name; do
                        # Extract NODE number from contig name
                        if [[ "$contig_name" =~ NODE_([0-9]+)_ ]]; then
                            echo "${{BASH_REMATCH[1]}}"
                        fi
                    done >> "$temp_binned_contigs"
                fi
            done < {input.bin_list}

            # Sort and unique the binned contigs list
            sort -u "$temp_binned_contigs" > "$temp_binned_contigs.sorted"

            # Parse GeNomad provirus predictions and keep only unbinned ones
            temp_unbinned_ids=$(mktemp)
            awk 'NR>1 {{
                if (match($2, /NODE_([0-9]+)_/, arr)) {{
                    print arr[1] "\t" $3 "\t" $4 "\tgenomad\t" NR-1 "\tnone"
                }}
            }}' "$provirus_tsv" | while IFS=$'\t' read -r contig start end tool pred_id bin; do
                # Check if contig is NOT in binned list
                if ! grep -q "^$contig$" "$temp_binned_contigs.sorted"; then
                    echo -e "$contig\t$start\t$end\t$tool\t$pred_id\t$bin" >> {output.prophage_bed}
                    echo "$pred_id" >> "$temp_unbinned_ids"
                fi
            done

            # Extract unbinned prophage sequences from FASTA
            if [ -f "$provirus_fna" ] && [ -s "$temp_unbinned_ids" ]; then
                # Convert prediction IDs to sequence names in FASTA
                # GeNomad names sequences like: final_filtered_contigs|provirus_1
                while read pred_id; do
                    grep -A1 "|provirus_$pred_id" "$provirus_fna" || true
                done < "$temp_unbinned_ids" >> {output.prophage_fasta}
            fi

            # Create prophage table with source column
            awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $4, $6, "unbinned"}}' {output.prophage_bed} >> {output.prophage_table}

            rm -f "$temp_binned_contigs" "$temp_binned_contigs.sorted" "$temp_unbinned_ids"
        fi 2> {log}
        """

rule extract_free_phages:
    input:
        genomad_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad_complete")
    conda: config["conda_envs"]["phage_all"]
    output:
        free_phage_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "free_phages.fasta"),
        free_phage_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "free_phage_table.tsv"),
        virus_summary = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "virus_summary.tsv")
    log:
        os.path.join(config["outdir"], "logs", "extract_free_phages", "{sample}.log")
    shell:
        """
        # Find genomad virus predictions (not proviruses)
        virus_fna=$(find {input.genomad_dir} -name "*_virus.fna" | head -1)
        virus_summary=$(find {input.genomad_dir} -name "*_virus_summary.tsv" | head -1)

        if [ -f "$virus_fna" ] && [ -s "$virus_fna" ]; then
            # Copy free phage sequences
            cp "$virus_fna" {output.free_phage_fasta} 2> {log}

            # Create table with source column
            echo -e "contig\tlength\ttopology\tcoordinates\tn_genes\tsource" > {output.free_phage_table}
            if [ -f "$virus_summary" ]; then
                # Extract relevant columns and add source
                awk 'NR>1 {{
                    if (match($1, /NODE_([0-9]+)_/, arr)) {{
                        print arr[1] "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\tfree_phage"
                    }}
                }}' "$virus_summary" >> {output.free_phage_table} 2>> {log}
                cp "$virus_summary" {output.virus_summary}
            else
                touch {output.virus_summary}
            fi
        else
            # No free phages found
            touch {output.free_phage_fasta}
            echo -e "contig\tlength\ttopology\tcoordinates\tn_genes\tsource" > {output.free_phage_table}
            touch {output.virus_summary}
            echo "No free phages found" >> {log}
        fi
        """

rule mask_prophage_regions:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "filt_4000_seqs_to_keep.fasta"),
        prophage_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.bed")
    conda: config["conda_envs"]["phage_all"]
    output:
        masked_contigs = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "masked_contigs.fasta"),
        mask_stats = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "masking_stats.txt")
    log:
        os.path.join(config["outdir"], "logs", "mask_prophages", "{sample}.log")
    shell:
        """
        # Check if there are any prophages to mask
        if [ -s {input.prophage_bed} ]; then
            # Count regions to mask
            n_regions=$(wc -l < {input.prophage_bed})
            echo "Masking $n_regions prophage regions" > {output.mask_stats}

            # Convert BED to proper format for bedtools
            # BED columns: contig start end, need to add NODE_ prefix
            awk '{{print "NODE_" $1 "\t" $2 "\t" $3}}' {input.prophage_bed} > {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/mask_regions.bed

            # Use bedtools maskfasta to mask regions with N's
            bedtools maskfasta -fi {input.contigs} -bed {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/mask_regions.bed -fo {output.masked_contigs} 2>> {log}

            # Report masking statistics
            echo "Original contigs: $(grep -c '^>' {input.contigs})" >> {output.mask_stats}
            echo "Masked contigs: $(grep -c '^>' {output.masked_contigs})" >> {output.mask_stats}
        else
            # No prophages to mask, just copy original contigs
            echo "No prophages found - no masking needed" > {output.mask_stats}
            cp {input.contigs} {output.masked_contigs}
        fi
        """

rule checkv_unbinned:
    input:
        prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.fasta"),
        free_phages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "free_phages.fasta"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_unbinned", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_unbinned", "{sample}_bmrk.txt")
    shell:
        """
        # Combine prophages and free phages for checkv
        cat {input.prophages} {input.free_phages} > {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/all_phages_for_checkv.fasta

        # Run checkv
        checkv end_to_end \
        {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/all_phages_for_checkv.fasta \
        {output} \
        -t {threads} \
        -d {input.db} 2> {log}
        """
