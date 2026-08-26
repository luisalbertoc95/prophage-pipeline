# Prophage and free phage prediction workflow for unbinned contigs
# Runs genomad, extracts free phages, masks prophage regions

import os

rule extract_unbinned_contigs:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        bin_list = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "bin_list.txt")
    conda: config["conda_envs"]["phage_all"]
    output:
        unbinned_contigs = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "unbinned_contigs.fasta")
    log:
        os.path.join(config["outdir"], "logs", "extract_unbinned_contigs", "{sample}.log")
    shell:
        """
        # Build set of binned contigs from bin FASTA files
        bin_dir={config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins
        temp_binned_contigs=$(mktemp)

        # Extract all contig IDs from bins
        while IFS= read -r bin_num; do
            bin_file="$bin_dir/bin.$bin_num.fa"
            if [ -f "$bin_file" ]; then
                grep "^>" "$bin_file" | sed 's/^>//' >> "$temp_binned_contigs"
            fi
        done < {input.bin_list}

        # Sort and unique the binned contigs list
        sort -u "$temp_binned_contigs" > "$temp_binned_contigs.sorted"

        # Extract unbinned contigs using seqkit
        # Create exclude pattern file (one pattern per line)
        seqkit grep -v -f "$temp_binned_contigs.sorted" {input.contigs} > {output.unbinned_contigs} 2> {log}

        rm -f "$temp_binned_contigs" "$temp_binned_contigs.sorted"
        """

rule genomad_unbinned:
    input:
        unbinned_contigs = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "unbinned_contigs.fasta"),
        db = config["genomad_database"]
    threads: 24
    conda: config["conda_envs"]["genomad"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad"))
    log:
        os.path.join(config["outdir"], "logs", "genomad_unbinned", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "genomad_unbinned", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output}
        # ROBUSTNESS (fork): genomad aborts on an empty/degenerate FASTA. When a sample has
        # no unbinned contigs to scan, leave the output dir empty -> identify_unbinned_genomad
        # (find|| true + touch) yields 0 unbinned prophages instead of crashing.
        if grep -q '^>' {input.unbinned_contigs} 2>/dev/null; then
            genomad end-to-end --cleanup --threads {threads} \
            {input.unbinned_contigs} {output} {input.db} 2> {log}
        else
            echo "No unbinned contigs; skipping genomad (0 unbinned prophages)" > {log}
        fi
        """

rule identify_unbinned_genomad:
    input:
        genomad_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad")
    conda: config["conda_envs"]["phage_all"]
    output:
        genomad_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad_prophages.bed"),
        genomad_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad_prophages.fasta")
    log:
        os.path.join(config["outdir"], "logs", "identify_unbinned_genomad", "{sample}.log")
    shell:
        """
        # Find GeNomad provirus predictions
        provirus_tsv=$(find {input.genomad_dir} -name "*_provirus.tsv" | head -1 || true)
        provirus_fna=$(find {input.genomad_dir} -name "*_provirus.fna" | head -1 || true)

        # Initialize output files
        touch {output.genomad_bed}
        touch {output.genomad_fasta}

        if [ -f "$provirus_tsv" ] && [ -s "$provirus_tsv" ]; then
            # Parse GeNomad provirus predictions
            awk 'NR>1 {{
                if (match($2, /NODE_([0-9]+)_/, arr)) {{
                    print arr[1] "\t" $3 "\t" $4 "\tgenomad\t" NR-1 "\tnone"
                }}
            }}' "$provirus_tsv" > {output.genomad_bed} 2> {log}

            # Copy prophage sequences
            if [ -f "$provirus_fna" ] && [ -s "$provirus_fna" ]; then
                cp "$provirus_fna" {output.genomad_fasta} 2>> {log}
            fi
        fi 2>> {log}
        """

rule merge_unbinned_prophages:
    input:
        genomad_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad_prophages.bed"),
        genomad_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad_prophages.fasta"),
        checkv_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_prophages.bed"),
        checkv_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_prophages.fasta")
    conda: config["conda_envs"]["phage_all"]
    output:
        prophage_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.bed"),
        prophage_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.fasta"),
        prophage_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophage_table.tsv"),
        checkv_unique_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_unique.bed")
    log:
        os.path.join(config["outdir"], "logs", "merge_unbinned_prophages", "{sample}.log")
    shell:
        """
        # Find CheckV predictions that DON'T overlap with GeNomad
        if [ -s {input.genomad_bed} ] && [ -s {input.checkv_bed} ]; then
            # Use bedtools to find non-overlapping CheckV predictions
            bedtools intersect -a {input.checkv_bed} -b {input.genomad_bed} -v > {output.checkv_unique_bed} 2>> {log}
        elif [ -s {input.checkv_bed} ]; then
            # No GeNomad predictions, all CheckV are unique
            cp {input.checkv_bed} {output.checkv_unique_bed}
        else
            # No CheckV predictions
            touch {output.checkv_unique_bed}
        fi

        # Merge: all GeNomad + unique CheckV
        cat {input.genomad_bed} {output.checkv_unique_bed} > {output.prophage_bed} 2>> {log}

        # Merge FASTA sequences: GeNomad + unique CheckV
        > {output.prophage_fasta}
        if [ -s {input.genomad_fasta} ]; then
            cat {input.genomad_fasta} >> {output.prophage_fasta}
        fi

        # Extract unique CheckV sequences if any
        if [ -s {output.checkv_unique_bed} ]; then
            # Get list of unique CheckV provirus IDs
            awk '{{print $5}}' {output.checkv_unique_bed} > {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/checkv_unique_ids.txt

            # Extract those sequences from CheckV FASTA
            if [ -s {input.checkv_fasta} ]; then
                seqkit grep -f {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/checkv_unique_ids.txt {input.checkv_fasta} >> {output.prophage_fasta} 2>> {log}
            fi
        fi

        # Create prophage table with source column
        echo -e "contig\tstart\tend\ttool\tbin\tsource" > {output.prophage_table}
        awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $4, $6, "unbinned"}}' {output.prophage_bed} >> {output.prophage_table}
        """

rule extract_free_phages:
    input:
        genomad_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad")
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
        virus_fna=$(find {input.genomad_dir} -name "*_virus.fna" | head -1 || true)
        virus_summary=$(find {input.genomad_dir} -name "*_virus_summary.tsv" | head -1 || true)

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
        contigs = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "unbinned_contigs.fasta"),
        prophage_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad_prophages.bed")
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

            # Build lookup table: contig_number -> full_contig_name
            # FASTA headers have sample prefix, e.g., SAMPLE_NODE_1215_length_95187_cov_12.9194
            grep "^>" {input.contigs} | sed 's/>//' | awk -F'_NODE_' '{{
                split($2, parts, "_")
                contig_num = parts[1]
                print contig_num "\t" $0
            }}' > {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/contig_lookup.tsv

            # Convert BED to proper format for bedtools using full contig names
            # BED input columns: contig_num start end ...
            # Output: full_contig_name start end
            awk 'NR==FNR {{lookup[$1]=$2; next}} $1 in lookup {{print lookup[$1] "\t" $2 "\t" $3}}' \
                {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/contig_lookup.tsv \
                {input.prophage_bed} > {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/mask_regions.bed

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

rule checkv_prophage_prediction_unbinned:
    input:
        unbinned_contigs = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "unbinned_contigs.fasta"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_prophage_prediction"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_prophage_prediction_unbinned", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_prophage_prediction_unbinned", "{sample}_bmrk.txt")
    shell:
        """
        # Run CheckV end_to_end on unbinned contigs to find additional prophages
        # Uses same input as GeNomad - deduplication happens via bedtools intersect
        checkv end_to_end \
        {input.unbinned_contigs} \
        {output} \
        -t {threads} \
        -d {input.db} 2> {log}
        """

rule extract_checkv_prophages_unbinned:
    input:
        checkv_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_prophage_prediction")
    conda: config["conda_envs"]["phage_all"]
    output:
        checkv_prophage_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_prophages.bed"),
        checkv_prophage_fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv_prophages.fasta")
    log:
        os.path.join(config["outdir"], "logs", "extract_checkv_prophages_unbinned", "{sample}.log")
    shell:
        """
        # Initialize output files
        touch {output.checkv_prophage_bed}
        touch {output.checkv_prophage_fasta}

        # Find CheckV provirus predictions
        provirus_fna=$(find {input.checkv_dir} -name "proviruses.fna" 2>/dev/null | head -1 || true)
        provirus_tsv=$(find {input.checkv_dir} -name "proviruses.tsv" 2>/dev/null | head -1 || true)

        if [ -f "$provirus_fna" ] && [ -s "$provirus_fna" ] && [ -f "$provirus_tsv" ] && [ -s "$provirus_tsv" ]; then
            # Copy CheckV prophage sequences
            cp "$provirus_fna" {output.checkv_prophage_fasta} 2> {log}

            # Parse CheckV provirus TSV to create BED file with actual coordinates
            # CheckV proviruses.tsv columns: contig_id, start, end, ...
            awk 'NR>1 {{
                if (match($1, /NODE_([0-9]+)_/, arr)) {{
                    # Extract provirus ID from contig_id (format: contig|provirus_X)
                    split($1, parts, "|")
                    provirus_id = parts[2]
                    print arr[1] "\t" $2 "\t" $3 "\tcheckv\t" provirus_id "\tnone"
                }}
            }}' "$provirus_tsv" > {output.checkv_prophage_bed} 2>> {log}
        else
            echo "No CheckV prophages found" >> {log}
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

        # ROBUSTNESS (fork): guard empty FASTA (low-biomass samples) -> header-only summary
        # instead of crashing checkv ("input file is empty").
        mkdir -p {output}
        if grep -q '^>' {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/all_phages_for_checkv.fasta 2>/dev/null; then
            checkv end_to_end \
            {config[outdir]}/{wildcards.sample}/phage_analysis/unbinned/all_phages_for_checkv.fasta \
            {output} \
            -t {threads} \
            -d {input.db} 2> {log}
        else
            printf 'contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n' > {output}/quality_summary.tsv
            echo "Input FASTA has no sequences; wrote header-only quality_summary.tsv" > {log}
        fi
        """
