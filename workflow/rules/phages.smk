# Shared phage-related rules and final output generation
# Database downloads and final prophage table/FASTA creation

import os

# Database download rules
rule genomad_db:
    output: directory(config["genomad_database"])
    conda: config["conda_envs"]["genomad"]
    shell:
        "genomad download-database ref"

# ORIGINAL APPROACH: Run GeNomad on complete assembly for full metagenomic context
# This rule is kept for comparison but not used in the per-MAG approach
rule genomad_complete_assembly_ORIGINAL:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        db = config["genomad_database"]
    threads: 24
    conda: config["conda_envs"]["genomad"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad_complete_ORIGINAL"))
    log:
        os.path.join(config["outdir"], "logs", "genomad_complete_assembly_ORIGINAL", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "genomad_complete_assembly_ORIGINAL", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output}
        genomad end-to-end --cleanup --threads {threads} \
        {input.contigs} {output} {input.db} 2> {log}
        """

rule download_bakta_db:
    output: directory(config["bakta_database"])
    conda: config["conda_envs"]["bakta"]
    shell:
        "bakta_db download --output {output} --type full"

rule checkv_db:
    output:
        directory(config["checkv_database"])
    conda: config["conda_envs"]["checkv"]
    shell:
        """
        checkv download_database {output}
        """

# Merge prophage tables from MAGs and unbinned contigs
rule merge_prophage_tables:
    input:
        mag_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "prophage_table.tsv"),
        unbinned_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophage_table.tsv")
    output:
        merged_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv")
    shell:
        """
        # Combine tables (both have header with: contig, start, end, tool, bin, source)
        head -1 {input.mag_table} > {output.merged_table}
        tail -n +2 {input.mag_table} >> {output.merged_table}
        tail -n +2 {input.unbinned_table} >> {output.merged_table}
        """

# Add taxonomy information to prophage table
rule add_taxonomy_to_prophage_table:
    input:
        basic_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv"),
        mmseqs_taxonomy = os.path.join(config["outdir"], "{sample}", "taxonomy", "mmseqs"),
        gtdbtk_taxonomy = os.path.join(config["outdir"], "{sample}", "taxonomy", "gtdbtk")
    conda: config["conda_envs"]["phage_all"]
    output:
        table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    log:
        os.path.join(config["outdir"], "logs", "add_taxonomy_to_prophage_table", "{sample}.log")
    script:
        "../scripts/add_taxonomy_to_prophage_table.R"

# Create final combined FASTA outputs
rule final_prophage_outputs:
    input:
        mag_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "all_prophages.fasta"),
        unbinned_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophages.fasta"),
        free_phages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "free_phages.fasta"),
        prophage_table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    output:
        all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "all_prophages_combined.fasta"),
        prophages_from_mags = os.path.join(config["outdir"], "{sample}", "phage_analysis", "prophages_from_mags.fasta"),
        prophages_from_unbinned = os.path.join(config["outdir"], "{sample}", "phage_analysis", "prophages_from_unbinned.fasta"),
        free_phages_final = os.path.join(config["outdir"], "{sample}", "phage_analysis", "free_phages.fasta")
    shell:
        """
        # Copy individual sources
        cp {input.mag_prophages} {output.prophages_from_mags}
        cp {input.unbinned_prophages} {output.prophages_from_unbinned}
        cp {input.free_phages} {output.free_phages_final}

        # Combine all prophages (MAG + unbinned, but not free phages)
        cat {input.mag_prophages} {input.unbinned_prophages} > {output.all_prophages}
        """

# CheckV quality assessment for all prophages (MAG + unbinned combined)
rule checkv_all_prophages:
    input:
        all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "all_prophages_combined.fasta"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv_all_prophages"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_all_prophages", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_all_prophages", "{sample}_bmrk.txt")
    shell:
        """
        # Run checkv on all prophages
        checkv end_to_end \
        {input.all_prophages} \
        {output} \
        -t {threads} \
        -d {input.db} 2> {log}
        """

# CheckV quality assessment for free phages only
rule checkv_free_phages:
    input:
        free_phages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "free_phages.fasta"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv_free_phages"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_free_phages", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_free_phages", "{sample}_bmrk.txt")
    shell:
        """
        # Run checkv on free phages
        checkv end_to_end \
        {input.free_phages} \
        {output} \
        -t {threads} \
        -d {input.db} 2> {log}
        """

# Optional comparison rule: Compare per-MAG vs complete assembly GeNomad approaches
rule compare_genomad_approaches:
    input:
        per_mag_mags = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "genomad_all.bed"),
        per_mag_unbinned_genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "genomad"),
        original_complete = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad_complete_ORIGINAL")
    conda: config["conda_envs"]["phage_all"]
    output:
        comparison_report = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad_comparison_report.txt"),
        per_mag_all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "per_mag_all_prophages.tsv"),
        original_all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "original_all_prophages.tsv")
    log:
        os.path.join(config["outdir"], "logs", "compare_genomad_approaches", "{sample}.log")
    shell:
        """
        # Collect per-MAG approach prophages
        echo -e "contig\tstart\tend\tsource" > {output.per_mag_all_prophages}

        # From MAGs
        if [ -s {input.per_mag_mags} ]; then
            awk '{{print $1 "\t" $2 "\t" $3 "\tMAG"}}' {input.per_mag_mags} >> {output.per_mag_all_prophages}
        fi

        # From unbinned
        unbinned_tsv=$(find {input.per_mag_unbinned_genomad} -name "*_provirus.tsv" 2>/dev/null | head -1)
        if [ -f "$unbinned_tsv" ] && [ -s "$unbinned_tsv" ]; then
            awk 'NR>1 {{
                if (match($2, /NODE_([0-9]+)_/, arr)) {{
                    print arr[1] "\t" $3 "\t" $4 "\tunbinned"
                }}
            }}' "$unbinned_tsv" >> {output.per_mag_all_prophages}
        fi

        # Collect original complete assembly approach prophages
        echo -e "contig\tstart\tend\tsource" > {output.original_all_prophages}
        original_tsv=$(find {input.original_complete} -name "*_provirus.tsv" 2>/dev/null | head -1)
        if [ -f "$original_tsv" ] && [ -s "$original_tsv" ]; then
            awk 'NR>1 {{
                if (match($2, /NODE_([0-9]+)_/, arr)) {{
                    print arr[1] "\t" $3 "\t" $4 "\tcomplete"
                }}
            }}' "$original_tsv" >> {output.original_all_prophages}
        fi

        # Create comparison report
        {{
            echo "GeNomad Approach Comparison Report"
            echo "==================================="
            echo ""
            echo "Sample: {wildcards.sample}"
            echo "Date: $(date)"
            echo ""
            echo "Per-MAG Approach (MAGs + Unbinned separately):"
            echo "  Total prophages: $(tail -n +2 {output.per_mag_all_prophages} | wc -l)"
            echo "  From MAGs: $(grep -c 'MAG' {output.per_mag_all_prophages} || echo 0)"
            echo "  From unbinned: $(grep -c 'unbinned' {output.per_mag_all_prophages} || echo 0)"
            echo ""
            echo "Complete Assembly Approach (all contigs at once):"
            echo "  Total prophages: $(tail -n +2 {output.original_all_prophages} | wc -l)"
            echo ""
            echo "Difference:"
            per_mag_count=$(tail -n +2 {output.per_mag_all_prophages} | wc -l)
            original_count=$(tail -n +2 {output.original_all_prophages} | wc -l)
            diff=$((original_count - per_mag_count))
            echo "  Complete assembly found $diff more prophages than per-MAG approach"
            echo ""
            echo "Files for detailed comparison:"
            echo "  Per-MAG: {output.per_mag_all_prophages}"
            echo "  Original: {output.original_all_prophages}"
        }} > {output.comparison_report} 2> {log}

        cat {output.comparison_report}
        """

# Final rule to complete all phage analysis
rule run_everything:
    input:
        coverm_stats = os.path.join(config["outdir"], "{sample}", "coverm", "{sample}_stats.txt"),
        checkm = os.path.join(config["outdir"], "{sample}", "binning", "checkm"),
        checkv_unbinned = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv"),
        checkv_all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv_all_prophages"),
        checkv_free_phages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv_free_phages"),
        all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "all_prophages_combined.fasta"),
        prophages_from_mags = os.path.join(config["outdir"], "{sample}", "phage_analysis", "prophages_from_mags.fasta"),
        prophages_from_unbinned = os.path.join(config["outdir"], "{sample}", "phage_analysis", "prophages_from_unbinned.fasta"),
        free_phages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "free_phages.fasta"),
        prophage_table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "touch {output}"
