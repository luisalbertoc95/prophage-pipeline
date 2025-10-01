# Shared phage-related rules and final output generation
# Database downloads and final prophage table/FASTA creation

import os

# Database download rules
rule genomad_db:
    output: directory(config["genomad_database"])
    conda: config["conda_envs"]["genomad"]
    shell:
        "genomad download-database ref"

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

# CheckM for MAG quality assessment (runs on original bins, not prophage-related)
rule checkm:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins")
    threads: 24
    conda: config["conda_envs"]["checkm"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "binning", "checkm"))
    log:
        os.path.join(config["outdir"], "logs", "checkm", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkm", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output}
        checkm lineage_wf -x fa \
        {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins/ \
        {output}/ -t {threads} --tab_table -f {output}/checkm_out.tsv 2> {log}
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

# Final rule to complete all phage analysis
rule run_everything:
    input:
        coverm_stats = os.path.join(config["outdir"], "{sample}", "coverm", "{sample}_stats.txt"),
        checkm = os.path.join(config["outdir"], "{sample}", "binning", "checkm"),
        checkv = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "checkv"),
        all_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "all_prophages_combined.fasta"),
        prophages_from_mags = os.path.join(config["outdir"], "{sample}", "phage_analysis", "prophages_from_mags.fasta"),
        prophages_from_unbinned = os.path.join(config["outdir"], "{sample}", "phage_analysis", "prophages_from_unbinned.fasta"),
        free_phages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "free_phages.fasta"),
        prophage_table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "touch {output}"
