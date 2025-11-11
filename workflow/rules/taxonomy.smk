
# Rule priority based on configuration
if config.get("taxonomy_scope", "prophage_only") == "all_contigs":
    ruleorder: mmseqs_taxonomy_all_contigs > mmseqs_taxonomy_prophage_only
else:
    ruleorder: mmseqs_taxonomy_prophage_only > mmseqs_taxonomy_all_contigs

rule mmseqs_taxonomy_all_contigs:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        masked_unbinned = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "masked_contigs.fasta")
    params:
        db = config["mmseqs_database"]
    threads: 24
    conda: config["conda_envs"]["mmseqs"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "taxonomy", "mmseqs"))
    log:
        os.path.join(config["outdir"], "logs", "mmseqs", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "mmseqs", "{sample}_bmrk.txt")
    shell:
        """
        set -ue
        mkdir -p {output}

        echo "Running taxonomy on masked unbinned contigs (comprehensive mode)" > {log}
        cp {input.masked_unbinned} {output}/contigs_for_taxonomy.fasta

        # Create temporary directory
        TMP_DIR=$(mktemp -d)

        # Run mmseqs easy-taxonomy (matches phage-analysis pipeline for better performance)
        mmseqs easy-taxonomy {output}/contigs_for_taxonomy.fasta {params.db} \
            {output}/contig {output}/tmp \
            --min-length 30 \
            -e 1e-15 \
            --search-type 2 \
            -s 4.0 \
            --shuffle 0 \
            --lca-mode 2 \
            -a \
            --tax-lineage 2 \
            --threads {threads} \
            --split-mode 0 \
            --orf-filter 1 \
            >> {log} 2>&1

        # Rename output to match expected filename
        if [ -f {output}/contig_lca.tsv ]; then
            mv {output}/contig_lca.tsv {output}/contig.taxonomy
        else
            echo "Warning: mmseqs did not produce expected output" >> {log}
            touch {output}/contig.taxonomy
        fi

        # Clean up temporary files
        rm -rf {output}/tmp $TMP_DIR
        """

rule mmseqs_taxonomy_prophage_only:
    input:
        masked_unbinned = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "masked_contigs.fasta"),
        prophage_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv")
    params:
        db = config["mmseqs_database"]
    threads: 24
    conda: config["conda_envs"]["mmseqs"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "taxonomy", "mmseqs"))
    log:
        os.path.join(config["outdir"], "logs", "mmseqs", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "mmseqs", "{sample}_bmrk.txt")
    shell:
        """
        set -ue
        mkdir -p {output}

        echo "Running taxonomy on masked prophage-containing contigs only (targeted mode)" > {log}

        # Extract list of unbinned prophage-containing contigs (bin == "none")
        awk 'NR>1 && $5=="none" {{print "NODE_" $1 "_"}}' {input.prophage_table} | sort -u > {output}/prophage_contigs.txt

        # Extract only unbinned prophage-containing contigs from masked contigs
        seqkit grep -r -f {output}/prophage_contigs.txt {input.masked_unbinned} > {output}/contigs_for_taxonomy.fasta 2>> {log}
        
        # Check if any contigs were extracted
        if [ ! -s {output}/contigs_for_taxonomy.fasta ]; then
            echo "No prophage-containing contigs found" >> {log}
            touch {output}/contig.taxonomy
            exit 0
        fi

        # Create temporary directory
        TMP_DIR=$(mktemp -d)

        # Run mmseqs easy-taxonomy (matches phage-analysis pipeline for better performance)
        mmseqs easy-taxonomy {output}/contigs_for_taxonomy.fasta {params.db} \
            {output}/contig {output}/tmp \
            --min-length 30 \
            -e 1e-15 \
            --search-type 2 \
            -s 4.0 \
            --shuffle 0 \
            --lca-mode 2 \
            -a \
            --tax-lineage 2 \
            --threads {threads} \
            --split-mode 0 \
            --orf-filter 1 \
            >> {log} 2>&1

        # Rename output to match expected filename
        if [ -f {output}/contig_lca.tsv ]; then
            mv {output}/contig_lca.tsv {output}/contig.taxonomy
        else
            echo "Warning: mmseqs did not produce expected output" >> {log}
            touch {output}/contig.taxonomy
        fi

        # Clean up temporary files
        rm -rf {output}/tmp $TMP_DIR
        """

rule gtdbtk_classify_bins:
    input:
        bins_done = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins"),
        mag_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "prophage_table.tsv")
    params:
        db = config["gtdbtk_database"]
    threads: 24
    conda: config["conda_envs"]["gtdbtk"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "taxonomy", "gtdbtk"))
    log:
        os.path.join(config["outdir"], "logs", "gtdbtk", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "gtdbtk", "{sample}_bmrk.txt")
    shell:
        """
        set -ue
        mkdir -p {output}/genomes
        
        # Copy only actual MAG bins (bin.*.fa) to genomes directory for GTDB-Tk
        find {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins -name "bin.*.fa" -exec cp {{}} {output}/genomes/ \\;
        
        # Only run GTDB-Tk if there are bin files to process
        if [ $(find {output}/genomes -name "*.fa" | wc -l) -gt 0 ]; then
            # Set GTDBTK_DATA_PATH environment variable
            export GTDBTK_DATA_PATH={params.db}
            
            # Run GTDB-Tk classify workflow on individual MAG bins
            gtdbtk classify_wf --genome_dir {output}/genomes --out_dir {output} \
            --cpus {threads} --extension fa --skip_ani_screen 2> {log}
        else
            echo "No MAG bins found for GTDB-Tk classification" > {log}
            # Create empty output files to satisfy Snakemake
            mkdir -p {output}/classify
            touch {output}/gtdbtk.bac120.summary.tsv
            touch {output}/gtdbtk.ar53.summary.tsv
        fi
        """

# Both taxonomy methods will run for hybrid approach
# GTDB-Tk for binned contigs, MMseqs for unbinned contigs
