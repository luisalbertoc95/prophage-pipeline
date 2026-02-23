
# =============================================================================
# MMseqs Taxonomy - Batched across all samples
# =============================================================================
# To reduce memory overhead, we combine contigs from all samples, run MMseqs
# once, then split results back to per-sample directories.
# This reduces N × 700GB jobs to 1 × 700GB job.
# =============================================================================

# Rule priority based on configuration
if config.get("taxonomy_scope", "prophage_only") == "all_contigs":
    ruleorder: combine_contigs_for_mmseqs_all > combine_contigs_for_mmseqs_prophage_only
else:
    ruleorder: combine_contigs_for_mmseqs_prophage_only > combine_contigs_for_mmseqs_all


# -----------------------------------------------------------------------------
# Step 1a: Prepare and combine contigs (prophage_only mode)
# -----------------------------------------------------------------------------
rule combine_contigs_for_mmseqs_prophage_only:
    input:
        masked_contigs = expand(
            os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "masked_contigs.fasta"),
            sample=SAMPLES
        ),
        prophage_tables = expand(
            os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv"),
            sample=SAMPLES
        )
    output:
        combined = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "combined_contigs.fasta"),
        sample_list = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "samples.txt")
    params:
        samples = SAMPLES,
        outdir = config["outdir"]
    log:
        os.path.join(config["outdir"], "logs", "mmseqs", "combine_contigs.log")
    conda: config["conda_envs"]["phage_all"]
    shell:
        """
        set -ue
        mkdir -p $(dirname {output.combined})

        echo "Combining prophage-containing contigs from all samples (prophage_only mode)" > {log}

        # Clear output files
        > {output.combined}
        > {output.sample_list}

        # Process each sample
        for sample in {params.samples}; do
            echo "$sample" >> {output.sample_list}

            masked_file="{params.outdir}/$sample/phage_analysis/unbinned/masked_contigs.fasta"
            prophage_table="{params.outdir}/$sample/phage_analysis/final_prophage_table.tsv"

            # Extract list of unbinned prophage-containing contigs (bin == "none")
            awk 'NR>1 && $5=="none" {{print "NODE_" $1 "_"}}' "$prophage_table" | sort -u > /tmp/${{sample}}_prophage_contigs.txt

            # Check if there are any prophage contigs
            if [ -s /tmp/${{sample}}_prophage_contigs.txt ]; then
                # Extract matching contigs and add sample prefix to headers
                seqkit grep -r -f /tmp/${{sample}}_prophage_contigs.txt "$masked_file" 2>> {log} | \
                    sed "s/^>/>${{sample}}__/" >> {output.combined}

                n_contigs=$(grep -c "^>" /tmp/${{sample}}_prophage_contigs.txt || echo 0)
                echo "  $sample: $n_contigs prophage-containing contigs" >> {log}
            else
                echo "  $sample: no prophage-containing contigs" >> {log}
            fi

            rm -f /tmp/${{sample}}_prophage_contigs.txt
        done

        total=$(grep -c "^>" {output.combined} || echo 0)
        echo "Total combined contigs: $total" >> {log}
        """


# -----------------------------------------------------------------------------
# Step 1b: Prepare and combine contigs (all_contigs mode)
# -----------------------------------------------------------------------------
rule combine_contigs_for_mmseqs_all:
    input:
        masked_contigs = expand(
            os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "masked_contigs.fasta"),
            sample=SAMPLES
        )
    output:
        combined = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "combined_contigs.fasta"),
        sample_list = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "samples.txt")
    params:
        samples = SAMPLES,
        outdir = config["outdir"]
    log:
        os.path.join(config["outdir"], "logs", "mmseqs", "combine_contigs.log")
    shell:
        """
        set -ue
        mkdir -p $(dirname {output.combined})

        echo "Combining all masked contigs from all samples (all_contigs mode)" > {log}

        # Clear output files
        > {output.combined}
        > {output.sample_list}

        # Process each sample
        for sample in {params.samples}; do
            echo "$sample" >> {output.sample_list}

            masked_file="{params.outdir}/$sample/phage_analysis/unbinned/masked_contigs.fasta"

            if [ -s "$masked_file" ]; then
                # Add sample prefix to headers and append
                sed "s/^>/>${{sample}}__/" "$masked_file" >> {output.combined}

                n_contigs=$(grep -c "^>" "$masked_file" || echo 0)
                echo "  $sample: $n_contigs contigs" >> {log}
            else
                echo "  $sample: no contigs" >> {log}
            fi
        done

        total=$(grep -c "^>" {output.combined} || echo 0)
        echo "Total combined contigs: $total" >> {log}
        """


# -----------------------------------------------------------------------------
# Step 2: Run MMseqs on combined contigs (single job for all samples)
# -----------------------------------------------------------------------------
rule mmseqs_taxonomy_combined:
    input:
        combined = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "combined_contigs.fasta")
    output:
        taxonomy = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "combined.taxonomy"),
        done = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "mmseqs.done")
    params:
        db = config["mmseqs_database"],
        outdir = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined")
    threads: 24
    conda: config["conda_envs"]["mmseqs"]
    log:
        os.path.join(config["outdir"], "logs", "mmseqs", "combined_taxonomy.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "mmseqs", "combined_taxonomy_bmrk.txt")
    shell:
        """
        set -ue

        echo "Running MMseqs taxonomy on combined contigs from all samples" > {log}

        # Check if there are any contigs to process
        if [ ! -s {input.combined} ]; then
            echo "No contigs to process" >> {log}
            touch {output.taxonomy}
            touch {output.done}
            exit 0
        fi

        # Run mmseqs easy-taxonomy
        mmseqs easy-taxonomy {input.combined} {params.db} \
            {params.outdir}/contig {params.outdir}/tmp \
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

        # Rename output to final filename
        if [ -f {params.outdir}/contig_lca.tsv ]; then
            mv {params.outdir}/contig_lca.tsv {output.taxonomy}
        else
            echo "Warning: mmseqs did not produce expected output" >> {log}
            touch {output.taxonomy}
        fi

        # Clean up temporary files
        rm -rf {params.outdir}/tmp

        touch {output.done}
        """


# -----------------------------------------------------------------------------
# Step 3: Split combined results back to per-sample directories
# -----------------------------------------------------------------------------
rule split_mmseqs_taxonomy_results:
    input:
        taxonomy = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "combined.taxonomy"),
        done = os.path.join(config["outdir"], "taxonomy", "mmseqs_combined", "mmseqs.done")
    output:
        directory(os.path.join(config["outdir"], "{sample}", "taxonomy", "mmseqs"))
    params:
        sample = "{sample}"
    log:
        os.path.join(config["outdir"], "logs", "mmseqs", "{sample}_split.log")
    shell:
        """
        set -ue
        mkdir -p {output}

        echo "Extracting taxonomy results for sample: {params.sample}" > {log}

        # Extract lines for this sample (header starts with SAMPLE__)
        # Remove the sample prefix from contig names in output
        grep "^{params.sample}__" {input.taxonomy} 2>/dev/null | \
            sed "s/^{params.sample}__//" > {output}/contig.taxonomy || true

        n_results=$(wc -l < {output}/contig.taxonomy || echo 0)
        echo "Extracted $n_results taxonomy assignments" >> {log}

        # Ensure file exists even if empty
        touch {output}/contig.taxonomy
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
