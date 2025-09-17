
# Function to determine taxonomy rule inputs based on scope
def get_mmseqs_input(wildcards):
    inputs = {
        "contigs": os.path.join(config["outdir"], wildcards.sample, "binning", "final_filtered_contigs.fasta")
    }
    # Only require prophage table for prophage_only mode
    if config.get("taxonomy_scope", "prophage_only") == "prophage_only":
        inputs["prophage_table"] = os.path.join(config["outdir"], wildcards.sample, "phage_analysis", "final_prophage_table.tsv")
    return inputs

rule mmseqs_taxonomy:
    input:
        unpack(get_mmseqs_input)
    params:
        db = config["mmseqs_database"],
        taxonomy_scope = config.get("taxonomy_scope", "prophage_only")
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
        
        if [ "{params.taxonomy_scope}" = "all_contigs" ]; then
            echo "Running taxonomy on ALL contigs (comprehensive mode)" > {log}
            # Use all contigs for taxonomy
            cp {input.contigs} {output}/contigs_for_taxonomy.fasta
        else
            echo "Running taxonomy on prophage-containing contigs only (targeted mode)" > {log}
            # Construct prophage table path manually (avoid {input.prophage_table} placeholder)
            sample_name=$(basename $(dirname $(dirname {output})))
            prophage_table_file="$(dirname $(dirname {output}))/phage_analysis/final_prophage_table.tsv"
            
            # Extract list of prophage-containing contigs
            awk 'NR>1 {{print "NODE_" $1 "_"}}' "$prophage_table_file" | sort -u > {output}/prophage_contigs.txt
            
            # Extract only prophage-containing contigs from the full contig set
            seqkit grep -r -f {output}/prophage_contigs.txt {input.contigs} > {output}/contigs_for_taxonomy.fasta 2>> {log}
        fi
        
        # Check if any contigs were extracted/available
        if [ ! -s {output}/contigs_for_taxonomy.fasta ]; then
            echo "No contigs found for taxonomy analysis" >> {log}
            touch {output}/contig.taxonomy
            exit 0
        fi
        
        # Convert contigs to mmseqs database format
        mmseqs createdb {output}/contigs_for_taxonomy.fasta {output}/queryDB 2>> {log}
        
        # Run mmseqs taxonomy against NR database
        mmseqs taxonomy {output}/queryDB {params.db} {output}/taxonomyResult {output}/tmp \
        --search-type 3 --tax-lineage 1 \
        --lca-ranks superkingdom,phylum,class,order,family,genus,species \
        --threads {threads} 2>> {log}
        
        # Convert results to TSV format
        mmseqs createtsv {output}/queryDB {output}/taxonomyResult {output}/contig.taxonomy 2>> {log}
        
        # Clean up temporary files
        rm -rf {output}/tmp {output}/queryDB* {output}/taxonomyResult*
        """

rule gtdbtk_classify_bins:
    input:
        bins_done = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins")
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
