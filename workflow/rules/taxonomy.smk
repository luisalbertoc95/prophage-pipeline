rule mmseqs_taxonomy:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
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
        
        # Convert input contigs to mmseqs database format
        mmseqs createdb {input} {output}/queryDB
        
        # Run mmseqs taxonomy against your NR database
        mmseqs taxonomy {output}/queryDB {params.db} {output}/taxonomyResult {output}/tmp \
        --search-type 3 --tax-lineage 1 \
        --lca-ranks superkingdom,phylum,class,order,family,genus,species \
        --threads {threads} 2> {log}
        
        # Convert results to TSV format
        mmseqs createtsv {output}/queryDB {output}/taxonomyResult {output}/contig.taxonomy 2>> {log}
        
        # Clean up temporary files
        rm -rf {output}/tmp {output}/queryDB* {output}/taxonomyResult*
        """

rule gtdbtk_classify_bins:
    input:
        bins_done = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins"),
        bins_dir = directory(os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}_DASTool_bins"))
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
        
        # Copy all bin files (excluding unbinned.fa) to genomes directory for GTDB-Tk
        find {input.bins_dir} -name "*.fa" ! -name "unbinned.fa" -exec cp {{}} {output}/genomes/ \;
        
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

# Conditional rule selection based on taxonomy method
if config["taxonomy_method"] == "gtdbtk":
    ruleorder: gtdbtk_classify_bins > mmseqs_taxonomy
else:
    ruleorder: mmseqs_taxonomy > gtdbtk_classify_bins
