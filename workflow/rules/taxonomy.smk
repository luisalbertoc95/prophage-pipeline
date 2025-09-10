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

rule gtdbtk_taxonomy:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
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
        
        # Copy input file to genome directory (GTDB-Tk expects genome files in a directory)
        cp {input} {output}/genomes/{wildcards.sample}.fasta
        
        # Set GTDBTK_DATA_PATH environment variable
        export GTDBTK_DATA_PATH={params.db}
        
        # Run GTDB-Tk classify workflow
        gtdbtk classify_wf --genome_dir {output}/genomes --out_dir {output} \
        --cpus {threads} --extension fasta --skip_ani_screen 2> {log}
        """

# Conditional rule selection based on taxonomy method
if config["taxonomy_method"] == "gtdbtk":
    ruleorder: gtdbtk_taxonomy > mmseqs_taxonomy
else:
    ruleorder: mmseqs_taxonomy > gtdbtk_taxonomy
