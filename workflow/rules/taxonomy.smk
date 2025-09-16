rule create_gtdb_mmseqs_db:
    input:
        gtdbtk_db = config["gtdbtk_database"]
    output:
        directory(config["gtdb_mmseqs_database"])
    conda: config["conda_envs"]["mmseqs"]
    threads: 8
    log:
        os.path.join(config["outdir"], "logs", "gtdb_mmseqs_db_creation.log")
    shell:
        """
        # Create GTDB-MMseqs database from GTDB-Tk database files
        mkdir -p {output}
        
        # Find GTDB genome files (they should be in the GTDB-Tk database directory)
        find {input.gtdbtk_db} -name "*.fna" -o -name "*.fa" -o -name "*.fasta" > {output}/genome_files.txt
        
        if [ ! -s {output}/genome_files.txt ]; then
            echo "ERROR: No GTDB genome files found in {input.gtdbtk_db}" > {log}
            echo "Please check GTDB-Tk database structure" >> {log}
            exit 1
        fi
        
        # Concatenate all GTDB genomes
        cat $(cat {output}/genome_files.txt) > {output}/gtdb_genomes.fna
        
        # Create MMseqs database
        mmseqs createdb {output}/gtdb_genomes.fna {output}/gtdb_mmseqs_db --threads {threads} 2> {log}
        
        # Create taxonomy database (this will need GTDB taxonomy mapping)
        # For now, create basic database - we'll enhance with taxonomy later
        echo "GTDB-MMseqs database created successfully" >> {log}
        """

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
