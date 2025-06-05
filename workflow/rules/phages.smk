rule genomad_db:
    output: directory(config["genomad_database"])
    conda: config["conda_envs"]["genomad"]
    shell:
        "genomad download-database ref"

rule genomad:
    input:
        contigs = os.path.join(config["outdir"], "binning", "{sample}", "final_filtered_contigs.fasta"),
        db = config["genomad_database"]
    threads: 24
    conda: config["conda_envs"]["genomad"]
    output:
        directory(os.path.join(config["outdir"], "phage_analysis", "{sample}", "genomad"))
    log:
        os.path.join(config["outdir"], "logs", "genomad", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "genomad", "{sample}_bmrk.txt")
    shell:
        """
        # Create the output directory
        mkdir -p {output}
        # Run genomad
        genomad end-to-end --cleanup --threads {threads} \
        --splits 16 \
        {input.contigs} {output} {input.db} 2> {log}
        """

rule download_bakta_db:
    output: directory(config["bakta_database"])
    conda: config["conda_envs"]["bakta"]
    shell:
        "bakta_db download --output {output} --type full"

rule bakta:
    input:
        contigs = os.path.join(config["outdir"], "binning", "{sample}", "final_filt_contigs_5000.fasta"),
        db = config["bakta_database"]
    threads: 24
    conda: config["conda_envs"]["bakta"]
    output: 
        directory(os.path.join(config["outdir"], "phage_analysis", "{sample}", "bakta"))
    log:
        os.path.join(config["outdir"], "logs", "bakta", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "bakta", "{sample}_bmrk.txt")
    shell:
        """
        bakta --db {input.db}/db --force --skip-plot --output {output} \
        --threads {threads} {input.contigs} 2> {log}
        """

rule phispy:
    input:
        os.path.join(config["outdir"], "phage_analysis", "{sample}", "bakta")
    threads: 24
    conda: config["conda_envs"]["phispy"]
    output:
        directory(os.path.join(config["outdir"], "phage_analysis", "{sample}", "phispy"))
    log:
        os.path.join(config["outdir"], "logs", "phispy", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "phispy", "{sample}_bmrk.txt")
    shell:
        "PhiSpy.py {input}/*.gbff -o {output} --output_choice 63 2> {log} || true"

rule phage_all:
    input:
        genomad = os.path.join(config["outdir"], "phage_analysis", "{sample}", "genomad"),
        phispy = os.path.join(config["outdir"], "phage_analysis", "{sample}", "phispy"),
        CAT = os.path.join(config["outdir"], "taxonomy", "{sample}", "CAT")
    conda: config["conda_envs"]["phage_all"]
    output:
        fasta = os.path.join(config["outdir"], "phage_analysis", "{sample}", "unique_phispy_prophage.fasta"),
        table = os.path.join(config["outdir"], "phage_analysis", "{sample}", "final_prophage_table.tsv"),
        table_with_taxonomy = os.path.join(config["outdir"], "phage_analysis", "{sample}", "final_prophage_table_with_host_taxonomy.tsv")
    script:
        "../scripts/merge_prophages.R"

rule final_prophage_output:
    input:
        genomad = os.path.join(config["outdir"], "phage_analysis", "{sample}", "genomad"),
        unique_phispy = os.path.join(config["outdir"], "phage_analysis", "{sample}", "unique_phispy_prophage.fasta")
    output:
        os.path.join(config["outdir"], "phage_analysis", "{sample}", "final_prophage.fasta")
    shell:
        """
        cp {input.genomad}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.fna \
        {config[outdir]}/phage_analysis/{wildcards.sample}/genomad_prophage.fasta

        cat {input.genomad}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.fna \
        {input.unique_phispy} > {output}
        """

rule checkv_db:
    output:
        directory(config["checkv_database"])
    conda: config["conda_envs"]["checkv"]
    shell:
        """
        checkv download_database {output}
        """    

rule checkv:
    input:
        fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage.fasta"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv"))
    log:
        os.path.join(config["outdir"], "logs", "checkv", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv", "{sample}_bmrk.txt")
    shell:
        """
        checkv end_to_end \
        {input.fasta} {output} \
        -t {threads} \
        -d {input.db}
        """

rule run_everything:
    input:
        coverm = os.path.join(config["outdir"], "coverm", "{sample}"),
        checkm = os.path.join(config["outdir"], "binning", "{sample}", "checkm"),
        checkv = os.path.join(config["outdir"], "phage_analysis", "{sample}", "checkv"),
        prophage = os.path.join(config["outdir"], "phage_analysis", "{sample}", "final_prophage.fasta")
    output:
        os.path.join(config["outdir"], "phage_analysis", "{sample}", "done")
    shell:
        "touch {config[outdir]}/phage_analysis/{wildcards.sample}/done"
