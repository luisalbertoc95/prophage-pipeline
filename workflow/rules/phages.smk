rule genomad_db:
    output: directory(config["genomad_database"])
    conda: config["conda_envs"]["genomad"]
    shell:
        "genomad download-database ref"

rule genomad:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        db = config["genomad_database"]
    threads: 24
    conda: config["conda_envs"]["genomad"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"))
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
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta"),
        db = config["bakta_database"]
    threads: 24
    conda: config["conda_envs"]["bakta"]
    output: 
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta"))
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
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta")
    threads: 24
    conda: config["conda_envs"]["phispy"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"))
    log:
        os.path.join(config["outdir"], "logs", "phispy", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "phispy", "{sample}_bmrk.txt")
    shell:
        "PhiSpy.py {input}/*.gbff -o {output} --output_choice 63 2> {log} || true"

rule download_pide_model:
    conda: config["conda_envs"]["pide"]
    output:
        model = config["pide_model"]
    log:
        os.path.join(config["outdir"], "logs", "pide_download.log")
    shell:
        """
        # Create directory and ensure we have write permissions
        mkdir -p $(dirname {output.model})
        cd $(dirname {output.model})
        
        # Remove any partial downloads
        rm -f PIDE.model.tar.gz*
        
        # Download with retry and timeout options
        wget --timeout=30 --tries=3 --retry-connrefused \
        https://zenodo.org/records/12759619/files/PIDE.model.tar.gz 2> {log}
        
        # Verify download succeeded before extracting
        if [ -f "PIDE.model.tar.gz" ]; then
            tar xzvf PIDE.model.tar.gz 2>> {log}
            rm PIDE.model.tar.gz
        else
            echo "Download failed - PIDE.model.tar.gz not found" >> {log}
            exit 1
        fi
        """

rule clone_pide:
    input:
        model = config["pide_model"]  # Wait for model download to complete first
    conda: config["conda_envs"]["pide"]
    output:
        pide_script = os.path.join(config["pide_repository"], "classification.py")
    log:
        os.path.join(config["outdir"], "logs", "pide_clone.log")
    shell:
        """
        mkdir -p $(dirname {config[pide_repository]})
        cd $(dirname {config[pide_repository]})
        git clone https://github.com/chyghy/PIDE.git 2> {log}
        """

rule pide:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        model = config["pide_model"],
        script = os.path.join(config["pide_repository"], "classification.py")
    resources:
        mem_mb=32000  # 32GB - generous allocation for ESM-2 model + safety margin
    threads: 24
    conda: config["conda_envs"]["pide"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "pide"))
    log:
        os.path.join(config["outdir"], "logs", "pide", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "pide", "{sample}_bmrk.txt")
    shell:
        """
        # Create output directory
        mkdir -p {output}
        
        # Run PIDE prophage detection
        cd {config[pide_repository]}
        python classification.py -o {output} {input.contigs} {input.model} 2> {log}
        """

rule prophage_tool_comparison:
    input:
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"),
        pide = os.path.join(config["outdir"], "{sample}", "phage_analysis", "pide")
    conda: config["conda_envs"]["phage_all"]
    output:
        raw_predictions = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "raw_predictions.tsv"),
        tool_stats = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "tool_statistics.tsv"),
        overlap_analysis = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "overlap_analysis.tsv"),
        agreement_summary = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "agreement_summary.tsv"),
        unique_predictions = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "unique_predictions.tsv"),
        unique_summary = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "unique_summary.tsv"),
        plots = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "comparison_plots.pdf")
    script:
        "../scripts/compare_prophage_tools.R"

rule phage_all:
    input:
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"),
        mmseqs = os.path.join(config["outdir"], "{sample}", "taxonomy", "mmseqs")
    conda: config["conda_envs"]["phage_all"]
    output:
        fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unique_phispy_prophage.fasta"),
        table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv"),
        table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    script:
        "../scripts/merge_prophages.R"

rule final_prophage_output:
    input:
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        unique_phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unique_phispy_prophage.fasta")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage.fasta")
    shell:
        """
        cp {input.genomad}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.fna \
        {config[outdir]}/{wildcards.sample}/phage_analysis/genomad_prophage.fasta

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

# Individual tool CheckV rules for tool comparison
rule checkv_genomad:
    input:
        genomad_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "tool_comparison", "checkv_genomad"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_genomad", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_genomad", "{sample}_bmrk.txt")
    shell:
        """
        # Use geNomad provirus sequences
        checkv end_to_end \
        {input.genomad_dir}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.fna \
        {output} -t {threads} -d {input.db} 2> {log}
        """

rule checkv_phispy:
    input:
        phispy_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "tool_comparison", "checkv_phispy"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_phispy", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_phispy", "{sample}_bmrk.txt")
    shell:
        """
        # Use PhiSpy prophage sequences
        checkv end_to_end \
        {input.phispy_dir}/phage.fasta \
        {output} -t {threads} -d {input.db} 2> {log}
        """

rule checkv_pide:
    input:
        pide_dir = os.path.join(config["outdir"], "{sample}", "phage_analysis", "pide"),
        db = config["checkv_database"]
    threads: 24
    conda: config["conda_envs"]["checkv"]
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "tool_comparison", "checkv_pide"))
    log:
        os.path.join(config["outdir"], "logs", "checkv_pide", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkv_pide", "{sample}_bmrk.txt")
    shell:
        """
        # Use PIDE prediction sequences - need to verify correct file name
        # This may need adjustment based on actual PIDE output format
        if [ -f {input.pide_dir}/predictions.fasta ]; then
            checkv end_to_end {input.pide_dir}/predictions.fasta {output} -t {threads} -d {input.db} 2> {log}
        elif [ -f {input.pide_dir}/prophage_sequences.fasta ]; then
            checkv end_to_end {input.pide_dir}/prophage_sequences.fasta {output} -t {threads} -d {input.db} 2> {log}
        else
            echo "Could not find PIDE sequence file - checking directory contents:" > {log}
            ls -la {input.pide_dir}/ >> {log}
            exit 1
        fi
        """

rule enhanced_tool_comparison:
    input:
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"),
        pide = os.path.join(config["outdir"], "{sample}", "phage_analysis", "pide"),
        checkv_genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "tool_comparison", "checkv_genomad"),
        checkv_phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "tool_comparison", "checkv_phispy"),
        checkv_pide = os.path.join(config["outdir"], "{sample}", "phage_analysis", "tool_comparison", "checkv_pide")
    conda: config["conda_envs"]["phage_all"]
    output:
        raw_predictions = os.path.join(config["outdir"], "{sample}", "phage_analysis", "enhanced_comparison", "raw_predictions_with_quality.tsv"),
        quality_stats = os.path.join(config["outdir"], "{sample}", "phage_analysis", "enhanced_comparison", "quality_statistics.tsv"),
        tool_performance = os.path.join(config["outdir"], "{sample}", "phage_analysis", "enhanced_comparison", "tool_performance.tsv"),
        high_confidence = os.path.join(config["outdir"], "{sample}", "phage_analysis", "enhanced_comparison", "high_confidence_predictions.tsv"),
        plots = os.path.join(config["outdir"], "{sample}", "phage_analysis", "enhanced_comparison", "enhanced_comparison_plots.pdf")
    script:
        "../scripts/compare_tools_with_quality.R"

rule run_prophage_comparison:
    input:
        comparison = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "comparison_plots.pdf")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison_done")
    shell:
        "touch {config[outdir]}/{wildcards.sample}/phage_analysis/comparison_done"

rule run_everything:
    input:
        coverm_stats = os.path.join(config["outdir"], "{sample}", "coverm", "{sample}_stats.txt"),
        checkm = os.path.join(config["outdir"], "{sample}", "binning", "checkm"),
        checkv = os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv"),
        prophage = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage.fasta"),
        comparison = os.path.join(config["outdir"], "{sample}", "phage_analysis", "comparison", "comparison_plots.pdf")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "touch {config[outdir]}/{wildcards.sample}/phage_analysis/done"
