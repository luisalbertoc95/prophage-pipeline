# Summary reporting rules for prophage pipeline

rule generate_summary_report:
    """
    Generate comprehensive HTML summary report aggregating results from all samples.
    Includes statistical analysis and interactive visualizations.
    """
    input:
        # Prophage tables
        prophage_tables = expand(os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv"), sample=SAMPLES),
        basic_prophage_tables = expand(os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv"), sample=SAMPLES),
        mag_tables = expand(os.path.join(config["outdir"], "{sample}", "phage_analysis", "mags", "prophage_table.tsv"), sample=SAMPLES),
        unbinned_tables = expand(os.path.join(config["outdir"], "{sample}", "phage_analysis", "unbinned", "prophage_table.tsv"), sample=SAMPLES),

        # Quality assessments
        checkv_all = expand(os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv_all_prophages", "quality_summary.tsv"), sample=SAMPLES),
        checkv_free = expand(os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv_free_phages", "quality_summary.tsv"), sample=SAMPLES),
        checkm = expand(os.path.join(config["outdir"], "{sample}", "binning", "checkm", "checkm_out.tsv"), sample=SAMPLES),

        # Assembly data
        assemblies = expand(os.path.join(config["outdir"], "{sample}", "assembly", "contigs.fasta"), sample=SAMPLES)
    params:
        samples = SAMPLES,
        outdir = config["outdir"]
    output:
        html = os.path.join(config["outdir"], "pipeline_summary_report.html")
    log:
        os.path.join(config["outdir"], "logs", "summary_report.log")
    conda:
        "../envs/reporting_env.yaml"
    script:
        "../scripts/generate_summary_report.py"
