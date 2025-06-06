# Summary and visualization rules for prophage pipeline results

rule create_summary_plots:
    """
    Generate summary plots using Python visualization script
    """
    input:
        expand(os.path.join(config["outdir"], "phage_analysis", "{sample}", "done"), sample=SAMPLES)
    conda: config["conda_envs"]["phage_all"]
    output:
        plots_dir = directory(os.path.join(config["outdir"], "summary_plots")),
        sample_summary = os.path.join(config["outdir"], "prophage_summary_by_sample.tsv"),
        tool_summary = os.path.join(config["outdir"], "tool_detection_summary.tsv")
    log:
        os.path.join(config["outdir"], "logs", "summary_plots.log")
    shell:
        """
        python {workflow.basedir}/scripts/create_summary_plots.py {config[outdir]} > {log} 2>&1
        """

rule create_summary_report:
    """
    Generate comprehensive HTML report using R script
    """
    input:
        expand(os.path.join(config["outdir"], "phage_analysis", "{sample}", "done"), sample=SAMPLES)
    conda: config["conda_envs"]["phage_all"]
    output:
        report = os.path.join(config["outdir"], "prophage_summary_report.html"),
        plots_dir = directory(os.path.join(config["outdir"], "summary_plots_R")),
        sample_summary = os.path.join(config["outdir"], "prophage_summary_by_sample_R.tsv"),
        tool_summary = os.path.join(config["outdir"], "tool_detection_summary_R.tsv")
    log:
        os.path.join(config["outdir"], "logs", "summary_report.log")
    script:
        "../scripts/create_summary_report.R"

rule all_summaries:
    """
    Generate all summary outputs (both Python plots and R report)
    """
    input:
        # Individual sample completion markers
        expand(os.path.join(config["outdir"], "phage_analysis", "{sample}", "done"), sample=SAMPLES),
        # Python visualization outputs
        os.path.join(config["outdir"], "summary_plots"),
        os.path.join(config["outdir"], "prophage_summary_by_sample.tsv"),
        os.path.join(config["outdir"], "tool_detection_summary.tsv"),
        # R report outputs
        os.path.join(config["outdir"], "prophage_summary_report.html"),
        os.path.join(config["outdir"], "summary_plots_R")