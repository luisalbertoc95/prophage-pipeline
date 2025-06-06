# Summary and visualization rules for prophage pipeline results

rule create_summary_report:
    """
    Generate comprehensive HTML report and visualizations using R script
    """
    input:
        expand(os.path.join(config["outdir"], "phage_analysis", "{sample}", "done"), sample=SAMPLES)
    conda: config["conda_envs"]["phage_all"]
    output:
        report = os.path.join(config["outdir"], "prophage_summary_report.html"),
        plots_dir = directory(os.path.join(config["outdir"], "summary_plots")),
        sample_summary = os.path.join(config["outdir"], "prophage_summary_by_sample.tsv"),
        tool_summary = os.path.join(config["outdir"], "tool_detection_summary.tsv")
    log:
        os.path.join(config["outdir"], "logs", "summary_report.log")
    script:
        "../scripts/create_summary_report.R"