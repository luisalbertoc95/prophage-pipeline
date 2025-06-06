# Comprehensive summary and visualization rules for prophage pipeline results

rule create_comprehensive_summary:
    """
    Generate comprehensive data warehouse tables, visualizations, and report
    """
    input:
        expand(os.path.join(config["outdir"], "phage_analysis", "{sample}", "done"), sample=SAMPLES)
    conda: config["conda_envs"]["phage_all"]
    output:
        summary_dir = directory(os.path.join(config["outdir"], "summary_results")),
        master_catalog = os.path.join(config["outdir"], "summary_results", "tables", "master_prophage_catalog.tsv"),
        sample_summary = os.path.join(config["outdir"], "summary_results", "tables", "sample_level_summary.tsv"),
        markdown_report = os.path.join(config["outdir"], "summary_results", "SUMMARY_REPORT.md")
    log:
        os.path.join(config["outdir"], "logs", "comprehensive_summary.log")
    script:
        "../scripts/create_comprehensive_summary.R"