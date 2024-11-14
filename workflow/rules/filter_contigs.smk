rule remove_short_contigs:
    input:
        hr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_hr.fastq.gz"),
        contigs = os.path.join(config["outdir"], "{sample}", "assembly", "contigs.fasta")
    threads: 12
    conda: "../envs/minimap_env.yaml"
    output:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "final_filtered_contigs.fasta"),
        contigs_5000bp = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta")
    log:
        os.path.join(config["outdir"], "logs", "binning_prep", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "binning_prep", "{sample}_bmrk.txt")
    shell:
        """
        cat {input.contigs} | seqkit seq -m 4000 > {output.contigs_filt}

        # Filter contigs for phispy input (5000bp filter)
        cat {output.contigs_filt} | seqkit seq -m 5000 > {output.contigs_5000bp}
        """
