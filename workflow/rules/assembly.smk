rule spades:
    input:
        hr1 = os.path.join(config["outdir"], "preprocessing", "{sample}", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "preprocessing", "{sample}", "{sample}_2_hr.fastq.gz"),
        #hr3 = os.path.join(config["reads"], "host_removed", "{sample}_unpaired_hr.fastq.gz")
    threads: 12
    conda: config["conda_envs"]["spades"]
    output:
        dir = directory(os.path.join(config["outdir"], "assembly", "{sample}")),
        contigs = os.path.join(config["outdir"], "assembly", "{sample}", "contigs.fasta")
    log:
        os.path.join(config["outdir"], "logs", "spades", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "spades", "{sample}_bmrk.txt")
    shell:
        """
        spades.py --meta --only-assembler \
        -m 400 -t {threads} \
        -1 {input.hr1} -2 {input.hr2} \
        -o {output.dir} 2> {log}
        """
