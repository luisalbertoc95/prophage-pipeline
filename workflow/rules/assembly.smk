rule megahit:
    input:
        hr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_hr.fastq.gz"),
    threads: 24
    conda: config["conda_envs"]["megahit"]
    output:
        dir = directory(os.path.join(config["outdir"], "{sample}", "assembly")),
        contigs = os.path.join(config["outdir"], "{sample}", "assembly", "contigs.fasta")
    log:
        os.path.join(config["outdir"], "logs", "megahit", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "megahit", "{sample}_bmrk.txt")
    shell:
        """
        # Run MEGAHIT assembly
        megahit -1 {input.hr1} -2 {input.hr2} \
        -t {threads} --min-contig-len 500 \
        --k-min 21 --k-max 141 --k-step 12 \
        --force -o {output.dir} 2> {log}
        
        # Convert MEGAHIT headers to SPAdes format for compatibility
        python3 -c "
import re
with open('{output.dir}/final.contigs.fa', 'r') as infile, open('{output.contigs}', 'w') as outfile:
    node_counter = 1
    for line in infile:
        if line.startswith('>'):
            # Extract length from MEGAHIT header: '>k141_1234 flag=1 multi=2.0000 len=4567'
            length_match = re.search(r'len=(\\d+)', line)
            cov_match = re.search(r'multi=([\\d.]+)', line)
            length = length_match.group(1) if length_match else '500'
            cov = cov_match.group(1) if cov_match else '1.0'
            # Create SPAdes-compatible header
            new_header = f'>NODE_{{node_counter}}_length_{{length}}_cov_{{cov}}\\n'
            outfile.write(new_header)
            node_counter += 1
        else:
            outfile.write(line)
" 2>> {log}
        """
