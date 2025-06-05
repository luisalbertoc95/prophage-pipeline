rule dastool:
    input: 
        contigs_filt = os.path.join(config["outdir"], "assembly", "{sample}", "contigs_filt_1000bp.fasta"),
        concoct_tsv = os.path.join(config["outdir"], "binning", "{sample}", "dastool", "concoct.contigs2bin.tsv"),
        maxbin_tsv = os.path.join(config["outdir"], "binning", "{sample}", "dastool", "maxbin.contigs2bin.tsv"),
        metabat_tsv = os.path.join(config["outdir"], "binning", "{sample}", "dastool", "metabat.contigs2bin.tsv")
    threads: 24
    conda: config["conda_envs"]["dastool"]
    output:
        os.path.join(config["outdir"], "binning", "{sample}", "done")
    log:
        os.path.join(config["outdir"], "logs", "dastool", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "dastool", "{sample}_bmrk.txt")
    shell:
        """
        touch {output}
        
        DAS_Tool --threads {threads} --write_bins --write_unbinned \
        -i {input.concoct_tsv},{input.maxbin_tsv},{input.metabat_tsv} \
        -l concoct,maxbin,metabat -c {input.contigs_filt} \
        -o {config[outdir]}/binning/{wildcards.sample}/dastool/{wildcards.sample} \
        2> {log} || true
        
        """

rule filter_unbinned:
    input:
        contigs_filt = os.path.join(config["outdir"], "assembly", "{sample}", "contigs_filt_1000bp.fasta"),
        dastool = os.path.join(config["outdir"], "binning", "{sample}", "done")
    conda: config["conda_envs"]["minimap"] 
    output:
        final_contigs = os.path.join(config["outdir"], "binning", "{sample}", "final_filtered_contigs.fasta"),
        unbinned_4000bp = os.path.join(config["outdir"], "binning", "{sample}","filt_4000_seqs_to_keep.fasta"),
        contigs_5000bp = os.path.join(config["outdir"], "binning", "{sample}", "final_filt_contigs_5000.fasta")
    shell:
        """
        # Extract the unbinned sequences >=4000bp
        cat {config[outdir]}/binning/{wildcards.sample}/dastool/{wildcards.sample}_DASTool_bins/unbinned.fa \
        | seqkit seq -m 4000 > {output.unbinned_4000bp}

        # Extract the IDs of the unbinned sequences <4000bp
        cat {config[outdir]}/binning/{wildcards.sample}/dastool/{wildcards.sample}_DASTool_bins/unbinned.fa \
        | seqkit seq -n -M 3999 \
        > {config[outdir]}/binning/{wildcards.sample}/filt_4000_seqs_to_discard.txt

        # From main contigs file, get all sequences except for those on this list
        seqkit grep -v -f {config[outdir]}/binning/{wildcards.sample}/filt_4000_seqs_to_discard.txt \
        {input.contigs_filt} -o {output.final_contigs}

        mv {config[outdir]}/binning/{wildcards.sample}/dastool/{wildcards.sample}_DASTool_bins/unbinned.fa \
        {config[outdir]}/binning/{wildcards.sample}/dastool

        # Filter contigs for phispy input (5000bp filter)
        cat {output.final_contigs} | seqkit seq -m 5000 > {output.contigs_5000bp}
        """  

rule separate_unbinned:
    input: 
        os.path.join(config["outdir"], "binning", "{sample}","filt_4000_seqs_to_keep.fasta")
    threads: 8
    output:
        os.path.join(config["outdir"], "binning", "{sample}", "dastool", "{sample}.bins")
    shell:
        """
        touch {output}
        cd {config[outdir]}/binning/{wildcards.sample}

        cat filt_4000_seqs_to_keep.fasta | awk '
        {{
            if (substr($0, 1, 1) == ">") {{ 
                filename = (substr($0, 2) ".fa") 
            }}
            print $0 >> filename
            close(filename)
        }}'

        mv NODE* dastool/{wildcards.sample}_DASTool_bins
        cd ../../../..
        """

rule checkm:
    input:
        os.path.join(config["outdir"], "binning", "{sample}", "dastool", "{sample}.bins")
    threads: 24
    conda: config["conda_envs"]["checkm"]
    output:
        checkm_dir = directory(os.path.join(config["outdir"], "binning", "{sample}", "checkm")),
        checkm_table = os.path.join(config["outdir"], "binning", "{sample}", "checkm", "checkm_out.tsv")
    log:
        os.path.join(config["outdir"], "logs", "checkm", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkm", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output.checkm_dir}
        checkm lineage_wf -x fa \
        {config[outdir]}/binning/{wildcards.sample}/dastool/{wildcards.sample}_DASTool_bins/ \
        {output.checkm_dir}/ -t {threads} --tab_table -f {output.checkm_table} 2> {log}
        """ 
