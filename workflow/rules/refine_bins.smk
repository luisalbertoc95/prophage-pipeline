rule dastool:
    input: 
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        concoct_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "concoct.contigs2bin.tsv"),
        maxbin_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "maxbin.contigs2bin.tsv"),
        metabat_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "metabat.contigs2bin.tsv")
    threads: 24
    conda: config["conda_envs"]["dastool"]
    output:
        os.path.join(config["outdir"], "{sample}", "binning", "done")
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
        -o {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample} \
        2> {log} || true
        
        """

rule filter_unbinned:
    input:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        dastool = os.path.join(config["outdir"], "{sample}", "binning", "done")
    conda: config["conda_envs"]["minimap"] 
    output:
        final_contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        unbinned_4000bp = os.path.join(config["outdir"], "{sample}", "binning","filt_4000_seqs_to_keep.fasta"),
        contigs_5000bp = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta")
    shell:
        """
        BINNING={config[outdir]}/{wildcards.sample}/binning
        BINS_DIR=$BINNING/dastool/{wildcards.sample}_DASTool_bins
        UNB=$BINS_DIR/unbinned.fa

        # ROBUSTNESS (fork): DAS_Tool omits unbinned.fa for low-biomass samples (no bins
        # written, or no contigs left unbinned), which crashes this rule under bash strict
        # mode. Reconstruct it directly: unbinned = contigs_filt minus contigs already in a
        # DAS_Tool bin. Yields all contigs when 0 bins, empty when all contigs binned.
        if [ ! -f "$UNB" ]; then
            mkdir -p "$BINS_DIR"
            BINFA=$(ls "$BINS_DIR"/*.fa 2>/dev/null | grep -v '/unbinned.fa$' || true)
            if [ -n "$BINFA" ]; then
                grep -h '^>' $BINFA | sed 's/^>//; s/[[:space:]].*//' \
                    > "$BINNING/binned_contig_ids.txt" || : > "$BINNING/binned_contig_ids.txt"
                if [ -s "$BINNING/binned_contig_ids.txt" ]; then
                    seqkit grep -v -f "$BINNING/binned_contig_ids.txt" {input.contigs_filt} > "$UNB" || : > "$UNB"
                else
                    cp {input.contigs_filt} "$UNB"
                fi
            else
                cp {input.contigs_filt} "$UNB"
            fi
        fi

        # Extract the unbinned sequences >=4000bp (empty-safe)
        seqkit seq -m 4000 "$UNB" > {output.unbinned_4000bp} || : > {output.unbinned_4000bp}

        # Extract the IDs of the unbinned sequences <4000bp (empty-safe)
        seqkit seq -n -M 3999 "$UNB" > "$BINNING/filt_4000_seqs_to_discard.txt" \
            || : > "$BINNING/filt_4000_seqs_to_discard.txt"

        # From main contigs file, get all sequences except for those on this list.
        # Empty discard list => nothing to remove => keep all contigs.
        if [ -s "$BINNING/filt_4000_seqs_to_discard.txt" ]; then
            seqkit grep -v -f "$BINNING/filt_4000_seqs_to_discard.txt" \
            {input.contigs_filt} -o {output.final_contigs}
        else
            cp {input.contigs_filt} {output.final_contigs}
        fi

        # Relocate unbinned.fa out of the bins dir (best-effort; already ensured above)
        mv "$UNB" $BINNING/dastool/ 2>/dev/null || true

        # Filter contigs for phispy input (5000bp filter, empty-safe)
        seqkit seq -m 5000 {output.final_contigs} > {output.contigs_5000bp} || : > {output.contigs_5000bp}
        """

rule separate_unbinned:
    input: 
        os.path.join(config["outdir"], "{sample}", "binning","filt_4000_seqs_to_keep.fasta")
    threads: 8
    output:
        os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins")
    shell:
        """
        touch {output}
        cd {config[outdir]}/{wildcards.sample}/binning

        cat filt_4000_seqs_to_keep.fasta | awk '
        {{
            if (substr($0, 1, 1) == ">") {{ 
                filename = (substr($0, 2) ".fa") 
            }}
            print $0 >> filename
            close(filename)
        }}'

        # Empty keep set (low-biomass sample) => no NODE files to move; don't fail.
        mv {wildcards.sample}_NODE* dastool/{wildcards.sample}_DASTool_bins 2>/dev/null || true
        cd ../../../..
        """

rule checkm:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins")
    threads: 24
    conda: config["conda_envs"]["checkm"]
    output:
        outdir = directory(os.path.join(config["outdir"], "{sample}", "binning", "checkm")),
        tsv = os.path.join(config["outdir"], "{sample}", "binning", "checkm", "checkm_out.tsv")
    log:
        os.path.join(config["outdir"], "logs", "checkm", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkm", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output.outdir}
        BINS_DIR={config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins
        # ROBUSTNESS (fork): checkm errors when there are no bins (0-MAG low-biomass
        # samples). Emit a header-only checkm_out.tsv so the sample completes with no MAGs.
        if ls "$BINS_DIR"/*.fa >/dev/null 2>&1; then
            checkm lineage_wf -x fa \
            "$BINS_DIR"/ \
            {output.outdir}/ -t {threads} --tab_table -f {output.tsv} 2> {log}
        else
            printf 'Bin Id\\tMarker lineage\\t# genomes\\t# markers\\t# marker sets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain heterogeneity\\n' > {output.tsv}
            echo "No DAS_Tool bins; wrote header-only checkm_out.tsv" > {log}
        fi
        """
