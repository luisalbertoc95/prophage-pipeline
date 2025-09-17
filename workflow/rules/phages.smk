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
        bakta --db {input.db}/db --force --skip-plot --keep-contig-headers --output {output} \
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

# Conditional input function for taxonomy data
def get_taxonomy_input(wildcards):
    if config["taxonomy_method"] == "gtdbtk":
        return os.path.join(config["outdir"], wildcards.sample, "taxonomy", "gtdbtk")
    else:
        return os.path.join(config["outdir"], wildcards.sample, "taxonomy", "mmseqs")

# Conditional input function for phage_all rule to avoid circular dependency
def get_phage_all_taxonomy_input(wildcards):
    # Only include taxonomy input if we're in all_contigs mode
    # In prophage_only mode, we'll run without taxonomy first, then add it later
    if config.get("taxonomy_scope", "prophage_only") == "all_contigs":
        return get_taxonomy_input(wildcards)
    else:
        # Return empty - we'll add taxonomy in a separate step
        return []

def get_phage_all_input(wildcards):
    inputs = {
        "genomad": os.path.join(config["outdir"], wildcards.sample, "phage_analysis", "genomad"),
        "phispy": os.path.join(config["outdir"], wildcards.sample, "phage_analysis", "phispy"),
        "mmseqs": get_taxonomy_input(wildcards)
    }
    
    return inputs

rule prophage_overlap_detection:
    input:
        unpack(get_phage_all_input)
    conda: config["conda_envs"]["phage_all"]
    output:
        merged_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "merged_prophages.bed"),
        phispy_unique_ids = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy_unique_ids.txt")
    log:
        os.path.join(config["outdir"], "logs", "prophage_overlap", "{sample}.log")
    shell:
        """
        # Create BED files from both tools
        # GeNomad: extract contig, start, end from TSV
        awk 'NR>1 {{
            # Extract NODE number from source_seq column (column 2)
            if (match($2, /NODE_([0-9]+)_/, arr)) {{
                print arr[1] "\t" $3 "\t" $4 "\tgenomad\t" NR-1
            }}
        }}' {input.genomad}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.tsv > {config[outdir]}/{wildcards.sample}/phage_analysis/genomad.bed 2> {log}
        
        # PhiSpy: extract contig, start, end from TSV  
        awk 'NR>1 {{
            # Extract NODE number from Contig column (column 2)
            if (match($2, /NODE_([0-9]+)_/, arr)) {{
                # Split prophage number (pp_X) to get index
                split($1, pp_parts, "_")
                print arr[1] "\t" $3 "\t" $4 "\tphispy\t" pp_parts[2]
            }}
        }}' {input.phispy}/prophage.tsv > {config[outdir]}/{wildcards.sample}/phage_analysis/phispy.bed 2>> {log}
        
        # Find PhiSpy predictions that DON'T overlap with geNomad (use -v flag)
        bedtools intersect -a {config[outdir]}/{wildcards.sample}/phage_analysis/phispy.bed \
                          -b {config[outdir]}/{wildcards.sample}/phage_analysis/genomad.bed \
                          -v > {config[outdir]}/{wildcards.sample}/phage_analysis/phispy_unique.bed 2>> {log}
        
        # Create merged BED: all geNomad + unique PhiSpy
        cat {config[outdir]}/{wildcards.sample}/phage_analysis/genomad.bed \
            {config[outdir]}/{wildcards.sample}/phage_analysis/phispy_unique.bed > {output.merged_bed} 2>> {log}
        
        # Extract PhiSpy unique IDs for FASTA extraction
        awk '{{print $5}}' {config[outdir]}/{wildcards.sample}/phage_analysis/phispy_unique.bed > {output.phispy_unique_ids} 2>> {log}
        """

rule create_basic_prophage_table:
    input:
        merged_bed = os.path.join(config["outdir"], "{sample}", "phage_analysis", "merged_prophages.bed"),
        phispy_unique_ids = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy_unique_ids.txt"),
        genomad = lambda wildcards: os.path.join(config["outdir"], wildcards.sample, "phage_analysis", "genomad"),
        phispy = lambda wildcards: os.path.join(config["outdir"], wildcards.sample, "phage_analysis", "phispy")
    conda: config["conda_envs"]["phage_all"]
    output:
        fasta = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unique_phispy_prophage.fasta"),
        table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv")
    log:
        os.path.join(config["outdir"], "logs", "create_basic_prophage_table", "{sample}.log")
    script:
        "../scripts/create_basic_prophage_table.R"

rule add_taxonomy_to_prophage_table:
    input:
        basic_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv"),
        taxonomy = get_taxonomy_input
    conda: config["conda_envs"]["phage_all"]
    output:
        table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    log:
        os.path.join(config["outdir"], "logs", "add_taxonomy_to_prophage_table", "{sample}.log")
    script:
        "../scripts/add_taxonomy_to_prophage_table.R"

rule final_prophage_output:
    input:
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        unique_phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "unique_phispy_prophage.fasta"),
        prophage_table = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table.tsv"),
        prophage_table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv"),
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
    output:
        final_prophage = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage.fasta"),
        contigs_with_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "contigs_with_prophages.fasta")
    shell:
        """
        # Create final prophage sequences (extracted prophages only)
        cp {input.genomad}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.fna \
        {config[outdir]}/{wildcards.sample}/phage_analysis/genomad_prophage.fasta

        cat {input.genomad}/final_filtered_contigs_find_proviruses/final_filtered_contigs_provirus.fna \
        {input.unique_phispy} > {output.final_prophage}

        # Create full contigs containing prophages
        awk 'NR>1 {{print "NODE_" $1 "_"}}' {input.prophage_table} | sort -u > {config[outdir]}/{wildcards.sample}/phage_analysis/prophage_contigs.txt
        
        seqkit grep -r -f {config[outdir]}/{wildcards.sample}/phage_analysis/prophage_contigs.txt \
        {input.contigs} > {output.contigs_with_prophages}
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

rule run_everything:
    input:
        coverm_stats = os.path.join(config["outdir"], "{sample}", "coverm", "{sample}_stats.txt"),
        checkm = os.path.join(config["outdir"], "{sample}", "binning", "checkm"),
        checkv = os.path.join(config["outdir"], "{sample}", "phage_analysis", "checkv"),
        prophage = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage.fasta"),
        contigs_with_prophages = os.path.join(config["outdir"], "{sample}", "phage_analysis", "contigs_with_prophages.fasta"),
        prophage_table_with_taxonomy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "final_prophage_table_with_host_taxonomy.tsv")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "touch {config[outdir]}/{wildcards.sample}/phage_analysis/done"
