"""Assembler rules"""


rule canu:
    input:
        get_fastq,
    output:
        contig="{sample}/canu/{sample}.contigs.fasta",
        done="{sample}/canu/canu.done",
    params:
        preset=config["canu"]["preset"],
        genome_size=config["canu"]["genome_size"],
        use_grid=config["canu"]["use_grid"],
        options=config["canu"]["options"],
    threads: config["canu"]["threads"]
    resources:
        **config["canu"]["resources"],
    wrapper:
        "main/wrappers/canu"


rule canu_correction:
    input:
        get_fastq,
    output:
        reads="{sample}/corrected_reads/{sample}.correctedReads.fasta.gz",
        done="{sample}/corrected_reads/canu-correct.done",
    params:
        step="-correct",
        preset=config["canu_correction"]["preset"],
        genome_size=config["canu_correction"]["genome_size"],
        use_grid=config["canu_correction"]["use_grid"],
        options=config["canu_correction"]["correction_options"],
    threads: config["canu_correction"]["threads"]
    resources:
        **config["canu_correction"]["resources"],
    wrapper:
        "main/wrappers/canu"


rule canu_trimming:
    input:
        "{sample}/corrected_reads/{sample}.correctedReads.fasta.gz",
    output:
        reads="{sample}/corrected_reads/{sample}.trimmedReads.fasta.gz",
        done="{sample}/corrected_reads/canu-trim.done",
    params:
        step="-trim",
        preset=config["canu_correction"]["preset"],
        genome_size=config["canu_correction"]["genome_size"],
        use_grid=config["canu_correction"]["use_grid"],
        options="-corrected " + config["canu_correction"]["trimming_options"],
    threads: config["canu_correction"]["threads"]
    resources:
        **config["canu_correction"]["resources"],
    wrapper:
        "main/wrappers/canu"


rule fasta_to_fastq:
    input:
        "{sample}/corrected_reads/{sample}.trimmedReads.fasta.gz",
    output:
        "{sample}/corrected_reads/{sample}.fastq",
    shell:
        """
        seqtk seq -F '#' {input} > {output}
        """


rule hifiasm:
    input:
        get_corrected_fastq,
    output:
        contig="{sample}/hifiasm/{sample}.contigs.fasta",
    params:
        config["hifiasm"]["options"],
    threads: config["hifiasm"]["threads"]
    resources:
        **config["hifiasm"]["resources"],
    shell:
        """
        mkdir -p {wildcards.sample}/hifiasm
        prefix={wildcards.sample}/hifiasm/{wildcards.sample}
        hifiasm -t{threads} {params} -o $prefix {input}\
            && awk '/^S/{{print ">"$2;print $3}}' ${{prefix}}.bp.p_ctg.gfa > {output}
        """


rule flye:
    input:
        get_corrected_fastq,
    output:
        contig="{sample}/flye/{sample}.contigs.fasta",
    params:
        preset=config["flye"]["preset"],
        options=config["flye"]["options"],
    threads: config["flye"]["threads"]
    resources:
        **config["flye"]["resources"],
    run:
        from os import path

        outdir = path.dirname(output["contig"])
        contigs = path.join(outdir, "assembly.fasta")
        shell(
            "flye {params.options} --{params.preset} {input} --out-dir {outdir} --threads {threads}"
            " && mv {contigs} {output}"
        )
