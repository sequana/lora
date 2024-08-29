"""Assembler rules"""


rule unicycler:
    input:
        get_fastq,
    output:
        fasta="{sample}/unicycler/{sample}.contigs.fasta",
        gfa="{sample}/unicycler/assembly.gfa",
    params:
        mode=config["unicycler"]["mode"],
        options=config["unicycler"]["options"],
        long_reads=True
    log:
        "{sample}/logs/unicycler.log",
    threads:
        config["unicycler"]["threads"]
    container:
        config["apptainers"]["unicycler"]
    resources:
        **config["unicycler"]["resources"],
    wrapper:
        f"{manager.wrappers}/wrappers/unicycler"


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
    threads:
        config["canu"]["threads"]
    resources:
        **config["canu"]["resources"],
    container:
        config['apptainers']['canu']
    wrapper:
        f"{manager.wrappers}/wrappers/canu"


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
    threads:
        config["canu_correction"]["threads"]
    resources:
        **config["canu_correction"]["resources"],
    container:
        config['apptainers']['canu']
    wrapper:
        f"{manager.wrappers}/wrappers/canu"


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
    threads:
        config["canu_correction"]["threads"]
    resources:
        **config["canu_correction"]["resources"],
    container:
        config['apptainers']['canu']
    wrapper:
        f"{manager.wrappers}/wrappers/canu"


rule fasta_to_fastq:
    input:
        "{sample}/corrected_reads/{sample}.trimmedReads.fasta.gz",
    output:
        "{sample}/corrected_reads/{sample}.fastq",
    container:
        config['apptainers']['seqtk']
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
    container:
        config['apptainers']['hifiasm']
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
        gfa="{sample}/flye/assembly_graph.gfa",
    params:
        preset=config["flye"]["preset"],
        options=config["flye"]["options"],
        genome_size=config["flye"]["genome_size"],
    container:
        config['apptainers']['flye']
    threads:
        config["flye"]["threads"]
    resources:
        **config["flye"]["resources"],
    shell:
        """
        outdir="$(dirname "{output.contig}")"

        flye --genome-size {params.genome_size} {params.options} --{params.preset} {input} --out-dir ${{outdir}} --threads {threads} \
            && mv ${{outdir}}/assembly.fasta {output.contig}
        """


def get_bandage_input(wildcards):
    if config["assembler"] == "flye":
        return rules.flye.output.gfa
    elif config["assembler"] == "unicycler":
        return rules.unicycler.output.gfa


rule bandage:
    input:
        gfa=get_bandage_input
    output:
        "{sample}/bandage/{sample}_graph.png"
    container:
        config["apptainers"]["bandage"]
    shell:
        """
        Bandage image {input.gfa} {output}
        """











