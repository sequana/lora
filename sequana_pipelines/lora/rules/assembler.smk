

rule canu:
    input:
        get_fastq
    output:
        contig = "{sample}/canu/{sample}.contigs.fasta",
        done = "{sample}/canu/canu.done"
    params:
        genome_size = config['canu']['genome_size'],
        use_grid = config['canu']['use_grid'],
        options = config['canu']['options']
    threads:
        config['canu']['threads']
    wrapper:
        "main/wrappers/canu"

rule canu_correction:
    input:
        get_fastq
    output:
        reads = "{sample}/corrected_reads/{sample}.correctedReads.fasta.gz",
        done = "{sample}/canu/canu.done"
    params:
        step = '-correct',
        genome_size = config['canu_correction']['genome_size'],
        use_grid = config['canu_correction']['use_grid'],
        options = config['canu_correction']['correction_options']
    threads:
        config['canu_correction']['threads']
    wrapper:
        "main/wrappers/canu"


rule canu_trimming:
    input:
        "{sample}/corrected_reads/{sample}.correctedReads.fasta.gz"
    output:
        "{sample}/corrected_reads/{sample}.trimmedReads.fasta.gz"
    params:
        step = '-trim',
        genome_size = config['canu_correction']['genome_size'],
        use_grid = config['canu_correction']['use_grid'],
        options = "-corrected " + config['canu_correction']['trimming_options'],
    threads:
        config['canu_correction']['threads']
    wrapper:
        "main/wrappers/canu"


rule fasta_to_fastq:
    input:
        "{sample}/corrected_reads/{sample}.trimmedReads.fasta.gz"
    output:
        "{sample}/corrected_reads/{sample}.fastq"
    shell:
        """
        seqtk seq -F '#' {input} > {output}
        """


rule hifiasm:
    input:
        get_hifi_fastq
    output:
        contig = "{sample}/hifiasm/{sample}.contigs.fasta"
    params:
        config['hifiasm']['options']
    threads:
        config['hifiasm']['threads']
    shell:
        """
        mkdir -p {wildcards.sample}/hifiasm
        prefix={wildcards.sample}/hifiasm/{wildcards.sample}
        hifiasm -t{threads} {params} -o $prefix {input}\
            && awk '/^S/{{print ">"$2;print $3}}' ${{prefix}}.bp.p_ctg.gfa > {output}
        """


rule circlator:
    input:
        contig = get_assembler_contigs,
        fastq = get_fastq
    output:
        "{sample}/circlator/{sample}.contigs.fasta"
    params:
        options = config['circlator']['options']
    threads:
        config['circlator']['threads']
    shell:
        """
        outdir={wildcards.sample}/circlator/
        circlator all {params.options} --threads {threads} --data_type pacbio-corrected {input.contig} {input.fastq}\
            ${{outdir}} && mv ${{outdir}}/06.fixstart.fasta {output}
        """