

CCS_MAX_CHUNKS = config['ccs']['max-chunks']

rule index_bam:
    input:
        manager.getrawdata()
    output:
        "{sample}/ccs/{sample}.bam"
    shell:
        """
        mkdir -p {wildcards.sample}/ccs
        ln -sfn {input} {wildcards.sample}/ccs/{wildcards.sample}.bam
        cd {wildcards.sample}/ccs && pbindex {wildcards.sample}.bam; cd -
        """


# intermediate files are needed to parallelize ccs computation ( ͡° ͜ʖ ͡°)
rule setup_ccs_chunk:
    input:
        "{sample}/ccs/{sample}.bam"
    output:
        temp(expand("{{sample}}/ccs/{{sample}}.{chunk}.txt", chunk=range(1, CCS_MAX_CHUNKS + 1)))
    shell:
        """
        touch {output}
        """


rule ccs:
    input:
        "{sample}/ccs/{sample}.{chunk}.txt",
        bam = "{sample}/ccs/{sample}.bam"
    output:
        bam = temp("{sample}/ccs/{sample}.ccs.{chunk}.bam"),
        pbi = temp("{sample}/ccs/{sample}.ccs.{chunk}.bam.pbi"),
        report = "{sample}/ccs/{sample}.{chunk}_report.txt"
    params:
        max_chunks = CCS_MAX_CHUNKS,
        min_rq = config['ccs']['min-rq'],
        min_passes = config['ccs']['min-passes'],
        options = config['ccs']['options']
    threads:
        config['ccs']['threads']
    shell:
        """
        ccs {input.bam} {output.bam} --chunk {wildcards.chunk}/{params.max_chunks} --min-rq {params.min_rq}\
            --min-passes {params.min_passes} --num-threads {threads} --report-file {output.report} {params.options}
        """


rule ccs_merge:
    input:
        expand("{{sample}}/ccs/{{sample}}.ccs.{chunk}.bam", chunk=range(1, CCS_MAX_CHUNKS + 1))
    output:
        "{sample}/ccs/{sample}.ccs.bam"
    params:
        config['samtools_merge']['options']
    threads:
        config['samtools_merge']['threads']
    shell:
        """
        samtools merge -@ $(({threads} - 1)) {params} {output} {input}
        """


rule bam_to_fastq:
    input:
        get_bam
    output:
        fastq = "{sample}/bam_to_fastq/{sample}.fastq"
    params:
        config['bam_to_fastq']['options']
    threads:
        config['bam_to_fastq']['threads']
    shell:
        """
        samtools bam2fq -@ $(({threads} - 1)) {params} {input} > {output.fastq}
        """