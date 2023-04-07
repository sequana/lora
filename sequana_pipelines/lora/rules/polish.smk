""" Assembly polishing rules"""


rule medaka_consensus:
    input:
        assembly=get_assembler_contigs,
        fastq=get_fastq
    output:
        "{sample}/medaka_consensus/{sample}.polish.fasta"
    params:
        model = config["medaka_consensus"]["model"],
        options = config["medaka_consensus"]["options"]
    threads:
        config["medaka_consensus"]["threads"]
    container:
        config["apptainers"]["medaka"]
    resources:
        **config["medaka_consensus"]["resources"]
    wrapper:
        "main/wrappers/medaka/consensus"


rule circlator:
    # if circlator requires data_type pacbio-corrected, uses ouput of canu 
    # otherwise, uses input fastq
    input:
        contig = get_medaka_consensus_contigs,
        fastq = get_corrected_fastq
    output:
        "{sample}/circlator/{sample}.circle.fasta"
    params:
        options = config['circlator']['options']
        data_type=config['circlator']['data_type']
    threads:
        config['circlator']['threads']
    resources:
        **config["circlator"]["resources"],
    container:
        config['apptainers']['circlator']
    shell:
        """
        outdir={wildcards.sample}/circlator/
        circlator all {params.options} --threads {threads} --force --data_type {params.data_type} {input.contig} {input.fastq}\
            ${{outdir}} && mv ${{outdir}}/06.fixstart.fasta {output}
        """


rule indexing_for_polypolish:
    input:
        contigs=get_circlator_contigs
    output:
        index="{sample}/logs/{sample}_indexing_polishing.done"
    container:
        config['apptainers']['bwa']
    shell:
        """
        bwa index {input.contigs} && touch {output.index}
        """


rule mapping_R1_for_polypolish:
    input:
        contigs=get_circlator_contigs,
        illumina=get_illumina_data,
        index="{sample}/logs/{sample}_indexing_polishing.done"
    output:
        aln="{sample}/preprocess_for_polypolish/{sample}.1.sam",
    resources:
        **config["polypolish_bwa"]["resources"],
    container:
        config['apptainers']['bwa']
    shell:
        """
        bwa mem -a {input.contigs} {input.illumina[0]} > {output.aln}
        """


rule mapping_R2_for_polypolish:
    input:
        contigs=get_circlator_contigs,
        illumina=get_illumina_data,
        index="{sample}/logs/{sample}_indexing_polishing.done"
    output:
        aln="{sample}/preprocess_for_polypolish/{sample}.2.sam"
    resources:
        **config["polypolish_bwa"]["resources"],
    container:
        config['apptainers']['bwa']
    shell:
        """
        bwa mem -a {input.contigs} {input.illumina[1]} > {output.aln}
        """


rule preprocess_for_polypolish:
    input:
        contigs=get_circlator_contigs,
        aln1="{sample}/preprocess_for_polypolish/{sample}.1.sam",
        aln2="{sample}/preprocess_for_polypolish/{sample}.2.sam"
    output:
        filter1="{sample}/preprocess_for_polypolish/{sample}.filter.1.sam",
        filter2="{sample}/preprocess_for_polypolish/{sample}.filter.2.sam"
    container:
        config['apptainers']['polypolish']
    shell:
        """
        polypolish_insert_filter.py \
           --in1 {input.aln1} \
           --in2 {input.aln2} \
           --out1 {output.filter1} \
           --out2 {output.filter2}
        """


rule polypolish:
    input:
       alignments=("{sample}/preprocess_for_polypolish/{sample}.filter.1.sam",
                   "{sample}/preprocess_for_polypolish/{sample}.filter.2.sam"),
       assembly=get_circlator_contigs
    output:
        fasta="{sample}/polypolish/{sample}.polish.fasta"
    params:
        options=config['polypolish'].get('options', "")
    log:
        "{sample}/logs/{sample}_polypolish.log"
    container:
        config['apptainers']['polypolish']
    resources:
        **config["polypolish"]["resources"],
    wrapper:
        "main/wrappers/polypolish"
