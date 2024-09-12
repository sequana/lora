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
        f"{manager.wrappers}/wrappers/medaka/consensus"


rule circlator:
    # if circlator requires data_type pacbio-corrected, uses ouput of canu 
    # otherwise, uses input fastq
    input:
        contig = get_medaka_consensus_contigs,
        fastq = get_corrected_fastq
    output:
        "{sample}/circlator/{sample}.circle.fasta"
    params:
        options = config['circlator']['options'],
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


rule polypolish_index:
    input:
        contigs=get_circlator_contigs
    output:
        index="{sample}/logs/{sample}_indexing_polishing.done"
    container:
        config['apptainers']['bwa']
    resources:
        **config["polypolish"]["resources"],
    shell:
        """
        bwa index {input.contigs} && touch {output.index}
        """


rule polypolish_R1_mapping:
    input:
        contigs=get_circlator_contigs,
        illumina=get_illumina_data,
        index="{sample}/logs/{sample}_indexing_polishing.done"
    output:
        aln="{sample}/preprocess_for_polypolish/{sample}.1.sam",
    resources:
        **config["polypolish"]["resources"],
    container:
        config['apptainers']['bwa']
    shell:
        """
        bwa mem -a {input.contigs} {input.illumina[0]} > {output.aln}
        """


if Ifact.paired:

    rule polypolish_R2_mapping:
        input:
            contigs=get_circlator_contigs,
            illumina=get_illumina_data,
            index="{sample}/logs/{sample}_indexing_polishing.done"
        output:
            aln="{sample}/preprocess_for_polypolish/{sample}.2.sam"
        resources:
            **config["polypolish"]["resources"],
        container:
            config['apptainers']['bwa']
        shell:
            """
            bwa mem -a {input.contigs} {input.illumina[1]} > {output.aln}
            """

    rule polypolish_filter:
        input:
            contigs=get_circlator_contigs,
            aln1="{sample}/preprocess_for_polypolish/{sample}.1.sam",
            aln2="{sample}/preprocess_for_polypolish/{sample}.2.sam"
        output:
            filter1="{sample}/preprocess_for_polypolish/{sample}.filter.1.sam",
            filter2="{sample}/preprocess_for_polypolish/{sample}.filter.2.sam"
        resources:
            **config["polypolish_filter"]["resources"],
        container:
            config['apptainers']['polypolish']
        shell:
            """
            polypolish filter \
               --in1 {input.aln1} \
               --in2 {input.aln2} \
               --out1 {output.filter1} \
               --out2 {output.filter2}
            """


def get_polypolish_alignment_files(wildcards):
    if Ifact.paired:
       return (
           f"{wildcards.sample}/preprocess_for_polypolish/{wildcards.sample}.filter.1.sam",
           f"{wildcards.sample}/preprocess_for_polypolish/{wildcards.sample}.filter.2.sam"
       )
    else:
       return (
           f"{wildcards.sample}/preprocess_for_polypolish/{wildcards.sample}.1.sam",
       )


rule polypolish:
    input:
        alignments=get_polypolish_alignment_files,
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
        f"{manager.wrappers}/wrappers/polypolish"



