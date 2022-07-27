""" Nanopore polishing rules"""


rule indexing_for_polypolish:
    input:
        contigs=get_contigs_before_polishing
    output:
        index="{sample}/logs/{sample}_indexing_polishing.done"
    shell:
        """
        bwa index {input.contigs} && touch {output.index}
        """


rule mapping_R1_for_polypolish:
    input:
        contigs=get_contigs_before_polishing,
        illumina=get_illumina_data,
        index="{sample}/logs/{sample}_indexing_polishing.done"
    output:
        aln="{sample}/preprocess_for_polypolish/{sample}.1.sam",
    shell:
        """
        bwa mem -a {input.contigs} {input.illumina[0]} > {output.aln}
        """


rule mapping_R2_for_polypolish:
    input:
        contigs=get_contigs_before_polishing,
        illumina=get_illumina_data,
        index="{sample}/logs/{sample}_indexing_polishing.done"
    output:
        aln="{sample}/preprocess_for_polypolish/{sample}.2.sam"
    shell:
        """
        bwa mem -a {input.contigs} {input.illumina[1]} > {output.aln}
        """


rule preprocess_for_polypolish:
    input:
        contigs=get_contigs_before_polishing,
        aln1="{sample}/preprocess_for_polypolish/{sample}.1.sam",
        aln2="{sample}/preprocess_for_polypolish/{sample}.2.sam"
    output:
        filter1="{sample}/preprocess_for_polypolish/{sample}.filter.1.sam",
        filter2="{sample}/preprocess_for_polypolish/{sample}.filter.2.sam"
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
                   "{sample}/preprocess_for_polypolish/{sample}.filter.2.sam")
       assembly=get_contigs_before_polishing
    output:
        fasta="{sample}/polypolish/{sample}.polish.fasta"
    params:
        options=config['polypolish'].get('options', "")
    log:
        "{sample}/logs/{sample}_polypolish.log"
    wrapper:
        "main/wrappers/polypolish"


