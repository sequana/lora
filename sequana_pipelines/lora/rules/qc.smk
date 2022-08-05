"""Rules to assess the quality of the assembly"""
from sequana_pipelines.lora import BLAST_KEY


rule seqkit_sort:
    input:
        get_final_contigs
    output:
        "{sample}/sorted_contigs/{sample}.fasta"
    threads:
        config["seqkit_sort"]["threads"]
    resources:
        **config["seqkit_sort"]["resources"],
    shell:
        """
        seqkit sort --threads {threads} --by-length --reverse {input} -o {output}
        """


rule minimap2_and_genomecov:
    input:
        fastq = get_fastq,
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        bam = "{sample}/minimap2/{sample}.bam",
        bed = "{sample}/minimap2/{sample}.bed"
    params:
        preset = config['minimap2']['preset'],
        options = config['minimap2']['options']
    threads:
        config['minimap2']['threads']
    resources:
        **config["minimap2"]["resources"],
    shell:
        """
        minimap2 {params.options} -t {threads} -ax {params.preset} {input.contigs} {input.fastq}\
            | samtools sort -@ $(({threads} - 1)) -o {output.bam}\
            && samtools index {output.bam}\
            && samtools depth -aa {output.bam} > {output.bed}
        """


rule sequana_coverage:
    input:
        bed = "{sample}/minimap2/{sample}.bed"
    output:
        html = "{sample}/sequana_coverage/multiqc_report.html"
    params:
        config['sequana_coverage']['options']
    resources:
        **config["sequana_coverage"]["resources"],
    shell:
        """
        sequana_coverage {params} -i {input.bed} -o --output-directory {wildcards.sample}/sequana_coverage
        """


rule quast:
    input:
        fastq = get_fastq,
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        "{sample}/quast/quast.done"
    params:
        preset = config['quast']['preset'],
        options = config['quast']['options']
    threads:
        config['quast']['threads']
    resources:
        **config["quast"]["resources"],
    shell:
        """
        quast.py {params.options} -t {threads} {input.contigs} --{params.preset} {input.fastq} -o {wildcards.sample}/quast\
            && touch {output}
        """


rule busco:
    input:
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        directory("{sample}/busco")
    log:
        "{sample}/logs/{sample}_busco.out"
    params:
        mode = "genome",
        lineage = config['busco']['lineage'],
        short_summary_filename = "short_summary_{sample}.txt",
        options = config['busco']['options']
    threads:
        config['busco']['threads']
    resources:
        **config["busco"]["resources"],
    wrapper:
        "main/wrappers/busco"


rule prokka:
    input:
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        "{sample}/prokka/{sample}.gbk"
    params:
        config['prokka']['options']
    threads:
        config['prokka']['threads']
    resources:
        **config["prokka"]["resources"],
    shell:
        """
        prokka {params} --force --cpus {threads} --outdir {wildcards.sample}/prokka --prefix {wildcards.sample} {input.contigs}
        """


rule seqkit_head:
    input:
        "{sample}/sorted_contigs/{sample}.fasta"
    output:
        "{sample}/subset_contigs/{sample}.subset.fasta"
    params:
        n_first = config["seqkit_head"]["n_first"]
    shell:
        """
        seqkit head -n {params.n_first} -o {output} {input}
        """


rule blast:
    input:
        contigs = "{sample}/subset_contigs/{sample}.subset.fasta"
    output:
        "{sample}/blast/{sample}.tsv"
    params:
        db = config['blast']['blastdb'],
        evalue = config['blast']['evalue'],
        outfmt = " ".join(BLAST_KEY),
        options = config['blast']['options']
    threads:
        config['blast']['threads']
    resources:
        **config["blast"]["resources"],
    shell:
        """
        export BLASTDB={params.db}
        blastn -query {input.contigs} -db nt -evalue {params.evalue} -out {output} -num_threads {threads} \
            {params.options} -outfmt "6 {params.outfmt}"
        """


rule multiqc:
    input:
        requested_output(manager)
    output:
        "multiqc/multiqc_report.html"
    shell:
        """
        multiqc -f . --outdir multiqc -m busco -m quast
        """
