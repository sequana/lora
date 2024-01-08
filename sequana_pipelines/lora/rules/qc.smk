"""Rules to assess the quality of the assembly"""
from sequana_pipelines.lora.src.enums import BLAST_KEY


rule seqkit_sort:
    input:
        get_final_contigs
    output:
        "{sample}/sorted_contigs/{sample}.fasta"
    threads:
        config["seqkit_sort"]["threads"]
    resources:
        **config["seqkit_sort"]["resources"],
    container:
        config['apptainers']['seqkit']
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
    container:
        config['apptainers']['minimap2']
    resources:
        **config["minimap2"]["resources"],
    shell:
        """
        minimap2 {params.options} -t {threads} -ax {params.preset} {input.contigs} {input.fastq} \
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
    container:
        config['apptainers']['sequana_coverage']
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
    container:
        config['apptainers']['quast']
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
    container:
        config['apptainers']['busco']
    resources:
        **config["busco"]["resources"],
    wrapper:
        f"{manager.wrappers}/wrappers/busco"


rule prokka:
    input:
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        "{sample}/prokka/{sample}.gbk"
    params:
        config['prokka']['options']
    threads:
        config['prokka']['threads']
    container:
        config['apptainers']['prokka']
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
    container:
        config['apptainers']['seqkit']
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
    container:
        config['apptainers']['blast']
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
    params:
        options=config["multiqc"]["options"],
        input_directory=config['multiqc']['input_directory'],
        config_file=config['multiqc']['config_file'],
        modules=config['multiqc']['modules']
    log:
        "multiqc/multiqc.log"
    resources:
        **config["multiqc"]["resources"]
    container:
        config["apptainers"]["multiqc"]
    wrapper:
       f"{manager.wrappers}/wrappers/multiqc"


rule checkm_marker:
    output:
        marker="marker"
    threads:
         1
    container:
        config['apptainers']['checkm']
    params:
        taxon_rank=config["checkm"]["taxon_rank"],
        taxon_name=config["checkm"]["taxon_name"]
    resources:
        **config["checkm"]["resources"],
    shell:
        """
        checkm taxon_set {params.taxon_rank} {params.taxon_name} {output.marker}
        """

rule checkm:
    input:
        fasta="{sample}/sorted_contigs/{sample}.fasta",
        marker=rules.checkm_marker.output.marker
    output:
        results="{sample}/checkm/results.txt",
        png="{sample}/checkm/{sample}.marker_pos_plot.png"
    params:
        taxon_rank=config["checkm"]["taxon_rank"],
        taxon_name=config["checkm"]["taxon_name"]
    threads:
        config['checkm']['threads']
    container:
        config['apptainers']['checkm']
    resources:
        **config["checkm"]["resources"],
    shell:
        """
        checkm analyze {input.marker} {wildcards.sample}/sorted_contigs/ {wildcards.sample}/checkm -x fasta -t {threads}
        checkm qa marker {wildcards.sample}/checkm/ -f {output.results}
        checkm marker_plot {wildcards.sample}/checkm/ {wildcards.sample}/sorted_contigs {wildcards.sample}/checkm/ -x fasta
        """

