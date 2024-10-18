"""Rules to assess the quality of the assembly"""
from sequana_pipelines.lora.src.enums import BLAST_KEY


rule seqkit_sort:
    input:
        get_final_contigs
    output:
        ctg="{sample}/sorted_contigs/{sample}.fasta",
        names="{sample}/sorted_contigs/{sample}.names.txt"
    threads:
        config["seqkit_sort"]["threads"]
    resources:
        **config["seqkit_sort"]["resources"],
    container:
        config['apptainers']['seqkit']
    shell:
        """
        seqkit sort --threads {threads} --by-length --reverse {input} -o {output.ctg}
        grep ">" {output.ctg} > {output.names}
        """


rule fasta2paf:
    input:
        rules.seqkit_sort.output.ctg
    output:
        paf="{sample}/graph/{sample}.paf"
    container:
        config["apptainers"]["minimap2"]
    shell:
        """
        minimap2 -x ava-pb -t4 {input} {input} > {output}
        """


rule paf2gfa:
    input:
        ctg=rules.seqkit_sort.output.ctg,
        paf=rules.fasta2paf.output.paf
    output:
        gfa="{sample}/graph/{sample}.gfa"
    container:
        config["apptainers"]["miniasm"]
    shell:
        """
        miniasm -f {input.ctg} {input.paf} > {output.gfa}
        """


rule minimap2:
    input:
        fastq = get_fastq,
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        bam = "{sample}/minimap2/{sample}.bam",
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
        """

rule bam2bed:
    input:
        bam=rules.minimap2.output.bam
    output:
        bed="{sample}/minimap2/{sample}.bed"
    container:
        config['apptainers']['mosdepth']
    params:
        mapq=0,
        second_mapq=35
    threads:
        config['bam2bed']['threads']
    resources:
        **config["bam2bed"]["resources"],
    shell:
        """
        mosdepth -t {threads} -Q {params.mapq}  -b1 {wildcards.sample}/minimap2/lenny1 {input.bam}
        mosdepth -t {threads} -Q {params.second_mapq} -b1 {wildcards.sample}/minimap2/lenny2 {input.bam}
        paste <(gunzip -c {wildcards.sample}/minimap2/lenny1.regions.bed.gz | cut -f 1,3,4) <(gunzip -c {wildcards.sample}/minimap2/lenny2.regions.bed.gz | cut -f 4 ) > {output.bed}
        rm -f {wildcards.sample}/minimap2/lenny?.*
        """


rule sequana_coverage:
    input:
        bed=rules.bam2bed.output.bed,
        contigs = "{sample}/sorted_contigs/{sample}.fasta"
    output:
        html="{sample}/sequana_coverage/multiqc_report.html"
    params:
        circular="-o " if config["sequana_coverage"]["circular"] else "",
        chunksize=config["sequana_coverage"]["chunksize"],
        double_threshold=config["sequana_coverage"]["double_threshold"],
        gc_window_size=config["sequana_coverage"]["gc_window_size"],
        high_threshold=config["sequana_coverage"]["high_threshold"],
        low_threshold=config["sequana_coverage"]["low_threshold"],
        mixture_models=config["sequana_coverage"]["mixture_models"],
        options=config["sequana_coverage"]["options"],
        window_size=config["sequana_coverage"]["window_size"],
        output_directory="{sample}/sequana_coverage"
    log:
        "logs/sequana_coverage/{sample}_sequana_coverage.log"
    resources:
        **config["sequana_coverage"]["resources"],
    container:
        config['apptainers']['sequana_coverage']
    shell:
        """
        sequana_coverage --input-file {input.bed} \
            -H {params.high_threshold} \
            -L {params.low_threshold} \
            --clustering-parameter {params.double_threshold} \
            --chunk-size {params.chunksize} \
            --window-gc {params.gc_window_size} \
            --mixture-models {params.mixture_models} \
            --output-directory {params.output_directory} \
            --window-median {params.window_size} \
            --reference-file {input.contigs} \
            {params.options}
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
        "{sample}/logs/busco.out"
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
    log:
        "{sample}/logs/prokka.out"
    container:
        config['apptainers']['prokka']
    resources:
        **config["prokka"]["resources"],
    shell:
        """
        prokka {params} --force --cpus {threads} --outdir {wildcards.sample}/prokka --prefix {wildcards.sample} {input.contigs} 2>&1 >{log}
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
        checkm taxon_set {params.taxon_rank} "{params.taxon_name}" {output.marker}
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

