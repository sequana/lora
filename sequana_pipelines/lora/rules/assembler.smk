"""Assembler rules"""


rule pecat_config:
    input:
        fastq=get_fastq
    output:
        config="{sample}/pecat/{sample}.txt",
        read_list="{sample}/pecat/fastq_list.txt"

    params:
        genome_size=config["pecat"]["genome_size"],
        min_required_length=config["pecat"]["min_required_length"],
    threads:
         config['pecat']['threads']
    container:
        config["apptainers"]["pecat"]
    shell:
        """
        pecat.pl config {output.config}

    	echo {input.fastq} > {output.read_list}

        # edit the config file
        sed -i -e "s/^project=/project={wildcards.sample}\/pecat/g" {output.config}
        sed -i -e "s/^genome_size=/genome_size={params.genome_size}/g" {output.config}
        sed -i -e "s/^reads=/reads={wildcards.sample}\/pecat\/fastq_list.txt/g" {output.config}
        sed -i -e "s/^prep_min_length=3000/prep_min_length={params.min_required_length}/g" {output.config}
        sed -i -e "s/^threads=4/threads={threads}/g" {output.config}
        # for haplotypes:
        sed -i -e "s/^asm2_assemble_options=/asm2_assemble_options= --max_trivial_length 10000 --contig_format dual,prialt/g" {output.config}

        """

rule pecat_correct:
    input:
        fastq=get_fastq,
        config=rules.pecat_config.output.config
    output:
        fasta="{sample}/pecat/1-correct/corrected_reads.fasta"
    log:
        "{sample}/logs/pecat.out"
    threads:
         config['pecat']['threads']
    resources:
        **config["pecat"]["resources"],
    container:
        config["apptainers"]["pecat"]
    shell:
        """
        pecat.pl correct {input.config} 2>&1 1>{log}
        """


rule pecat_assemble:
    input:
        config=rules.pecat_config.output.config,
        fasta=rules.pecat_correct.output.fasta
    output:
        primary="{sample}/pecat/3-assemble/primary.fasta",
        alternate="{sample}/pecat/3-assemble/alternate.fasta",
    log:
        "{sample}/logs/pecat.out"
    threads:
         config['pecat']['threads']
    resources:
        **config["pecat"]["resources"],
    container:
        config["apptainers"]["pecat"]
    shell:
        """
        pecat.pl assemble {input.config} 2>&1 1>{log}
        """

rule pecat_unzip:
    input:
        config=rules.pecat_config.output.config,
        fasta=rules.pecat_assemble.output.primary
    output:
        primary="{sample}/pecat/5-assemble/primary.fasta",
        alternate="{sample}/pecat/5-assemble/alternate.fasta",
        polished_primary="{sample}/pecat/6-polish/racon/primary.fasta",
        polished_alternate="{sample}/pecat/6-polish/racon/alternate.fasta",
    log:
        "{sample}/logs/pecat.out"
    threads:
         config['pecat']['threads']
    resources:
        **config["pecat"]["resources"],
    container:
        config["apptainers"]["pecat"]
    shell:
        """
        pecat.pl unzip {input.config} 2>&1 1>{log}
        """




    

rule necat_config:
    input:
        fastq=get_fastq,
    output:
        config="{sample}/necat/{sample}.txt",
        read_list="{sample}/necat/fastq_list.txt"
    params:
        genome_size=config["necat"]["genome_size"],
        min_required_length=config["necat"]["min_required_length"],
    threads:
         config['necat']['threads']
    shell:
        """
         # create config 
        necat.pl config {output.config}

	echo {input.fastq} >> {output.read_list}

        # edit the config file
        sed -i -e "s/^PROJECT=/PROJECT={wildcards.sample}\/necat/g" {output.config} 
        sed -i -e "s/^GENOME_SIZE=/GENOME_SIZE={params.genome_size}/g" {output.config} 
	# we must reuse wildcards (not output.read_list) because we need to escape the \ signs
        sed -i -e "s/^ONT_READ_LIST=/ONT_READ_LIST={wildcards.sample}\/necat\/fastq_list.txt/g" {output.config} 
        sed -i -e "s/^MIN_READ_LENGTH=3000/MIN_READ_LENGTH={params.min_required_length}/g" {output.config} 
        sed -i -e "s/^THREADS=4/THREADS={threads}/g" {output.config} 
        """


rule necat_correct:
    input:
        fastq=get_fastq,
        config=rules.necat_config.output.config
    output:	
        cns ="{sample}/necat/1-consensus/cns_final.fasta.gz"
    log:
        "{sample}/logs/necat.out"
    threads:
         config['necat']['threads']
    resources:
        **config["necat"]["resources"],
    #container:
    #    config["apptainers"]["necat"] 
    shell:
        """
        necat.pl correct {input.config} 2>&1 1>{log}
        """


rule necat_assemble:
    input:
        cns = rules.necat_correct.output.cns,
        config=rules.necat_config.output.config
    output:
        fasta="{sample}/necat/4-fsa/polished_contigs.fasta"
    log:
        "{sample}/logs/necat.out"
    threads:
         config['necat']['threads']
    resources:
        **config["necat"]["resources"],
    #container:
    #    config["apptainers"]["necat"] 
    shell:
        """
        necat.pl assemble {input.config} 2>&1 1>{log}
        """


rule necat_bridge:
    input:
        fasta=rules.necat_assemble.output.fasta,
        config=rules.necat_config.output.config
    output:
        fasta="{sample}/necat/6-bridge_contigs/polished_contigs.fasta"
    log:
        "{sample}/logs/necat.out"
    threads:
         config['necat']['threads']
    resources:
        **config["necat"]["resources"],
    #container:
    #    config["apptainers"]["necat"] 
    shell:
        """
        necat.pl bridge {input.config} 2>&1 1>{log}
        """


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
        "{sample}/logs/unicycler.out",
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
    log:
        "{sample}/logs/flye.out",
    container:
        config['apptainers']['flye']
    threads:
        config["flye"]["threads"]
    resources:
        **config["flye"]["resources"],
    shell:
        """
        outdir="$(dirname "{output.contig}")"

        flye --genome-size {params.genome_size} {params.options} --{params.preset} {input} --out-dir ${{outdir}} --threads {threads} 2>&1 1>{log}    && mv ${{outdir}}/assembly.fasta {output.contig}
        """


def get_bandage_input(wildcards):
    if config["assembler"] == "flye":
        return rules.flye.output.gfa
    elif config["assembler"] == "unicycler":
        return rules.unicycler.output.gfa
    elif config["assembler"] == "canu":
        return rules.paf2gfa.output.gfa



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











