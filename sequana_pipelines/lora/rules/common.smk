

def requested_output(manager):
    """ Resolve all needed output knowing the user config.
    """
    output_list = [expand("{sample}/quast/quast.done", sample=manager.samples)]
    if config['circlator']['do']:
        output_list += [expand("{sample}/circlator/{sample}.contigs.fasta", sample=manager.samples)]
    if config['sequana_coverage']['do']:
        output_list += [expand("{sample}/sequana_coverage/multiqc_report.html", sample=manager.samples)]
    if config['busco']['do']:
        output_list += [expand("{sample}/busco/", sample=manager.samples)]
    if config['prokka']['do']:
        output_list += [expand("{sample}/prokka/{sample}.gbk", sample=manager.samples)]
    if config['blast']['do']:
        output_list += [expand("{sample}/blast/{sample}.tsv", sample=manager.samples)]
    return output_list

    
def get_bam(wildcards):
    """ Use bam from ccs generation or raw bam.
    """
    if config['ccs']['do']:
        return f'{wildcards.sample}/ccs/{wildcards.sample}.ccs.bam'
    return manager.samples[wildcards.sample]


def get_fastq(wildcards):
    """ Convert bam to fastq format.
    """
    filename = manager.samples[wildcards.sample]
    if filename.endswith('.bam'):
        return f"{wildcards.sample}/bam_to_fastq/{wildcards.sample}.fastq"
    return filename

def get_hifi_fastq(wildcards):
    """  Use hifi data or corrected reads for hifiasm that needs high quality reads.
    """
    if config['input_hifi']:
        return get_fastq(wildcards)
    return f"{wildcards.sample}/corrected_reads/{wildcards.sample}.fastq"


def get_assembler_contigs(wildcards):
    """ Get contigs with canu or hifiasm.
    """
    if config['hifiasm']['do']:
        return f"{wildcards.sample}/hifiasm/{wildcards.sample}.contigs.fasta"
    return f'{wildcards.sample}/canu/{wildcards.sample}.contigs.fasta'


def get_final_contigs(wildcards):
    """ Get contigs did with assembler or circlator
    """
    if config['circlator']['do']:
        return f'{wildcards.sample}/circlator/{wildcards.sample}.contigs.fasta'
    return get_assembler_contigs(wildcards)