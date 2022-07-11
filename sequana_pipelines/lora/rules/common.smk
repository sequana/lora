"""Functions that define the pipeline path"""
from sequana_pipelines.lora import exceptions


assembler_output = {
    "canu": "{sample}/canu/{sample}.contigs.fasta",
    "hifiasm": "{sample}/hifiasm/{sample}.contigs.fasta",
    "flye": "{sample}/flye/{sample}.contigs.fasta"
}

polishing_output = {
    "medaka": "{sample}/medaka/{sample}.polish.fasta",
    "polypolish": "{sample}/polypolish/{sample}.polish.fasta"
}


def requested_output(manager):
    """Resolve all needed output knowing the user config."""
    output_list = [expand("{sample}/quast/quast.done", sample=manager.samples)]
    if config["sequana_coverage"]["do"]:
        output_list += [expand("{sample}/sequana_coverage/multiqc_report.html", sample=manager.samples)]
    if config["busco"]["do"]:
        output_list += [expand("{sample}/busco/", sample=manager.samples)]
    if config["prokka"]["do"]:
        output_list += [expand("{sample}/prokka/{sample}.gbk", sample=manager.samples)]
    if config["blast"]["do"]:
        output_list += [expand("{sample}/blast/{sample}.tsv", sample=manager.samples)]
    if config["circlator"]["do"]:
        output_list += [expand("{sample}/circlator/{sample}.contigs.fasta", sample=manager.samples)]
    return output_list


def get_raw_data(wildcards):
    """Get raw data."""
    return manager.samples[wildcards.sample]


def aggregate_flowcell(wildcards):
    checkpoint_output = checkpoints.index_bam.get(**wildcards).output[0]
    return expand(
        "{sample}/ccs/{sample}_{flowcell}.ccs.bam",
        sample=wildcards.sample,
        flowcell=glob_wildcards(os.path.join(checkpoint_output, "{sample}_{flowcell}.bam")).flowcell,
    )


def get_bam(wildcards):
    """Use bam from ccs generation or raw bam."""
    if config["ccs"]["do"]:
        return f"{wildcards.sample}/ccs/{wildcards.sample}.ccs.bam"
    return manager.samples[wildcards.sample][0]


def get_fastq(wildcards):
    """Convert bam to fastq format."""
    filename = manager.samples[wildcards.sample][0]
    if filename.endswith(".bam"):
        return f"{wildcards.sample}/bam_to_fastq/{wildcards.sample}.fastq"
    return filename


def get_corrected_fastq(wildcards):
    """Use corrected fastq or not"""
    if config["canu_correction"]["do"]:
        return f"{wildcards.sample}/corrected_reads/{wildcards.sample}.fastq"
    return get_fastq(wildcards)


def get_assembler_contigs(wildcards):
    """Get contigs with implemented assembler."""
    return assembler_output[config["assembler"]].format(sample=wildcards.sample)


def get_polished_contigs(wildcard):
    """Get contigs with implemented polisher."""
    if config["polishing"]["do"]:
        return polishing_output[config["polishing"]["tool"]].format(sample=wildcards.sample)
    return get_assembler_contigs(wildcards)


def get_final_contigs(wildcards):
    """Get contigs did with assembler or circlator"""
    if config["circlator"]["do"]:
        return f"{wildcards.sample}/circlator/{wildcards.sample}.contigs.fasta"
    return get_polished_contigs(wildcards)
