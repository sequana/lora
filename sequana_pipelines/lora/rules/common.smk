"""Functions that define the pipeline path"""
from sequana_pipelines.lora.src import exceptions


assembler_output = {
    "canu": "{sample}/canu/{sample}.contigs.fasta",
    "hifiasm": "{sample}/hifiasm/{sample}.contigs.fasta",
    "flye": "{sample}/flye/{sample}.contigs.fasta",
    "unicycler": "{sample}/unicycler/{sample}.contigs.fasta"
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
        output_list += [expand("{sample}/circlator/{sample}.circle.fasta", sample=manager.samples)]
    if config["polypolish"]["do"]:
        output_list += [expand("{sample}/polypolish/{sample}.polish.fasta", sample=manager.samples)]
    if config["checkm"]["do"]:
        output_list += [expand("{sample}/checkm/{sample}.marker_pos_plot.png", sample=manager.samples)]
    if config["assembler"] in ["flye", "unicycler"]:
        output_list += [expand("{sample}/bandage/{sample}_graph.png", sample=manager.samples)]

    return output_list


def get_raw_data(wildcards):
    """Get raw data."""
    return manager.samples[wildcards.sample]


# utility function to get dictionary of sample names and their input Illumina files

if config["polypolish"]["do"]:
    from sequana_pipetools.snaketools import FastQFactory
    input_pattern = config["polypolish"]["input_pattern"]
    input_readtag = config["polypolish"]["input_readtag"]
    input_directory = config["polypolish"]["input_directory"]
    Ifact = FastQFactory(os.sep.join([input_directory, input_pattern]), read_tag=input_readtag)
else:
    # this is maybe not required. create a dummy struct.
    Ifact = lambda: None
    setattr(illumina_factory, 'paired', False)  #

    
def _get_illumina_data():
    if config["polypolish"]["do"]:
        if Ifact.paired:
            samples = {
                tag: [
                      Ifact.get_file1(tag), 
                      Ifact.get_file2(tag)] for tag in Ifact.tags
            }
        else:
            samples = {
                tag: [Ifact.get_file1(tag)] for tag in Ifact.tags
            }
        return samples
    return {}
illumina_data = _get_illumina_data()


def get_illumina_data(wildcards):
    """Get raw illumina data."""
    if wildcards.sample not in illumina_data:
        logger.error(f"The sample name {wildcards.sample} is discordant (not found in the nanopore sample names)")
        sys.exit(1)
    return illumina_data[wildcards.sample]


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


def get_raw_fastq(wildcards):
    filename = manager.samples[wildcards.sample][0]
    if filename.endswith(".bam"):
        return f"{wildcards.sample}/bam_to_fastq/{wildcards.sample}.fastq"
    return filename


def get_fastq(wildcards):
    """Convert bam to fastq format."""
    if config["fastp"]["do"]:
        return f"{wildcards.sample}/fastp/{wildcards.sample}.fastq.gz"
    else:
        return get_raw_fastq(wildcards)


def get_corrected_fastq(wildcards):
    """Use corrected fastq or not"""
    if config["canu_correction"]["do"]:
        return f"{wildcards.sample}/corrected_reads/{wildcards.sample}.fastq"
    return get_fastq(wildcards)


def get_assembler_contigs(wildcards):
    """Get contigs with implemented assembler."""
    return assembler_output[config["assembler"]].format(sample=wildcards.sample)


def get_medaka_consensus_contigs(wildcards):
    """Get contigs polished with medaka."""
    if config["medaka_consensus"]["do"]:
        return f"{wildcards.sample}/medaka_consensus/{wildcards.sample}.polish.fasta"
    return get_assembler_contigs(wildcards)


def get_circlator_contigs(wildcards):
    """Get circularize contig or raw contig."""
    if config["circlator"]["do"]:
        return f"{wildcards.sample}/circlator/{wildcards.sample}.circle.fasta"
    return get_medaka_consensus_contigs(wildcards)


def get_final_contigs(wildcards):
    """Get contigs after ith assembler or circlator"""
    if config["polypolish"]["do"]:
        return f"{wildcards.sample}/polypolish/{wildcards.sample}.polish.fasta"
    return get_circlator_contigs(wildcards)
