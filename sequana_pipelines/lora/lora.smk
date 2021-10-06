""" LORA pipeline v1.0
"""
import csv
import os

from sequana_pipetools.snaketools import PipelineManagerDirectory, FileFactory

from sequana_pipelines.lora import exceptions

shell.executable('bash')

configfile: "config.yml"


manager = PipelineManagerDirectory('lora', config, schema="schema.yml")
csv_filename = config.get('input_csv')
input_directory = config.get('input_directory')
input_pattern = config.get('input_pattern', '*.bam')
if csv_filename:
    # fill samples raw data using input csv
    with open(csv_filename) as csv_file:
        csv_reader = csv.reader(csv_file, skipinitialspace=True)
        manager.samples = {sample: files for sample, *files in csv_reader}
elif input_directory and os.path.isdir(input_directory):
    # use input directory and pattern
    ff = FileFactory(os.path.join(input_directory, input_pattern))
    manager.samples = {sample: [file] for sample, file in zip(ff.filenames, ff.realpaths)}
else:
    raise exceptions.LoraException("Please add a valid input_csv or input_directory")



rule lora:
    input:
        "multiqc/multiqc_report.html"


include: "rules/common.smk"
include: "rules/ccs.smk"
include: "rules/assembler.smk"
include: "rules/qc.smk"


onsuccess:
    from sequana_pipelines.lora import create_report
    from sequana_pipetools.snaketools import OnSuccessCleaner

    # Create LORA summary report
    create_report(
        "lora_report.html",
        manager.samples,
        busco_done=config['busco']['do'],
        blast_done=config['blast']['do'],
        sequana_done=config['sequana_coverage']['do']
    )

    print("Please open the report lora_report.html or multiqc/multiqc_report.html.")
    shell("rm -f ./samples/*/*.done")
    shell("rm -f ./samples/*/*.log")
    shell("chmod -R g+w .")

    sc = OnSuccessCleaner()
    toremove = config["onsuccess"]["toclean"]
    sc.files_to_remove.append(toremove)
    sc.add_makefile()
    print("Once done, please clean up the directory using\n'make clean'")

onerror:
    from sequana_pipetools.errors import PipeError
    p = PipeError("LORA")
    p.status()
