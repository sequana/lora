#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Website:       https://github.com/sequana/lora
#  Documentation: http://sequana.readthedocs.io
#  Documentation: https://github.com/sequana/lora/README.rst
##############################################################################
""" LORA (LOng Read Assembly) pipeline"""
import csv
import os

from sequana_pipetools.snaketools import PipelineManagerDirectory, FileFactory, modules
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


localrules: lora, rulegraph

rule lora:
    input:
        ".sequana/rulegraph.svg",
        "multiqc/multiqc_report.html"


include: "rules/common.smk"
include: "rules/ccs.smk"
include: "rules/assembler.smk"
include: "rules/polish.smk"
include: "rules/qc.smk"
include: "rules/utils.smk"

onsuccess:
    from sequana_pipelines.lora import create_reports
    from sequana import logger

    logger.setLevel("INFO")
    manager.teardown()

    # Create LORA report and summary
    create_reports("summary.html", "lora.html", manager.samples, config)

    print("Please open the report summary.html or multiqc/multiqc_report.html.")
    shell("rm -f ./samples/*/*.done")
    shell("rm -f ./samples/*/*.log")
    shell("chmod -R g+w .")
    shell("rm -rf rulegraph")

onerror:
    from sequana_pipetools.errors import PipeError
    p = PipeError("lora")
    p.status()
