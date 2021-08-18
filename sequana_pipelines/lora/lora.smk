""" LORA pipeline v0.0
"""

from sequana_pipetools.snaketools import PipelineManager

shell.executable('bash')

configfile: "config.yml"

manager = PipelineManager('lora', config, fastq=False, pattern=config['input_pattern'], schema="schema.yml")


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
