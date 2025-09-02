


rule ragtag_scaffold:
    input:
        ctg=get_final_contigs,
        ref=config['reference_file']
    output:
        scaffold="{sample}/ragtag_scaffold/ragtag.scaffold.fasta",
    threads:
        config["ragtag_scaffold"]["threads"]
    log:
        "{sample}/ragtag_scaffold/sequana.log",
    resources:
        **config["ragtag_scaffold"]["resources"],
    container:
        config['apptainers']['ragtag']
    shell:
        """
        ragtag_scaffold.py {input.ctg} {input.ref} -w -o {wildcards.sample}/ragtag_scaffold 2>{log}
        """

