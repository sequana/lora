"""Utils rules"""

rule rulegraph:
    input:
        workflow.snakefile
    output:
        svg = ".sequana/rulegraph.svg",
    params:
        configname = "config.yml",
        required_local_files = ["schema.yml", "rules"]
    wrapper:
        "main/wrappers/rulegraph"
