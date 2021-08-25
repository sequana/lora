import base64
import glob
import json
import os
from contextlib import ExitStack
from collections import defaultdict

from pandas import read_csv
from jinja2 import Environment, PackageLoader

from .enums import BLAST_KEY, BUSCO_KEY, QUAST_KEY, SET_QUAST_KEY


def create_report(output_name, samples, lora_dir=".", busco_done=True, blast_done=True, sequana_done=True):
    """Create lora summary report.
    :params list samples: List of sample names.
    :params str lora_dir: Directory where Lora pipeline was executed.
    """
    # get quast information
    summary = {"header": QUAST_KEY.copy()}
    summary["results"] = get_quast_information(samples, lora_dir)

    # get busco information
    if busco_done:
        summary["header"] += BUSCO_KEY
        summary["results"] = get_busco_information(summary["results"], lora_dir)

    # get blast and sequana information per contigs
    analysis = defaultdict(lambda: defaultdict(dict))
    if blast_done:
        analysis = get_blast_result(analysis, samples, lora_dir)
    if sequana_done:
        analysis = get_sequana_coverage(analysis, samples, lora_dir)

    # create html
    env = Environment(loader=PackageLoader("sequana_pipelines.lora.src", "templates"))
    template = env.get_template("lora.html")
    report_output = template.render(summary=summary, analysis=analysis)
    with open(output_name, "w") as fout:
        print(report_output, file=fout)


def get_quast_information(samples, lora_dir):
    """Get quast information."""
    quast_report = f"{lora_dir}/{{}}/quast/report.tsv"
    with ExitStack() as stack:
        # open all report file
        files = [(sample, stack.enter_context(open(quast_report.format(sample)))) for sample in samples]
        # get lines of interest from all samples
        iter_quast_results = (
            (sample, (value for key, value in (line.rstrip().split("\t") for line in filin) if key in SET_QUAST_KEY))
            for sample, filin in files
        )
        quast_results = {
            sample: [result for result in _iter_value_to_float(results)] for sample, results in iter_quast_results
        }
    return quast_results


def get_busco_information(summary, lora_dir):
    """Get busco information."""
    busco_report = f"{lora_dir}/{{0}}/busco/short_summary_{{0}}.txt"
    with ExitStack() as stack:
        files = [(sample, stack.enter_context(open(busco_report.format(sample)))) for sample in summary.keys()]
        for sample, filin in files:
            iter_busco = (next(iter(line.strip().split())) for line in filin if any(key in line for key in BUSCO_KEY))
            busco_subresult = [result for result in _iter_value_to_float(iter_busco)]
            summary[sample] += busco_subresult
    return summary


def get_blast_result(analysis_dict, samples, lora_dir):
    """Get blast results."""
    blast_report = f"{lora_dir}/{{0}}/blast/{{0}}.tsv"
    for sample in samples:
        df = read_csv(blast_report.format(sample), sep="\t", names=BLAST_KEY, index_col=0)
        # blast results are grouped by seqId and stitle. The best bitscore for each couple is kept.
        top_alignments = df.groupby([df.index, "stitle"], sort=False).max("bitscore").reset_index(level=-1)
        for contig in top_alignments.index.unique():
            try:
                analysis_dict[sample][contig]["blast"] = (
                    top_alignments.loc[contig].reset_index(drop=True).head(5).to_dict("records")
                )
            except TypeError:
                # There are only one alignment
                analysis_dict[sample][contig]["blast"] = [top_alignments.loc[contig].to_dict()]
            except KeyError:
                # No alignment for this contig
                pass
    return analysis_dict


def get_sequana_coverage(analysis_dict, samples, lora_dir):
    """Get sequana results."""
    sequana_report = f"{lora_dir}/{{0}}/sequana_coverage/*/sequana_summary_coverage.json"
    # iter over sequana json file
    iter_json = (
        (sample, sequana_file) for sample in samples for sequana_file in glob.iglob(sequana_report.format(sample))
    )
    with ExitStack() as stack:
        files = [(sample, stack.enter_context(open(sequana_file))) for sample, sequana_file in iter_json]
        for sample, json_file in files:
            contig_results = json.load(json_file)
            contig = contig_results["data"]["chrom_name"]
            # get general mapping information
            analysis_dict[sample][contig]["coverage"] = {
                "length": contig_results["data"]["length"],
                "DOC": f"{contig_results['data']['DOC']:.3f}",
                "BOC": f"{contig_results['data']['BOC']:.3f}",
            }
            # add base64 image
            image_name = os.path.join(os.path.dirname(json_file.name), "coverage.png")
            with open(image_name, "rb") as f:
                analysis_dict[sample][contig]["cov_image"] = base64.b64encode(f.read()).decode("utf-8")
    return analysis_dict


def _iter_value_to_float(values):
    """Generator to transform float value."""
    for value in values:
        try:
            yield int(value)
        except ValueError:
            try:
                yield f"{float(value):.3f}"
            except ValueError:
                yield value
