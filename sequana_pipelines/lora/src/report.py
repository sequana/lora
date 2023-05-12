import base64
import glob
import json
import os
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import DefaultDict, Dict, Iterator, List, Tuple, Union

from jinja2 import Environment, PackageLoader
from pandas import read_csv

from sequana_pipelines.lora import version

from .enums import BLAST_KEY, BUSCO_KEY, QUAST_KEY, SET_QUAST_KEY
from .utils import get_tools_versions


def create_reports(summary_name: str, lora_name: str, samples: List[str], config: Dict, lora_dir=Path(".")) -> None:
    """Create LORA summary report and the main LORA report.
    :param str summary_name: LORA summary report name.
    :param str lora_name: LORA report name.
    :param list samples: List of sample names.
    :param dict config: Config used by the pipeline.
    :param Path lora_dir: Directory where Lora pipeline was executed.
    """
    # get quast information
    summary_header = QUAST_KEY.copy()
    summary_results = get_quast_information(samples, lora_dir)

    # get busco information
    if config["busco"]["do"]:
        summary_header += BUSCO_KEY
        for sample, busco_result in get_busco_information(samples, lora_dir):
            summary_results[sample] += busco_result

    # get blast and sequana information per contigs
    analysis = defaultdict(lambda: defaultdict(dict))
    if config["blast"]["do"]:
        fill_blast_result(analysis, samples, lora_dir)
    if config["sequana_coverage"]["do"]:
        fill_sequana_coverage(analysis, samples, lora_dir)

    # create summary report
    create_summary(summary_name, lora_name, summary_results, config, lora_dir)

    # create lora report
    env = Environment(loader=PackageLoader("sequana_pipelines.lora.src", "templates"))
    template = env.get_template("lora.html")
    report_output = template.render(
        summary_header=summary_header, summary_results=summary_results, analysis=analysis, version=version
    )
    with open(lora_name, "w") as fout:
        print(report_output, file=fout)


def create_summary(summary_name: str, lora_name: str, quast_info: Dict, config: Dict, lora_dir: Path) -> None:
    """Create the summary lora report.
    :param str summary_name: LORA summary report name.
    :param str lora_name: LORA report name.
    :param dict quast_info: Dict with quast informations by samples.
    :param dict config: Config used by the pipeline.
    :param Path lora_dir: Directory where LORA pipeline was executed.
    """
    lora_dir = Path(lora_dir)
    env = Environment(loader=PackageLoader("sequana_pipelines.lora.src", "templates"))
    template = env.get_template("summary.html")
    report_output = template.render(
        lora_dir=lora_dir,
        lora_report=lora_name,
        coverage_done=config["sequana_coverage"]["do"],
        samples=((sample, quast_result[0]) for sample, quast_result in quast_info.items()),
        version=version,
        rulegraph=(lora_dir / ".sequana" / "rulegraph.svg").read_text(),
        dependencies=get_tools_versions(config),
    )
    with open(summary_name, "w") as fout:
        print(report_output, file=fout)


def get_quast_information(samples: List[str], lora_dir: Path) -> Dict:
    """Get quast information."""
    quast_report = f"{lora_dir}/{{}}/quast/report.tsv"
    quast_results = {}
    with ExitStack() as stack:
        # open all report file
        files = [(sample, stack.enter_context(open(quast_report.format(sample)))) for sample in samples]
        # get lines of interest from all samples
        iter_quast_results = (
            (sample, (value for key, value in (line.rstrip().split("\t") for line in filin) if key in SET_QUAST_KEY))
            for sample, filin in files
        )
        quast_results = {sample: list(_iter_value_to_float(results)) for sample, results in iter_quast_results}
    return quast_results


def get_busco_information(samples: List[str], lora_dir: Path) -> Iterator[Tuple[str, List]]:
    """Get busco information."""
    busco_report = f"{lora_dir}/{{0}}/busco/short_summary_{{0}}.txt"
    with ExitStack() as stack:
        files = [(sample, stack.enter_context(open(busco_report.format(sample)))) for sample in samples]
        for sample, filin in files:
            iter_busco = (line.strip().split()[0] for line in filin if any(key in line for key in BUSCO_KEY))
            yield sample, list(_iter_value_to_float(iter_busco))


def fill_blast_result(
    analysis_dict: DefaultDict[str, DefaultDict[str, Dict]], samples: List[str], lora_dir: Path
) -> None:
    """Get blast results."""
    blast_report = f"{lora_dir}/{{0}}/blast/{{0}}.tsv"
    for sample in samples:
        df = read_csv(blast_report.format(sample), sep="\t", names=BLAST_KEY)
        # blast results are grouped by seqId and stitle. The best bitscore for each couple is kept.
        top_alignments = df.loc[df.groupby(["qseqid", "stitle"], sort=False)["bitscore"].idxmax()].set_index("qseqid")
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


def fill_sequana_coverage(
    analysis_dict: DefaultDict[str, DefaultDict[str, Dict]], samples: List[str], lora_dir: Path
) -> None:
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


def _iter_value_to_float(values) -> Iterator[Union[int, str]]:
    """Generator to transform float value."""
    for value in values:
        try:
            yield int(value)
        except ValueError:
            try:
                yield f"{float(value):.3f}"
            except ValueError:
                yield value
