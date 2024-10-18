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

    # make sure analysis exsits so that lora.html is populated
    analysis = defaultdict(lambda: defaultdict(dict))
    for sample in samples:
        analysis[sample] = defaultdict(lambda: defaultdict(dict))

    # get blast/checkm/sequana information per contigs
    if config["blast"]["do"]:
        fill_blast_result(analysis, samples, lora_dir)
    if config["sequana_coverage"]["do"]:
        fill_sequana_coverage(analysis, samples, lora_dir)
    if config["checkm"]["do"]:
        fill_checkm_result(analysis, samples, lora_dir)
        for sample in samples:
            analysis[sample]["genus"] = config["checkm"]["taxon_rank"]
            analysis[sample]["rank"] = config["checkm"]["taxon_name"]
    for sample in samples:
        if os.path.exists(f"{sample}/bandage/{sample}_graph.png"):
            analysis[sample]["bandage_image"] = f"{sample}/bandage/{sample}_graph.png"

    # we want to sort the results by contig size so let us store the list of contigs
    # sorted alphabetically.
    for sample in samples:
        # for older reports
        if os.path.exists(lora_dir / sample / f"sorted_contigs/{sample}.names.txt"):
            with open(lora_dir / sample / f"sorted_contigs/{sample}.names.txt", "r") as fin:
                data = fin.readlines()
                ctg_names = [x.strip(">\n").split()[0] for x in data if len(x.strip("\n"))]
            analysis[sample]["contig_order"] = ctg_names
        else:
            analysis[sample]["contig_order"] = sorted(list(analysis[sample].keys()))

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
            (sample, (value for key, value in (line.rstrip().split("\t") for line in filin) if key in QUAST_KEY))
            for sample, filin in files
        )
        quast_results = {sample: list(_iter_value_to_float(results)) for sample, results in iter_quast_results}

    # quast sometimes report mapping rate > 100%
    # use QUAST_KEY here and above to make sure we retrive the correct key position
    index_mapped = QUAST_KEY.index("Mapped (%)")
    for sample in samples:
        m = float(quast_results[sample][index_mapped])
        quast_results[sample][index_mapped] = min(m, 100.0)
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
        df = read_csv(blast_report.format(sample), sep="\t", names=BLAST_KEY, low_memory=False)
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
            info = {
                "length": contig_results["data"]["length"],
                "DOC": f"{contig_results['data']['DOC']:.3f}",
                "BOC": f"{contig_results['data']['BOC']:.3f}",
            }
            if contig in analysis_dict[sample]:
                analysis_dict[sample][contig]["coverage"] = info
            else:
                analysis_dict[sample][contig] = {}
                analysis_dict[sample][contig]["coverage"] = info
            # add base64 image
            image_name = os.path.join(os.path.dirname(json_file.name), "coverage.png")
            with open(image_name, "rb") as f:
                analysis_dict[sample][contig]["cov_image"] = base64.b64encode(f.read()).decode("utf-8")


def fill_checkm_result(
    analysis_dict: DefaultDict[str, DefaultDict[str, Dict]], samples: List[str], lora_dir: Path
) -> None:
    """Get checkm information."""
    checkm_report = f"{lora_dir}/{{}}/checkm/results.txt"
    checkm_image = f"{lora_dir}/{{}}/checkm/{{}}.marker_pos_plot.png"

    with ExitStack() as stack:
        # open all report file
        files = [(sample, stack.enter_context(open(checkm_report.format(sample)))) for sample in samples]
        # get lines of interest from all samples

        def _filter(line, sample):
            # the file from checkm is not CSV or tabulated. one pattern is at least 2 spaces between
            # items. we also get rid of last \n
            items = [x.strip() for x in line.split("  ") if len(x.strip())]
            if items[0] == sample:
                return items
            # return None otherwise on purpose

        for sample, filin in zip(samples, files):
            values = [
                value for value in (_filter(line, sample) for line in filin[1].readlines() if _filter(line, sample))
            ]
            values = values[0]
            info = {"Completeness": values[11], "Contamination": values[12], "Heterogenity": values[13]}
            analysis_dict[sample]["checkm"] = info

            # add base64 image
            image_name = checkm_image.format(sample, sample)
            with open(image_name, "rb") as f:
                analysis_dict[sample]["checkm_image"] = base64.b64encode(f.read()).decode("utf-8")


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
