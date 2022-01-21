import sys
import os
import argparse

from sequana_pipetools.options import SlurmOptions, SnakemakeOptions, InputOptions, GeneralOptions, before_pipeline
from sequana_pipetools.misc import Colors
from sequana_pipetools.info import sequana_epilog, sequana_prolog
from sequana_pipetools import SequanaManager

col = Colors()

NAME = "lora"


class Options(argparse.ArgumentParser):
    def __init__(self, prog=NAME, epilog=None):
        usage = col.purple(sequana_prolog.format(**{"name": NAME}))
        super(Options, self).__init__(
            usage=usage,
            prog=prog,
            description="",
            epilog=epilog,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

        # add a new group of options to the parser
        so = SlurmOptions()
        so.add_options(self)

        # add a snakemake group of options to the parser
        so = SnakemakeOptions(working_directory=NAME)
        so.add_options(self)

        so = InputOptions(input_pattern="*.bam", add_input_readtag=False)
        so.add_options(self)

        so = GeneralOptions()
        so.add_options(self)

        pipeline_group = self.add_argument_group("pipeline_general")
        pipeline_group.add_argument(
            "--input-csv",
            dest="input_csv",
            help="Simple CSV file with the samples names and files. LORA will generate ccs and merge your files."
            " If you do not want to do ccs, you can put only one file for each samples.",
        )
        pipeline_group.add_argument(
            "--assembler",
            dest="assembler",
            default="canu",
            choices=["canu", "hifiasm", "flye"],
            help="An assembler in canu, hifiasm, flye",
        )
        pipeline_group.add_argument(
            "--mode",
            dest="mode",
            default="bacteria",
            choices=["eukaryotes", "bacteria"],
            help="If bacteria, blast, circlator, busco, prokka, sequana_coverage are ON, else these options are OFF",
        )
        pipeline_group.add_argument(
            "--do-correction",
            dest="do_correction",
            action="store_true",
            help="Run canu correction before hifiasm or flye.",
        )
        pipeline_group.add_argument(
            "--do-circlator",
            dest="do_circlator",
            action="store_true",
            help="Run circlator after assembler."
        )
        pipeline_group.add_argument("--blastdb", dest="blastdb", help="Path to your blast database")
        pipeline_group.add_argument("--lineage", dest="lineage", 
            help="""Lineage or path to lineage file for BUSCO. Note that we support only version 5 of the BUSCO lineage.""")

        pipeline_group = self.add_argument_group("ccs")
        pipeline_group.add_argument("--ccs-min-passes", default=3, type=int,
            help="""mini number of passes required to build the CCS. Set to 3 for HIFI quality""")
        pipeline_group.add_argument("--ccs-min-rq",  default=0.7, type=int,
            help="minimum quality required to build the CCS. Set to 0.99 for HIFI quality")



    def parse_args(self, *args):
        args_list = list(*args)
        if "--from-project" in args_list:
            if len(args_list) > 2:
                msg = (
                    "WARNING [sequana]: With --from-project option, "
                    + "pipeline and data-related options will be ignored."
                )
                print(col.error(msg))
            for action in self._actions:
                if action.required is True:
                    action.required = False
        options = super(Options, self).parse_args(*args)
        return options


def main(args=None):
    if not args:
        args = sys.argv

    # whatever needs to be called by all pipeline before the options parsing
    before_pipeline(NAME)

    # option parsing including common epilog
    options = Options(NAME, epilog=sequana_epilog).parse_args(args[1:])

    # the real stuff is here
    manager = SequanaManager(options, NAME)

    # create the beginning of the command and the working directory
    manager.setup()

    if options.from_project is None:
        # fill the config file with input parameters
        cfg = manager.config.config
        cfg.input_directory = os.path.abspath(options.input_directory)
        cfg.input_pattern = options.input_pattern
        cfg.input_csv = os.path.abspath(options.input_csv) if options.input_csv else ""

        if options.mode == "bacteria":
            cfg.circlator['do'] = True
            cfg.blast['do'] = True
            cfg.busco['do'] = True
            cfg.prokka['do'] = True
            cfg.sequana_coverage['do'] = True
        elif options.mode == "eukaryotes":
            cfg.circlator['do'] = False
            cfg.blast['do'] = False
            cfg.busco['do'] = True
            cfg.busco['options'] += " --metaeuk "
            cfg.prokka['do'] = False
            cfg.sequana_coverage['do'] = True
        else:
            cfg.circlator['do'] = options.do_circlator
            if options.blastdb:
                cfg.blast["blastdb"] = options.blastdb
            if options.lineage:
                cfg.busco["lineage"] = options.lineage
            


        cfg.canu_correction['do'] = options.do_correction
        cfg.assembler = options.assembler

        cfg.ccs['min-rq'] = options.ccs_min_rq
        cfg.ccs['min-passes'] = options.ccs_min_passes



    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown(check_input_files=False)


if __name__ == "__main__":
    main()
