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
            choices=["canu", "hifiasm"],
            help="An assembler in canu, hifiasm",
        )
        pipeline_group.add_argument(
            "--hifi", dest="is_hifi", action="store_true", help="Run assembler with adapted parameters for hifi data"
        )
        pipeline_group.add_argument("--blastdb", dest="blastdb", help="Path to your blast database")
        pipeline_group.add_argument("--lineage", dest="lineage", help="Lineage or path to lineage file for BUSCO")

        pipeline_group = self.add_argument_group("pipeline")

        pipeline_group.add_argument("--TODO", dest="TODO", default=4, type=int)

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
        cfg.input_csv = os.path.abspath(options.input_csv) if options.input_csv else ""
        cfg.is_hifi = options.is_hifi
        cfg.assembler = options.assembler
        if options.blastdb:
            cfg.blast["blastdb"] = options.blastdb
        if options.lineage:
            cfg.busco["lineage"] = options.lineage

    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown(check_input_files=False)


if __name__ == "__main__":
    main()
