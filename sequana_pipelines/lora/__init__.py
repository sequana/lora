import pkg_resources

version = pkg_resources.require("sequana_lora")[0].version

from .src import create_report, BLAST_KEY, exceptions
