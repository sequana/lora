import pkg_resources

version = pkg_resources.require("lora")[0].version

from .src import create_reports, BLAST_KEY, exceptions
