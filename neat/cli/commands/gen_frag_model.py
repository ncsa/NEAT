"""
Command line interface for NEAT's compute fragment length function
"""

import argparse

from ...gen_frag_model import compute_fraglen_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a model of the fragment length distribution, derived from real data
    """
    name = "gen-frag-model"
    description = "Generate fragment length model from a BAM or SAM file."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """
        parser.add_argument('-i',
                            type=str,
                            metavar="input_file",
                            dest="input_file",
                            required=True,
                            default=None,
                            help="Bam or sam input file.")

        parser.add_argument('--min_reads',
                            type=int,
                            metavar="minimum_number_reads",
                            required=False,
                            default=2,
                            help="Minimum number of reads for a fragment length to consider it in the model. The "
                                 "default is 2, to handle smaller datasets. Set to 0 to turn off filtering. "
                                 "For a larger dataset, try 100 and adjust from there.")

        parser.add_argument('--overwrite',
                            action='store_true',
                            required=False,
                            default=False,
                            help="Set this flag to overwrite existing output.")

        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command

        :param arguments: The namespace with arguments and their values.
        """
        compute_fraglen_runner(arguments.input_file, arguments.min_reads, arguments.output, arguments.overwrite)
