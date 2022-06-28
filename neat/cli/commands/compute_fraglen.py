"""
Command line interface for NEAT's compute fragment length function
"""

import argparse

from ...compute_fraglen import compute_fraglen_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a model of the fragment length distribution, derived from real data
    """
    name = "compute_fraglen"
    description = "Generate fragment length model from a BAM or SAM file."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """
        parser.add_argument('-i',
                            type=str,
                            metavar="input_file",
                            required=True,
                            default=None,
                            help="Bam or sam input file.")

        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command

        :param arguments: The namespace with arguments and their values.
        """
        compute_fraglen_runner(arguments.input_file, arguments.output)
