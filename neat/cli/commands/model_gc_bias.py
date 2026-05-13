"""
Command line interface for NEAT's compute GC bias function
"""

import argparse

from ...model_gc_bias import compute_gc_bias_runner
from .base import BaseCommand
from .options import output_group

class Command(BaseCommand):
    """
    Class that generates a model of GC bias, derived from real data
    """
    name = "model-gc-bias"
    description = "Generate GC bias model from a BAM file and a reference FASTA."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser
        """
        parser.add_argument('-i',
                            type=str,
                            metavar="input_bam",
                            dest="input_bam",
                            required=True,
                            help="BAM input file.")

        parser.add_argument('-r',
                            type=str,
                            metavar="reference",
                            dest="reference",
                            required=True,
                            help="Reference FASTA file.")

        parser.add_argument('-w',
                            type=int,
                            metavar="window_size",
                            dest="window_size",
                            required=False,
                            default=100,
                            help="Window size for GC bias calculation. Default is 100.")

        parser.add_argument('--overwrite',
                            action='store_true',
                            required=False,
                            default=False,
                            help="Set this flag to overwrite existing output.")
        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command
        """
        compute_gc_bias_runner(
            arguments.input_bam,
            arguments.reference,
            arguments.output_dir,
            arguments.prefix,
            arguments.window_size,
            arguments.overwrite
        )
