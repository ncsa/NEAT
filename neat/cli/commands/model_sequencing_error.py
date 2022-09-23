"""
Command line interface for NEAT's sequencing error model creation function.
"""

import argparse

from ...model_sequencing_error import model_seq_err_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a model of the sequencing error, derived from real data
    """
    name = "model-seq-err"
    description = "Generate sequencing error model from a FASTQ, BAM, or SAM file."

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
                            help="fastq(.gz) or sam/bam file")

        parser.add_argument('-q',
                            type=int,
                            metavar="quality_offset",
                            dest="quality_offset",
                            required=False,
                            default=33,
                            help="quality score offset")

        parser.add_argument('-Q',
                            type=int,
                            metavar="quality_scores",
                            dest="quality_scores",
                            nargs='+',
                            required=False,
                            default=[0, 12, 24, 36],
                            help="Quality score bins. Enter as a list for binned scores, "
                                 "or enter a single maximum score for a full range")

        parser.add_argument('-n',
                            type=int,
                            metavar="max_num",
                            dest="max_num",
                            required=False,
                            default=-1,
                            help="Max number of reads to process [all].")

        parser.add_argument('-s',
                            type=int,
                            metavar="num_iterations",
                            dest="num_iterations",
                            required=False,
                            default=1000000,
                            help="Number of simulation iterations (higher number should improve accuracy).")

        parser.add_argument('--plot',
                            required=False,
                            action='store_true',
                            default=False,
                            help="Output optional plots (filename based on input name)")

        parser.add_argument('--overwrite',
                            required=False,
                            action='store_true',
                            default=False,
                            help="Overwrite previous output files. "
                                 "Default is to throw an error if the file already exists.")

        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command

        :param arguments: The namespace with arguments and their values.
        """
        model_seq_err_runner(
            arguments.input_file,
            arguments.quality_offset,
            arguments.quality_scores,
            arguments.max_num,
            arguments.num_iterations,
            arguments.plot,
            arguments.overwrite,
            arguments.output
        )
