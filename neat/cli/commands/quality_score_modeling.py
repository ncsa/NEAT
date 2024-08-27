"""
Command line interface for NEAT's quality score error modeling function.
"""

import argparse

from ...quality_score_modeling import quality_score_model_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a model of the quality score errors, derived from a real BAM file
    """
    name = "model-qual-score"
    description = "Simulate quality scores across a read with a Markov model, given a BAM or SAM file."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """
        parser.add_argument('-b',
                            type=str,
                            metavar="FILE",
                            dest="input_file",
                            nargs='+',
                            required=True,
                            help="bam file to process. To enter a sam input, please convert to a bam file"
                                 "using an external tool, such as samtools.")

        parser.add_argument('-c',
                            type=str,
                            metavar="FILE",
                            dest="output_csv",
                            nargs='+',
                            required=False,
                            help="csv file that is optional to output.")

        parser.add_argument('-p',
                            type=str,
                            metavar="FILE",
                            dest="output_pickle",
                            nargs='+',
                            required=True,
                            help="pickle file to store Markov model quality score predictions across"
                                 "all bases along the read, given user-specified read length.")

        parser.add_argument('--overwrite',
                            required=False,
                            action='store_true',
                            default=False,
                            help="Overwrite previous output files. "
                                 "Default is to throw an error if the file_list already exists.")

        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command

        :param arguments: The namespace with arguments and their values.
        """
        quality_score_model_runner(
            arguments.input_file,
            arguments.output_csv,
            arguments.output_pickle,
            arguments.overwrite,
        )