"""
Command line interface for NEAT's sequencing error model creation function.
"""

import argparse
import numpy as np

from ...model_sequencing_error import model_seq_err_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a model of the sequencing error, derived from real data
    """
    name = "model-seq-err"
    description = "Generate sequencing error model from a FASTQ, BAM, or SAM file_list."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """
        parser.add_argument('-i',
                            type=str,
                            metavar="FILE",
                            dest="input_files",
                            nargs='+',
                            required=True,
                            help="fastq(.gz) to process. Entering 1 fastqs will set the modeler to "
                                 "single-ended mode. Use -i2 for paired ended models."
                                 "To use a sam/bam input first convert to fastq using an external tool, "
                                 "such as samtools.")

        parser.add_argument('-q',
                            type=int,
                            metavar="OFFSET",
                            dest="quality_offset",
                            required=False,
                            default=33,
                            help="quality score offset [33]")

        parser.add_argument('-Q',
                            type=int,
                            metavar="QUAL_SCORE",
                            dest="quality_scores",
                            nargs='+',
                            required=False,
                            default=42,
                            help="Quality score max. The default 42. The lowest possible score is 1. To used binned"
                                 "scoring, enter a space separated list of scores, e.g., -Q 2 15 23 37")

        parser.add_argument('-m',
                            type=int,
                            metavar="MAX",
                            dest="max_num",
                            required=False,
                            default=np.inf,
                            help="Max number of reads to process [default: all].")

        # parser.add_argument('--pileup',
        #                     type=str,
        #                     metavar="FILE",
        #                     required=False,
        #                     help="Pileup statistics file from running samtools mpileup. Not yet implemented.")
        #
        # parser.add_argument('--plot',
        #                     required=False,
        #                     action='store_true',
        #                     default=False,
        #                     help="Output optional plots (filename based on input name). Not yet implemented.")

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
        model_seq_err_runner(
            arguments.input_files,
            arguments.quality_offset,
            arguments.quality_scores,
            arguments.max_num,
            # arguments.pileup,
            # arguments.plot,
            arguments.overwrite,
            arguments.output_dir,
            arguments.prefix
        )
