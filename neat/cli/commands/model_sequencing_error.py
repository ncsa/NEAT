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
                            dest="input_file",
                            required=True,
                            help="fastq(.gz) to process. Entering 1 fastqs will set the modeler to "
                                 "single-ended mode. Use -i2 for paired ended models."
                                 "To use a sam/bam input first convert to fastq using an external tool, "
                                 "such as samtools.")

        parser.add_argument('-i2',
                            type=str,
                            metavar="FILE2",
                            dest="input_file2",
                            help="fastq(.gz) file_list to process. Entering a value for -i2 will set the modeler to "
                                 "paired-ended mode, modeling each strand separately.")

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
                            default=[2, 11, 25, 37],
                            help="Quality score bins. Enter as a list for binned scores, "
                                 "or enter a single maximum score for a full range (i.e., entering 42 will give error"
                                 "scores from 1-42). The default is binned quality scores: [2, 11, 24, 37]. Note that"
                                 "using quality score bins on an unbinned fastq will result in a binned model, at the"
                                 "cost of some inaccuracy.")

        parser.add_argument('-m',
                            type=int,
                            metavar="MAX",
                            dest="max_num",
                            required=False,
                            default=np.inf,
                            help="Max number of reads to process [all].")

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
            arguments.input_file,
            arguments.input_file2,
            arguments.quality_offset,
            arguments.quality_scores,
            arguments.max_num,
            # arguments.pileup,
            # arguments.plot,
            arguments.overwrite,
            arguments.output
        )
