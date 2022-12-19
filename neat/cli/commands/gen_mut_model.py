"""
Command line interface for NEAT's compute mutation function
"""

import argparse

from ...gen_mut_model import compute_mut_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a model of the mutation distribution, derived from real data
    """
    name = "gen-mut-model"
    description = "Generate mutation model from a pickle or BED file and user input."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """

        parser.add_argument('reference',
                            type=str,
                            metavar='reference.fa',
                            help="Reference file for organism in fasta format")

        parser.add_argument('mutations',
                            type=str,
                            metavar='mutation.vcf',
                            help="Mutation file for organism in VCF format")

        parser.add_argument('-b',
                            '--bed',
                            type=str,
                            required=False,
                            default=None,
                            help="Bed file with regions to use in the model")

        parser.add_argument('--outcounts',
                            type=str,
                            required=False,
                            default=None,
                            help="Path to trinucleotide counts file for reference. Note, this "
                                 "file will not be used if a bed file is also used.")

        parser.add_argument('--show-trinuc',
                            action='store_true',
                            required=False,
                            default=False,
                            help='Shows trinucleotide counts, for reference')

        parser.add_argument('--save_trinuc',
                            action='store_true',
                            required=False,
                            default=False,
                            help='Saves trinucleotide counts to a file called <out>.counts')

        parser.add_argument('--human_sample',
                            action='store_true',
                            required=False,
                            default=False,
                            help='Only use numbered chroms, X, Y, and MT. '
                                 'Omit this flag to include all chroms in reference.')

        parser.add_argument('--skip_common',
                            action='store_true',
                            required=False,
                            default=False,
                            help="Includes a list of common variants, if you want to visualize "
                                 "common variants with plot_mut_model.py.")

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
        compute_mut_runner(arguments.reference, arguments.mutations, arguments.bed, arguments.outcounts,
                           arguments.show_trinuc, arguments.save_trinuc, arguments.human_sample,
                           arguments.skip_common, arguments.output, arguments.overwrite)
