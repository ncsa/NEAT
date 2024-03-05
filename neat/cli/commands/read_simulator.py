"""
Command line interface for NEAT's generate reads function
"""

import argparse

from ...read_simulator import read_simulator_runner
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates a Dataset of simulated NGS reads. NEAT first generates a set of variants to insert, then
        generates fragments to sample from and adds sampling errors to the reads.

    Optional outputs include a golden vcf showing all inserted true variants (i.e., not the simulated errors),
        a golden BAM showing a fast alignment against the region where the read came from,
        a fasta file containing the inserted variants only (no errors).
    """

    name = "read-simulator"
    description = "Simulate NGS reads dataset (See README for complete description of the config input)."

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to.
        """
        parser.add_argument(
            "-c", "--config",
            metavar="config",
            type=str,
            required=False,
            help="Path (including filename) to the configuration file for this run."
        )

        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command.

        :param arguments: The namespace with arguments and their values.
        """
        read_simulator_runner(arguments.config, arguments.output)
