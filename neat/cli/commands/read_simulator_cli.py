"""
Command line interface for NEAT's generate reads function, with
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
<<<<<<< HEAD
        a fasta file containing the inserted variants only (no errors).
    """

    name = "read_simulator_cli"
    description = "Simulate NGS reads dataset (See README for complete description of the config input)."
=======
        a fasta file_list containing the inserted variants only (no errors).
    """

    name = "read-simulator-cli"
    description = "Simulate NGS reads dataset using command line inputs (which get turned into a config)."
>>>>>>> 5d3b24e267938502bbd68617dafca289f2135e02

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to.
        """

<<<<<<< HEAD

        output_group.add_to_parser(parser)
=======
        pass
        # output_group.add_to_parser(parser)
>>>>>>> 5d3b24e267938502bbd68617dafca289f2135e02

    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command.

        :param arguments: The namespace with arguments and their values.
        """

<<<<<<< HEAD
        read_simulator_runner(arguments.config, arguments.output)
=======
        pass
        # read_simulator_runner(arguments.config, arguments.output)
>>>>>>> 5d3b24e267938502bbd68617dafca289f2135e02
