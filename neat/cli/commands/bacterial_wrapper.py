"""
Command line interface for NEAT's bacterial wrapper function
"""

import argparse
import subprocess
import os
from pathlib import Path

from ...bacterial_wrapper import bacterial_wrapper
from ...bacterial_wrapper import stitch_all_outputs
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Class that generates wrapped bacterial models
    """
    name = "bacterial-wrapper"
    description = "Generate wrapped bacterial model reads"

    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """

        parser.add_argument('reference',
                            type=str,
                            metavar='reference.fa',
                            help="Reference file for organism in fasta format.")

        parser.add_argument('bacteria_name',
                            type=str,
                            metavar='bacteria_name',
                            help="Name of the bacteria.")

        parser.add_argument(
            "-c", "--config",
            metavar="config",
            type=str,
            required=True,
            help="Path (including filename) to the configuration file for the reference run."
        )

        output_group.add_to_parser(parser)
        
    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command

        :param arguments: The namespace with arguments and their values.
        """
        bacterial_wrapper(arguments.reference, arguments.bacteria_name, arguments.config, arguments.output_dir)

        output_path = Path(arguments.output_dir)
        file_list = os.listdir(output_path / "Regular")  # same file names for both Regular and Wrapped folders

        output_files = []
        
        for file in file_list:
            reg_file_path = output_path / "Regular" / file
            wrap_file_path = output_path / "Wrapped" / file            
            
            if ("vcf" in file):
                subprocess.run(["gzip", "-d", reg_file_path])
                subprocess.run(["gzip", "-d", wrap_file_path])

                file = file[:-3]
                reg_file_path = output_path / "Regular" / file
                wrap_file_path = output_path / "Wrapped" / file   
            
            output_files.append(reg_file_path)
            output_files.append(wrap_file_path)

        stitch_all_outputs(output_files, arguments.output_dir)