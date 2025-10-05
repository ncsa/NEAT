"""
Definitions of shared subcommand options.
"""

__all__ = ["output_group"]

import os

from .base import Group

output_group = Group("output", is_mutually_exclusive=True, required=True)
output_group.add_argument(
    "-o",
    "--output_dir",
    dest="output_dir",
    type=str,
    help="Path to the output directory. Will create if not present.",
    default=os.getcwd()
)

output_group.add_argument(
    "-p",
    "--prefix",
    dest="prefix",
    type=str,
    help="Prefix to use to name files",
    default="neat_sim"
)