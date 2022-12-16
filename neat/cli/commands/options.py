"""
Definitions of shared subcommand options.
"""

__all__ = ["output_group"]

from .base import Group

output_group = Group("output", is_mutually_exclusive=True, required=True)
output_group.add_argument(
    "-o",
    "--output",
    dest="output",
    type=str,
<<<<<<< HEAD
    help="Path (including filename prefix) to the output file",
=======
    help="Path (including filename prefix) to the output file_list",
>>>>>>> 5d3b24e267938502bbd68617dafca289f2135e02
    default=None
)
