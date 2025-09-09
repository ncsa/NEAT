"""
Command line interface for parallelized wrapper of NEAT.
"""

import argparse
from pathlib import Path
from typing import List

from .base import BaseCommand
from ...parallel_read_simulator.parallelize import parallelize_main as pipeline_main


class Command(BaseCommand):
    """
    Split the reference, run read simulator, and stitch outputs together.
    """
    name = "parallel"
    description = (
        "Split the reference, run read-simulator in parallel, and stitch outputs together."
    )

    def add_arguments(self, parser: argparse.ArgumentParser) -> None:
        """
        Register CLI arguments for the parallel read simulator.
        """

        parser.add_argument(
            "-c",
            "--config",
            type=Path,
            required=True,
            help="NEAT YAML/YML config containing the 'reference:' field",
        )

        parser.add_argument(
            "--outdir",
            type=Path,
            required=False,
            default=None,
            help="Top-level directory for splits and stitched results (optional)",
        )

        # Splitting options
        split = parser.add_argument_group("splitting options")
        split.add_argument(
            "--by",
            choices=["contig", "size"],
            default=None,
            help="Split mode",
        )
        split.add_argument(
            "--size",
            type=int,
            default=None,
            help="Target chunk size when --by size",
        )
        split.add_argument(
            "--cleanup-splits",
            action=argparse.BooleanOptionalAction,
            default=None,
            help="Delete the 'splits' directory after stitching completes",
        )
        split.add_argument(
            "--reuse-splits",
            action=argparse.BooleanOptionalAction,
            default=None,
            help="Skip splitting and reuse existing YAML/FASTA files in 'splits'",
        )

        # Simulation options
        sim = parser.add_argument_group("simulation options")
        sim.add_argument(
            "--jobs",
            type=int,
            default=None,
            help="Maximum number of parallel NEAT jobs",
        )
        sim.add_argument(
            "--neat-cmd",
            default=None,
            help="Command used to launch the read simulator (e.g. 'neat read-simulator')",
        )

        # Stitching options
        stitch = parser.add_argument_group("stitching options")
        stitch.add_argument(
            "--samtools",
            default=None,
            help="Path to samtools executable used by stitch_outputs.py",
        )
        stitch.add_argument(
            "--final-prefix",
            type=Path,
            default=None,
            help="Prefix (no extension) for stitched outputs",
        )

        # Optional YAML/JSON describing parallel settings
        parser.add_argument(
            "--parallel-config",
            type=Path,
            help="Optional YAML/JSON file with parallelization settings (jobs, by, size, etc.)",
        )

    def execute(self, arguments: argparse.Namespace) -> None:
        # Optionally overlay values from a parallel-config file
        if arguments.parallel_config and arguments.parallel_config.is_file():
            import json, yaml
            ext = arguments.parallel_config.suffix.lower()
            with open(arguments.parallel_config, "r") as fh:
                overrides = yaml.safe_load(fh) if ext in (".yml", ".yaml") else json.load(fh)
            for k, v in overrides.items():
                if hasattr(arguments, k):
                    setattr(arguments, k, v)

        argv: List[str] = [str(arguments.config)]

        # Only forward flags the user actually set
        if arguments.outdir is not None:
            argv += ["--outdir", str(arguments.outdir)]
        if arguments.by is not None:
            argv += ["--by", arguments.by]
        if arguments.size is not None and arguments.by == "size":
            argv += ["--size", str(arguments.size)]

        # Handle booleans
        if arguments.cleanup_splits is not None:
            argv += ["--cleanup-splits"] if arguments.cleanup_splits else ["--no-cleanup-splits"]
        if arguments.reuse_splits is not None:
            argv += ["--reuse-splits"] if arguments.reuse_splits else ["--no-reuse-splits"]

        # Other parameters
        if arguments.jobs is not None:
            argv += ["--jobs", str(arguments.jobs)]
        if arguments.neat_cmd is not None:
            argv += ["--neat-cmd", arguments.neat_cmd]
        if arguments.samtools is not None:
            argv += ["--samtools", arguments.samtools]
        if arguments.final_prefix is not None:
            argv += ["--final-prefix", str(arguments.final_prefix)]

        pipeline_main(argv)
