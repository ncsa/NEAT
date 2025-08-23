"""
Command line interface for parallelized wrapper of NEAT.
"""

import argparse
from pathlib import Path
from typing import List

from .base import BaseCommand
from ...read_simulator.utils.parallelize import main as pipeline_main


class Command(BaseCommand):
    """
    Split the reference, run read-simulator, and stitch outputs together.
    """
    name = "parallel"
    description = "Split the reference, run read-simulator, and stitch outputs together."

    def add_arguments(self, parser: argparse.ArgumentParser) -> None:
        """
        Add the command's arguments to its parser

        :param parser: The parser to add arguments to
        """
        parser.add_argument("-c",
                            "--config",
                            type=Path,
                            required=True,
                            help="NEAT YAML/YML config containing the 'reference:' field")
        parser.add_argument("--outdir",
                            type=Path,
                            required=True,
                            help="Top-level directory to hold splits, per-chunk outputs, and final stitched results")

        split = parser.add_argument_group("splitting options")
        split.add_argument("--by", choices=["contig", "size"], default="contig", help="Split mode")
        split.add_argument("--size", type=int, default=1000000, help="Target chunk size when --by size")
        split.add_argument("--cleanup-splits", action="store_true", help="Delete the 'splits' directory after stitching completes")
        split.add_argument("--reuse-splits", action="store_true", help="Skip splitting and reuse existing YAML/FASTA files in 'splits'")

        sim = parser.add_argument_group("simulation options")
        sim.add_argument("--jobs", type=int, default=2, help="Maximum number of parallel NEAT jobs")
        sim.add_argument("--neat-cmd", default="neat read-simulator", help="Command used to launch the read simulator")

        stitch = parser.add_argument_group("stitching options")
        stitch.add_argument("--samtools", default="samtools", help="Path to samtools executable used by stitch_outputs.py")
        stitch.add_argument("--final-prefix", type=Path, default=Path("stitched/final"), help="Prefix (no extension) for stitched outputs")

        parser.add_argument("--parallel-config",
                            type=Path,
                            help="Optional YAML/JSON file with parallelization settings (jobs, by, size, etc.)")

    def execute(self, arguments: argparse.Namespace) -> None:
        # Optionally overlay values from a parallel-config file
        if arguments.parallel_config and arguments.parallel_config.is_file():
            import json
            import yaml
            ext = arguments.parallel_config.suffix.lower()
            with open(arguments.parallel_config, "r") as fh:
                overrides = yaml.safe_load(fh) if ext in (".yml", ".yaml") else json.load(fh)
            for k, v in overrides.items():
                if hasattr(arguments, k):
                    setattr(arguments, k, v)

        # Build argv for your pipeline's main()
        argv: List[str] = [
            str(arguments.config),
            "--outdir", str(arguments.outdir),
            "--by", arguments.by,
        ]

        if arguments.by == "size":
            argv += ["--size", str(arguments.size)]
        if arguments.cleanup_splits:
            argv += ["--cleanup-splits"]
        if arguments.reuse_splits:
            argv += ["--reuse-splits"]

        argv += [
            "--jobs", str(arguments.jobs),
            "--neat-cmd", arguments.neat_cmd,
            "--samtools", arguments.samtools,
            "--final-prefix", str(arguments.final_prefix),
        ]

        pipeline_main(argv)
